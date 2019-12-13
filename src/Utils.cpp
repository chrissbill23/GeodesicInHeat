
#include <igl/octree.h>
#include <igl/knn.h>
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/fit_plane.h>
#include <igl/centroid.h>
#include <igl/adjacency_list.h>
#include <igl/dijkstra.h>

#include <chrono>
#include <list>
#include <set>
#include <vector>

//#include <igl/copyleft/cgal/point_areas.h>


using namespace std;
using namespace Eigen;

int factorial(int x){
    if(x <= 1)
      return 1;
    return x * factorial(x-1);
}
double gaussianweight(const RowVectorXd& x0, const RowVectorXd& x1){
   double dist = (x0 - x1).norm();
   return pow(4*M_PI*0.25*exp(dist/(4*0.5)),-1);
}
double wendlandweight(const RowVectorXd& x0, const RowVectorXd& x1){
   double d = (x0 - x1).norm();
   double h = d*1.5;
   h = h <= 0 ? 1. : h;
   return pow(1-d/h,4)*(4*d/h+1);
}
void computeTriMeshFromCentroid(const MatrixXd &V, const RowVectorXi &F, const RowVectorXd& c, MatrixXd &V2, MatrixXi &F2){
    F2.setZero(F.cols(), 3);
    V2.setZero(F.cols()+1, 3);
    V2.row(0) = c;
    for(size_t i = 0; i < F.cols()-1; ++i){
       V2.row(i+1) = V.row(F(i));
       V2.row(i+2) = V.row(F(i+1));
       F2(i,1) = i+1;
       F2(i,2) = i+2;
    }
       F2(F.cols()-1,1) = F.cols();
       F2(F.cols()-1,2) = 1;
}
void computeMeanCentroidOfFace(const MatrixXd &V, const RowVectorXi &F, RowVectorXd& c){
    c.setZero(3);
    for(size_t i = 0; i < F.cols(); ++i){
        c += V.row(F(i));
    }
    c/=double(F.cols());
}

VectorXd computeAreaPolygon(const MatrixXd &V, const MatrixXi &F){
  VectorXd area(F.rows()); area.setZero();
  if( F.cols() > 2){
	  if(F.cols() == 3 || F.cols() == 4){
	    igl::doublearea(V,F,area);
	  } else {
	      for(size_t i = 0; i < F.rows(); ++i){
		 RowVectorXd c;
		 computeMeanCentroidOfFace(V, F.row(i), c);
                 MatrixXd V2; MatrixXi F2;
                 computeTriMeshFromCentroid(V, F.row(i),c, V2, F2);
                 VectorXd tmp;
                 igl::doublearea(V2,F2,tmp);
                 area(i) = tmp.sum();
	      }
	  }
  }
  return area*0.5;
}

void normalPlanePointClouds(const MatrixXd &V, const MatrixXi &I, MatrixXd &N, MatrixXd &P){
     N = MatrixXd::Zero(V.rows(),3);
     P = MatrixXd::Zero(V.rows(),3);
     for(size_t i = 0; i < V.rows(); ++i){
        MatrixXd Vtemp(I.cols(),V.cols());
        RowVector3d Ntmp, Ptmp;
        for(size_t j = 0; j < I.cols(); ++j)
            Vtemp.row(j) = V.row(I(i,j));
        igl::fit_plane(Vtemp, Ntmp, Ptmp);
        N.row(i) = Ntmp;
        P.row(i) = Ptmp;
     }
}
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A,MatrixXi &F, 
                    SparseMatrix<double> &L, int k, MatrixXd& N, MatrixXd& P, SparseMatrix<double>& M, SparseMatrix<double>& M_inv){
  cout << " LAPLACIAN COMPUTATION USING OCTREE"<<endl;
  MatrixXi I;
  vector<vector<int > > O_I;
  MatrixXi O_CH;
  MatrixXd O_CN;
  VectorXd O_W;
  igl::octree(V,O_I,O_CH,O_CN,O_W);
  igl::knn(V,k,O_I,O_CH,O_CN,O_W,I);
  
  normalPlanePointClouds(V, I,N,P);
  //igl::copyleft::cgal::point_areas(V,I,N,Areas);

  int rows = V.rows();
  VectorXd Areas(rows);

   for(size_t i = 0; i < rows; ++ i){
      Areas(i) = computeAreaPolygon(V, I.row(i))(0);
   }

  D = MatrixXd::Zero(rows,rows);
  A = MatrixXi::Zero(rows,rows);
  L = SparseMatrix<double>(rows,rows);
  F = MatrixXi::Zero(rows,I.cols());

  M = SparseMatrix<double>(rows,rows);
  M_inv = SparseMatrix<double>(rows,rows);
  for(size_t i = 0; i < rows; ++ i){
      for(size_t j = 1; j < I.cols(); ++ j){
        size_t y = I(i,j), y0 = I(i,j-1);
        F(i,j) = y;
        F(i,j-1) = y0;
	double dist = (V.row(y0) - V.row(y)).norm();
        D(y0,y) = dist;
        D(y,y0) = dist;
        if(A(y0,y) == 0){
           A(y0,y) = 1;
           A(y,y0) = 1;
        }
      }
      if(I.cols() > 1 && A(i,I(i,I.cols()-1)) == 0 ){
        size_t y = I(i,I.cols()-1);
	double dist = (V.row(i) - V.row(y)).norm();
        D(i,y) = dist;
        D(y,i) = dist;
        A(i,y) = 1;
        A(y,i) = 1;
      }
      double sum = 0;
      for(size_t k = i+1; k < rows; ++k){
         double a = Areas(i) * Areas(k) * gaussianweight(V.row(i), V.row(k));
         sum += a;
         L.coeffRef(i,k) = a; L.coeffRef(k,i) = a;
         L.coeffRef(k,k) -= a;
      }
      L.coeffRef(i,i) -= sum;
      M.coeffRef(i,i) = Areas(i);
      M_inv.coeffRef(i,i) = Areas(i) == 0 ? 0 : 1./ (Areas(i));
  }
  L = M_inv * L;
}
void laplacianTriMesh(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double> &L){
  cout << " LAPLACIAN COMPUTATION USING MESH"<<endl;
  igl::cotmatrix(V,F,L);
}



void massMatrix(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double>& M, SparseMatrix<double>& M_inv, int type){
     if(F.cols() == 3 || F.cols() == 4){
        SparseMatrix<double> tmp;
        igl::massmatrix(V,F,(igl::MassMatrixType)type, tmp);
        M = tmp;
     } else{
       M = SparseMatrix<double>(V.rows(), V.rows());
       auto areas = computeAreaPolygon(V, F);
       for(size_t i = 0; i < V.rows(); ++i){
          double tot = 0;
          for(size_t j = 0; j < F.rows(); ++j){
             for(size_t k = 0; k < F.cols(); ++k){
                 if(i == F(j,k)){
                    tot += areas(j);
                    break;
                 }
             }
          }
          tot /= 3;
          M.insert(i,i) = tot;
       }
    }
    M_inv = SparseMatrix<double>(V.rows(), V.rows());
    for(size_t i = 0; i < V.rows(); ++i){
        M_inv.insert(i,i) = 1./M.coeffRef(i,i);
    }
}
double cotangent(const Vector3d& x0, const Vector3d& x1){
       return x0.dot(x1)/x0.cross(x1).norm();
}
double cotangent(const Vector3d& x0, const Vector3d& x1, const Vector3d& x2){
       return x0.norm()/x1.dot(x2);
}
void buildCoefMatrix(const RowVector3d &v, MatrixXd &A, RowVectorXd &b){
     int x = 0, y = 1, z = 2;
     A.setZero(20, 20);
     A.row(0) << 1 , v(x), v(y) , v(z) , v(x)*v(x) , v(x)*v(y) , v(y)*v(y) , v(y)*v(z) , v(z)*v(z) , v(x)*v(z) , v(x)*v(x)*v(x),
        v(x)*v(x)*v(y) , v(x)*v(y)*v(z) , v(y)*v(y)*v(z) , v(y)*v(y)*v(y) , v(y)*v(z)*v(z) , v(x)*v(z)*v(z) , v(x)*v(x)*v(z) , v(x)*v(y)*v(y) , v(z)*v(z)*v(z);
     b = A.row(0);
     for(size_t i = 1; i < 20; ++i){
         A.row(i) = A.row(0)*A(0,i);
     } 
}
void leastSquareInterp3D3Degree(const MatrixXd &V, const MatrixXi &neights, const VectorXd &f, MatrixXd &c){
   c.setZero(V.rows(), 20); 
   for(size_t i = 0; i < V.rows(); ++i){
      RowVectorXd b = RowVectorXd::Zero(20);
      MatrixXd A; A.setZero(20, 20); 
      for(size_t j = 0; j < neights.cols(); ++j){
           MatrixXd tmp; RowVectorXd btmp;
           buildCoefMatrix(V.row(neights(i,j)), tmp, btmp);
           A += tmp * wendlandweight(V.row(i), V.row(neights(i,j))); b += btmp*f(neights(i,j));
      }
      LDLT<MatrixXd> solver(A);
      c.row(i) = solver.solve(b.transpose());
    }
}
void gradientLeastSquared3D3Degree(const MatrixXd &v, const MatrixXd &C, MatrixXd& G){
     if(C.cols() == 20){
        G.setZero(v.rows(), 3); int x = 0, y = 1, z =2;
        for(size_t i = 0; i < v.rows(); ++i){
           G(i,x) =  C(i,1) + v(i,x)*2.*C(i,4) + v(i,y)*C(i,5)  + v(i,z)*C(i,9) + v(i,x)*v(i,x)*3.*C(i,10) + v(i,x)*2.*v(i,y)*C(i,11) + v(i,y)*v(i,z)*C(i,12) + v(i,z)*v(i,z)*C(i,16) + v(i,x)*2.*v(i,z)*C(i,17) + v(i,y)*v(i,y)*C(i,18);
           G(i,y) =   C(i,2) + v(i,x)*C(i,5) + v(i,y)*2.*C(i,6) + v(i,z)*C(i,7) +  v(i,x)*v(i,x)*C(i,11) + v(i,x)*v(i,z)*C(i,12) + v(i,y)*2.*v(i,z)*C(i,13) + v(i,y)*3.*v(i,y)*C(i,14) + v(i,z)*v(i,z)*C(i,15)  +  v(i,x)*2.*v(i,y)*C(i,18);
           G(i,z) = C(i,3) + v(i,y)*C(i,7) + v(i,z)*2.*C(i,8) + v(i,x)*C(i,9) + v(i,x)*v(i,y)*C(i,12) + v(i,y)*v(i,y)*C(i,13)  + v(i,y)*2.*v(i,z)*C(i,15) + v(i,x)*2.*v(i,z)*C(i,16) + v(i,x)*v(i,x)*C(i,17)  + v(i,z)*3.*v(i,z)*C(i,19);
        }
     }
}
double evaluateLeastSquared3D3Degree(const RowVector3d &v, const RowVectorXd &C){
      int x = 0, y = 1, z = 2;
      return C(0) + v(x)*C(1) + v(y)*C(2)  + v(z)*C(3)  + v(x)*v(x)*C(4)  + v(x)*v(y)*C(5)  + v(y)*v(y)*C(6)  + v(y)*v(z)*C(7)  + v(z)*v(z)*C(8)  + v(x)*v(z)*C(9)  + v(x)*v(x)*v(x)*C(10)  +
        v(x)*v(x)*v(y)*C(11)  + v(x)*v(y)*v(z)*C(12)  + v(y)*v(y)*v(z)*C(13)  + v(y)*v(y)*v(y)*C(14)  + v(y)*v(z)*v(z)*C(15)  + v(x)*v(z)*v(z)*C(16)  + v(x)*v(x)*v(z)*C(17)  + v(x)*v(y)*v(y)*C(18)  + v(z)*v(z)*v(z)*C(19) ;
}



int dijkstra(const list<int> &sources, const MatrixXd &targets,const MatrixXi& F, VectorXd &min_distance){
    chrono::time_point<chrono::high_resolution_clock> start =  chrono::high_resolution_clock::now();
    if(sources.size() > 0){
        vector<vector<int> > A; vector<double > D;
        igl::adjacency_list(F,A);
        min_distance.setConstant(targets.rows(), numeric_limits<VectorXd::Scalar>::infinity());
        list <int>::const_iterator it;
        VectorXi previous;
        for(it = sources.begin(); it != sources.end(); ++it){
            VectorXi visited = VectorXi::Zero(targets.rows());
            int index = *it; D.clear();
           for(size_t j = 0; j < targets.rows(); ++j)
               D.push_back((targets.row(index)-targets.row(j)).norm());
           int count = targets.rows();
           visited(index) = 1;
           min_distance(index) = 0.;
           while(visited.sum() < targets.rows() && count > 0){
               set<int > V;
               for(size_t i = 0; i < targets.rows(); ++i)
                   if(visited(i) == 0)
                       V.insert(i);
               VectorXd min_distancetmp;
               int min = igl::dijkstra(index,V,A, D,min_distancetmp,previous);
               visited(min) = 1;
               double val = min_distancetmp(min);
               if(val == 0.)
                   min_distancetmp(min) = D[min];
               if(min_distance(min) == numeric_limits<VectorXd::Scalar>::infinity())
                  min_distance(min) = min_distancetmp(min);
               else {
                  min_distance(min) += min_distancetmp(min);
               }
               --count;
            }
        }
        min_distance /= sources.size(); 
    }
    auto stop = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
}

