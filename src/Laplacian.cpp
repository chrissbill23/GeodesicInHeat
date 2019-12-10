
#include <igl/octree.h>
#include <igl/knn.h>
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/fit_plane.h>
#include <igl/centroid.h>
//#include <igl/copyleft/cgal/point_areas.h>


using namespace Eigen;
using namespace std;

int factorial(int x){
    if(x <= 1)
      return 1;
    return x * factorial(x-1);
}
double gaussianweight(const RowVectorXd& x0, const RowVectorXd& x1){
   double dist = (x0 - x1).norm();
   return pow(4*M_PI*0.25*exp(dist/(4*0.5)),-1);
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
         L.coeffRef(i,k) = L.coeffRef(k,i) = a;
         L.coeffRef(k,k) -= a;
      }
      L.coeffRef(i,i) -= sum;
      M.coeffRef(i,i) = Areas(i);
      M_inv.coeffRef(i,i) = Areas(i) == 0 ? 0 : 1./ (Areas(i));
  }
  L = M_inv * L;
}
void laplacianTriMesh(const MatrixXd &V, const MatrixXi &F, MatrixXd &D, MatrixXi &A, SparseMatrix<double> &L){
  cout << " LAPLACIAN COMPUTATION USING MESH"<<endl;
  int rows = F.rows();
  D = MatrixXd::Zero(rows,rows);
  A = MatrixXi::Zero(rows,rows);
  igl::cotmatrix(V,F,L);
  for(size_t i = 0; i < rows; ++ i){
      int v1 = F(i,0), v2 = F(i,1), v3 = F(i,2);
      const RowVectorXd& ve1 = V.row(v1);
      const RowVectorXd& ve2 = V.row(v2);
      const RowVectorXd& ve3 = V.row(v3);

      double dist1 = (ve1 - ve2).norm();
      double dist2 = (ve1 - ve3).norm();
      double dist3 = (ve3 - ve2).norm();

      D(v1,v2) = D(v2,v1) = dist1; 
      D(v1,v3) = D(v3,v1) = dist2; 
      D(v3,v2) = D(v2,v3) = dist3;
      
      A(v1,v2) = A(v2,v1) = 1; 
      A(v1,v3) = A(v3,v1) = 1; 
      A(v3,v2) = A(v2,v3) = 1;
  }
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
void leastSquareInterp(const MatrixXd &V, const VectorXd &f, VectorXd &c, int m){
       int d = V.cols();
       int components = factorial(m+d)/(factorial(m)*factorial(d)); 
}



