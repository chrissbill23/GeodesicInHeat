
#include <igl/octree.h>
#include <igl/knn.h>
#include <iostream>
#include <igl/cotmatrix.h>

using namespace Eigen;
using namespace std;

/*void distances(const MatrixXd &V, MatrixXd &D, VectorXd& mind){
  int rows = V.rows();
  D = MatrixXd::Zero(rows,rows);
  mind.setConstant(rows, 1000000000.0);
  for(int i = 0; i < rows; ++ i){
      for(int j = i+1; j < rows; ++ j){
	double dist = (V.row(i) - V.row(j)).norm();
        D(i,j) = dist;
        D(j,i) = dist;
        if(dist < mind(i))
            mind(i) = dist;
        if(dist < mind(j))
            mind(j) = dist;
    
      }
  }

}
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A, MatrixXd &L, double maxthreshold){
  cout << " NAIVE LAPLACIAN COMPUTATION "<<endl;
  int rows = V.rows();
  L = MatrixXd::Zero(rows,rows);
  A = MatrixXi::Zero(rows,rows);
  VectorXd mind;
  distances(V,D,mind);
  for(int i = 0; i < rows; ++ i){
      for(int j = i+1; j < rows; ++ j){
	double dist = (V.row(i) - V.row(j)).norm();
        D(i,j) = dist;
        D(j,i) = dist;
        if(dist < mind(i))
            mind(i) = dist;
        if(dist < mind(j))
            mind(j) = dist;
      }
  }

  for(int i = 0; i < rows; ++i){
      for(int j = i+1; j < rows; ++j){
	  if(D(i,j) <= mind(i)*maxthreshold || D(i,j) <= mind(j)*maxthreshold){
             A(i,j) = 1;
             A(j,i) = 1;
             L(i,i) += 1;
             L(j,j) += 1;
             L(i,j) = -1;
             L(j,i) = -1;
          }
      }
  }
}*/
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A,MatrixXi &F, MatrixXd &L, int k){
  cout << " LAPLACIAN COMPUTATION USING OCTREE"<<endl;
  MatrixXi I;
  vector<vector<int > > O_I;
  MatrixXi O_CH;
  MatrixXd O_CN;
  VectorXd O_W;
  igl::octree(V,O_I,O_CH,O_CN,O_W);
  igl::knn(V,k,O_I,O_CH,O_CN,O_W,I);

  int rows = V.rows();
  D = MatrixXd::Zero(rows,rows);
  A = MatrixXi::Zero(rows,rows);
  L = MatrixXd::Zero(rows,rows);
  F = MatrixXi::Zero(rows,I.cols());
  for(int i = 0; i < rows; ++ i){
      for(int j = 1; j < I.cols(); ++ j){
        int y = I(i,j), y0 = I(i,j-1);
        F(i,j) = y;
        F(i,j-1) = y0;
	double dist = (V.row(y0) - V.row(y)).norm();
        D(y0,y) = dist;
        D(y,y0) = dist;
        if(A(y0,y) == 0){
	   L(y0,y0) += 1;
	   L(y,y) += 1;
	   L(y0,y) = -1;
	   L(y,y0) = -1;
           A(y0,y) = 1;
           A(y,y0) = 1;
        }
      }
      if(I.cols() > 1 && A(i,I(i,I.cols()-1)) == 0 ){
        int y = I(i,I.cols()-1);
	double dist = (V.row(i) - V.row(y)).norm();
        D(i,y) = dist;
        D(y,i) = dist;
        A(i,y) = 1;
        A(y,i) = 1;
        L(i,i) += 1;
        L(y,y) += 1;
        L(i,y) = -1;
        L(y,i) = -1;
      }
  }
}
void laplacianTriMesh(const MatrixXd &V, const MatrixXi &F, MatrixXd &D, MatrixXi &A, MatrixXd &L){
  cout << " LAPLACIAN COMPUTATION USING MESH"<<endl;
  int rows = F.rows();
  D = MatrixXd::Zero(rows,rows);
  A = MatrixXi::Zero(rows,rows);
  SparseMatrix<double> L2;
  igl::cotmatrix(V,F,L2);
  L = L2;
  for(int i = 0; i < rows; ++ i){
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


