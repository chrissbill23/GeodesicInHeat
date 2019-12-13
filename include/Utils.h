
#ifndef UTILS_H
#define UTILS_H

#include <list>
#include <set>

void distances(const MatrixXd &V, MatrixXd &D, MatrixXd &A,VectorXd& mind);
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A,MatrixXi &F, SparseMatrix<double> &L, int k, 
                     MatrixXd& N, MatrixXd& P, SparseMatrix<double>& M, SparseMatrix<double>& M_inv);
void laplacianTriMesh(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double> &L);
void computeMeanCentroidOfFace(const MatrixXd &V, const RowVectorXi &F, RowVectorXd& c);
VectorXd computeAreaPolygon(const MatrixXd &V, const MatrixXi &F);
void massMatrix(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double>& M, SparseMatrix<double>& M_inv, int type = 3);
double cotangent(const Vector3d& x0, const Vector3d& x1);
double cotangent(const Vector3d& x0, const Vector3d& x1, const Vector3d& x2);
void leastSquareInterp3D3Degree(const MatrixXd &V, const MatrixXi &neights, const VectorXd &f, MatrixXd &c);
int factorial(int x);
double gaussianweight(const RowVectorXd& x0, const RowVectorXd& x1);
void gradientLeastSquared3D3Degree(const MatrixXd &v,  const MatrixXd &C, MatrixXd& G);
double evaluateLeastSquared3D3Degree(const RowVector3d &v, const RowVectorXd &C);
double wendlandweight(const RowVectorXd& x0, const RowVectorXd& x1);
int dijkstra(const list<int> &sources, const MatrixXd &targets,const MatrixXi& F, VectorXd &min_distance);
#endif



