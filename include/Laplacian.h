
#ifndef LAPLACIAN_H
#define LAPLACIAN_H

void distances(const MatrixXd &V, MatrixXd &D, MatrixXd &A,VectorXd& mind);
//void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A, MatrixXd &L, double maxthreshold);
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A,MatrixXi &F, SparseMatrix<double> &L, int k, 
                     MatrixXd& N, MatrixXd& P, SparseMatrix<double>& M, SparseMatrix<double>& M_inv);
void laplacianTriMesh(const MatrixXd &V, const MatrixXi &F, MatrixXd &D, MatrixXi &A, SparseMatrix<double> &L);
void computeMeanCentroidOfFace(const MatrixXd &V, const RowVectorXi &F, RowVectorXd& c);
VectorXd computeAreaPolygon(const MatrixXd &V, const MatrixXi &F);
void massMatrix(const MatrixXd &V, const MatrixXi &F, SparseMatrix<double>& M, SparseMatrix<double>& M_inv, int type = 3);
double cotangent(const Vector3d& x0, const Vector3d& x1);
double cotangent(const Vector3d& x0, const Vector3d& x1, const Vector3d& x2);
void leastSquareInterp(const MatrixXd &V, const VectorXd &f, VectorXd &c, int m);
int factorial(int x);
double gaussianweight(const RowVectorXd& x0, const RowVectorXd& x1);
#endif



