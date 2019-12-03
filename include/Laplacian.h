
#ifndef LAPLACIAN_H
#define LAPLACIAN_H

void distances(const MatrixXd &V, MatrixXd &D, MatrixXd &A,VectorXd& mind);
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A, MatrixXd &L, double maxthreshold);
void laplacianClouds(const MatrixXd &V, MatrixXd &D, MatrixXi &A,MatrixXi &F, MatrixXd &L, int k);
void laplacianTriMesh(const MatrixXd &V, const MatrixXi &F, MatrixXd &D, MatrixXi &A, MatrixXd &L);
#endif



