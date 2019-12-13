
#include "PointCloudHeat.h"
#include "Utils.h"

void PointCloudHeat::init_attr(){
  Faces = new MatrixXi;
  laplacianClouds(*Vertices, D, A, *Faces, L, 3,normals,tangents, M, M_inv);
  smooth *= 0.005;
  laplacianFactor = LLT<MatrixXd> (L);
}
void PointCloudHeat::computeTangentNormals(){}

void PointCloudHeat::solveVectorField(){
  MatrixXd& V = *Vertices; MatrixXi& F = *Faces;
     poisson.setZero(V.rows());
     MatrixXd coef;
     leastSquareInterp3D3Degree(V,F,heat,coef);
     MatrixXd G;
     gradientLeastSquared3D3Degree(V, coef, G);
     MatrixXd Di(V.rows(), V.rows());
     for(size_t i = 0; i < V.rows(); ++i){
         for(size_t j = 0; j < V.rows(); ++j){
             Di(i,j) = wendlandweight(V.row(i),V.row(j))*abs(evaluateLeastSquared3D3Degree(V.row(j), coef.row(i)) - heat(j));
             /*if(A(i,j) == 1){
                G.row(j).normalize();
                poisson(i) += G(i,0)*heat(j) + G(i,1)*heat(j)+ G(i,2)*heat(j);

             }*/
         }
         G.row(i).normalize();
         poisson(i) += G(i,0)*heat(i) + G(i,1)*heat(i) + G(i,2)*heat(i);
     }
     poisson *= -1;
     poisson = Di.transpose() * M * poisson;
}
void PointCloudHeat::computeTimeStep(){
     
  MatrixXd& V = *Vertices; MatrixXi& F = *Faces;
     time = 0.;
     MatrixXi comp = MatrixXi::Zero(A.rows(), A.cols());
     int totEdges = 0;
     for(size_t i = 0; i < F.rows(); ++i){
        for(size_t j = 0; j < F.cols()-1; ++j) {
            int v0 = F(i,j), v1 = F(i,j+1);
            if(comp(v0,v1) == 0){
                 time += D(v0,v1);
                 comp(v0,v1) = 1;
                 comp(v1,v0) = 1;
                 ++totEdges;
            }
        }
     }
     if(totEdges > 0)
        time /= totEdges;
     time*=smooth;
}
void PointCloudHeat::solvePoisson(){
     
  MatrixXd& V = *Vertices; MatrixXi& F = *Faces;
     poisson = laplacianFactor.solve(poisson);
     double mean = 0.;
     for(size_t i = 0; i < sources.rows(); ++i){
        if(sources(i) == 1)
           mean += poisson(i);
     }
     mean /= sources.norm();
     poisson = poisson - VectorXd::Constant(sources.rows(),mean);
     if(poisson.mean() < 0)
        poisson *= -1;
}
