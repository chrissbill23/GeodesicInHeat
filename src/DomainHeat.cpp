
#include "DomainHeat.h"
#include "Laplacian.h"
#include <iostream>
using namespace std;

void DomainHeat::init() {
  init_attr();
  massMatrix(V, F, M);
  computeTimeStep();
} 
void DomainHeat::computeTimeStep(){
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
}


