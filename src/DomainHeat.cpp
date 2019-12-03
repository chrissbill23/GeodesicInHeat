
#include "DomainHeat.h"
#include <iostream>
using namespace std;

void DomainHeat::init() {
  init_attr(); 
  computeTimeStep();
} 
void DomainHeat::computeTimeStep(){
     time = 0.;
     MatrixXi comp = MatrixXi::Zero(A.rows(), A.cols());
     int totEdges = 0;
     for(size_t i = 0; i < A.rows(); ++i){
        for(size_t j = 0; j < A.cols(); ++j) {
            if(i!=j && A(i,j) == 1 && comp(i,j) == 0){
                 time += D(i,j);
                 comp(i,j) = 1;
                 comp(j,i) = 1;
                 ++totEdges;
            }
        }
     }
     if(totEdges > 0)
        time /= totEdges;
}


