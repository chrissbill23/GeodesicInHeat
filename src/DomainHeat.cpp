
#include "DomainHeat.h"
#include <iostream>
using namespace std;

void DomainHeat::init() {
  cout<<"******** FACTORIZATION ********"<<endl;
  startTimer();
  init_attr();
  sources.setConstant(V.rows(), 1);
  computeTimeStep();
  cout << " FACTORIZATION OF THE HEAT FLOW OPERATOR"<<endl;
  heatFlowOperatorFactor = LDLT<MatrixXd>(heatFlowOperator());
  computeTangentNormals();
  stopTimer();
} 

void DomainHeat::compute(int index, double smoothtime){
     if(index < 0 || index >= V.rows())
        return;
     if(smoothtime > -1.){
        smooth = smoothtime;
        cout << "time smoothness has changed. Factor revaluation..."<<endl;
        computeTimeStep();
        heatFlowOperatorFactor = LDLT<MatrixXd>(heatFlowOperator());
     }
     cout<<"******** COMPUTATION OF GEODESIC DISTANCE FROM SOURCE "<<index<<" AND TIME t = "<<time<<" ********"<<endl;
     startTimer();
     sources = VectorXd::Zero(V.rows());
     sources(index) = 1;
     solveHeat(); 
     solveVectorField();
     solvePoisson();
     stopTimer();
}
VectorXi DomainHeat::boundaryIndexes() const{
       return VectorXi();
}
void DomainHeat::solveHeat(){
     heat = heatFlowOperatorFactor.solve(sources);;
}

/*
void DomainHeat::computeTangentNormals(){
  // compute the normals and tangents using PCA
  tangents = normals = MatrixXd::Zero(V.rows(),3);
  for(int i = 0; i<V.rows(); ++i){
        int k = A.row(i).sum();
	MatrixXd knns(k+1,3);
	for(int j = 0; j < A.cols(); ++ j){
           if(A(i,j) == 1)
	   knns.row(j) = V.row(j);
	}
	knns.row(k) = V.row(i);
        MatrixXd yi = knns.rowwise() - knns.colwise().mean();
        MatrixXd cov;
	if(k > 1)
 	    cov = (yi.adjoint() * yi) / double(knns.rows() - 1);
	else 
	   cov = (yi.adjoint() * yi);

       EigenSolver<MatrixXd> s;
       s.compute(cov);
       int index=0;
       VectorXd evals = s.pseudoEigenvalueMatrix().diagonal();
       evals.minCoeff(&index);
       VectorXd normal = s.pseudoEigenvectors().col(index).transpose();
       double d = -normal(0)*V(i,0)-normal(1)*V(i,1)-normal(2)*V(i,2);
       VectorXd tangent(3);
       tangent(0) = -d/normal(0);
       tangent(1) = -d/normal(1);
       tangent(2) = -d/normal(2);
       
       normals.row(i) = normal;
       tangents.row(i) = tangent;
   }
}*/
void DomainHeat::startTimer(){
     start = std::chrono::high_resolution_clock::now();
}
void DomainHeat::stopTimer(){
     auto stop = std::chrono::high_resolution_clock::now();
     cout<<"------->END - Time elapsed = "<<std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count()<<" milliseconds <-------"<<endl;
}


