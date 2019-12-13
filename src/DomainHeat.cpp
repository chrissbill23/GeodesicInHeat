
#include "DomainHeat.h"
#include <iostream>

using namespace std;

void DomainHeat::init() {
  cout<<"******** FACTORIZATION ********"<<endl;
  startTimer();
  init_attr();
  sources = VectorXd::Constant(Vertices->rows(), 1);
  computeTimeStep();
  cout << " FACTORIZATION OF THE HEAT FLOW OPERATOR"<<endl;
  heatFlowOperatorFactor = LDLT<MatrixXd>(heatFlowOperator());
  computeTangentNormals();
  stopTimer();
} 

int DomainHeat::compute(int index, double smoothtime, bool erase){
     if(index < 0 || index >= Vertices->rows())
        return 0.0;
     if(smoothtime > 0 && smoothtime != smooth){
        smooth = smoothtime;
        cout << "time smoothness has changed. Factor revaluation..."<<endl;
        computeTimeStep();
        heatFlowOperatorFactor = LDLT<MatrixXd>(heatFlowOperator());
     }
     cout<<"******** COMPUTATION OF GEODESIC DISTANCE FROM SOURCE "<<index<<" AND TIME t = "<<time<<" ********"<<endl;
     startTimer();
     if(erase)
        sources = VectorXd::Zero(Vertices->rows());
     sources(index) = 1;
     solveHeat(); 
     solveVectorField();
     solvePoisson();
     return stopTimer();
}
int DomainHeat::reload(double smoothtime){
    startTimer();
    if(smoothtime > 0 && smoothtime != smooth){
        smooth = smoothtime;
        cout << "time smoothness has changed. Factor revaluation..."<<endl;
        computeTimeStep();
        heatFlowOperatorFactor = LDLT<MatrixXd>(heatFlowOperator());
        solveHeat(); 
        solveVectorField();
        solvePoisson();
     }
     return stopTimer();
}


void DomainHeat::solveHeat(){
     heat = heatFlowOperatorFactor.solve(sources);;
}
void DomainHeat::startTimer(){
     start = std::chrono::high_resolution_clock::now();
}
int DomainHeat::stopTimer(){
     auto stop = std::chrono::high_resolution_clock::now();
     auto time = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
     cout<<"------->END - Time elapsed = "<<time<<" milliseconds <-------"<<endl;
     return time;
}


