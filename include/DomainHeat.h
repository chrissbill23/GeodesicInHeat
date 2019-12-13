
#ifndef DOMAINHEAT_H
#define DOMAINHEAT_H

#include <igl/min_quad_with_fixed.h>
#include <chrono>
#include <list>
#include <thread>

using namespace Eigen;
using namespace std;

class DomainHeat {
      private:
          chrono::time_point<chrono::high_resolution_clock> start;
      protected:
          MatrixXd* Vertices = nullptr;
          MatrixXi* Faces = nullptr;
          VectorXd sources;
          SparseMatrix<double> L;
          SparseMatrix<double> M;
          SparseMatrix<double> M_inv;
          double time;
          double smooth = 1.;
          VectorXd heat;
          MatrixXd normals;
          MatrixXd tangents;
          VectorXd poisson;
          LDLT<MatrixXd> heatFlowOperatorFactor;

          virtual void init_attr() = 0;
          virtual void computeTangentNormals() = 0;
          virtual void solveHeat();
          virtual void solveVectorField() = 0;
          virtual void solvePoisson() = 0;
          virtual void computeTimeStep() = 0;

          void startTimer();
          int stopTimer();
      public:
          int compute(int, double = 1., bool = true);
          void init();
          int reload(double);


          inline VectorXd getDistance()const{return poisson;}
          inline double getTime()const{ return time;}
          inline const SparseMatrix<double>& getLaplacian()const{ return L;}
          inline SparseMatrix<double> laplacianOperator()const{ return M_inv*L;}
          inline SparseMatrix<double> heatFlowOperator()const{ return M - time * L;}
          
};
#endif

