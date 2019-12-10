
#ifndef DOMAINHEAT_H
#define DOMAINHEAT_H

#include <igl/min_quad_with_fixed.h>
#include <chrono>

using namespace Eigen;

class DomainHeat {
      private:
          std::chrono::time_point<std::chrono::high_resolution_clock> start;
      protected:
          MatrixXd V;
          MatrixXi F;
          VectorXd sources;
          SparseMatrix<double> L;
          SparseMatrix<double> M;
          SparseMatrix<double> M_inv;
          MatrixXi A;
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
          void stopTimer();
      public:
          void compute(int, double = -1.);
          void init();

          virtual VectorXi boundaryIndexes()const;


          inline VectorXd getDistance()const{return poisson;}
          inline double getTime()const{ return time;}
          inline const SparseMatrix<double>& getLaplacian()const{ return L;}
          inline const MatrixXi& getFaces()const{ return F;}
          inline const MatrixXi& getAdjacency()const{ return A;}
          inline SparseMatrix<double> laplacianOperator()const{ return M_inv*L;}
          inline SparseMatrix<double> heatFlowOperator()const{ return M - time * L;}
          
};
#endif

