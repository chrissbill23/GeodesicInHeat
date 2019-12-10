
#ifndef TRIANGLEMESHHEAT_H
#define TRIANGLEMESHHEAT_H

#include "DomainHeat.h"
#include "HalfedgeDS.cpp"

class TriangleMeshHeat: public DomainHeat {
      public:
          enum TypeHeat{PURE_NEUMAN, PURE_DIRICHLET, BOUNDARY_COND_NEUMAN, BOUNDARY_COND_DIRICHLET, BOUNDARY_COND_BOTH};
      private:
          const HalfedgeDS& he;
          VectorXi boundaryIndexes;
          LLT<MatrixXd> laplacianFactor;
          VectorXd heatNeuman;
          VectorXd heatDirichlet;
          TypeHeat boundarycondition = PURE_NEUMAN;
          void solveHeatNeuman(bool = false);
          void solveHeatDirichlet(bool = false);
          void findBoundaries();
      protected:
         void init_attr();
         void computeTimeStep();
         void computeTangentNormals();
         void solveHeat();
         void solvePoisson();
      public:
          inline TriangleMeshHeat(const MatrixXd& Ve, const MatrixXi& Fa, const HalfedgeDS& h): he(h){ V=Ve; F = Fa; }
          void solveVectorField();
};
#endif
