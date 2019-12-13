
#ifndef TRIANGLEMESHHEAT_H
#define TRIANGLEMESHHEAT_H

#include "DomainHeat.h"
#include "HalfedgeDS.cpp"

class TriangleMeshHeat: public DomainHeat {
      public:
          enum TypeHeat{PURE_NEUMAN, PURE_DIRICHLET, BOUNDARY_COND_NEUMAN, BOUNDARY_COND_DIRICHLET, BOUNDARY_COND_BOTH};
      private:
          HalfedgeDS& he;
          list<int> boundaryIndexes;
          LLT<MatrixXd> laplacianFactor;
          VectorXd heatNeuman;
          VectorXd heatDirichlet;
          TypeHeat boundarycondition;
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
          inline TriangleMeshHeat(MatrixXd& Ve, MatrixXi& Fa, HalfedgeDS& h, TypeHeat bc = PURE_NEUMAN): he(h),boundarycondition(bc){ Vertices=&Ve; Faces = &Fa; }
          void solveVectorField();
          inline void setHalfEdge(HalfedgeDS& h){he = h;}
};
#endif
