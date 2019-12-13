
#ifndef POINTCLOUDHEAT_H
#define POINTCLOUDHEAT_H

#include "DomainHeat.h"

class PointCloudHeat: public DomainHeat {
      private:
        MatrixXd D;
        MatrixXd C;
        MatrixXi A;
        LLT<MatrixXd> laplacianFactor;
      protected:
         void init_attr();
         void computeTangentNormals();
         void computeTimeStep();
         void solvePoisson();
      public:
          inline PointCloudHeat(MatrixXd& Ve){Vertices=&Ve;}
          void solveVectorField();
};
#endif
