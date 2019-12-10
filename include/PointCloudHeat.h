
#ifndef POINTCLOUDHEAT_H
#define POINTCLOUDHEAT_H

#include "DomainHeat.h"

class PointCloudHeat: public DomainHeat {
      private:
        MatrixXd D;
      protected:
         void init_attr();
         void computeTangentNormals();
         void computeTimeStep();
         void solvePoisson();
      public:
          inline PointCloudHeat(const MatrixXd& Ve){V=Ve;}
          void solveVectorField();
};
#endif
