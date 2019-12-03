
#ifndef POINTCLOUDHEAT_H
#define POINTCLOUDHEAT_H

#include "DomainHeat.h"

class PointCloudHeat: public DomainHeat {
      protected:
         void init_attr();
      public:
          inline PointCloudHeat(const MatrixXd& Ve){V=Ve;}
          void compute(int);
};
#endif
