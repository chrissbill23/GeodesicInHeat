
#ifndef POINTCLOUD_HEAT
#define POINTCLOUD_HEAT

#include "DomainHeat.h"

class PointCloudHeat: public DomainHeat {
      public:
          inline PointCloudHeat(HalfedgeDS h) : DomainHeat(h){}
          void compute(int);
};
#endif
