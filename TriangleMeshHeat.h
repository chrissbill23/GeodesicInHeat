
#ifndef TRIANGLEMESHHEAT_H
#define TRIANGLEMESHHEAT_H

#include "DomainHeat.h"

class TriangleMeshHeat: public DomainHeat {
      public:
          inline TriangleMeshHeat(HalfedgeDS h) : DomainHeat(h){}
          void compute(int);
};
#endif
