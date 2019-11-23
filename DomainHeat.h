
#ifndef DOMAINHEAT_H
#define DOMAINHEAT_H

#include "HalfedgeDS.cpp"

using namespace Eigen;
using namespace std;

class DomainHeat {
      protected:
          HalfedgeDS he;
          MatrixXd flow;

      public:
          inline DomainHeat(HalfedgeDS h = HalfedgeDS(0,0)) : he(h){}
          virtual void compute(int) = 0;
};
#endif

