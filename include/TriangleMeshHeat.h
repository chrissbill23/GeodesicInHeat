
#ifndef TRIANGLEMESHHEAT_H
#define TRIANGLEMESHHEAT_H

#include "DomainHeat.h"

class TriangleMeshHeat: public DomainHeat {
      protected:
         void init_attr();
      public:
          inline TriangleMeshHeat(const MatrixXd& Ve, const MatrixXi& Fa){ V=Ve; F = Fa; }
          void compute(int);
};
#endif
