
#ifndef DOMAINHEAT_H
#define DOMAINHEAT_H

#include <igl/opengl/glfw/Viewer.h>
using namespace Eigen;

class DomainHeat {
      private:
          MatrixXd M;
      protected:
          MatrixXd V;
          MatrixXi F;
          MatrixXd flow;
          MatrixXd L;
          MatrixXd D;
          MatrixXi A;
          double time;
          virtual void init_attr() = 0;
      public:
          virtual void compute(int) = 0;
          virtual void computeTimeStep();


          void init();
          inline double getTime(){ return time;}
          inline const MatrixXd& getLaplacian(){ return L;}
          inline const MatrixXd& getFlow(){ return flow;}
          inline const MatrixXi& getFaces(){ return F;}
          inline const MatrixXd& getDistances(){ return D;}
          inline const MatrixXi& getAdjacency(){ return A;}
          inline MatrixXd laplacianOperator(){ return M.inverse()*L;}
          
};
#endif

