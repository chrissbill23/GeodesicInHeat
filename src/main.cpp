#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>

#include "HalfedgeBuilder.cpp"
#include "TriangleMeshHeat.h"
#include "PointCloudHeat.h"
#include "Laplacian.h"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;


int main(int argc, char *argv[])
{

   MatrixXd V;
   MatrixXi F;
   if (argc < 2)
   {
	igl::readOFF("../data/homer.off", V, F);
   } else {
            igl::readOFF(argv[1], V, F);
          }
   cout<<"Tot vertices: "<<V.rows()<<" Tot faces: "<<F.rows()<<endl;
   igl::opengl::glfw::Viewer viewer; // create the 3d viewer
   viewer.data().set_mesh(V, F);
   viewer.core(0).align_camera_center(V, F);
   viewer.launch(); // run the viewer
   MatrixXd D; MatrixXd A; MatrixXd L; MatrixXi F2;
   DomainHeat* t = new TriangleMeshHeat(V,F);
   t->init();

}
