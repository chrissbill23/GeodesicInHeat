
#include "TriangleMeshHeat.h"
#include "Laplacian.h"

void TriangleMeshHeat::compute(int) { 
  //TODO

}
void TriangleMeshHeat::init_attr(){
  laplacianTriMesh(V, F,D, A, L);
}

