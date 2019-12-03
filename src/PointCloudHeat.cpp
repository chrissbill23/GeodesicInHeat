
#include "PointCloudHeat.h"
#include "Laplacian.h"

void PointCloudHeat::compute(int) { 
  //TODO

}
void PointCloudHeat::init_attr(){
  laplacianClouds(V, D, A, F, L, 3);
}
