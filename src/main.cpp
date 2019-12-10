#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/unproject_onto_mesh.h>

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
   } 
   else 
   {
     igl::readOFF(argv[1], V, F);
   }
Eigen::MatrixXd C = Eigen::MatrixXd::Constant(V.rows(),3,1);
   HalfedgeBuilder builder;  //

   HalfedgeDS he=builder.createMeshWithFaces(V.rows(), F);  //
   cout<<"Tot vertices: "<<V.rows()<<" Tot faces: "<<F.rows()<<endl;

   DomainHeat* t = new TriangleMeshHeat(V,F,he);
   //DomainHeat* t = new PointCloudHeat(V);
   t->init();


  igl::opengl::glfw::Viewer viewer;


  bool down_on_mesh = false;
  const auto update = [&]()->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view,
      viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
    {
      // 3d position of hit
      const Eigen::RowVector3d m3 =
        V.row(F(fid,0))*bc(0) + V.row(F(fid,1))*bc(1) + V.row(F(fid,2))*bc(2);
      int cid = 0;
      Eigen::Vector3d(
          (V.row(F(fid,0))-m3).squaredNorm(),
          (V.row(F(fid,1))-m3).squaredNorm(),
          (V.row(F(fid,2))-m3).squaredNorm()).minCoeff(&cid);
      const int vid = F(fid,cid);
      C.row(vid)<<1,0,0;
      t->compute(vid);
      Eigen::VectorXd D = t->getDistance();
      //std::cout <<std::endl<<vid<<std::endl<<D<<std::endl<<std::endl;
      viewer.data().set_colors((D/D.maxCoeff()).replicate(1,3));
      return true;
    }
    return false;
  };
  viewer.callback_mouse_down =
    [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    if(update())
    {
      down_on_mesh = true;
      return true;
    }
    return false;
  };
  viewer.callback_mouse_up =
    [&down_on_mesh](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    down_on_mesh = false;
    return false;
  };

  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.data().show_lines = false;
  viewer.launch_init(true,false);

  viewer.data().meshgl.init();
  igl::opengl::destroy_shader_program(viewer.data().meshgl.shader_mesh);

  {
    std::string mesh_vertex_shader_string =
R"(#version 150
uniform mat4 view;
uniform mat4 proj;
uniform mat4 normal_matrix;
in vec3 position;
in vec3 normal;
out vec3 position_eye;
out vec3 normal_eye;
in vec4 Ka;
in vec4 Kd;
in vec4 Ks;
in vec2 texcoord;
out vec2 texcoordi;
out vec4 Kai;
out vec4 Kdi;
out vec4 Ksi;
void main()
{
  position_eye = vec3 (view * vec4 (position, 1.0));
  normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
  normal_eye = normalize(normal_eye);
  gl_Position = proj * vec4 (position_eye, 1.0); //proj * view * vec4(position, 1.0);
  Kai = Ka;
  Kdi = Kd;
  Ksi = Ks;
  texcoordi = texcoord;
})";

    std::string mesh_fragment_shader_string =
R"(#version 150
uniform mat4 view;
uniform mat4 proj;
uniform vec4 fixed_color;
in vec3 position_eye;
in vec3 normal_eye;
uniform vec3 light_position_eye;
vec3 Ls = vec3 (1, 1, 1);
vec3 Ld = vec3 (1, 1, 1);
vec3 La = vec3 (1, 1, 1);
in vec4 Ksi;
in vec4 Kdi;
in vec4 Kai;
in vec2 texcoordi;
uniform sampler2D tex;
uniform float specular_exponent;
uniform float lighting_factor;
uniform float texture_factor;
out vec4 outColor;
void main()
{
vec3 Ia = La * vec3(Kai);    // ambient intensity
float ni = 30.0; // number of intervals
float t = 1.0-round(ni*Kdi.r)/ni; // quantize and reverse
vec3 Kdiq = clamp(vec3(2.*t,2.*t-1.,6.*t-5.),0,1); // heat map
vec3 vector_to_light_eye = light_position_eye - position_eye;
vec3 direction_to_light_eye = normalize (vector_to_light_eye);
float dot_prod = dot (direction_to_light_eye, normalize(normal_eye));
float clamped_dot_prod = max (dot_prod, 0.0);
vec3 Id = Ld * Kdiq * clamped_dot_prod;    // Diffuse intensity
vec3 reflection_eye = reflect (-direction_to_light_eye, normalize(normal_eye));
vec3 surface_to_viewer_eye = normalize (-position_eye);
float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);
dot_prod_specular = float(abs(dot_prod)==dot_prod) * max (dot_prod_specular, 0.0);
float specular_factor = pow (dot_prod_specular, specular_exponent);
vec3 Kfi = 0.5*vec3(Ksi);
vec3 Lf = Ls;
float fresnel_exponent = 2*specular_exponent;
float fresnel_factor = 0;
{
  float NE = max( 0., dot( normalize(normal_eye), surface_to_viewer_eye));
  fresnel_factor = pow (max(sqrt(1. - NE*NE),0.0), fresnel_exponent);
}
vec3 Is = Ls * vec3(Ksi) * specular_factor;    // specular intensity
vec3 If = Lf * vec3(Kfi) * fresnel_factor;     // fresnel intensity
vec4 color = vec4(lighting_factor * (If + Is + Id) + Ia +
  (1.0-lighting_factor) * Kdiq,(Kai.a+Ksi.a+Kdi.a)/3);
outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;
if (fixed_color != vec4(0.0)) outColor = fixed_color;
})";

    igl::opengl::create_shader_program(
      mesh_vertex_shader_string,
      mesh_fragment_shader_string,
      {},
      viewer.data().meshgl.shader_mesh);
  }

  viewer.launch_rendering(true);
  viewer.launch_shut();

}
