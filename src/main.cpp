#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/unproject_onto_mesh.h>
#include "HalfedgeBuilder.cpp"
#include "TriangleMeshHeat.h"
#include "PointCloudHeat.h"
#include "Utils.h"
#include <UIMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

using namespace Eigen; 
using namespace std;

thread th;
thread th2, th3;
int main(int argc, char *argv[])
{

   MatrixXd V;
   MatrixXi F;
   list<int> sources;
   Eigen::VectorXd D;
   Eigen::VectorXd Dij;
   int timegeo = 0;
   int timedji = 0;
   double smoothtime = 1.;
   bool initterminated = false;
   string filename = argc < 2 ? "../data/sphere.off" : argv[1];
   igl::readOFF(filename, V, F);

   Eigen::MatrixXd C = Eigen::MatrixXd::Constant(V.rows(),3,1);
   HalfedgeBuilder builder;  //

   HalfedgeDS he=builder.createMeshWithFaces(V.rows(), F);  //
   cout<<"Tot vertices: "<<V.rows()<<" Tot faces: "<<F.rows()<<endl;
   TriangleMeshHeat tri(V,F,he , TriangleMeshHeat::BOUNDARY_COND_BOTH);
   DomainHeat* t = &tri;

  igl::opengl::glfw::Viewer viewer;

   
   UIMenu menu;

   
  bool down_on_mesh = false;
  bool multipoints = false;
  bool dijicomplete = false;
  bool computediji = false;
  bool clickable = true;

  th2 = thread([&](){
          initterminated = false;
          clickable = false;
          t->init();
          initterminated = true;
          clickable = true;
      });
  const auto update = [&]()->bool
  {
    int fid;
    Eigen::Vector3f bc;
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view,
      viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
    {
      if(clickable){
         computediji = false; dijicomplete = false; clickable = false;
         const Eigen::RowVector3d m3 =
         V.row(F(fid,0))*bc(0) + V.row(F(fid,1))*bc(1) + V.row(F(fid,2))*bc(2);
         int cid = 0;
         Eigen::Vector3d(
           (V.row(F(fid,0))-m3).squaredNorm(),
           (V.row(F(fid,1))-m3).squaredNorm(),
           (V.row(F(fid,2))-m3).squaredNorm()).minCoeff(&cid);
         const int vid = F(fid,cid);
         C.row(vid)<<1,0,0;

         timegeo = t->compute(vid, smoothtime, !multipoints);
         if(multipoints == false){
            sources.clear();
          }
          sources.push_back(vid);
          D = t->getDistance();
          viewer.data().set_colors((D/D.maxCoeff()).replicate(1,3));
          clickable = true;
      }
      return true;
    }
    return false;
  };

//MENU  ------ BEGIN
   
  viewer.plugins.push_back(&menu);

  menu.callback_draw_custom_window = [&]()
  {
    double h = 160, w = 200.;
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(w, h), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "SOURCES", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );
    
     list <int> :: const_iterator it; bool exists = false;
     for(it = sources.begin(); it != sources.end(); ++it){
        size_t i = *it;
         if(i >= 0 && i < V.rows()){
           ImGui::Text("Vertex %d : (%.2f, %.2f, %.2f)", i, V(i,0), V(i,1), V(i,2));
           exists = true;
         }
     }

    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), h + 20), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(w, h), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "MIN DISTANCES FROM SOURCES", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );
    if(exists){
        for(size_t i = 0; i < D.rows(); ++i){
           ImGui::Text("Vertex %d : %.9f", i, D(i));
        }
     }

    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), h*2 + 20), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(w, h), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "DIJKSTRA MIN DISTANCES", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );
    if(!computediji){
       if (ImGui::Button("Compute Dijkstra", ImVec2(-1,0)))
      {
        if(th.joinable())
           th.join();
        computediji = true;
        clickable = false;
        th = thread([&](){
           timedji = dijkstra(sources, V,F, Dij);
           dijicomplete = true;
           clickable = true;
        });
      }
    }
    else if(exists && dijicomplete){
        for(size_t i = 0; i < Dij.rows(); ++i){
           ImGui::Text("Vertex %d : %.9f", i, Dij(i));
        }
     }

    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), h*3 + 20), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(w, h/2), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "TIME PERFORMANCES", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );
    if(exists){
        ImGui::Text("Heat Geodesic: %d msec", timegeo);
        if(dijicomplete)
           ImGui::Text("Dijkstra: %d msec", timedji);
     }

    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(0.f * menu.menu_scaling(), h*3.65 ), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(w*1.1, h*0.6), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "SETTING", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );
    ImGui::PushItemWidth(-80);
    if (ImGui::InputDouble("Scale time step", &smoothtime, 0.00000005f, 0, "%.8f")){
      if(sources.size() > 0 && clickable){
        if(th3.joinable())
            th3.join();
        th3 = thread([&](){
        clickable = false;
        timegeo = t->reload(smoothtime);
        D = t->getDistance();
        viewer.data().set_colors((D/D.maxCoeff()).replicate(1,3));
        clickable = true;});
      }
    }  
    ImGui::PopItemWidth();

    ImGui::Text("Current time step: %.16f ", t->getTime());
    ImGui::Checkbox("Multi sources", &multipoints);
    ImGui::End();



    
    if(!initterminated){
          ImGui::SetNextWindowPos(ImVec2(180.f * 3 * menu.menu_scaling(), h*2 + 20), ImGuiSetCond_FirstUseEver);
          ImGui::SetNextWindowSize(ImVec2(w*1.25, h/2), ImGuiSetCond_FirstUseEver);
          ImGui::Begin(
                 "INITIALIZING HEAT OBJECT......", nullptr,
                  ImGuiWindowFlags_NoSavedSettings
           );
           ImGui::Text("See terminal for more information on the init.");


           ImGui::End();
    }
  };
//MENU  ------ END


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

    viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer& , unsigned int key, int mod)->bool
  {
    switch(key)
    {
    default:
      return false;
    case 'm': multipoints = true; break;
    case 'M': multipoints = true; break;
    case 's': multipoints = false; break;
    case 'S': multipoints = false; break;
    }
    return true;
  };

  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.data().show_lines = false;
  viewer.launch_init(true,false);

  const auto updatemeshgl = [&]()->void{
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
  };
  updatemeshgl();
  menu.callback_post_load = [&]()
    {
       viewer.erase_mesh(0);
       V = viewer.data().V;
       F = viewer.data().F;
       sources.clear();
       timegeo = 0;
       timedji = 0;
       multipoints = false;
       computediji = false;

      C = Eigen::MatrixXd::Constant(V.rows(),3,1);
      he=builder.createMeshWithFaces(V.rows(), F);
      cout<<"Tot vertices: "<<V.rows()<<" Tot faces: "<<F.rows()<<endl;
      viewer.data().set_colors(C);
      viewer.data().show_lines = false;
      if(th2.joinable())
         th2.join();
      th2 = thread([&](){
          initterminated = false;
          clickable = false;
          t->init();
          initterminated = true;
          clickable = true;
      });

      updatemeshgl();
      
      
      return false;
   };
  //viewer.launch();
  viewer.launch_rendering(true);
  viewer.launch_shut();

}
