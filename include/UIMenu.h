
#ifndef UIMENU_H
#define UIMENU_H

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

class UIMenu: public igl::opengl::glfw::imgui::ImGuiMenu {
      public:
          std::function<bool(void)> callback_post_load;
          IGL_INLINE virtual bool post_load() override
          {
             if (callback_post_load) { callback_post_load(); }
             return false;
          }
};
#endif
