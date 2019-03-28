add_executable(MUUM main.cpp)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} ) #${X11_INCLUDE_DIR})

target_link_libraries(MUUM ${OPENGL_LIBRARIES} ${GLUT_LIBRARY}) #${X11_LIBRARY})

set_property(TARGET MUUM PROPERTY CXX_STANDARD 14)
install(TARGETS MUUM DESTINATION ${MUUM_INSTALL_DIR})

