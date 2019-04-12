cmake_minimum_required(VERSION 3.1)
project(Mus_World)

# Bump up warning levels appropriately for each compiler

#when compiled at Mac

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -march=native")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -DDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -fno-exceptions")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -fno-exceptions")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -fno-exceptions")

set(Mus_World_INSTALL_DIR "${Mus_World_SOURCE_DIR}/bin")


#when compiled at Win
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
	if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /Zi")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX /arch:AVX2")
	set(MUUM_INSTALL_DIR "${KYLE_SOURCE_DIR}/bin")
endif()

add_executable(Mus_World main.cpp)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} ) #${X11_INCLUDE_DIR})
target_link_libraries(Mus_World ${OPENGL_LIBRARIES} ${GLUT_LIBRARY}) #${X11_LIBRARY})

set_property(TARGET Mus_World PROPERTY CXX_STANDARD 14)
install(TARGETS Mus_World DESTINATION ${Mus_World_INSTALL_DIR})
