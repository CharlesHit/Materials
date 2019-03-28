#pragma once
#include "matrix.h"

#define X                1
#define Y                2
#define Z                3

#define INFINITE_PLANE   0
#define PLANE            1
#define SPHERE           2

#define EPSILON          0.000001
#define N_OBJECTS        1024
#define MAX_INTENSITY    255.0

#define Ex               5.0
#define Ey               5.0
#define Ez               2.5

#define Gx               0.0
#define Gy               0.0
#define Gz              -1.0

#define UPx              0.0
#define UPy              0.0
#define UPz              1.0

#define Lx               10.0
#define Ly               5.0
#define Lz               5.0

#define M_PI 3.14159265

const double Near = 1.0;
const double Far = 25.0;

const double THETA = 45.0;
const int windowW = 300;
const int windowH = 300;
const double ASPECT = (double)windowW / (double)windowH;

vec Eye = vec(Ex, Ey, Ez, 1);
vec Gaze = vec(Gx, Gy, Gz, 1);
vec Up = vec(UPx, UPy, UPz, 1);
vec Light = vec(Lx, Ly, Lz, 1);