//
// Created by Charles on 2019-04-11.
//

#ifndef MUS_WORLD_MAIN_H
#define MUS_WORLD_MAIN_H

#include "GLUT/glut.h"

#include <math.h>
#include <limits>
#include <vector>

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

#define Near             1.0
#define Far              25.0

#define THETA            45.0

const int H = 768;
const int W = 1024;
const double ASPECT = H/W;

const long nChars = W * H * 3;

void Quit() {
	exit(0);
}

void OnKeyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 'q':
			Quit();
			break;
	}
}


#endif MUS_WORLD_MAIN_H
