//
// Created by Charles on 2019-04-11.
//

#ifndef KYLE_MAIN2_H
#define KYLE_MAIN2_H

#include "GLUT/glut.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include "matrix_ori.h"
#include <vector>

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

#define Near             1.0
#define Far              25.0

#define THETA            45.0

#define H                800
#define W                1600
const double ASPECT = H/W;

#define M_PI 3.14159265

const long nChars = W * H * 3;

void QuitX() {
	exit(0);
}

void OnKeyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 'q':
			QuitX();
			break;
	}
}


#endif //KYLE_MAIN2_H
