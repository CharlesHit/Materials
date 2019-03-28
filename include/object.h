#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

class Camera {
public:
	vec UP;
	vec E;
	vec G;
	vec n;//n = GE = E - G
	vec u, v;


	Camera::Camera ( )
	{	
		E = Eye;
		G = Gaze;
		n = E - G; n.normalized();
		UP = Up;
		u = UP; u.cross(n).normalized();
		v = n; v.cross(u).normalized();
	}
};

struct color_t {
	double r, g, b;

	color_t() {}
	color_t(double _r, double _g, double _b) {
		r = _r;
		g = _g;
		b = _b;
	}
	color_t(const color_t& initializer) {
		r = initializer.r;
		g = initializer.g;
		b = initializer.b;
	}

	color_t operator*(double scalar) {
		return color_t(r * scalar, g * scalar, b * scalar);
	}

	color_t operator+(const color_t& col) {
		return color_t(r + col.r, g + col.g, b + col.b);
	}
};

struct object_t {
	int type;
	double(*intersection)(dmatrix_t *, dmatrix_t *);
	vec(*normal)(dmatrix_t *);
	dmatrix_t M, Minv;
	color_t specular_color, diffuse_color, ambient_color;
	double density, reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f;
};

struct light_t {
	dmatrix_t position;
	color_t color;
	color_t intensity;
};