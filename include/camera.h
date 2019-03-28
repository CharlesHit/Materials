#pragma once
#include "matrix.h"
#include "data.h"

struct camera_t {
	dmatrix_t UP;
	dmatrix_t E;
	dmatrix_t G;
	dmatrix_t u, v, n;

	camera_t(camera_t *Camera, window_t *Window);
};