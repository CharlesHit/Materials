//
// Created by Charles on 2019-04-11.
//

#ifndef MUS_WORLD_WINDOWS_H
#define MUS_WORLD_WINDOWS_H

#include "matrix.h"

class Windows
{
public:
	int width;
	int height;
	double fov;
	int maxDepth;
	Vec3 backgroundColor;
	double bias;
	Windows();
	Windows(int w, int h, Vec3 backcolor, int maxRecur);

};

class Light
{
public:
	Light(const Vec3 &p, const Vec3 &i) : position(p), intensity(i) {}
	Vec3 position;
	Vec3 intensity;
};


#endif //MUS_WORLD_WINDOWS_H
