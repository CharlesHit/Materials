//
// Created by Charles on 2019-04-11.
//

#ifndef MUS_WORLD_OBJECT_H
#define MUS_WORLD_OBJECT_H

#include "matrix.h"
#include "windows.h"

enum MaterialType { LAMBERT, GRASS, REFLECTION, METAL, STEEL };

class Object
{
public:
	Object();
	virtual ~Object() = default;
	//get the intersection
	virtual bool intersect(const Vec3 &, const Vec3 &, double &, uint32_t &, Vec2 &) const = 0;
	//get the normal vector
	virtual void getSurfaceProperties(const Vec3 &, const Vec3 &, const uint32_t &, const Vec2 &, Vec3 &, Vec2 &) const = 0;
	virtual Vec3 evalDiffuseColor(const Vec2 &) const;
	MaterialType materialType;
	double ior;
	double Kd, Ks;
	Vec3 diffuseColor;
	double specularExponent;
};

bool solveQuadratic(const double &a, const double &b, const double &c, double &x0, double &x1);

class Sphere : public Object{
public:
	Vec3 center;
	double radius, radius2;
	Sphere(const Vec3 &c, const double &r);
	bool intersect(const Vec3 &orig, const Vec3 &dir, double &tnear, uint32_t &index, Vec2 &uv) const override;
	void getSurfaceProperties(const Vec3 &P, const Vec3 &I, const uint32_t &index, const Vec2 &uv, Vec3 &N, Vec2 &st) const override;
};

class Cone : public Object{
public:
	Vec3 center;
	double radius, height;
	Cone(const Vec3 &center, const double &radius,  const double &height);
	bool intersect(const Vec3 &orig, const Vec3 &dir, double &tnear, uint32_t &index, Vec2 &uv) const override;
	void getSurfaceProperties(const Vec3 &P, const Vec3 &I, const uint32_t &index, const Vec2 &uv, Vec3 &N, Vec2 &st) const override;
};

bool rayTriangleIntersect(
		const Vec3 &v0, const Vec3 &v1, const Vec3 &v2,
		const Vec3 &orig, const Vec3 &dir,
		double &tnear, double &u, double &v);

class Plane : public Object
{
public:
	Plane(Vec3 color_1, Vec3 color_2, double depth_mesh, double near_mesh, double far_mesh, double left_mesh, double right_mesh);

	bool intersect(const Vec3 &orig, const Vec3 &dir, double &tnear, uint32_t &index, Vec2 &uv) const override;
	void getSurfaceProperties(const Vec3 &P, const Vec3 &I, const uint32_t &index, const Vec2 &uv, Vec3 &N, Vec2 &st) const override;
	Vec3 evalDiffuseColor(const Vec2 &st) const override;

	std::unique_ptr<Vec3[]> vertices;
	uint32_t numTriangles;
	std::unique_ptr<uint32_t[]> vertexIndex;
	std::unique_ptr<Vec2[]> stCoordinates;
	Vec3 color0;
	Vec3 color1;
};

bool trace(
		const Vec3 &orig, const Vec3 &dir,
		const std::vector<std::unique_ptr<Object>> &objects,
		double &tNear, uint32_t &index, Vec2 &uv, Object **hitObject);

Vec3 castRay(
		const Vec3 &orig, const Vec3 &dir,
		const std::vector<std::unique_ptr<Object>> &objects,
		const std::vector<std::unique_ptr<Light>> &lights,
		const Windows &options,
		uint32_t depth,
		bool test = false);

#endif MUS_WORLD_OBJECT_H
