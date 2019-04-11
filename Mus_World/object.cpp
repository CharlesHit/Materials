//
// Created by Mu.H on 2019-04-01.
//

#include "object.h"

Object::Object () :
		materialType(DIFFUSE_AND_GLOSSY),
		ior(1.3), Kd(0.8), Ks(0.2), diffuseColor(Vec3(0.2)), specularExponent(25) {}

Vec3 Object::evalDiffuseColor ( const Vec2 & ) const
{ return diffuseColor; }

bool solveQuadratic ( const double &a, const double &b, const double &c, double &x0, double &x1 )
{
	double discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0) x0 = x1 = - 0.5 * b / a;
	else {
		double q = (b > 0) ?
		           -0.5 * (b + sqrt(discr)) :
		           -0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}
	if (x0 > x1) std::swap(x0, x1);
	return true;
}

bool
rayTriangleIntersect ( const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec3 &orig, const Vec3 &dir, double &tnear,
                       double &u, double &v )
{
	Vec3 edge1 = v1 - v0;
	Vec3 edge2 = v2 - v0;
	Vec3 pvec = crossProduct(dir, edge2);
	double det = dotProduct(edge1, pvec);
	if (det == 0 || det < 0) return false;

	Vec3 tvec = orig - v0;
	u = dotProduct(tvec, pvec);
	if (u < 0 || u > det) return false;

	Vec3 qvec = crossProduct(tvec, edge1);
	v = dotProduct(dir, qvec);
	if (v < 0 || u + v > det) return false;

	double invDet = 1 / det;

	tnear = dotProduct(edge2, qvec) * invDet;
	u *= invDet;
	v *= invDet;

	return true;
}

bool trace ( const Vec3 &orig, const Vec3 &dir, const std::vector<std::unique_ptr<Object>> &objects, double &tNear,
             uint32_t &index, Vec2 &uv, Object ** hitObject )
{
	*hitObject = nullptr;
	for (const auto & object : objects) {
		double tNearK = kInfinity;
		uint32_t indexK;
		Vec2 uvK;
		if (object->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear) {
			*hitObject = object.get();
			tNear = tNearK;
			index = indexK;
			uv = uvK;
		}
	}

	return (*hitObject != nullptr);
}

Vec3 castRay ( const Vec3 &orig, const Vec3 &dir, const std::vector<std::unique_ptr<Object>> &objects,
               const std::vector<std::unique_ptr<Light>> &lights, const Windows &options, uint32_t depth, bool test )
{
	if (depth > options.maxDepth) {
		return options.backgroundColor;
	}

	Vec3 hitColor = options.backgroundColor;
	double tnear = kInfinity;
	Vec2 uv;
	uint32_t index = 0;
	Object *hitObject = nullptr;
	if (trace(orig, dir, objects, tnear, index, uv, &hitObject)) {
		Vec3 hitPoint = orig + dir * tnear;
		Vec3 N; // normal
		Vec2 st; // st coordinates
		hitObject->getSurfaceProperties(hitPoint, dir, index, uv, N, st);
		switch (hitObject->materialType) {
			case REFLECTION_AND_REFRACTION:
			{
				Vec3 reflectionDirection = normalize(reflect(dir, N));
				Vec3 refractionDirection = normalize(refract(dir, N, hitObject->ior));
				Vec3 reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?
				                         hitPoint - N * options.bias :
				                         hitPoint + N * options.bias;
				Vec3 refractionRayOrig = (dotProduct(refractionDirection, N) < 0) ?
				                         hitPoint - N * options.bias :
				                         hitPoint + N * options.bias;
				Vec3 reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1, true);
				Vec3 refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1, true);
				double kr;
				fresnel(dir, N, hitObject->ior, kr);
				hitColor = reflectionColor * kr + refractionColor * (1 - kr);
				break;
			}
			case REFLECTION:
			{
				double kr;
				fresnel(dir, N, hitObject->ior, kr);
				Vec3 reflectionDirection = reflect(dir, N);
				Vec3 reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?
				                         hitPoint + N * options.bias :
				                         hitPoint - N * options.bias;
				hitColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1) * kr;
				break;
			}
			default:
			{
				// [comment]
				// We use the Phong illumation model int the default case. The phong model
				// is composed of a diffuse and a specular reflection component.
				// [/comment]
				Vec3 lightAmt = 0, specularColor = 0;
				Vec3 shadowPointOrig = (dotProduct(dir, N) < 0) ?
				                       hitPoint + N * options.bias :
				                       hitPoint - N * options.bias;
				// [comment]
				// Loop over all lights in the scene and sum their contribution up
				// We also apply the lambert cosine law here though we haven't explained yet what this means.
				// [/comment]
				for (const auto & light : lights) {
					Vec3 lightDir = light->position - hitPoint;
					// square of the distance between hitPoint and the light
					double lightDistance2 = dotProduct(lightDir, lightDir);
					lightDir = normalize(lightDir);
					double LdotN = std::max(0.0, dotProduct(lightDir, N));
					Object *shadowHitObject = nullptr;
					double tNearShadow = kInfinity;
					// is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
					bool inShadow = trace(shadowPointOrig, lightDir, objects, tNearShadow, index, uv, &shadowHitObject) &&
					                tNearShadow * tNearShadow < lightDistance2;
					lightAmt += (1 - inShadow) * light->intensity * LdotN;
					Vec3 reflectionDirection = reflect(-lightDir, N);
					specularColor += powf(std::max(0.0, -dotProduct(reflectionDirection, dir)), hitObject->specularExponent) * light->intensity;
				}
				hitColor = lightAmt * hitObject->evalDiffuseColor(st) * hitObject->Kd + specularColor * hitObject->Ks;
				break;
			}
		}
	}

	return hitColor;
}

Sphere::Sphere ( const Vec3 &c, const double &r ) : center(c), radius(r), radius2(r * r) {}

bool Sphere::intersect ( const Vec3 &orig, const Vec3 &dir, double &tnear, uint32_t &index, Vec2 &uv ) const
{
	// analytic solution
	Vec3 L = orig - center;
	double a = dotProduct(dir, dir);
	double b = 2 * dotProduct(dir, L);
	double c = dotProduct(L, L) - radius2;
	double t0, t1;
	if (!solveQuadratic(a, b, c, t0, t1)) return false;
	if (t0 < 0) t0 = t1;
	if (t0 < 0) return false;
	tnear = t0;

	return true;
}

void Sphere::getSurfaceProperties ( const Vec3 &P, const Vec3 &I, const uint32_t &index, const Vec2 &uv, Vec3 &N,
                                    Vec2 &st ) const
{ N = normalize(P - center); }

bool rayTriangleIntersect(
		const Vec3 &v0, const Vec3 &v1, const Vec3 &v2,
		const Vec3 &orig, const Vec3 &dir,
		float &tnear, float &u, float &v)
{
	Vec3 edge1 = v1 - v0;
	Vec3 edge2 = v2 - v0;
	Vec3 pvec = crossProduct(dir, edge2);
	float det = dotProduct(edge1, pvec);
	if (det == 0 || det < 0) return false;

	Vec3 tvec = orig - v0;
	u = dotProduct(tvec, pvec);
	if (u < 0 || u > det) return false;

	Vec3 qvec = crossProduct(tvec, edge1);
	v = dotProduct(dir, qvec);
	if (v < 0 || u + v > det) return false;

	float invDet = 1 / det;

	tnear = dotProduct(edge2, qvec) * invDet;
	u *= invDet;
	v *= invDet;

	return true;
}

MeshTriangle::MeshTriangle ( const Vec3 * verts, const uint32_t * vertsIndex, const uint32_t &numTris, const Vec2 * st )
{
	uint32_t maxIndex = 0;
	for (uint32_t i = 0; i < numTris * 3; ++i)
		if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
	maxIndex += 1;
	vertices = std::unique_ptr<Vec3[]>(new Vec3[maxIndex]);
	memcpy(vertices.get(), verts, sizeof(Vec3) * maxIndex);
	vertexIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
	memcpy(vertexIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
	numTriangles = numTris;
	stCoordinates = std::unique_ptr<Vec2[]>(new Vec2[maxIndex]);
	memcpy(stCoordinates.get(), st, sizeof(Vec2) * maxIndex);
}

bool MeshTriangle::intersect ( const Vec3 &orig, const Vec3 &dir, double &tnear, uint32_t &index, Vec2 &uv ) const
{
	bool intersect = false;
	for (uint32_t k = 0; k < numTriangles; ++k) {
		const Vec3 & v0 = vertices[vertexIndex[k * 3]];
		const Vec3 & v1 = vertices[vertexIndex[k * 3 + 1]];
		const Vec3 & v2 = vertices[vertexIndex[k * 3 + 2]];
		double t, u, v;
		if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
			tnear = t;
			uv.x = u;
			uv.y = v;
			index = k;
			intersect |= true;
		}
	}

	return intersect;
}

void MeshTriangle::getSurfaceProperties ( const Vec3 &P, const Vec3 &I, const uint32_t &index, const Vec2 &uv, Vec3 &N,
                                          Vec2 &st ) const
{
	const Vec3 &v0 = vertices[vertexIndex[index * 3]];
	const Vec3 &v1 = vertices[vertexIndex[index * 3 + 1]];
	const Vec3 &v2 = vertices[vertexIndex[index * 3 + 2]];
	Vec3 e0 = normalize(v1 - v0);
	Vec3 e1 = normalize(v2 - v1);
	N = normalize(crossProduct(e0, e1));
	const Vec2 &st0 = stCoordinates[vertexIndex[index * 3]];
	const Vec2 &st1 = stCoordinates[vertexIndex[index * 3 + 1]];
	const Vec2 &st2 = stCoordinates[vertexIndex[index * 3 + 2]];
	st = st0 * (1 - uv.x - uv.y) + st1 * uv.x + st2 * uv.y;
}

Vec3 MeshTriangle::evalDiffuseColor ( const Vec2 &st ) const
{
	double scale = 5;
	double pattern = (fmodf(st.x * scale, 1) > 0.5) ^ (fmodf(st.y * scale, 1) > 0.5);
	return mix(Vec3(0.815, 0.235, 0.031), Vec3(0.937, 0.937, 0.231), pattern);
}
