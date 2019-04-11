/*            PURPOSE : Simple framework for ray-tracing

		PREREQUISITES : matrix.h
 */
#include "GLUT/glut.h"
//#include "GL/freeglut.h"
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
#define ASPECT           1.5

#define H                300

#define M_PI 3.14159265

struct Window {
	int width, height, maxDepth;
	vec backgroundColor;
	double bias = 0.00001;
	double fov = 90;
};

struct camera_t {
	mat UP;
	mat E;
	mat G;
	mat u, v, n;
};

class Color:public vec{
public:

	double r = m[1][1];
	double g = m[2][1];
	double b = m[3][1];

	Color () {
		this->l = 3;
		this->c = 1;
		m[1][1] = 0;
		m[2][1] = 0;
		m[3][1] = 0;
		r = m[1][1];
		g = m[2][1];
		b = m[3][1];
	}

	Color(double _r, double _g, double _b) {
		this->l = 3;
		this->c = 1;
		m[1][1] = _r;
		m[2][1] = _g;
		m[3][1] = _b;
		r = m[1][1];
		g = m[2][1];
		b = m[3][1];
	}

	Color (Color const &color) {
		this->l = 3;
		this->c = 1;
		m[1][1] = color.r;
		m[2][1] = color.g;
		m[3][1] = color.b;
		r = m[1][1];
		g = m[2][1];
		b = m[3][1];
	}

	vec toVec()
	{
		vec v = vec(r, g, b);
		v.l = 3; v.c = 1;
		return v;
	}

//Color's operator
//vec + vec, vec - vec
	inline Color operator+(const Color &b )
	{
		return Color { m[1][1] + b.m[1][1], m[2][1] + b.m[2][1],
		             m[3][1] + b.m[3][1] };
	}
	inline Color& operator+=(const Color &b )
	{
		m[1][1] += b.m[1][1];
		m[2][1] += b.m[2][1];
		m[3][1] += b.m[3][1];
		return *this;
	}
	inline Color operator-( const Color &b )
	{
		return Color { m[1][1] - b.m[1][1], m[2][1] - b.m[2][1],
		             m[3][1] - b.m[3][1] };
	}
//vec * num, num* vec, vec*vec (dot product)
	inline Color operator*(double s )
	{
		return Color { m[1][1] * s, m[2][1] * s, m[3][1] * s };
	}
	//notice here
	friend Color operator*( double s, Color &v )
	{
		return v * s;
	}
	inline double operator*(const Color &v )
	{
		return this->dot ( v );
	}
//Color / num
	inline Color operator/(double s )
	{
		return (*this) * ( 1 / s );
	}
	inline Color operator/(const Color &b )
	{
		return Color { m[1][1] / b.m[1][1], m[2][1] / b.m[2][1],
		             m[3][1] / b.m[3][1] };
	}
	friend Color operator/( const Color &a, const Color &b )
	{
		return Color { a.m[1][1] / b.m[1][1], a.m[2][1] / b.m[2][1],
		             a.m[3][1] / b.m[3][1] };
	}
	friend Color operator/( double s, const Color &v )
	{
		return Color { s / v.m[1][1], s / v.m[2][1], s / v.m[3][1] };
	}
// Color = -Color
	inline Color operator-()
	{
		return Color { -m[1][1], -m[2][1], -m[3][1] };
	}

//equal
	friend bool operator==( const Color &a, const Color &b){
		if(a.l == b.l && a.l == 3 && a.c == b.c && a.c == 1)
		{
			return a.m[1][1] == b.m[1][1] && a.m[1][2] == b.m[1][2] && a.m[1][3] == b.m[1][3] ;
		}
		else return false;
	}
};

//Color's info
inline ostream &operator<<( ostream &os, const Color &a )
{
	if(0 < a.m[1][1] && a.m[1][1] < 1)
		os << "Color:\n\tr = " << a.m[1][1] * 255 << "\n\tg = " << a.m[2][1] * 255
	   << "\n\tb = " << a.m[3][1] * 255;
	else
		os << "Color:\n\tr = " << a.m[1][1] * 1 << "\n\tg = " << a.m[2][1] * 1
	        << "\n\tb = " << a.m[3][1] * 1;
	return os;
}

class Obj {
public:
	int type;
	//double(*intersection)(mat *, mat *);
	virtual bool intersect(const vec &, const vec &, double &, int &, vec2 &) const = 0;
	virtual void getSurfaceProperties(const vec &, const vec &, const int &, const vec2 &, vec &, vec2 &) const = 0;
	virtual vec evalDiffuseColor(const vec2 &) const { return diffuse_color; }
	double specular_exp = 25;
	mat *(*normal)(mat *);
	mat M, Minv;//Minv = inverse matrix
	Color specular_color, diffuse_color, ambient_color;
	double density, reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f, ior;

	Obj() :type(0),
			ior(1.3), density(0.8), reflectivity(0.2), diffuse_color(Color(0.2, 0.2, 0.2)), specular_coeff(25) {}
	virtual ~Obj() {}
};

class Light {
public:
	vec position;
	vec intensity;
	vec color;
	Light (){ }
	Light(const vec &p):position(p){intensity = vec(1.0, 1.0, 1.0);}
	Light(const vec &p, const vec &i):position(p),intensity(i){intensity = vec(1.0, 1.0, 1.0); color = vec(1.0, 1.0, 1.0); }
	Light(const vec &p, const vec &i, const vec &c):position(p),intensity(i), color(c){ }
};

void OnDisplay();
void OnKeyboard(unsigned char key, int x, int y);
void Render ();
Color shade ( Light * light, Obj * object, mat * e, mat * d, Color color, Color background,
                int level);

const int nChars = ((int)(H * ASPECT)) * H * 3;

std::vector<char> frame(nChars, 255);
Color foregroundColor;//It's the most interesting design I've ever met!!!

void initGLUT(int argc, char** argv, Window& window) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(window.width, window.height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Assignment 3");

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_FLAT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glutDisplayFunc(OnDisplay);
	glutKeyboardFunc(OnKeyboard);
}

void SetBackColor ( Color u) {
	foregroundColor = u;
}

void SetPixel ( Window &window, int i, int j) {
	if (i >= window.width || j >= window.height)
		return;

	unsigned int index = 3 * (j * window.width + i);
	frame[index] = (int)(255 * (foregroundColor.r));
	frame[index + 1] = (int)(255 * (foregroundColor.g));
	frame[index + 2] = (int)(255 * (foregroundColor.b));
}

void QuitX() {
	exit(0);
}

// Allocates and creates a rotation matrix
mat *rotate(double Vx, double Vy, double Vz, double angle)

{
	mat *I, *J, *V;

	I = (mat *)malloc(sizeof(mat));
	J = (mat *)malloc(sizeof(mat));
	V = (mat *)malloc(sizeof(mat));

	mat_alloc(I, 3, 3);
	mat_alloc(J, 3, 3);
	mat_alloc(V, 3, 1);

	I = dmat_identity(I);
	J = dmat_init(J, 0.0);

	(*V).m[1][1] = Vx;
	(*V).m[2][1] = Vy;
	(*V).m[3][1] = Vz;

	V = dmat_normalize(V);

	(*J).m[2][3] = -(*V).m[1][1];
	(*J).m[3][2] = (*V).m[1][1];

	(*J).m[1][3] = (*V).m[2][1];
	(*J).m[3][1] = -(*V).m[2][1];

	(*J).m[1][2] = -(*V).m[3][1];
	(*J).m[2][1] = (*V).m[3][1];

	mat* ret = to_homogeneous(dmat_add(I, dmat_add(dmat_scalar_mult(J, sin(angle)), dmat_scalar_mult(dmat_mult(J, J), 1.0 - cos(angle)))), 1.0);
	delete_dmatrix(I);
	delete_dmatrix(J);
	delete_dmatrix(V);
	return ret;
}

// Allocates and creates a translation matrix
mat *translate(double Tx, double Ty, double Tz)

{
	mat *T;

	T = (mat *)malloc(sizeof(mat));
	mat_alloc(T, 4, 4);

	T = dmat_identity(T);

	(*T).m[1][4] = Tx;
	(*T).m[2][4] = Ty;
	(*T).m[3][4] = Tz;

	return T;
}

// Allocates and creates a scale matrix
mat *scale(double Sx, double Sy, double Sz)

{
	mat *S;

	S = (mat *)malloc(sizeof(mat));
	mat_alloc(S, 4, 4);

	S = dmat_identity(S);

	(*S).m[1][1] = Sx;
	(*S).m[2][2] = Sy;
	(*S).m[3][3] = Sz;

	return S;
}

Window *build_window(Window *window, int height, double aspect, int maxDepth) {

	window->height = height;
	window->width = (int)(aspect * height);
	window->maxDepth = maxDepth;

	window->fov = 90;
	//the write ivory as back
	window->backgroundColor = vec(250.0 / 255.0, 255.0 / 255.0, 240.0 / 255.0);
	window->maxDepth = 3;
	window->bias = 0.00001;

	return(window);
}

camera_t *build_camera(camera_t *Camera, Window *Window) {

	mat_alloc(&Camera->E, 4, 1);

	Camera->E.m[X][1] = Ex;
	Camera->E.m[Y][1] = Ey;
	Camera->E.m[Z][1] = Ez;
	Camera->E.m[4][1] = 1.0;

	mat_alloc(&Camera->G, 4, 1);

	Camera->G.m[X][1] = Gx;
	Camera->G.m[Y][1] = Gy;
	Camera->G.m[Z][1] = Gz;
	Camera->G.m[4][1] = 1.0;

	mat_alloc(&Camera->n, 4, 1);
	Camera->n = *dmat_normalize(dmat_sub(&Camera->E, &Camera->G));
	Camera->n.l = 3;

	mat_alloc(&Camera->UP, 4, 1);

	Camera->UP.l = 3;

	Camera->UP.m[X][1] = UPx;
	Camera->UP.m[Y][1] = UPy;
	Camera->UP.m[Z][1] = UPz;
	Camera->UP.m[4][1] = 1.0;

	mat_alloc(&Camera->u, 4, 1);

	Camera->u = *dmat_normalize(dcross_product(&Camera->UP, &Camera->n));
	Camera->v = *dmat_normalize(dcross_product(&Camera->n, &Camera->u));

	return(Camera);
}

mat *intersection_coordinates(mat *e, mat *direction, double t) {

	mat *intersection;

	intersection = (mat *)malloc(sizeof(mat));
	mat_alloc(intersection, 4, 1);

	intersection->m[X][1] = e->m[X][1] + direction->m[X][1] * t;
	intersection->m[Y][1] = e->m[Y][1] + direction->m[Y][1] * t;
	intersection->m[Z][1] = e->m[Z][1] + direction->m[Z][1] * t;
	intersection->m[4][1] = 1.0;

	return intersection;
}

double infinite_plane_intersection(mat *e, mat *d) {

	double t;

	if (fabs(d->m[Z][1]) < EPSILON) {
		t = -1.0;
	}
	else {
		t = -e->m[Z][1] / d->m[Z][1];
		if (t <= 0.0) {
			t = -1.0;
		}
		else {
			t = -1.0*e->m[Z][1] / d->m[Z][1];
		}
	}
	return t;
}

double plane_intersection(mat *e, mat *d) {

	double t;
	mat *intersection;

	if (fabs(d->m[Z][1]) < EPSILON) {
		t = -1.0;
	}
	else {
		t = -e->m[Z][1] / d->m[Z][1];
		if (t <= 0.0) {
			t = -1.0;
		}
		else {
			intersection = intersection_coordinates(e, d, t);
			if ((fabs(intersection->m[X][1]) > 1.0) || (fabs(intersection->m[Y][1]) > 1.0)) {
				t = -1.0;
			}
			delete_dmatrix(intersection);
		}
	}
	return t;
}

double solve_quadratic(const double &a, const double &b, const double &c, double &t1, double &t2) {

	double discriminant;

	discriminant = b * b - 4 * a * c;
	if (discriminant < 0.0) {
		return -1.0;
	}
	else {
		if (discriminant < EPSILON) {
			t1 = t2 = - 0.5 * b / a;
		}
		else {
			double q = (b > 0) ?
			          -0.5 * (b + sqrt(discriminant)) :
			          -0.5 * (b - sqrt(discriminant));
			t1 = q / a;
			t2 = c / q;
		}if (t1 > t2)std::swap(t1, t2);
		return 1.0;
	}
}

double sphere_intersection(mat *e, mat *d) {

	double a = ddot_product(d, d);
	vec ee = toVec(*e, 3);
	vec dd = toVec(*d, 3);
	double b = ee * dd;
	double c = ee * ee - 1.0;
	double t1, t2;
	return solve_quadratic(a, b, c, t1, t2);
}

mat *sphere_normal(mat *intersection) {

	mat *normal;

	normal = (mat *)malloc(sizeof(mat));
	mat_alloc(normal, 4, 1);

	normal->m[X][1] = intersection->m[X][1];
	normal->m[Y][1] = intersection->m[Y][1];
	normal->m[Z][1] = intersection->m[Z][1];
	normal->m[4][1] = 0.0;

	return dmat_normalize(normal);
}

mat *plane_normal(mat *intersection) {

	mat *normal;

	normal = (mat *)malloc(sizeof(mat));
	mat_alloc(normal, 4, 1);

	normal->m[X][1] = 0.0;
	normal->m[Y][1] = 0.0;
	normal->m[Z][1] = 1.0;
	normal->m[4][1] = 0.0;

	return dmat_normalize(normal);
}

int find_min_hit_time(double t0[N_OBJECTS]) {

	double min_t = std::numeric_limits<double>::max();
	int position = -1;

	for (int i = 0; i < 1024; i++) {
		if (t0[i] != -1.0) {
			if (t0[i] < min_t) {
				min_t = t0[i];
				position = i;
			}
		}
	}
	return position;
}

mat *ray_direction(camera_t *Camera, Window *Window, double height, double width, double i, double j) {

	int k;
	mat *d;

	d = (mat *)malloc(sizeof(mat));
	mat_alloc(d, 3, 1);

	for (k = 1; k <= 3; k++) {
		d->m[k][1] = -1.0*Near*Camera->n.m[k][1] + width * (2.0*i / Window->width - 1.0)*Camera->u.m[k][1] + height * (2.0*j / Window->height - 1.0)*Camera->v.m[k][1];
	}

	mat* ret = to_homogeneous(d, 0.0);
	delete_dmatrix(d);
	return ret;
}

mat *vector_to_light_source(mat *intersection, mat *light_position) {

	return dmat_normalize(dmat_sub(light_position, intersection));
}

mat *vector_to_center_of_projection(mat *intersection, mat *e) {

	return dmat_normalize(dmat_sub(e, intersection));
}

mat *vector_to_specular_reflection(mat *N, mat *S) {

	return dmat_normalize(dmat_add(dmat_scalar_mult(S, -1.0), dmat_scalar_mult(N, 2.0*ddot_product(N, S))));
}

//int shadowed(mat *e, mat *d) {
//
//	int h, k;
//	double t0[N_OBJECTS];
//
//	for (k = 0; k < nobjects; k++) {
//		t0[k] = (object[k].intersection)(dmat_mult(&object[k].Minv, e), dmat_normalize(dmat_mult(&object[k].Minv, d)));
//	}
//	h = find_min_hit_time(t0);
//	return h != -1;
//}

//----------------------------------------------------
//
//class version
//
//----------------------------------------------------
class Sphere:public Obj
{
public:
	vec center;
	double radius, radius2;

	Sphere(const vec &c, const double &r) : center(c), radius(r), radius2(r * r) {}
	bool intersect(const vec &orig, const vec &dir, double &tnear, int &index, vec2 &uv) const
	{
		// analytic solution
		vec L = orig - center;
		double a = dir * dir;
		double b = 2 * (dir * L);
		double c = (L * L) - radius2;
		double t0, t1;
		double num = solve_quadratic(a, b, c, t0, t1);
		if (num < 0) return false;
		if (t0 < 0) t0 = t1;
		if (t0 < 0) return false;
		tnear = t0;

		return true;
	}
	void getSurfaceProperties(const vec &P, const vec &I, const int &index, const vec2 &uv, vec &N, vec2 &st) const
	{ N = (P - center).normalized(); }
};

bool rayTriangleIntersect(
		const vec &v0, const vec &v1, const vec &v2,
		const vec &orig, const vec &dir,
		double &tnear, double &u, double &v)
{
	vec edge1 = v1 - v0;
	vec edge2 = v2 - v0;
	vec pvec = dir.cross(edge2);
	double det = (edge1*pvec);
	if (det == 0 || det < 0) return false;

	vec tvec = orig - v0;
	u = (tvec * pvec);
	if (u < 0 || u > det) return false;

	vec qvec = tvec.cross(edge1);
	v = (dir * qvec);
	if (v < 0 || u + v > det) return false;

	double invDet = 1 / det;

	tnear = (edge2 * qvec) * invDet;
	u *= invDet;
	v *= invDet;

	return true;
}

class MeshTriangle : public Obj
{
public:
	MeshTriangle(
			const vec *verts,
			const uint32_t *vertsIndex,
			const uint32_t &numTris,
			const vec2 *st)
	{
		uint32_t maxIndex = 0;
		for (uint32_t i = 0; i < numTris * 3; ++i)
			if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
		maxIndex += 1;
		vertices = std::unique_ptr<vec[]>(new vec[maxIndex]);
		memcpy(vertices.get(), verts, sizeof(vec) * maxIndex);
		vertexIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
		memcpy(vertexIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
		numTriangles = numTris;
		stCoordinates = std::unique_ptr<vec2[]>(new vec2[maxIndex]);
		memcpy(stCoordinates.get(), st, sizeof(vec2) * maxIndex);
	}

	bool intersect(const vec &orig, const vec &dir, double &tnear, uint32_t &index, vec2 &uv) const
	{
		bool intersect = false;
		for (int k = 0; k < numTriangles; ++k) {
			const vec & v0 = vertices[vertexIndex[k * 3]];
			const vec & v1 = vertices[vertexIndex[k * 3 + 1]];
			const vec & v2 = vertices[vertexIndex[k * 3 + 2]];
			double t, u, v;
			if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
				tnear = t;
				uv.m[1][1] = u;
				uv.m[2][1] = v;
				index = k;
				intersect |= true;
			}
		}

		return intersect;
	}

	void getSurfaceProperties(const vec &P, const vec &I, const int &index, const vec2 &uv, vec &N, vec2 &st) const
	{
		const vec &v0 = vertices[vertexIndex[index * 3]];
		const vec &v1 = vertices[vertexIndex[index * 3 + 1]];
		const vec &v2 = vertices[vertexIndex[index * 3 + 2]];
		vec e0 = (v1 - v0).normalized();
		vec e1 = (v2 - v1).normalized();
		N = (e0.cross(e1)).normalized();
		const vec2 &st0 = stCoordinates[vertexIndex[index * 3]];
		const vec2 &st1 = stCoordinates[vertexIndex[index * 3 + 1]];
		const vec2 &st2 = stCoordinates[vertexIndex[index * 3 + 2]];
		st = st0 * (1 - uv[1] - uv[2]) + st1 * uv[1] + st2 * uv[2];
	}

	vec evalDiffuseColor(const vec2 &st) const
	{
		double scale = 5;
		double pattern = (fmodf(st.m[1][1] * scale, 1) > 0.5) ^ (fmodf(st.m[2][1] * scale, 1) > 0.5);
		return mix(vec(0.815, 0.235, 0.031), vec(0.937, 0.937, 0.231), pattern);
	}

	std::unique_ptr<vec[]> vertices;
	uint32_t numTriangles;
	std::unique_ptr<uint32_t[]> vertexIndex;
	std::unique_ptr<vec2[]> stCoordinates;
};

/*
Obj build_object(int object_type, mat *M, Color ambient_color, Color diffuse_color, Color specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) {

	Obj object;

	object.type = object_type;

	mat_alloc(&object.M, 4, 4);
	object.M = *dmat_duplicate(M);

	mat_alloc(&object.Minv, 4, 4);
	object.Minv = *dmat_inverse(&object.M);

	object.specular_color = specular_color;
	object.diffuse_color = diffuse_color;
	object.ambient_color = ambient_color;

	object.specular_coeff = specular_coeff;
	object.diffuse_coeff = diffuse_coeff;
	object.ambient_coeff = ambient_coeff;

	object.f = f;
	object.reflectivity = reflectivity;

	switch (object_type) {

	case SPHERE:

		object.intersection = &sphere_intersection;
		object.normal = &sphere_normal;
		break;

	case PLANE:
		object.intersection = &plane_intersection;
		object.normal = &plane_normal;
		break;

	case INFINITE_PLANE:

		object.intersection = &infinite_plane_intersection;
		object.normal = &plane_normal;
		break;

	default:
		break;

	}
	nobjects++;
	return(object);
}
*/

bool trace(
		const vec &orig, const vec &dir,
		const std::vector<std::unique_ptr<Obj>> &objects,
		double &tNear, uint32_t &index, vec &uv, Obj **hitObject)
{
	*hitObject = nullptr;
	for (uint32_t k = 0; k < objects.size(); ++k) {
		double tNearK = kInfinity;
		int indexK;
		vec2 uvK;
		if (objects[k]->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear) {
			*hitObject = objects[k].get();
			tNear = tNearK;
			index = indexK;
			uv = uvK.toVec3();
		}
	}

	return (*hitObject != nullptr);
}

vec castRay(
		const vec &orig, const vec &dir,
		const std::vector<std::unique_ptr<Obj>> &objects,
		const std::vector<std::unique_ptr<Light>> &lights,
		const Window &options,
		uint32_t depth,
		bool test = false)
{
	if (depth > options.maxDepth) {
		return options.backgroundColor;
	}

	vec hitColor = options.backgroundColor;
	double tnear = kInfinity;
	vec2 uv;
	uint32_t index = 0;
	Obj *hitObject = nullptr;
	if (trace(orig, dir, objects, tnear, index, uv, &hitObject)) {
		vec hitPoint = orig + dir * tnear;
		vec N; // normal
		vec2 st; // st coordinates
		hitObject->getSurfaceProperties(hitPoint, dir, index, uv, N, st);
		vec tmp = hitPoint;
		switch (hitObject->type) {
			case 0:
			{
				vec reflectionDirection = (reflect(dir, N).normalized());
				vec refractionDirection = (refract(dir, N, hitObject->ior).normalized());
				vec reflectionRayOrig = ((reflectionDirection * N) < 0) ?
				                          hitPoint - N * options.bias :
				                          hitPoint + N * options.bias;
				vec refractionRayOrig = ((refractionDirection * N) < 0) ?
				                          hitPoint - N * options.bias :
				                          hitPoint + N * options.bias;
				vec reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1, 1);
				vec refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1, 1);
				double kr;
				fresnel(dir, N, hitObject->ior, kr);
				hitColor = reflectionColor * kr + refractionColor * (1 - kr);
				break;
			}
			case 1:
			{
				double kr;
				fresnel(dir, N, hitObject->ior, kr);
				vec reflectionDirection = reflect(dir, N);
				vec reflectionRayOrig = ((reflectionDirection * N) < 0) ?
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
				vec lightAmt, specularColor;
				vec shadowPointOrig = ((dir * N) < 0) ?
				                        hitPoint + N * options.bias :
				                        hitPoint - N * options.bias;
				// [comment]
				// Loop over all lights in the scene and sum their contribution up
				// We also apply the lambert cosine law here though we haven't explained yet what this means.
				// [/comment]
				for (uint32_t i = 0; i < lights.size(); ++i) {
					vec lightDir = lights[i]->position - hitPoint;
					// square of the distance between hitPoint and the light
					double lightDistance2 = (lightDir * lightDir);
					lightDir = (lightDir.normalized());
					double LdotN = std::max(0.0, (lightDir * N));
					Obj *shadowHitObject = nullptr;
					double tNearShadow = kInfinity;
					// is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
					bool inShadow = trace(shadowPointOrig, lightDir, objects, tNearShadow, index, uv, &shadowHitObject) &&
					                tNearShadow * tNearShadow < lightDistance2;
					lightAmt += (1 - inShadow) * lights[i]->intensity * LdotN;
					vec reflectionDirection = reflect(-lightDir, N);
					specularColor += powf(std::max(0.0, -(reflectionDirection * dir)), hitObject->specular_exp) * lights[i]->intensity;
				}
				hitColor = lightAmt.times(hitObject->evalDiffuseColor(st)) * hitObject->density + specularColor * hitObject->specular_coeff;
				break;
			}
		}
	}

	return hitColor;
}

camera_t Camera;
Window win;
Light light;
Color background;
std::vector<std::unique_ptr<Obj>> objects;
std::vector<std::unique_ptr<Light>> lights;

void render( )
		//const std::vector<Light> &lights)
{
	vec *framebuffer = new vec[win.width * win.height];
	vec *pix = framebuffer;
	float scale = tan((win.fov * 0.5*M_PI/180));
	float imageAspectRatio = win.width / (float)win.height;
	vec orig(0, 0, 0);
	for (int j = 0; j < win.height; ++j) {
		for (int i = 0; i < win.width; ++i) {
			// generate primary ray direction
			float x = (2 * (i + 0.5) / (float)win.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)win.height) * scale;
			vec dir = (vec(x, y, -1).normalized());
			*(pix++) = castRay(orig, dir, objects, lights, win, 0);
		}
	}

	for (int i = 0; i < win.width; i++) {
		for (int j = 0; j < win.height; j++)
		{
			unsigned int index = 3 * ( j * win.width + i );
			frame[index] = ( int ) ( 255 * framebuffer[i].m[1][1] );
			frame[index + 1] = ( int ) ( 255 * framebuffer[i].m[2][1] );
			frame[index + 2] = ( int ) ( 255 * framebuffer[i].m[3][1] );
		}
	}

	delete [] framebuffer;
}

void OnDisplay() {
	frame.resize(nChars, 255);
	render();
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels((int)(H * ASPECT), H, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte*)&frame[0]);
	glutSwapBuffers();
	glFlush();
}

int main(int argc, char** argv) {
	/* Set the background color */
	background = Color(250.0/255.0, 255.0/255.0, 240.0/255.0);

	/* Set up light position, intensity, and color */
	vec light_position = vec(Lx, Ly, Lz);

	Color light_intensity(1.0, 1.0, 1.0);
	Color light_color(1.0, 1.0, 1.0);
	light = Light(light_position, light_color, light_intensity);

	/* Build display window and synthetic camera */
	win = *build_window(&win, H, ASPECT, 3);
	Camera = *build_camera(&Camera, &win);

	/* Build a sphere */
	mat M = *translate(0.0, 0.0, 0.0);

	Color specular_color = Color(1, 1, 1);
	Color diffuse_color = Color(0.0, 0.0, 1.0);
	Color ambient_color = Color(0, 0, 1);
	double specular_coeff = 0.4;
	double diffuse_coeff = 0.4;
	double ambient_coeff = 0.2;
	double f = 10.0;
	double reflectivity = 0.0;

	//	Sphere(const vec &c, const double &r) : center(c), radius(r), radius2(r * r) {}
	Sphere *sph1 = new Sphere(vec(0, 0, 0), 3);
	sph1->type = 0;
	sph1->diffuse_color = Color(0.6, 0, 0.8);

//	vec verts[4] = {{-5,-3,-6}, {5,-3,-6}, {5,-3,-16}, {-5,-3,-16}};
//	int vertIndex[6] = {0, 1, 3, 1, 2, 3};
//	vec2 st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
//	MeshTriangle *mesh = new MeshTriangle(verts, vertIndex, 2, st);
//	mesh->type = 0;

	lights.push_back(std::unique_ptr<Light>(new Light(vec(-20, 70, 20))));
	lights.push_back(std::unique_ptr<Light>(new Light(vec(30, 50, -12))));
	objects.push_back(std::unique_ptr<Sphere>(sph1));
//	objects.push_back(std::unique_ptr<MeshTriangle>(mesh));

	initGLUT(argc, argv, win);
	glutMainLoop();

	return 0;
}

bool trace(const mat* orig, const mat *dir, double tNear);

void OnKeyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'q':
		QuitX();
		break;
	}
}
