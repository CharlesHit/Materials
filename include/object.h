#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

using namespace std;

struct window_t
{
public:
	int width, height;
	inline window_t (int w, int h ):width(w), height(h)
	{
		width = w;
		height = h;
	}
};

class Ray:public vec{
public:
	vec o, d;
	double t_min, t_max;
	double active;

	/*
	 * Create a new group of active rays
	 * Note: only the sign bit of the active mask will be set, as this is all that's used by blendv
	 * and movemask
	 */
	Ray ( vec o, vec d, double t_min_ = 0, double t_max_ = INFINITY )
		: o ( o ), d ( d ), t_min (( t_min_ ) ), t_max ( ( t_max_ ) ),
		active ( ( -0.f ) )
	{
	}
	vec at ( int t ) const
	{
		return o + t * d;
	}
};

class Camera {
public:
	Camera () { }

	vec E;
	vec u, v, n;//n = E - G		
	vec top_left, screen_du, screen_dv;
	double aspect, fovy; //based on window's size.
	Camera (vec E, vec G, vec UP, double aspect, double fovy):E(E),fovy(60), aspect(1)
	{	
		//TODO: E-G or G-E?
		n = E - G; n = n.normalized();
		u = UP.cross(n).normalized();
		v = n.cross(u).normalized();

		float dim_y = 2.f * sin ( ( fovy / 2.f ) * ( M_PI ) / 180.f );
		float dim_x = dim_y * aspect;
		top_left = n - 0.5f * dim_x * u - 0.5f * dim_y * v;
		screen_du = u * dim_x;
		screen_dv = v * dim_y;
	}

	void generate_rays ( Ray &rays, const vec &samples ) const
	{
		rays.o = this->E;
		rays.d = top_left;
		const auto u_step = samples.m[1][1] * screen_du;
		const auto v_step = samples.m[2][1] * screen_dv;
		rays.d = rays.d + u_step + v_step;
		rays.d.normalized ( );
		rays.t_min = 0;
		rays.t_max = INFINITY;
	}
};

unsigned char frame[windowH*windowW * 3];
int nobjects = 0;

class Colour:public vec {
public:
	double r, g, b; //this stupid setting is for compatibility.
	double* ptr_r = this->m[1].data() + 1;
	double* ptr_g = this->m[2].data() + 1;
	double* ptr_b = this->m[3].data() + 1;

	inline Colour() {}
	inline Colour(double _r, double _g, double _b): r(_r), g(_g), b(_b) {
		*ptr_r = r;
		*ptr_g = g;
		*ptr_b = b;
	}
	inline Colour to_sRGB() const {
		const double a = 0.055f;
		const double b = 1.f / 2.4f;
		Colour srgb;
		for (int i = 0; i < 3; ++i)
			if (this->m[i][1] <= 0.0031308f)srgb.m[i][1] = 12.92f * this->m[i][1];
			else srgb.m[i][1] = (1.f + a) * std::pow( this->m[i][1], b) - a;
		return srgb;
	}
	inline bool isBlack() const {
		return *ptr_r == 0 && *ptr_g == 0 && *ptr_b == 0;
	}
	inline double luminance() const {
		return 0.2126f * *ptr_r + 0.7152f * *ptr_g + 0.0722f * *ptr_b;
	}
};

std::ostream& operator<<(std::ostream &os, const Colour &c){
	os << "Colour: [r = " << *c.ptr_r << ", g = " << *c.ptr_r
		<< ", b = " << *c.ptr_r << "] ";
	return os;
}


struct light_t {
	vec position;
	Colour color;
	Colour intensity;
};

// Create/allocate a light
light_t *build_light ( light_t *light, vec position, Colour color, Colour intensity )
{
	light->position = vec(position);
	light->color.r = color.r;
	light->color.g = color.g;
	light->color.b = color.b;
	light->intensity.r = intensity.r;
	light->intensity.g = intensity.g;
	light->intensity.b = intensity.b;
	return light;
}

vec intersection_coordinates ( vec e, vec direction, double t )
{

	vec intersection = vec ( );
	intersection.m[X][1] = e.m[X][1] + direction.m[X][1] * t;
	intersection.m[Y][1] = e.m[Y][1] + direction.m[Y][1] * t;
	intersection.m[Z][1] = e.m[Z][1] + direction.m[Z][1] * t;
	intersection.m[4][1] = 1.0;

	return intersection;
}

double infinite_plane_intersection ( vec e, vec d )
{

	double t;

	if ( fabs ( d.m[Z][1] ) < EPSILON ) {
		t = -1.0;
	} else {
		t = -e.m[Z][1] / d.m[Z][1];
		if ( t <= 0.0 ) {
			t = -1.0;
		} else {
			t = -1.0*e.m[Z][1] / d.m[Z][1];
		}
	}
	return t;
}

double plane_intersection ( vec e, vec d )
{

	double t;
	vec intersection;

	if ( fabs ( d.m[Z][1] ) < EPSILON ) {
		t = -1.0;
	} else {
		t = -e.m[Z][1] / d.m[Z][1];
		if ( t <= 0.0 ) {
			t = -1.0;
		} else {
			intersection = intersection_coordinates ( e, d, t );
			if ( ( fabs ( intersection.m[X][1] ) > 1.0 ) || ( fabs ( intersection.m[Y][1] ) > 1.0 ) ) {
				t = -1.0;
			}
		}
	}
	return t;
}

double solve_quadratic ( double a, double b, double c )
{

	double discriminant, t1, t2, min;

	discriminant = b * b - a * c;
	if ( discriminant < 0.0 ) {
		return -1.0;
	} else {
		if ( discriminant < EPSILON ) {
			return -b / a;
		} else {
			t1 = -b / a - sqrtl ( discriminant ) / a;
			t2 = -b / a + sqrtl ( discriminant ) / a;

			if ( t1 < t2 ) {
				min = t1;
			} else {
				min = t2;
			}

			if ( min > EPSILON ) {
				return min;
			} else {
				return -1.0;
			}
		}
	}
}

double sphere_intersection ( vec e, vec d )
{

	double a = d*d;
	e.set ( 4, 1, 0 );
	d.set ( 4, 1, 0 );
	double b = e * d;
	double c = e * e - 1.0;

	return solve_quadratic ( a, b, c );
}

vec sphere_normal ( dmatrix_t intersection )
{

	vec normal = vec ( intersection.m[X][1], intersection.m[Y][1], intersection.m[Z][1], 0 );

	return normal.normalized ( );
}

vec plane_normal ( dmatrix_t intersection )
{

	vec normal = vec ( 0, 0, 1, 0 );

	return normal.normalized ( );
}

int find_min_hit_time ( double t0[N_OBJECTS] )
{

	double min_t = std::numeric_limits<double>::max ( );
	int position = -1;

	for ( int i = 0; i < nobjects; i++ ) {
		if ( t0[i] != -1.0 ) {
			if ( t0[i] < min_t ) {
				min_t = t0[i];
				position = i;
			}
		}
	}
	return position;
}

dmatrix_t *ray_direction ( Camera camera, window_t *Window, double height, double width, double i, double j )
{

	int k;
	dmatrix_t *d;

	d = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( d, 3, 1 );

	for ( k = 1; k <= 3; k++ ) {
		d->m[k][1] = -1.0*Near*camera.n.m[k][1] + width * ( 2.0*i / Window->width - 1.0 )*camera.u.m[k][1] + height * ( 2.0*j / Window->height - 1.0 )*camera.v.m[k][1];
	}

	dmatrix_t* ret = to_homogeneous ( d, 0.0 );
	delete_dmatrix ( d );
	return ret;
}

dmatrix_t *vector_to_light_source ( dmatrix_t *intersection, dmatrix_t *light_position )
{

	return dmat_normalize ( dmat_sub ( light_position, intersection ) );
}

dmatrix_t *vector_to_center_of_projection ( dmatrix_t *intersection, dmatrix_t *e )
{

	return dmat_normalize ( dmat_sub ( e, intersection ) );
}

dmatrix_t vector_to_specular_reflection ( vec N, vec S )
{

	return ((-S) + N.scalared(2.0* ( N.dot(S) ) ) ).normalized().toMat();
}

class object_t
{
public:
	int type;
	dmatrix_t M, Minv;
	Colour specular_color, diffuse_color, ambient_color;
	double density, reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f;

	double ( *intersection )( vec, vec );
	vec ( *normal )( dmatrix_t );

	inline object_t ( ){ }

	object_t
	( int object_type,
	  dmatrix_t m,
	  Colour ambient,
	  Colour diffuse,
	  Colour specular,
	  double ambient_coeff,
	  double diffuse_coeff,
	  double specular_coeff,
	  double f,
	  double reflectivity ) :type ( object_type ), M ( m ), specular_color ( specular ), diffuse_color ( diffuse ), ambient_color ( ambient ), specular_coeff ( specular_coeff ), diffuse_coeff ( diffuse_coeff ), ambient_coeff ( ambient_coeff ), f ( f ), reflectivity ( reflectivity )
	{
		Minv = dmat_inverse ( m );//Todo: is this really useful?

		switch ( object_type ) {

			case SPHERE:

				intersection = &sphere_intersection;
				normal = sphere_normal;
				break;

			case PLANE:
				intersection = &plane_intersection;
				normal = plane_normal;
				break;

			case INFINITE_PLANE:

				intersection = &infinite_plane_intersection;
				normal = plane_normal;
				break;

			default:
				break;

		}
		nobjects++;
	}
};

Colour foregroundColor;
object_t object[N_OBJECTS];

Colour shade ( light_t *light, object_t *object, vec e, dmatrix_t *d, Colour color, Colour background, int level )
{

	color.r = 0.0;
	color.g = 0.0;
	color.b = 0.0;

	return color;
}
int find_min_hit_time ( double t0[N_OBJECTS] );
int shadowed ( vec e, vec d )
{

	int h, k;
	double t0[N_OBJECTS];

	for ( k = 0; k < nobjects; k++ ) {
		t0[k] = ( object[k].intersection )( dmat_mult (&object[k].Minv, e ),  dmat_mult ( &object[k].Minv, d ).normalized() );
	}
	h = find_min_hit_time ( t0 );
	return h != -1;
}

object_t *build_object ( int object_type, dmatrix_t *M, Colour ambient_color, Colour diffuse_color, Colour specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity )
{

	object_t *object;

	object = (object_t*) malloc ( sizeof ( *object ) );
	object->type = object_type;

	dmat_alloc ( &object->M, 4, 4 );
	object->M.set ( *M );

	dmat_alloc ( &object->Minv, 4, 4 );
	object->Minv = dmat_inverse ( object->M );

	object->specular_color = Colour ( specular_color );
	object->diffuse_color = Colour ( diffuse_color );
	object->ambient_color = Colour ( ambient_color );

	object->specular_coeff = specular_coeff;
	object->diffuse_coeff = diffuse_coeff;
	object->ambient_coeff = ambient_coeff;

	object->f = f;
	object->reflectivity = reflectivity;

	switch ( object_type ) {

		case SPHERE:

			object->intersection = &sphere_intersection;
			object->normal = sphere_normal;
			break;

		case PLANE:
			object->intersection = &plane_intersection;
			object->normal = plane_normal;
			break;

		case INFINITE_PLANE:

			object->intersection = &infinite_plane_intersection;
			object->normal = plane_normal;
			break;

		default:
			break;

	}
	nobjects++;
	return( object );
}