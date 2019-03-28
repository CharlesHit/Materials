#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

struct window_t
{
	int width, height;
	inline window_t ( )
	{
		width = windowW;
		height = windowH;
	}
};

class Camera {
public:
	vec UP;
	vec E;
	vec G;
	vec n;//n = GE = E - G
	vec u, v;


	Camera ( )
	{	
		E = Eye;
		G = Gaze;
		n = E - G; n.normalized();
		UP = Up;
		u = UP; u.cross(n).normalized();
		v = n; v.cross(u).normalized();
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
};

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