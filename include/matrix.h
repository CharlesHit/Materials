
/*
 * matrix.h, ver 2.0. 
 *
 * author: Charles.H.
 * date: 19-3-26
 */

#pragma once
#include <cassert>
#include <cmath>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#define quote(x) #x

class dmatrix_t
{
public:
	vector<vector<double>> m;
	int l, c;
	
	dmatrix_t():l(0), c(0)
	{
		m.resize(2, vector<double>(2, 0));
		l = 1;
		c = 1;
	}

	dmatrix_t ( int sizex, int sizey ):l(sizex), c(sizey)
	{
		m.resize(sizex+1, vector<double>(sizey+1, 0));
		l = sizex;
		c = sizey;
	}

	void set ( dmatrix_t A )
	{
		m = A.m;
		l = A.l;
		c = A.c;
	}

	void set ( int x, int y, double value )
	{
		m[x+1][y+1] = value;
	}

	~dmatrix_t ( )
	{
		//thanks vector! No need to free
	}

	void print ( )
	{
		int i, j;
		for ( i = 1; i <= l; i++ ) {
			for ( j = 1; j <= c; j++ ) {
				printf ( "%7.4f ", m[i][j] );
			}
			printf ( "\n" );
		}
	}

	void scalar_mult ( int scalar )
	{
		for ( int i = 1; i <= l; i++ )
			for ( int j = 1; j <= c; j++ )
				m[i][j] = m[i][j] * scalar;
	}

	//I used wanna to design two way: B.identify(A), but finally it doesn't be better than this one:
	dmatrix_t identify ( )
	{
		dmatrix_t temp;
		temp.m = this->m;
		temp.c = this->c;
		temp.l = this->l;
		for ( int i = 0; i <= l; i++ )
			for ( int j = 0; j <= c; j++ )
				if(i == j)temp.m[i][j] = 1;
				else temp.m[i][j] = 0;
		return temp;
	}
};

// a 3*1 or 4*1 vector. In 4*1-case, the last one actually don't attend
// calculations.
class vec
{
public:
	vector<vector<double>> m;
	int l, c;

	inline vec ( ):l(3), c(1){
		m.resize(5, vector<double>(2, 0));
		m[4][1] = 1;
	}
	inline vec ( double x, double y, double z ):l(3), c(1)
	{
		m.resize(5, vector<double>(2, 0));
		m[4][1] = 1;
		m[1][1] = x;
		m[2][1] = y;
		m[3][1] = z;
	}
	inline vec ( double x, double y, double z, double buffer ):l(4), c(1)
	{
		m.resize(5, vector<double>(2, 0));
		m[1][1] = x;
		m[2][1] = y;
		m[3][1] = z;
		m[4][1] = buffer;
	}
	inline vec ( const vec &v )
	{
		m = v.m;
		l = v.l;
		c = v.c;
	}

	inline ~vec ()
	{
		//thanks vector again!
	}

	inline double dot ( const vec &v ) const
	{
		return m[1][1] * v.m[1][1] + m[2][1] * v.m[2][1]  + m[3][1] * v.m[3][1];
	}
	inline vec cross ( const vec &v ) const
	{
		return vec { m[2][1] * v.m[3][1] - m[3][1] * v.m[2][1],
				   m[3][1] * v.m[1][1] - m[1][1] * v.m[3][1],
				   m[1][1] * v.m[2][1] - m[2][1] * v.m[1][1] };
	}
	double length_sqr ( )
	{
		return m[1][1] * m[1][1] + m[2][1] * m[2][1] + m[3][1] * m[3][1];
	}
	double length ( )
	{
		return sqrt ( length_sqr ( ) );
	}
	inline vec normalized ( )
	{
		double s = 1.f / length ( );
		return vec { m[1][1] * s, m[2][1] * s, m[3][1] * s };
	}
	inline void set ( int x, int y, double value )
	{
		m[x+1][y+1] = value;
	}					

	void print ( )
	{
		cout << "x:"<< m[1][1] << " y:" << m[2][1] << " z:" << m[3][1];
		if(l==4) cout<< " buffer" << m[4][1] <<endl;
		else if (l ==3) cout <<endl;
	}
};

inline vec operator+( const vec &a, const vec &b )
{
	return vec { a.m[1][1] + b.m[1][1], a.m[2][1] + b.m[2][1],
			   a.m[3][1] + b.m[3][1] };
}
inline vec operator-( const vec &a, const vec &b )
{
	return vec { a.m[1][1] - b.m[1][1], a.m[2][1] - b.m[2][1],
			   a.m[3][1] - b.m[3][1] };
}
inline vec operator*( const vec &v, float s )
{
	return vec { v.m[1][1] * s, v.m[2][1] * s, v.m[3][1] * s };
}
inline vec operator*( float s, const vec &v )
{
	return v * s;
}
inline vec operator/( const vec &v, float s )
{
	return v * ( 1 / s );
}
inline vec operator/( float s, const vec &v )
{
	return vec { s / v.m[1][1], s / v.m[2][1], s / v.m[3][1] };
}
inline vec operator/( const vec &a, const vec &b )
{
	return vec { a.m[1][1] / b.m[1][1], a.m[2][1] / b.m[2][1],
			   a.m[3][1] / b.m[3][1] };
}
inline vec operator-( const vec &v )
{
	return vec { -v.m[1][1], -v.m[2][1], -v.m[3][1] };
}
inline std::ostream &operator<<( std::ostream &os, const vec &v )
{
	os << "vec:\n\tx = " << v.m[1][1] << "\n\ty = " << v.m[2][1]
		<< "\n\tz = " << v.m[3][1];
	return os;
}

// vec3/4 actually is same as 3. It won't have anyother opreation in this
// project.
typedef vec vec3;
typedef vec vec4;

void error ( const char error_text[] )
{
	fprintf ( stderr, "Run-time error...\n" );
	fprintf ( stderr, "%s\n", error_text );
	exit ( 1 );
}

double **dmatrix ( int nrl, int nrh, int ncl, int nch )

{
	int i;
	double **m;

	m = (double **) malloc ( (unsigned) ( nrh - nrl + 1 ) * sizeof ( double ) );
	if ( !m ) {
		error ( "MATRIX.H: allocation failure" );
	}
	m -= nrl;

	for ( i = nrl; i <= nrh; i++ ) {
		m[i] = (double *) malloc ( (unsigned) ( nch - ncl + 1 ) * sizeof ( double ) );
		if ( !m[i] ) {
			error ( "MATRIX.H: allocation failure" );
		}
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch )

{
	int i;

	for ( i = nrh; i >= nrl; i-- ) {
		free ( (char *) ( m[i] + ncl ) );
	}
	free ( (char *) ( m + nrl ) );
}

void delete_dmatrix ( dmatrix_t *A )

{
	// free_dmatrix ( A->m, 1, A->l, 1, A->c );
	// free ( A );
}

void dmat_alloc ( dmatrix_t *A, int l, int c )

{
	// ( *A ).m = dmatrix ( 1, l, 1, c );
	// ( *A ).l = l;
	// ( *A ).c = c;
}

dmatrix_t *dmat_init ( dmatrix_t *A, double a )

{
	int i, j;

	for ( i = 1; i <= ( *A ).l; i++ ) {
		for ( j = 1; j <= ( *A ).c; j++ ) {
			( *A ).m[i][j] = a;
		}
	}
	return A;
}

dmatrix_t *dmat_identity ( dmatrix_t *A )

{
	int i, j;

	for ( i = 1; i <= ( *A ).l; i++ ) {
		for ( j = 1; j <= ( *A ).c; j++ ) {
			if ( i == j ) {
				( *A ).m[i][j] = 1.0;
			} else {
				( *A ).m[i][j] = 0.0;
			}
		}
	}
	return A;
}

dmatrix_t *dmat_scalar_mult(dmatrix_t *A, double a)

{
	dmatrix_t *B;
	int i, j;

	B = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(B, (*A).l, (*A).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = (*A).m[i][j] * a;
		}
	}
	return B;
}

dmatrix_t *dmat_mult ( dmatrix_t *A, dmatrix_t *B )

{
	dmatrix_t *C;
	double s;
	int i, j, k;

	if ( ( *A ).c != ( *B ).l ) {
		error ( "MATRIX.H: incompatible matrix sizes" );
	}
	C = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( C, ( *A ).l, ( *B ).c );

	for ( i = 1; i <= ( *C ).l; i++ ) {
		for ( j = 1; j <= ( *C ).c; j++ ) {
			for ( s = 0.0, k = 1; k <= ( *A ).c; k++ ) {
				s += ( *A ).m[i][k] * ( *B ).m[k][j];
			}
			( *C ).m[i][j] = s;
		}
	}
	return C;
}

dmatrix_t *dmat_add ( dmatrix_t *A, dmatrix_t *B )
{
	dmatrix_t *C;
	int i, j;

	if ( ( *A ).l != ( *B ).l || ( *A ).c != ( *B ).c ) {
		error ( "MATRIX.H: incompatible matrix sizes" );
	}
	C = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( C, ( *A ).l, ( *A ).c );

	for ( i = 1; i <= ( *C ).l; i++ ) {
		for ( j = 1; j <= ( *C ).c; j++ ) {
			( *C ).m[i][j] = ( *A ).m[i][j] + ( *B ).m[i][j];
		}
	}
	return C;
}

dmatrix_t *dmat_sub ( dmatrix_t *A, dmatrix_t *B )
{
	dmatrix_t *C;
	int i, j;

	if ( ( *A ).l != ( *B ).l || ( *A ).c != ( *B ).c ) {
		error ( "MATRIX.H: incompatible matrix sizes" );
	}
	C = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( C, ( *A ).l, ( *A ).c );

	for ( i = 1; i <= ( *C ).l; i++ ) {
		for ( j = 1; j <= ( *C ).c; j++ ) {
			( *C ).m[i][j] = ( *A ).m[i][j] - ( *B ).m[i][j];
		}
	}
	return C;
}

dmatrix_t *dmat_transpose ( dmatrix_t *A )

{
	dmatrix_t *B;
	int i, j;

	B = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( B, ( *A ).c, ( *A ).l );

	for ( i = 1; i <= ( *A ).l; i++ ) {
		for ( j = 1; j <= ( *A ).c; j++ ) {
			( *B ).m[j][i] = ( *A ).m[i][j];
		}
	}
	return B;
}

double ddot_product ( dmatrix_t *A, dmatrix_t *B )

{
	dmatrix_t *C;

	C = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );

	if ( ( *A ).c == ( *B ).c && ( *A ).l == 1 && ( *B ).l == 1 ) {
		C = dmat_mult ( A, dmat_transpose ( B ) );
	} else if ( ( *A ).c == ( *B ).l && ( *A ).l == 1 && ( *B ).c == 1 ) {
		C = dmat_mult ( A, B );
	} else if ( ( *A ).l == ( *B ).c && ( *A ).c == 1 && ( *B ).l == 1 ) {
		C = dmat_mult ( B, A );
	} else if ( ( *A ).l == ( *B ).l && ( *A ).c == 1 && ( *B ).c == 1 ) {
		C = dmat_mult ( dmat_transpose ( A ), B );
	} else
		error ( "MATRIX.H: Incompatible matrix sizes" );
	return ( *C ).m[1][1];
}

dmatrix_t *cross_product_matrix ( dmatrix_t *A, dmatrix_t *B )

{
	int i, j;
	dmatrix_t *C;

	C = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( C, ( *A ).l, ( *A ).l );

	for ( i = 1; i <= ( *A ).l; i++ ) {
		for ( j = 1; j <= ( *A ).l; j++ ) {
			if ( i == 1 ) {
				( *C ).m[i][j] = 1.0;
			} else if ( i == 2 ) {
				( *C ).m[i][j] = ( *A ).m[j][1];
			} else {
				( *C ).m[i][j] = ( *B ).m[j][1];
			}
		}
	}
	return C;
}

dmatrix_t *sub_matrix ( dmatrix_t *A, int r, int c )

{
	int i, j, k, l;
	dmatrix_t *B;

	if ( r < 1 || r > ( *A ).l || c < 1 || c > ( *A ).c || ( *A ).c < 2 || ( *A ).l < 2 ) {
		error ( "MATRIX.H: erroneous indices" );
	}
	B = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( B, ( *A ).l - 1, ( *A ).c - 1 );

	for ( i = 1, k = 1; i <= ( *A ).l; i++ ) {
		for ( j = 1, l = 1; j <= ( *A ).c; j++ ) {
			if ( j != c && i != r ) {
				( *B ).m[k][l] = ( *A ).m[i][j];
			}
			if ( j != c )
				l++;
		}
		if ( i != r )
			k++;
	}
	return B;
}

double determinant ( dmatrix_t *A )

{
	int i;
	double det;

	if ( ( *A ).l < 1 || ( *A ).c < 1 ) {
		error ( "MATRIX.H: erroneous matrix size" );
	} else if ( ( *A ).l != ( *A ).c ) {
		error ( "MATRIX.H: not a square matrix" );
	} else if ( ( *A ).l == 1 ) {
		det = ( *A ).m[1][1];
	} else {
		det = 0.0;
		for ( i = 1; i <= ( *A ).c; i++ ) {
			det +=
				pow ( -1.0, i + 1.0 ) * ( *A ).m[1][i] * determinant ( sub_matrix ( A, 1, i ) );
		}
	}
	return det;
}

dmatrix_t *cofactor ( dmatrix_t *A )

{
	int i, j;
	dmatrix_t *B;

	B = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	dmat_alloc ( B, ( *A ).l, ( *A ).c );

	for ( i = 1; i <= ( *A ).l; i++ ) {
		for ( j = 1; j <= ( *A ).c; j++ ) {
			( *B ).m[i][j] = pow ( -1.0, i + j ) * determinant ( sub_matrix ( A, i, j ) );
		}
	}
	return B;
}

dmatrix_t *dmat_inverse ( dmatrix_t *A )

{
	( A->scalar_mult ( ( dmat_transpose ( cofactor ( A ) ), 1.0 / determinant ( A ) ) ) );
	return A;
}

dmatrix_t *to_homogeneous ( dmatrix_t *A, double l )

{
	int i, j;
	dmatrix_t *B;

	if ( ( *A ).l <= 0 || ( *A ).c <= 0 ) {
		error ( "MATRIX.H: erroneous matrix size" );
	}

	B = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );

	if ( ( *A ).c == 1 ) {
		dmat_alloc ( B, ( *A ).l + 1, 1 );
		for ( i = 1; i < ( *B ).l; i++ ) {
			( *B ).m[i][1] = ( *A ).m[i][1];
		}
		( *B ).m[( *B ).l][1] = l;
	} else if ( ( *A ).l == 1 ) {
		dmat_alloc ( B, 1, ( *A ).c + 1 );
		for ( i = 1; i < ( *B ).c; i++ ) {
			( *B ).m[1][i] = ( *A ).m[1][i];
		}
		( *B ).m[1][( *B ).c] = l;
	} else {
		dmat_alloc ( B, ( *A ).l + 1, ( *A ).c + 1 );
		B = dmat_init ( B, 0.0 );
		for ( i = 1; i < ( *B ).l; i++ ) {
			for ( j = 1; j < ( *B ).c; j++ ) {
				( *B ).m[i][j] = ( *A ).m[i][j];
			}
		}
		( *B ).m[( *B ).l][( *B ).c] = l;
	}
	return B;
}

double dmat_norm(dmatrix_t *A)

{
	double s;
	int i;

	if ((*A).l != 1 && (*A).c != 1) {
		s = 0;
		error("MATRIX.H: incompatible matrix sizes");
	}
	else {
		if ((*A).l == 1) {
			for (s = 0.0, i = 1; i <= (*A).c; i++) {
				s += pow((*A).m[1][i], 2.0);
			}
		}
		else {
			for (s = 0.0, i = 1; i <= (*A).l; i++) {
				s += pow((*A).m[i][1], 2.0);
			}
		}
	}
	return sqrt(s);
}

dmatrix_t *dmat_normalize(dmatrix_t *A)

{
	return dmat_scalar_mult(A, 1.0 / dmat_norm(A));
}

dmatrix_t *from_homogeneous ( dmatrix_t *A )

{
	int i, j;
	dmatrix_t *B;

	if ( ( *A ).l < 1 || ( *A ).c < 1 ) {
		error ( "MATRIX.H: erroneous matrix size" );
	}

	B = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );

	if ( ( *A ).c == 1 ) {
		dmat_alloc ( B, ( *A ).l - 1, 1 );
		for ( i = 1; i < ( *A ).l; i++ ) {
			( *B ).m[i][1] = ( *A ).m[i][1];
		}
	} else if ( ( *A ).l == 1 ) {
		dmat_alloc ( B, 1, ( *A ).c - 1 );
		for ( i = 1; i < ( *A ).c; i++ ) {
			( *B ).m[1][i] = ( *A ).m[1][i];
		}
	} else {
		dmat_alloc ( B, ( *A ).l - 1, ( *A ).c - 1 );
		for ( i = 1; i < ( *A ).l; i++ ) {
			for ( j = 1; j < ( *A ).c; j++ ) {
				( *B ).m[i][j] = ( *A ).m[i][j];
			}
		}
	}
	return B;
}