#pragma once

/*
 * matrix.h, ver 1.0 & 1.5
 *
 * author: Mu.H. & Dr.S.B.
 * date: 19-3-26
 */

#include <cassert>
#include <cmath>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define quote(x) #x
//#define dmatrix_t mat

using namespace std;

//-----------------------------------------------------------------------
// matrix.h - ver 1.5. copyright@ Mu He
// 19-3-26 ~ 19-3-28
//------------------------------------------------------------------------
void error ( const char error_text[] )
{
	fprintf ( stderr, "Run-time error...\n" );
	fprintf ( stderr, "%s\n", error_text );
	exit ( 1 );
}

//classical matrix.
class dmatrix_t
{
public:
	vector<vector<double>> m;
	int l, c;
	
	dmatrix_t():l(0), c(0)
	{
		m.resize(5, vector<double>(5, 0));
		l = 4;
		c = 4;
	}

	dmatrix_t ( int sizex, int sizey, double value):l(sizex), c(sizey)
	{
		m.resize(sizex+1, vector<double>(sizey+1, value));
		l = sizex;
		c = sizey;
	}

	/*dmatrix_t ( vec* A ) :l((*A).l), c( ( *A ).c), m( ( *A ).m){ }*/

	inline void set ( dmatrix_t A )
	{
		this->m = A.m;
		this->l = A.l;
		this->c = A.c;
	}

	inline void set ( int x, int y, double value )
	{
		this->m[x][y] = value;
	}

	inline void copy ( dmatrix_t A, int destination_start_line, int destination_start_col, int source_start_line, int source_start_col )
	{
		int d_l = destination_start_line;
		int d_c = destination_start_col;
		int s_l = source_start_line;
		int s_c = source_start_col;

		int min_l = ( (this->l - d_l) < (A.l - s_l) ? ( this->l - d_l ) : ( A.l - s_l ) );
		int min_c = ( ( this->c - d_c ) < ( A.c - s_c ) ? ( this->c - d_c ) : ( A.c - s_c ) );
		for ( int cc = 0; cc <= min_c; cc++ )
			for ( int ll = 0; ll <= min_l; ll++ )
				this->m[ll+d_l][cc+d_c] = A.m[ll+s_l][cc+s_c];
	}

	inline void print ( )
	{
		int i, j;
		for ( i = 1; i <= l; i++ ) {
			for ( j = 1; j <= c; j++ ) {
				printf ( "%7.4f ", this->m[i][j] );
			}
			printf ( "\n" );
		}
		cout << "//--------------------------------" << endl;
	}

	~dmatrix_t ( )
	{
		//thanks vector! No need to free
	}

	//I used wanna to design two way: B.identify(A), but finally it doesn't be better than this one:
	inline dmatrix_t identify ( )
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

	inline dmatrix_t translated ( double x, double y, double z )
	{
		dmatrix_t temp = dmatrix_t ( 5, 5, 0 );
		temp = temp.identify ( );
		temp.set ( 1, 4, x );
		temp.set ( 2, 4, y );
		temp.set ( 3, 4, z );
		return temp;
	}

	inline dmatrix_t scalared ( int scalar )
	{
		return scalared ( scalar, scalar, scalar );
	}

	inline dmatrix_t scalared ( double x, double y, double z )
	{
		dmatrix_t temp = dmatrix_t ( 5, 5, 0 );
		temp = temp.identify ( );
		temp.set ( 1, 1, x );
		temp.set ( 2, 2, y );
		temp.set ( 3, 3, z );
		return temp;
	}

	inline dmatrix_t rotated ( double x, double y, double z, double angle, int direction)
	{
		dmatrix_t I = dmatrix_t ( 5, 5, 0 );
		I = I.identify ( );
		dmatrix_t J = dmatrix_t ( 5, 5, 0 );
		//vec V = vec ( x, y, z );
		//V.normalized ( );

		return I;

		//TODO:
	}

	inline dmatrix_t transpose ( )
	{
		dmatrix_t B = dmatrix_t ( c, l, 0 );
		for (int i = 1; i <= l; i++ ) 
			for (int j = 1; j <= c; j++ ) 
				B.m[j][i] = m[i][j];
		return B;
	}

	inline dmatrix_t part (int line, int col )
	{
		dmatrix_t B = dmatrix_t ( c, l, 0 );

		int i, j, a, b;

		if (line < 1 || line > l || col < 1 || col > c || c < 2 || l < 2 ) error ( "MATRIX.H: mat.part: erroneous indices" );
	
		for ( i = 1, a = 1; i <= l; i++ ) {
			for ( j = 1, b = 1; j <= c; j++ ) {
				if ( j != col && i != line )B.m[a][b] = m[i][j];
				if ( j != col )b++;
			}
			if ( i != line )a++;
		}
		B.l -= 1;
		B.c -= 1;
		return B;
	}

	inline dmatrix_t inverse ( )
	{
		dmatrix_t temp;
		temp.m = this->m;
		temp.c = c;
		temp.l = l;
	
		dmatrix_t cofactor = dmatrix_t ( l, c, 0 );

		for ( int i = 1; i <= l; i++ )
			for ( int j = 1; j <= c; j++ )
				cofactor.m[i][j] = pow ( -1.0, i + j ) * temp.part ( i, j ).det ( ) * ( 1.0 / this->det ( ) );

		return cofactor.transpose ( );
	}

	inline double det ( )
	{
		double det = 0;

		if ( l < 1 || c < 1 || l != c) error ( "MATRIX.H: mat.det(): matrix size must be n*n" );
		else if (l == 1 )det = m[1][1];
		else
			for ( int i = 1; i <= c; i++ )
				det += pow ( -1.0, i + 1.0 ) * m[1][i] * this->part ( 1, i ).det();
		return det;
	}

	//Source:
	//// Allocates and creates a rotation matrix
	//dmatrix_t *rotate ( double Vx, double Vy, double Vz, double angle )
	//{
	//	dmatrix_t *I, *J;
	//	vec V = vec ( Vx, Vy, Vz );

	//	I = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	//	J = (dmatrix_t *) malloc ( sizeof ( dmatrix_t ) );
	//	dmat_alloc ( I, 3, 3 );
	//	dmat_alloc ( J, 3, 3 );

	//	I = dmat_identity ( I );
	//	J = dmat_init ( J, 0.0 );

	//	V.normalized ( );

	//	( *J ).m[2][3] = -V.m[1][1];
	//	( *J ).m[3][2] = V.m[1][1];

	//	( *J ).m[1][3] = V.m[2][1];
	//	( *J ).m[3][1] = -V.m[2][1];

	//	( *J ).m[1][2] = -V.m[3][1];
	//	( *J ).m[2][1] = V.m[3][1];

	//	dmatrix_t* ret = to_homogeneous ( dmat_add ( I, dmat_add ( dmat_scalar_mult ( J, sin ( angle ) ), dmat_scalar_mult ( dmat_mult ( J, J ), 1.0 - cos ( angle ) ) ) ), 1.0 );
	//	return ret;
	//}

};

// a 3*1 or 4*1 vector. In 4*1-case, the last one actually don't attend any calculations.
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

	inline void scalared ( int scalar )
	{
		for ( int i = 1; i <= l; i++ )
			for ( int j = 1; j <= c; j++ )
				this->m[i][j] = this->m[i][j] * scalar;
	}

	inline vec scalared ( double x, double y, double z )
	{
		vec temp = vec ( );
		temp.m[1][1] *= x;
		temp.m[2][1] *= y;
		temp.m[3][1] *= z;
		return temp;
	}

	inline vec scalared ( double index )
	{
		return scalared ( index, index, index );
	}

	inline dmatrix_t toMat ( )
	{
		dmatrix_t mat = dmatrix_t(this->l,this->c,0);
		mat.m = this->m;
		return mat;
	}

	void print ( )
	{
		cout << "x:"<< m[1][1] << " y:" << m[2][1] << " z:" << m[3][1];
		if(l==4) cout<< " buffer" << m[4][1] <<endl;
		else if (l ==3) cout <<endl;
	}
};

//vec's operator

//vec + vec, vec - vec
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
//vec * num, num* vec, vec*vec (dot product)
inline vec operator*( const vec &v, double s )
{
	return vec { v.m[1][1] * s, v.m[2][1] * s, v.m[3][1] * s };
}
inline vec operator*( double s, const vec &v )
{
	return v * s;
}
inline double operator*( const vec &u, const vec &v )
{
	return u.dot ( v );
}
//vec / num
inline vec operator/( const vec &v, double s )
{
	return v * ( 1 / s );
}
inline vec operator/( double s, const vec &v )
{
	return vec { s / v.m[1][1], s / v.m[2][1], s / v.m[3][1] };
}
inline vec operator/( const vec &a, const vec &b )
{
	return vec { a.m[1][1] / b.m[1][1], a.m[2][1] / b.m[2][1],
			   a.m[3][1] / b.m[3][1] };
}
// vec = -vec
inline vec operator-( const vec &v )
{
	return vec { -v.m[1][1], -v.m[2][1], -v.m[3][1] };
}
//vec's info
inline std::ostream &operator<<( std::ostream &os, const vec &v )
{
	os << "vec:\n\tx = " << v.m[1][1] << "\n\ty = " << v.m[2][1]
		<< "\n\tz = " << v.m[3][1];
	return os;
}

//mat's operator

//mat + num & num + mat
inline dmatrix_t operator+( const dmatrix_t &v, double index )
//this three + are all mat + vec: add vec as a column in the right.
{
	dmatrix_t temp;
	temp.m = v.m;
	temp.c = v.c;
	temp.l = v.l;
	for ( int i = 0; i <= v.l; i++ )
		for ( int j = 0; j <= v.c; j++ )
			temp.m[i][j] += index;
	return temp;
}
inline dmatrix_t operator+( double index, const dmatrix_t &v )//this + is mat + double: add a number to all the entries.
{
	return v + index;
}
//mat + vec and vec + mac
inline dmatrix_t operator+( const dmatrix_t &a, const vec &v )
{
	int newLine = ( a.l > v.l ? a.l : v.l );
	dmatrix_t temp = dmatrix_t ( newLine, a.c + 1, 0 );
	temp.l = newLine;
	temp.c = a.c + 1;

	temp.copy ( a, 1, 1, 1, 1 );

	for ( int j = 1; j <= v.l; j++ )
		temp.m[j][temp.c] = v.m[j][1];
	return temp;
}
inline dmatrix_t operator+(const vec &v, const dmatrix_t &a)
{
	int newLine = ( a.l > v.l ? a.l : v.l );
	dmatrix_t temp = dmatrix_t ( newLine, a.c + 1, 0 );
	temp.l = newLine;
	temp.c = a.c + 1;

	temp.copy ( a, 1,2,1,1);

	for ( int j = 1; j <= v.l; j++ )
		temp.m[j][1] = v.m[j][1];
	return temp;
}
inline void operator+=(dmatrix_t &a, const vec &v )
{
	a = a + v;
}
inline void operator+=(const vec &v, dmatrix_t &a)
{
	a = a + v;
}
//mat + mat
inline dmatrix_t operator+( const dmatrix_t &A, const dmatrix_t &B )
{
	int i, j;
	if ( A.l != B.l || A.c != B.c ) {
		error ( "MATRIX.H: matA + matB: matrix sizes are unequal;" );
	}
	dmatrix_t C = dmatrix_t ( A.l, A.c, 0 );

	for ( i = 1; i <= C.l; i++ )
		for ( j = 1; j <= C.c; j++ )
			C.m[i][j] = A.m[i][j] + B.m[i][j];
	return C;
}
inline void operator+=( dmatrix_t &A, const dmatrix_t &B )
{
	A = A + B;
}
//mat = -mat
inline dmatrix_t operator-( const dmatrix_t &v )
{
	dmatrix_t temp;
	temp.m = v.m;
	temp.c = v.c;
	temp.l = v.l;
	for ( int i = 0; i <= v.l; i++ )
		for ( int j = 0; j <= v.c; j++ )
			temp.m[i][j] *= -1;
	return temp;
}
//mat - mat
inline dmatrix_t operator-( const dmatrix_t &A, const dmatrix_t &B )
{
	return A + ( -B );
}
inline void operator-=(dmatrix_t &A, const dmatrix_t &B )
{
	A = A + ( -B );
}
// mat = a*mat, mat = mat/a
inline dmatrix_t operator*( const dmatrix_t &v, double index )
{
	dmatrix_t temp;
	temp.m = v.m;
	temp.c = v.c;
	temp.l = v.l;
	for ( int i = 0; i <= v.l; i++ )
		for ( int j = 0; j <= v.c; j++ )
			temp.m[i][j] *= index;
	return temp;
}
inline dmatrix_t operator*( double s, const dmatrix_t &v )
{
	return v * s;
}
inline dmatrix_t operator/( const dmatrix_t &v, double s )
{
	return v * ( 1 / s );
}
//info of a mat
inline std::ostream &operator<<( std::ostream &os, const dmatrix_t &v )
{
	;
}

// vec3/4 actually is same as 3. It won't have anyother opreation in this
// project.
typedef vec vec3;
typedef vec vec4;

//-----------------------------------------------------------------------
// matrix.h - ver 1.0. copyright@Dr. Steven S. Beauchemin
//------------------------------------------------------------------------
																			 
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

vec dmat_mult ( dmatrix_t *A, vec B )

{
	vec C;
	double s;
	int i, j, k;

	if ( ( *A ).c != ( B ).l ) {
		error ( "MATRIX.H: incompatible matrix sizes" );
	}
	C = vec ( ( *A ).l + 1, B.c + 1, 0 );

	for ( i = 1; i < C.l; i++ ) {
		for ( j = 1; j < C.c; j++ ) {
			for ( s = 0.0, k = 1; k <= ( *A ).c; k++ ) {
				s += ( *A ).m[i][k] * B.m[k][j];
			}
			C.m[i][j] = s;
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

dmatrix_t dmat_transpose ( dmatrix_t A )

{
	dmatrix_t B;
	int i, j;

	B = dmatrix_t ( A.c, A.l, 0 );

	for ( i = 1; i <= A.l; i++ ) {
		for ( j = 1; j <= A.c; j++ ) {
			B.m[j][i] = A.m[i][j];
		}
	}
	return B;
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

dmatrix_t sub_matrix ( dmatrix_t A, int r, int c )

{
	int i, j, k, l;
	dmatrix_t B = dmatrix_t ( A.l, A.c, 0 );

	if ( r < 1 || r > A.l || c < 1 || c > A.c || A.c < 2 || A.l < 2 ) {
		error ( "MATRIX.H: erroneous indices" );
	}

	for ( i = 1, k = 1; i <= A.l; i++ ) {
		for ( j = 1, l = 1; j <= A.c; j++ ) {
			if ( j != c && i != r ) {
				B.m[k][l] = A.m[i][j];
			}
			if ( j != c )
				l++;
		}
		if ( i != r )
			k++;
	}
	B.l -= 1;
	B.c -= 1;
	return B;
}

double determinant ( dmatrix_t A )

{
	int i;
	double det;

	if ( A.l < 1 || A.c < 1 ) {
		error ( "MATRIX.H: erroneous matrix size" );
	} else if ( A.l != A.c ) {
		error ( "MATRIX.H: not a square matrix" );
	} else if ( A.l == 1 ) {
		det = A.m[1][1];
	} else {
		det = 0.0;
		for ( i = 1; i <= A.c; i++ ) {
			det +=
				pow ( -1.0, i + 1.0 ) * A.m[1][i] * determinant ( sub_matrix ( A, 1, i ) );
		}
	}
	return det;
}

dmatrix_t cofactor ( dmatrix_t A )
{
	int i, j;
	dmatrix_t B = dmatrix_t (A.l, A.c , 0);

	for ( i = 1; i <= A.l; i++ ) 
		for ( j = 1; j <= A.c; j++ ) 
			B.m[i][j] = pow ( -1.0, i + j ) * A.part ( i, j ).det ( );
	return B;
}

dmatrix_t dmat_inverse ( dmatrix_t A )
{
	return cofactor ( A ).transpose( ) * (1.0 / A.det( ));
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

dmatrix_t from_homogeneous ( dmatrix_t *A )

{
	int i, j;
	dmatrix_t B;

	if ( ( *A ).l < 1 || ( *A ).c < 1 ) {
		error ( "MATRIX.H: erroneous matrix size" );
	}

	if ( ( *A ).c == 1 ) {
		B = dmatrix_t(( *A ).l, 1, 0);
		for ( i = 1; i < ( *A ).l; i++ ) {
			B.m[i][1] = ( *A ).m[i][1];
		}
	} else if ( ( *A ).l == 1 ) {
		B = dmatrix_t ( 1, ( *A ).c, 0 );
		for ( i = 1; i < ( *A ).c; i++ ) {
			B.m[1][i] = ( *A ).m[1][i];
		}
	} else {
		B = dmatrix_t ( ( *A ).l, ( *A ).c, 0 );
		for ( i = 1; i < ( *A ).l; i++ ) {
			for ( j = 1; j < ( *A ).c; j++ ) {
				B.m[i][j] = ( *A ).m[i][j];
			}
		}
	}
	return B;
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