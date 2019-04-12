#pragma once

/*
 * matrix.h, ver 1.0 & 1.5 & 2.0
 *
 * Mainly add vector's operation on it.
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
#include <iostream>
#include <vector>

using namespace std;

const double kInfinity = std::numeric_limits<double>::max();

struct mat {
public:
	double **m;
	int l, c;
};

void mat_alloc ( mat * A, int l, int c);

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

	vec times(const vec& b)
	{
		return vec(m[1][1]*b.m[1][1], m[2][1]*b.m[2][1], m[3][1]*b.m[3][1]);
	}


	inline mat toMat (int dim)
	{
		mat *B = (mat *)malloc(sizeof(mat));
		if(dim==3)
		{
			mat_alloc(B, 3, 1);
			B->m[1][1] = m[1][1];
			B->m[2][1] = m[2][1];
			B->m[3][1] = m[3][1];
		}
		if(dim==4)
		{
			mat_alloc(B, 4, 1);
			B->m[1][1] = m[1][1];
			B->m[2][1] = m[2][1];
			B->m[3][1] = m[3][1];
			B->m[4][1] = 1;
		}
		return *B;
	}

	void print ( )
	{
		cout << "x:"<< m[1][1] << " y:" << m[2][1] << " z:" << m[3][1];
		if(l==4) cout<< " buffer" << m[4][1] <<endl;
		else if (l ==3) cout <<endl;
	}

	//vector as an array
	inline const double& operator[](int i) const;

	//vec's operator

//vec + vec, vec - vec
	inline vec operator+(const vec &b )
	{
		return vec { m[1][1] + b.m[1][1], m[2][1] + b.m[2][1],
		             m[3][1] + b.m[3][1] };
	}

	friend vec operator+(const vec &a, const vec &b )
	{
		return vec { a.m[1][1] + b.m[1][1],  a.m[2][1] + b.m[2][1],
		             a.m[3][1] + b.m[3][1] };
	}
	inline vec& operator+=(const vec &b )
	{
		m[1][1] += b.m[1][1];
		m[2][1] += b.m[2][1];
		m[3][1] += b.m[3][1];
		return *this;
	}
	inline vec operator-( vec &b )
	{
		return vec { m[1][1] - b.m[1][1], m[2][1] - b.m[2][1],
		             m[3][1] - b.m[3][1] };
	}
	friend vec operator-( const vec &a, const vec &b )
	{
		return vec { a.m[1][1] - b.m[1][1], a.m[2][1] - b.m[2][1],
		             a.m[3][1] - b.m[3][1] };
	}
//vec * num, num* vec, vec*vec (dot product)
//	inline vec operator*(double s )
//	{
//		return vec { m[1][1] * s, m[2][1] * s, m[3][1] * s };
//	}
	friend vec operator*(const vec a, double s )
	{
		return vec { a.m[1][1] * s, a.m[2][1] * s, a.m[3][1] * s };
	}
	//notice here
	friend vec operator*( double s, const vec &v )
	{
		return v * s;
	}
	friend vec operator*( double s, vec &v )
	{
		return v * s;
	}

	friend double operator*( const vec &v, const vec &w  )
	{
		return v.dot ( w );
	}

//vec / num
	inline vec operator/(double s )
	{
		return (*this) * ( 1 / s );
	}
	inline vec operator/(const vec &b )
	{
		return vec { m[1][1] / b.m[1][1], m[2][1] / b.m[2][1],
		             m[3][1] / b.m[3][1] };
	}
	friend vec operator/( const vec &a, const vec &b )
	{
		return vec { a.m[1][1] / b.m[1][1], a.m[2][1] / b.m[2][1],
		             a.m[3][1] / b.m[3][1] };
	}
	friend vec operator/( double s, const vec &v )
	{
		return vec { s / v.m[1][1], s / v.m[2][1], s / v.m[3][1] };
	}
// vec = -vec
	inline vec operator-()
	{
		return vec { -m[1][1], -m[2][1], -m[3][1] };
	}

	friend vec operator-(const vec &a)
	{
		return vec { -a.m[1][1], -a.m[2][1], -a.m[3][1] };
	}

//equal
	bool operator==( const vec &b){
		return m[1][1] == b.m[1][1] && m[2][1] == b.m[2][1] && m[3][1] == b.m[3][1] && m[4][1] == b.m[4][1] ;
	}
//vec's info
	inline ostream &operator<<( ostream &os )
	{
		os << "vec:\n\tx = " << m[1][1] << "\n\ty = " << m[2][1]
		   << "\n\tz = " << m[3][1];
		return os;
	}
};

const double &vec::operator[] ( int i ) const
{
	switch (i){
		case 1:
			return this->m[1][1];
			break;
		case 2:
			return this->m[2][1];
			break;
		case 3:
			return this->m[3][1];
			break;
		case 4:
			return this->m[4][1];
			break;
		default:
			return this->m[1][1];
			break;
	}
}

vec toVec(mat mat, int dimension)
{
	vec temp = vec();
	if(dimension == 3)
	{
		temp = vec(mat.m[1][1],mat.m[2][1],mat.m[3][1]);
		return temp;
	}
	else
	{
		temp = vec(mat.m[1][1], mat.m[2][1], mat.m[3][1], 1);
		return temp;
	}
}

class vec_2:public vec
{
public:
	vector<vector<double>> m;
	int l, c;

	inline vec_2 ( ):l(2), c(1){
		m.resize(3, vector<double>(2, 0));
	}
	inline vec_2 ( double x, double y):l(2), c(1)
	{
		m.resize(3, vector<double>(2, 0));
		m[1][1] = x;
		m[2][1] = y;
	}
	inline vec_2 ( const vec &v )
	{
		m.resize(3, vector<double>(2, 0));
		m[1][1] = v.m[1][1];
		m[2][1] = v.m[2][1];
		l = v.l;
		c = v.c;
	}
	inline vec_2 ( const vec_2 &v )
	{
		m = v.m;
		l = v.l;
		c = v.c;
	}

	inline vec toVec3()
	{
		vec a = vec(m[1][1], m[2][1], 0);
		return a;
	}

	inline vec toVec4()
	{
		vec a = vec(m[1][1], m[2][1], 0, 1);
		return a;
	}

	inline ~vec_2 ()
	{
		//thanks vector again!
	}

	//vec + vec, vec - vec
	inline vec_2 operator+(const vec &b )
	{
		return vec_2 { m[1][1] + b.m[1][1], m[2][1] + b.m[2][1] };
	}
	inline vec_2& operator+=(const vec &b )
	{
		m[1][1] += b.m[1][1];
		m[2][1] += b.m[2][1];
		return *this;
	}
	inline vec_2 operator-( vec &b )
	{
		return vec_2 { m[1][1] - b.m[1][1], m[2][1] - b.m[2][1] };
	}
	friend vec_2 operator-( const vec_2 &a, const vec_2 &b )
	{
		return vec_2 { a.m[1][1] - b.m[1][1], a.m[2][1] - b.m[2][1]};
	}
//vec * num, num* vec, vec*vec (dot product)
	inline vec_2 operator*(double s )
	{
		return vec_2 { m[1][1] * s, m[2][1] * s };
	}
	//notice here
	friend vec_2 operator*( double s, vec_2 &v )
	{
		return v * s;
	}
	friend double operator*( const vec_2 &v, const vec_2 &w  )
	{
		return v.m[1][1]*w.m[1][1] + v.m[2][1]*w.m[2][1] ;
	}

//vec / num
	inline vec_2 operator/(double s )
	{
		return (*this) * ( 1 / s );
	}
	inline vec_2 operator/(const vec_2 &b )
	{
		return vec_2 { m[1][1] / b.m[1][1], m[2][1] / b.m[2][1]};
	}
	friend vec_2 operator/( const vec_2 &a, const vec_2 &b )
	{
		return vec_2 { a.m[1][1] / b.m[1][1], a.m[2][1] / b.m[2][1] };
	}
	friend vec_2 operator/( double s, const vec_2 &v )
	{
		return vec_2 { s / v.m[1][1], s / v.m[2][1] };
	}
// vec = -vec
	inline vec_2 operator-()
	{
		return vec_2 { -m[1][1], -m[2][1] };
	}

//equal
	bool operator==( const vec_2 &b){
		return m[1][1] == b.m[1][1] && m[2][1] == b.m[2][1] ;
	}
//vec's info
	inline ostream &operator<<( ostream &os )
	{
		os << "vec_2:\n\tx = " << m[1][1] << "\n\ty = " << m[2][1] << endl;
		return os;
	}
};

class Vec3{
public:
	Vec3() : x(0), y(0), z(0) {}
	Vec3(double xx) : x(xx), y(xx), z(xx) {}
	Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
	Vec3 operator * (const double &r) const { return Vec3(x * r, y * r, z * r); }
	Vec3 operator * (const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
	Vec3 operator - (const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
	Vec3 operator + (const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
	Vec3 operator - () const { return Vec3(-x, -y, -z); }
	Vec3& operator += (const Vec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
	friend Vec3 operator * (const double &r, const Vec3 &v)
	{ return Vec3(v.x * r, v.y * r, v.z * r); }
	friend std::ostream & operator << (std::ostream &os, const Vec3 &v)
	{ return os << v.x << ", " << v.y << ", " << v.z; }
	double x, y, z;
};

class Vec2
{
public:
	Vec2() : x(0), y(0) {}
	Vec2(double xx) : x(xx), y(xx) {}
	Vec2(double xx, double yy) : x(xx), y(yy) {}
	Vec2 operator * (const double &r) const { return Vec2(x * r, y * r); }
	Vec2 operator + (const Vec2 &v) const { return Vec2(x + v.x, y + v.y); }
	double x, y;
};

Vec3 normalize(const Vec3 &v)
{
	double mag2 = v.x * v.x + v.y * v.y + v.z * v.z;
	if (mag2 > 0) {
		double invMag = 1 / sqrtf(mag2);
		return Vec3(v.x * invMag, v.y * invMag, v.z * invMag);
	}

	return v;
}

inline
double dotProduct(const Vec3 &a, const Vec3 &b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

Vec3 crossProduct(const Vec3 &a, const Vec3 &b)
{
	return Vec3(
			a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x
	);
}

inline
double clamp(const double &lo, const double &hi, const double &v)
{ return std::max(lo, std::min(hi, v)); }

inline
double deg2rad(const double &deg)
{ return deg * M_PI / 180; }

inline
Vec3 mix(const Vec3 &a, const Vec3& b, const double &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

Vec3 reflect(const Vec3 &I, const Vec3 &N)
{
	return I - 2 * dotProduct(I, N) * N;
}

Vec3 refract(const Vec3 &I, const Vec3 &N, const double &ior)
{
	double cosi = clamp(-1, 1, dotProduct(I, N));
	double etai = 1, etat = ior;
	Vec3 n = N;
	if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
	double eta = etai / etat;
	double k = 1 - eta * eta * (1 - cosi * cosi);
	return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}

void fresnel(const Vec3 &I, const Vec3 &N, const double &ior, double &kr)
{
	float cosi = clamp(-1, 1, dotProduct(I, N));
	float etai = 1, etat = ior;
	if (cosi > 0) {  std::swap(etai, etat); }
	// Compute sini using Snell's law
	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	// Total internal reflection
	if (sint >= 1) {
		kr = 1;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		kr = (Rs * Rs + Rp * Rp) / 2;
	}
}


void error(const char error_text[])
{
	fprintf(stderr, "Run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	exit(1);
}


void write_dmatrix(mat *M)

{
	int i, j;

	for (i = 1; i <= (*M).l; i++) {
		for (j = 1; j <= (*M).c; j++) {
			printf("%7.4f ", (*M).m[i][j]);
		}
		printf("\n");
	}
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)

{
	int i;
	double **m;

	m = (double **)malloc((unsigned)(nrh - nrl + 1) * sizeof(double));
	if (!m) {
		error("MATRIX.H: allocation failure");
	}
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (double *)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
		if (!m[i]) {
			error("MATRIX.H: allocation failure");
		}
		m[i] -= ncl;
	}
	return m;
}


void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)

{
	int i;

	for (i = nrh; i >= nrl; i--) {
		free((char *)(m[i] + ncl));
	}
	free((char *)(m + nrl));
}


void delete_dmatrix(mat *A)

{
	free_dmatrix(A->m, 1, A->l, 1, A->c);
	free(A);
}


void mat_alloc ( mat * A, int l, int c)

{
	(*A).m = dmatrix(1, l, 1, c);
	(*A).l = l;
	(*A).c = c;
}


mat *dmat_duplicate(mat *A)

{
	int i, j;
	mat *B;
	B = (mat *)malloc(sizeof(mat));
	mat_alloc(B, ( *A ).l, ( *A ).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = (*A).m[i][j];
		}
	}
	return B;
}


mat *dmat_init(mat *A, double a)

{
	int i, j;

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*A).m[i][j] = a;
		}
	}
	return A;
}


mat *dmat_identity(mat *A)

{
	int i, j;

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			if (i == j) {
				(*A).m[i][j] = 1.0;
			}
			else {
				(*A).m[i][j] = 0.0;
			}
		}
	}
	return A;
}


mat *dmat_scalar_mult(mat *A, double a)

{
	mat *B;
	int i, j;

	B = (mat *)malloc(sizeof(mat));
	mat_alloc(B, ( *A ).l, ( *A ).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = (*A).m[i][j] * a;
		}
	}
	return B;
}


mat *dmat_mult(mat *A, mat *B)

{
	mat *C;
	double s;
	int i, j, k;

	if ((*A).c != (*B).l) {
		error("MATRIX.H: incompatible matrix sizes");
	}
	C = (mat *)malloc(sizeof(mat));
	mat_alloc(C, ( *A ).l, ( *B ).c);

	for (i = 1; i <= (*C).l; i++) {
		for (j = 1; j <= (*C).c; j++) {
			for (s = 0.0, k = 1; k <= (*A).c; k++) {
				s += (*A).m[i][k] * (*B).m[k][j];
			}
			(*C).m[i][j] = s;
		}
	}
	return C;
}


mat *dmat_add(mat *A, mat *B)
{
	mat *C;
	int i, j;

	if ((*A).l != (*B).l || (*A).c != (*B).c) {
		error("MATRIX.H: incompatible matrix sizes");
	}
	C = (mat *)malloc(sizeof(mat));
	mat_alloc(C, ( *A ).l, ( *A ).c);

	for (i = 1; i <= (*C).l; i++) {
		for (j = 1; j <= (*C).c; j++) {
			(*C).m[i][j] = (*A).m[i][j] + (*B).m[i][j];
		}
	}
	return C;
}


mat *dmat_sub(mat *A, mat *B)
{
	mat *C;
	int i, j;

	if ((*A).l != (*B).l || (*A).c != (*B).c) {
		error("MATRIX.H: incompatible matrix sizes");
	}
	C = (mat *)malloc(sizeof(mat));
	mat_alloc(C, ( *A ).l, ( *A ).c);

	for (i = 1; i <= (*C).l; i++) {
		for (j = 1; j <= (*C).c; j++) {
			(*C).m[i][j] = (*A).m[i][j] - (*B).m[i][j];
		}
	}
	return C;
}


double dmat_norm(mat *A)

{
	double s;
	int i;

	if ((*A).l != 1 && (*A).c != 1) {
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


mat *dmat_normalize(mat *A)

{
	return dmat_scalar_mult(A, 1.0 / dmat_norm(A));
}


mat *dmat_transpose(mat *A)

{
	mat *B;
	int i, j;

	B = (mat *)malloc(sizeof(mat));
	mat_alloc(B, ( *A ).c, ( *A ).l);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[j][i] = (*A).m[i][j];
		}
	}
	return B;
}


double ddot_product(mat *A, mat *B)

{
	mat *C;

	C = (mat *)malloc(sizeof(mat));

	if ((*A).c == (*B).c && (*A).l == 1 && (*B).l == 1) {
		C = dmat_mult(A, dmat_transpose(B));
	}
	else if ((*A).c == (*B).l && (*A).l == 1 && (*B).c == 1) {
		C = dmat_mult(A, B);
	}
	else if ((*A).l == (*B).c && (*A).c == 1 && (*B).l == 1) {
		C = dmat_mult(B, A);
	}
	else if ((*A).l == (*B).l && (*A).c == 1 && (*B).c == 1) {
		C = dmat_mult(dmat_transpose(A), B);
	}
	else error("MATRIX.H: Incompatible matrix sizes");
	return (*C).m[1][1];
}


mat *cross_product_matrix(mat *A, mat *B)

{
	int i, j;
	mat *C;

	C = (mat *)malloc(sizeof(mat));
	mat_alloc(C, ( *A ).l, ( *A ).l);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).l; j++) {
			if (i == 1) {
				(*C).m[i][j] = 1.0;
			}
			else if (i == 2) {
				(*C).m[i][j] = (*A).m[j][1];
			}
			else {
				(*C).m[i][j] = (*B).m[j][1];
			}
		}
	}
	return C;
}


mat *sub_matrix(mat *A, int r, int c)

{
	int i, j, k, l;
	mat *B;

	if (r < 1 || r >(*A).l || c < 1 || c >(*A).c || (*A).c < 2 || (*A).l < 2) {
		error("MATRIX.H: erroneous indices");
	}
	B = (mat *)malloc(sizeof(mat));
	mat_alloc(B, ( *A ).l - 1, ( *A ).c - 1);

	for (i = 1, k = 1; i <= (*A).l; i++) {
		for (j = 1, l = 1; j <= (*A).c; j++) {
			if (j != c && i != r) {
				(*B).m[k][l] = (*A).m[i][j];
			}
			if (j != c) l++;
		}
		if (i != r) k++;
	}
	return B;
}


double determinant(mat *A)

{
	int i;
	double det;

	if ((*A).l < 1 || (*A).c < 1) {
		error("MATRIX.H: erroneous matrix size");
	}
	else if ((*A).l != (*A).c) {
		error("MATRIX.H: not a square matrix");
	}
	else if ((*A).l == 1) {
		det = (*A).m[1][1];
	}
	else {
		det = 0.0;
		for (i = 1; i <= (*A).c; i++) {
			det += pow(-1.0, i + 1.0)*(*A).m[1][i] * determinant(sub_matrix(A, 1, i));
		}
	}
	return det;
}


mat *cofactor(mat *A)

{
	int i, j;
	mat *B;

	B = (mat *)malloc(sizeof(mat));
	mat_alloc(B, ( *A ).l, ( *A ).c);

	for (i = 1; i <= (*A).l; i++) {
		for (j = 1; j <= (*A).c; j++) {
			(*B).m[i][j] = pow(-1.0, i + j)*determinant(sub_matrix(A, i, j));
		}
	}
	return B;
}


mat *dmat_inverse(mat *A)

{
	return dmat_scalar_mult(dmat_transpose(cofactor(A)), 1.0 / determinant(A));
}


mat *dcross_product(mat *A, mat *B)

{
	int i;
	mat *C;

	if ((*A).l != 3 || (*A).c != 1 || (*B).l != 3 || (*B).c != 1) {
		error("MATRIX.H: Incompatible matrix sizes");
	}

	C = (mat *)malloc(sizeof(mat));
	mat_alloc(C, ( *A ).l, ( *A ).c);

	for (i = 1; i <= (*A).l; i++) {
		(*C).m[i][1] = pow(-1.0, i + 1)*determinant(sub_matrix(cross_product_matrix(A, B), 1, i));
	}
	return C;
}


mat *to_homogeneous(mat *A, double l)

{
	int i, j;
	mat *B;

	if ((*A).l <= 0 || (*A).c <= 0) {
		error("MATRIX.H: erroneous matrix size");
	}

	B = (mat *)malloc(sizeof(mat));

	if ((*A).c == 1) {
		mat_alloc(B, ( *A ).l + 1, 1);
		for (i = 1; i < (*B).l; i++) {
			(*B).m[i][1] = (*A).m[i][1];
		}
		(*B).m[(*B).l][1] = l;
	}
	else if ((*A).l == 1) {
		mat_alloc(B, 1, ( *A ).c + 1);
		for (i = 1; i < (*B).c; i++) {
			(*B).m[1][i] = (*A).m[1][i];
		}
		(*B).m[1][(*B).c] = l;
	}
	else {
		mat_alloc(B, ( *A ).l + 1, ( *A ).c + 1);
		B = dmat_init(B, 0.0);
		for (i = 1; i < (*B).l; i++) {
			for (j = 1; j < (*B).c; j++) {
				(*B).m[i][j] = (*A).m[i][j];
			}
		}
		(*B).m[(*B).l][(*B).c] = l;
	}
	return B;
}


mat *from_homogeneous(mat *A)

{
	int i, j;
	mat *B;

	if ((*A).l < 1 || (*A).c < 1) {
		error("MATRIX.H: erroneous matrix size");
	}

	B = (mat *)malloc(sizeof(mat));

	if ((*A).c == 1) {
		mat_alloc(B, ( *A ).l - 1, 1);
		for (i = 1; i < (*A).l; i++) {
			(*B).m[i][1] = (*A).m[i][1];
		}
	}
	else if ((*A).l == 1) {
		mat_alloc(B, 1, ( *A ).c - 1);
		for (i = 1; i < (*A).c; i++) {
			(*B).m[1][i] = (*A).m[1][i];
		}
	}
	else {
		mat_alloc(B, ( *A ).l - 1, ( *A ).c - 1);
		for (i = 1; i < (*A).l; i++) {
			for (j = 1; j < (*A).c; j++) {
				(*B).m[i][j] = (*A).m[i][j];
			}
		}
	}
	return B;
}


class Mat
{
public:
	vector<vector<double>> m;
	int l, c;

	Mat():l(0), c(0)
	{
		m.resize(5, vector<double>(5, 0));
		l = 4;
		c = 4;
	}

	Mat ( int sizex, int sizey, double value):l(sizex), c(sizey)
	{
		m.resize(sizex+1, vector<double>(sizey+1, value));
		l = sizex;
		c = sizey;
	}

	/*Mat ( vec* A ) :l((*A).l), c( ( *A ).c), m( ( *A ).m){ }*/

	inline void set ( Mat A )
	{
		this->m = A.m;
		this->l = A.l;
		this->c = A.c;
	}

	inline void set ( int x, int y, double value )
	{
		this->m[x][y] = value;
	}

	inline void copy ( Mat A, int destination_start_line, int destination_start_col, int source_start_line, int source_start_col )
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

	~Mat ( )
	{
		//thanks vector! No need to free
	}

	//I used wanna to design two way: B.identify(A), but finally it doesn't be better than this one:
	inline Mat identify ( )
	{
		Mat temp;
		temp.m = this->m;
		temp.c = this->c;
		temp.l = this->l;
		for ( int i = 0; i <= l; i++ )
			for ( int j = 0; j <= c; j++ )
				if(i == j)temp.m[i][j] = 1;
				else temp.m[i][j] = 0;
		return temp;
	}

	inline Mat translated ( double x, double y, double z )
	{
		Mat temp = Mat ( 5, 5, 0 );
		temp = temp.identify ( );
		temp.set ( 1, 4, x );
		temp.set ( 2, 4, y );
		temp.set ( 3, 4, z );
		return temp;
	}

	inline Mat scalared ( int scalar )
	{
		return scalared ( scalar, scalar, scalar );
	}

	inline Mat scalared ( double x, double y, double z )
	{
		Mat temp = Mat ( 5, 5, 0 );
		temp = temp.identify ( );
		temp.set ( 1, 1, x );
		temp.set ( 2, 2, y );
		temp.set ( 3, 3, z );
		return temp;
	}

	inline Mat rotated ( double x, double y, double z, double angle, int direction)
	{
		Mat I = Mat ( 5, 5, 0 );
		I = I.identify ( );
		Mat J = Mat ( 5, 5, 0 );
		//vec V = vec ( x, y, z );
		//V.normalized ( );

		return I;

		//TODO:
	}

	inline Mat transpose ( )
	{
		Mat B = Mat ( c, l, 0 );
		for (int i = 1; i <= l; i++ )
			for (int j = 1; j <= c; j++ )
				B.m[j][i] = m[i][j];
		return B;
	}

	inline Mat part (int line, int col )
	{
		Mat B = Mat ( c, l, 0 );

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

	inline Mat inverse ( )
	{
		Mat temp;
		temp.m = this->m;
		temp.c = c;
		temp.l = l;

		Mat cofactor = Mat ( l, c, 0 );

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
	//Mat *rotate ( double Vx, double Vy, double Vz, double angle )
	//{
	//	Mat *I, *J;
	//	vec V = vec ( Vx, Vy, Vz );

	//	I = (Mat *) malloc ( sizeof ( Mat ) );
	//	J = (Mat *) malloc ( sizeof ( Mat ) );
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

	//	Mat* ret = to_homogeneous ( dmat_add ( I, dmat_add ( dmat_scalar_mult ( J, sin ( angle ) ), dmat_scalar_mult ( dmat_mult ( J, J ), 1.0 - cos ( angle ) ) ) ), 1.0 );
	//	return ret;
	//}

};

//mat's operator

//mat + num & num + mat
inline Mat operator+( const Mat &v, double index )
//this three + are all mat + vec: add vec as a column in the right.
{
	Mat temp;
	temp.m = v.m;
	temp.c = v.c;
	temp.l = v.l;
	for ( int i = 0; i <= v.l; i++ )
		for ( int j = 0; j <= v.c; j++ )
			temp.m[i][j] += index;
	return temp;
}
inline Mat operator+( double index, const Mat &v )//this + is mat + double: add a number to all the entries.
{
	return v + index;
}
//mat + vec and vec + mac
inline Mat operator+( const Mat &a, const vec &v )
{
	int newLine = ( a.l > v.l ? a.l : v.l );
	Mat temp = Mat ( newLine, a.c + 1, 0 );
	temp.l = newLine;
	temp.c = a.c + 1;

	temp.copy ( a, 1, 1, 1, 1 );

	for ( int j = 1; j <= v.l; j++ )
		temp.m[j][temp.c] = v.m[j][1];
	return temp;
}
inline Mat operator+(const vec &v, const Mat &a)
{
	int newLine = ( a.l > v.l ? a.l : v.l );
	Mat temp = Mat ( newLine, a.c + 1, 0 );
	temp.l = newLine;
	temp.c = a.c + 1;

	temp.copy ( a, 1,2,1,1);

	for ( int j = 1; j <= v.l; j++ )
		temp.m[j][1] = v.m[j][1];
	return temp;
}
inline void operator+=(Mat &a, const vec &v )
{
	a = a + v;
}
inline void operator+=(const vec &v, Mat &a)
{
	a = a + v;
}
//mat + mat
inline Mat operator+( const Mat &A, const Mat &B )
{
	int i, j;
	if ( A.l != B.l || A.c != B.c ) {
		error ( "MATRIX.H: matA + matB: matrix sizes are unequal;" );
	}
	Mat C = Mat ( A.l, A.c, 0 );

	for ( i = 1; i <= C.l; i++ )
		for ( j = 1; j <= C.c; j++ )
			C.m[i][j] = A.m[i][j] + B.m[i][j];
	return C;
}
inline void operator+=( Mat &A, const Mat &B )
{
	A = A + B;
}
//mat = -mat
inline Mat operator-( const Mat &v )
{
	Mat temp;
	temp.m = v.m;
	temp.c = v.c;
	temp.l = v.l;
	for ( int i = 0; i <= v.l; i++ )
		for ( int j = 0; j <= v.c; j++ )
			temp.m[i][j] *= -1;
	return temp;
}
//mat - mat
inline Mat operator-( const Mat &A, const Mat &B )
{
	return A + ( -B );
}
inline void operator-=(Mat &A, const Mat &B )
{
	A = A + ( -B );
}
// mat = a*mat, mat = mat/a
inline Mat operator*( const Mat &v, double index )
{
	Mat temp;
	temp.m = v.m;
	temp.c = v.c;
	temp.l = v.l;
	for ( int i = 0; i <= v.l; i++ )
		for ( int j = 0; j <= v.c; j++ )
			temp.m[i][j] *= index;
	return temp;
}
inline Mat operator*( double s, const Mat &v )
{
	return v * s;
}
inline Mat operator/( const Mat &v, double s )
{
	return v * ( 1 / s );
}
//equal
inline bool operator==(const Mat &a, const Mat &b){
	if(a.c!=b.c || a.l!=b.l)return false;
	for(int i = 1; i <= a.c; i++)
		for(int j = 1; j <= a.l; j++)
			if (a.m[j][i] != b.m[j][i])
				return false;
	return true;
}
//info of a mat
inline std::ostream &operator<<( std::ostream &os, const Mat &v )
{
	;
}

class Camera : public Vec3
{
public:
	double x, y, z;
	Camera(double xx, double yy, double zz):x(xx),y(yy),z(zz){ }
	friend Camera operator*( const Mat &v, const Camera &c )
	{
		double xx = c.x*v.m[1][1]+ c.x*v.m[1][2]+ c.x*v.m[1][3];
		double yy = c.y*v.m[2][1]+ c.y*v.m[2][2]+ c.y*v.m[2][3];
		double zz = c.z*v.m[3][1]+ c.z*v.m[3][2]+ c.z*v.m[3][3];
		return {xx,yy,zz};
	}

	Vec3 toVec()
	{
		return {x,y,z};
	}
};
