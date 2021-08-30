#ifndef CMATH_H
#define CMATH_H
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include  <math.h>
#include <string.h>
typedef void fatalfunc_t(const char*); /* fatal callback function type */
class CMath
{
public:
	CMath();
    ~CMath();
public:
/**
@brief  allocate memory of matrix
@param
@����ֵ
*/
static double* mat(int n, int m);
/**
@brief  * allocate memory of integer matrix

@param
@����ֵ
*/
static int* imat(int n, int m);
/**
@brief  generate new zero matrix

@param
@����ֵ
*/
static double* zeros(int n, int m);
/**
@brief  generate new identity matrix

@param
@����ֵ
*/
static double* eye(int n);
/**
@brief  inner product of vectors

@param
@����ֵ
*/
static double dot(const double* a, const double* b, int n);
/**
@brief  euclid norm of vector

@param
@����ֵ
*/
static double norm(const double* a, int n);
/**
@brief  outer product of 3d vectors

@param
@����ֵ
*/
static void cross3(const double* a, const double* b, double* c);
/**
@brief  normalize 3d vector

@param
@����ֵ
*/
static int  normv3(const double* a, double* b);
/**
@brief copy matrix 

@param
@����ֵ
*/
static void matcpy(double* A, const double* B, int n, int m);
/**
@brief multiply matrix 

@param
@����ֵ
*/
static void matmul(const char* tr, int n, int k, int m, double alpha,
    const double* A, const double* B, double beta, double* C);
/**
@brief inverse of matrix 

@param
@����ֵ
*/
static int  matinv(double* A, int n);
/**
@brief inverse of matrix

@param
@����ֵ
*/
static int  solve(const char* tr, const double* A, const double* Y, int n,
    int m, double* X);
/**
@brief least square estimation 

@param
@����ֵ
*/
static int  lsq(const double* A, const double* y, int n, int m, double* x,
    double* Q);
/**
@brief kalman filter

@param 
@����ֵ
*/
static int  filter(double* x, double* P, const double* H, const double* v,
    const double* R, int n, int m);
/**
@brief combine forward and backward filters by fixed-interval smoother as follows

@param
@����ֵ
*/
static int  smoother(const double* xf, const double* Qf, const double* xb,
    const double* Qb, int n, double* xs, double* Qs);
/**
@brief /* LD factorization (Q=L'*diag(D)*L)

@param
@����ֵ
*/
static int LD(int n, const double* Q, double* L, double* D);
/**
@brief ���������ָ���ļ�

@param
@����ֵ
*/
static void matfprint(const double A[], int n, int m, int p, int q, FILE* fp);

/**
@brief ���������ָ���ļ�

@param
@����ֵ
*/
static void matprint(const double A[], int n, int m, int p, int q);
};



#endif // !CMATH_H


