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
@返回值
*/
static double* mat(int n, int m);
/**
@brief  * allocate memory of integer matrix

@param
@返回值
*/
static int* imat(int n, int m);
/**
@brief  generate new zero matrix

@param
@返回值
*/
static double* zeros(int n, int m);
/**
@brief  generate new identity matrix

@param
@返回值
*/
static double* eye(int n);
/**
@brief  inner product of vectors

@param
@返回值
*/
static double dot(const double* a, const double* b, int n);
/**
@brief  euclid norm of vector

@param
@返回值
*/
static double norm(const double* a, int n);
/**
@brief  outer product of 3d vectors

@param
@返回值
*/
static void cross3(const double* a, const double* b, double* c);
/**
@brief  normalize 3d vector

@param
@返回值
*/
static int  normv3(const double* a, double* b);
/**
@brief copy matrix 

@param
@返回值
*/
static void matcpy(double* A, const double* B, int n, int m);
/**
@brief multiply matrix 

@param
@返回值
*/
static void matmul(const char* tr, int n, int k, int m, double alpha,
    const double* A, const double* B, double beta, double* C);
/**
@brief inverse of matrix 

@param
@返回值
*/
static int  matinv(double* A, int n);
/**
@brief inverse of matrix

@param
@返回值
*/
static int  solve(const char* tr, const double* A, const double* Y, int n,
    int m, double* X);
/**
@brief least square estimation 

@param
@返回值
*/
static int  lsq(const double* A, const double* y, int n, int m, double* x,
    double* Q);
/**
@brief kalman filter

@param 
@返回值
*/
static int  filter(double* x, double* P, const double* H, const double* v,
    const double* R, int n, int m);
/**
@brief combine forward and backward filters by fixed-interval smoother as follows

@param
@返回值
*/
static int  smoother(const double* xf, const double* Qf, const double* xb,
    const double* Qb, int n, double* xs, double* Qs);
/**
@brief /* LD factorization (Q=L'*diag(D)*L)

@param
@返回值
*/
static int LD(int n, const double* Q, double* L, double* D);
/**
@brief 矩阵输出到指定文件

@param
@返回值
*/
static void matfprint(const double A[], int n, int m, int p, int q, FILE* fp);

/**
@brief 矩阵输出到指定文件

@param
@返回值
*/
static void matprint(const double A[], int n, int m, int p, int q);
};



#endif // !CMATH_H


