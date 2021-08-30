#include <direct.h>
#include "CCoordTrans.h"
CCoordTrans::CCoordTrans()
{

}

CCoordTrans::~CCoordTrans()
{

}
/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
void CCoordTrans::xyz2enu(const double* pos, double* E)
{
    double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);

    E[0] = -sinl;      E[3] = cosl;       E[6] = 0.0;
    E[1] = -sinp * cosl; E[4] = -sinp * sinl; E[7] = cosp;
    E[2] = cosp * cosl;  E[5] = cosp * sinl;  E[8] = sinp;
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void CCoordTrans::ecef2pos(const double* r, double* pos)
{
    double e2 = FE_WGS84 * (2.0 - FE_WGS84), r2 = CMath::dot(r, r, 2), z, zk, v = RE_WGS84, sinp;

    for (z = r[2], zk = 0.0; fabs(z - zk) >= 1E-4;) {
        zk = z;
        sinp = z / sqrt(r2 + z * z);
        v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);
        z = r[2] + v * e2 * sinp;
    }
    pos[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (r[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
    pos[1] = r2 > 1E-12 ? atan2(r[1], r[0]) : 0.0;
    pos[2] = sqrt(r2 + z * z) - v;
}

/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *r        O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void CCoordTrans::pos2ecef(const double* pos, double* r)
{
    double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);
    double e2 = FE_WGS84 * (2.0 - FE_WGS84), v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);

    r[0] = (v + pos[2]) * cosp * cosl;
    r[1] = (v + pos[2]) * cosp * sinl;
    r[2] = (v * (1.0 - e2) + pos[2]) * sinp;
}

/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
void CCoordTrans::ecef2enu(const double* pos, const double* r, double* e)
{
    double E[9];

    CCoordTrans::xyz2enu(pos, E);
    CMath::matmul("NN", 3, 1, 3, 1.0, E, r, 0.0, e);
}
/* transform local vector to ecef coordinate -----------------------------------
* transform local tangental coordinate vector to ecef
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *e        I   vector in local tangental coordinate {e,n,u}
*          double *r        O   vector in ecef coordinate {x,y,z}
* return : none
*-----------------------------------------------------------------------------*/
void CCoordTrans::enu2ecef(const double* pos, const double* e, double* r)
{
    double E[9];

    CCoordTrans::xyz2enu(pos, E);
    CMath::matmul("TN", 3, 1, 3, 1.0, E, e, 0.0, r);
}

/* transform covariance to local tangental coordinate --------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *P        I   covariance in ecef coordinate
*          double *Q        O   covariance in local tangental coordinate
* return : none
*-----------------------------------------------------------------------------*/
void CCoordTrans::covenu(const double* pos, const double* P, double* Q)
{
    double E[9], EP[9];

    CCoordTrans::xyz2enu(pos, E);
    CMath::matmul("NN", 3, 3, 3, 1.0, E, P, 0.0, EP);
    CMath::matmul("NT", 3, 3, 3, 1.0, EP, E, 0.0, Q);
}
/* transform local enu coordinate covariance to xyz-ecef -----------------------
* transform local enu covariance to xyz-ecef coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *Q        I   covariance in local enu coordinate
*          double *P        O   covariance in xyz-ecef coordinate
* return : none
*-----------------------------------------------------------------------------*/
void CCoordTrans::covecef(const double* pos, const double* Q, double* P)
{
    double E[9], EQ[9];

    CCoordTrans::xyz2enu(pos, E);
    CMath::matmul("TN", 3, 3, 3, 1.0, E, Q, 0.0, EQ);
    CMath::matmul("NN", 3, 3, 3, 1.0, EQ, E, 0.0, P);
}

/* astronomical arguments: f={l,l',F,D,OMG} (rad) ----------------------------*/
static void ast_args(double t, double* f)
{
    static const double fc[][5] = { /* coefficients for iau 1980 nutation */
        { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470},
        { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149},
        {  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417},
        { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169},
        { 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939}
    };
    double tt[4];
    int i, j;

    for (tt[0] = t, i = 1; i < 4; i++) tt[i] = tt[i - 1] * t;
    for (i = 0; i < 5; i++) {
        f[i] = fc[i][0] * 3600.0;
        for (j = 0; j < 4; j++) f[i] += fc[i][j + 1] * tt[j];
        f[i] = fmod(f[i] * AS2R, 2.0 * PI);
    }
}

/* iau 1980 nutation ---------------------------------------------------------*/
static void nut_iau1980(double t, const double* f, double* dpsi, double* deps)
{
    static const double nut[106][10] = {
        {   0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9},
        {   0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1},
        {   0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5},
        {   0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5},
        {   0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1},
        {   1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0},
        {   0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6},
        {   0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0},
        {   1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1},
        {   0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3},
        {  -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0},
        {  -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0},
        {   1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0},
        {   0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0},
        {  -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0},
        {  -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0},
        {   1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0},
        {  -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0},
        {   0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0},
        {   2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0},
        {   2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0},
        {   1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0},
        {   0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0},
        {  -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0},
        {   0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0},
        {   0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0},
        {  -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0},
        {   0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0},
        {   1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0},
        {   0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0},
        {   2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0},
        {  -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0},
        {   1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0},
        {   0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0},
        {   0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0},
        {   1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0},
        {  -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0},
        {   0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0},
        {   2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0},
        {   1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0},
        {   0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0},
        {   0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0},
        {   2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0},
        {   1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0},
        {   1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0},
        {   0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0},
        {   0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0},
        {   1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0},
        {   2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0},
        {   0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0},
        {   1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0},
        {   1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0},
        {  -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0},
        {   0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0},
        {   1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0},
        {   3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0},
        {   1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0},
        {  -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0},
        {   1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0},
        {  -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0},
        {   0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0},
        {  -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0},
        {   2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0},
        {   3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0},
        {   1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0},
        {   0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0},
        {   1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0},
        {   1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0},
        {   0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0},
        {   0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0},
        {   1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0},
        {   2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0},
        {  -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0},
        {   2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0},
        {   0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0},
        {   0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0},
        {   1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0},
        {   3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0},
        {  -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0},
        {  -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0},
        {   0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0},
        {  -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0},
        {   2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0},
        {   2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0},
        {   2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0},
        {   2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0},
        {   1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0},
        {  -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0},
        {  -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0}
    };
    double ang;
    int i, j;

    *dpsi = *deps = 0.0;

    for (i = 0; i < 106; i++) {
        ang = 0.0;
        for (j = 0; j < 5; j++) ang += nut[i][j] * f[j];
        *dpsi += (nut[i][6] + nut[i][7] * t) * sin(ang);
        *deps += (nut[i][8] + nut[i][9] * t) * cos(ang);
    }
    *dpsi *= 1E-4 * AS2R; /* 0.1 mas -> rad */
    *deps *= 1E-4 * AS2R;
}
/* eci to ecef transformation matrix -------------------------------------------
* compute eci to ecef transformation matrix
* args   : gtime_t tutc     I   time in utc
*          double *erpv     I   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*          double *U        O   eci to ecef transformation matrix (3 x 3)
*          double *gmst     IO  greenwich mean sidereal time (rad)
*                               (NULL: no output)
* return : none
* note   : see ref [3] chap 5
*          not thread-safe
*-----------------------------------------------------------------------------*/
void CCoordTrans::eci2ecef(gtime_t tutc, const double* erpv, double* U, double* gmst)
{
    const double ep2000[] = { 2000,1,1,12,0,0 };
    static gtime_t tutc_;
    static double U_[9], gmst_;
    gtime_t tgps;
    double eps, ze, th, z, t, t2, t3, dpsi, deps, gast, f[5];
    double R1[9], R2[9], R3[9], R[9], W[9], N[9], P[9], NP[9];
    int i;


    if (fabs(CTime::timediff(tutc, tutc_)) < 0.01) { /* read cache */
        for (i = 0; i < 9; i++) U[i] = U_[i];
        if (gmst) *gmst = gmst_;
        return;
    }
    tutc_ = tutc;

    /* terrestrial time */
    tgps = CTime::utc2gpst(tutc_);
    t = (CTime::timediff(tgps, CTime::epoch2time(ep2000)) + 19.0 + 32.184) / 86400.0 / 36525.0;
    t2 = t * t; t3 = t2 * t;

    /* astronomical arguments */
    ast_args(t, f);

    /* iau 1976 precession */
    ze = (2306.2181 * t + 0.30188 * t2 + 0.017998 * t3) * AS2R;
    th = (2004.3109 * t - 0.42665 * t2 - 0.041833 * t3) * AS2R;
    z = (2306.2181 * t + 1.09468 * t2 + 0.018203 * t3) * AS2R;
    eps = (84381.448 - 46.8150 * t - 0.00059 * t2 + 0.001813 * t3) * AS2R;
    Rz(-z, R1); Ry(th, R2); Rz(-ze, R3);
    CMath::matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
    CMath::matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

    /* iau 1980 nutation */
    nut_iau1980(t, f, &dpsi, &deps);
    Rx(-eps - deps, R1); Rz(-dpsi, R2); Rx(eps, R3);
    CMath::matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, R);
    CMath::matmul("NN", 3, 3, 3, 1.0, R, R3, 0.0, N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

    /* greenwich aparent sidereal time (rad) */
    gmst_ = CTime::utc2gmst(tutc_, erpv[2]);
    gast = gmst_ + dpsi * cos(eps);
    gast += (0.00264 * sin(f[4]) + 0.000063 * sin(2.0 * f[4])) * AS2R;

    /* eci to ecef transformation matrix */
    Ry(-erpv[0], R1); Rx(-erpv[1], R2); Rz(gast, R3);
    CMath::matmul("NN", 3, 3, 3, 1.0, R1, R2, 0.0, W);
    CMath::matmul("NN", 3, 3, 3, 1.0, W, R3, 0.0, R); /* W=Ry(-xp)*Rx(-yp) */
    CMath::matmul("NN", 3, 3, 3, 1.0, N, P, 0.0, NP);
    CMath::matmul("NN", 3, 3, 3, 1.0, R, NP, 0.0, U_); /* U=W*Rz(gast)*N*P */

    for (i = 0; i < 9; i++) U[i] = U_[i];
    if (gmst) *gmst = gmst_;

   
}

/* convert degree to deg-min-sec -----------------------------------------------
* convert degree to degree-minute-second
* args   : double deg       I   degree
*          double *dms      O   degree-minute-second {deg,min,sec}
*          int    ndec      I   number of decimals of second
* return : none
*-----------------------------------------------------------------------------*/
void CCoordTrans::deg2dms(double deg, double* dms, int ndec)
{
    double sign = deg < 0.0 ? -1.0 : 1.0, a = fabs(deg);
    double unit = pow(0.1, ndec);
    dms[0] = floor(a); a = (a - dms[0]) * 60.0;
    dms[1] = floor(a); a = (a - dms[1]) * 60.0;
    dms[2] = floor(a / unit + 0.5) * unit;
    if (dms[2] >= 60.0) {
        dms[2] = 0.0;
        dms[1] += 1.0;
        if (dms[1] >= 60.0) {
            dms[1] = 0.0;
            dms[0] += 1.0;
        }
    }
    dms[0] *= sign;
}

/*convert degree - minute - second to degree
* args   : double* dms      I   degree - minute - second{ deg,min,sec }
*return : degree
* ---------------------------------------------------------------------------- - */
double CCoordTrans::dms2deg(const double* dms)
{
    double sign = dms[0] < 0.0 ? -1.0 : 1.0;
    return sign * (fabs(dms[0]) + dms[1] / 60.0 + dms[2] / 3600.0);
}