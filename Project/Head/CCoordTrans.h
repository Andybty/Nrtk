#ifndef CCOORDTRANS_H
#define CCOORDTRANS_H
#include "dpi_types_basic.h"
#include "CMath.h"
#include "CTime.h"
#define PI          3.1415926535897932                 /* pi */
#define D2R         (PI/180.0)                         /* deg to rad */
#define AS2R        (D2R/3600.0)                       /* arc sec to radian */
#define RE_WGS84    6378137.0                          /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563)                /* earth flattening (WGS84) */
/* coordinate rotation matrix------------------------------------------------*/
#define Rx(t,X) do { \
    (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
    (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
    (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
    (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

class CCoordTrans
{
public:
	CCoordTrans();
	~CCoordTrans();
public:
	/**
	@brief  ecef to local coordinate transfromation matrix
	@param
	@返回值
	*/
	static void xyz2enu(const double* pos, double* E);
	/**
	@brief  transform ecef to geodetic postion
	@param
	@返回值
	*/
	static void ecef2pos(const double* r, double* pos);
	/**
	@brief  transform geodetic to ecef position 
	@param
	@返回值
	*/
	static void pos2ecef(const double* pos, double* r);
	/**
	@brief transform ecef vector to local tangental coordinate
	@param
	@返回值
	*/
	static void ecef2enu(const double* pos, const double* r, double* e);
	/**
	@brief  transform local vector to ecef coordinate
	@param
	@返回值
	*/
	static void enu2ecef(const double* pos, const double* e, double* r);
	/**
	@brief transform covariance to local tangental coordinate 
	@param
	@返回值
	*/
	static void covenu(const double* pos, const double* P, double* Q);
	/**
	@brief transform local enu coordinate covariance to xyz-ecef
	@param
	@返回值
	*/
	static void covecef(const double* pos, const double* Q, double* P);	
	/**
	@brief eci to ecef transformation matrix
	@param
	@返回值
	*/
	static void eci2ecef(gtime_t tutc, const double* erpv, double* U, double* gmst);
	/**
	@brief  convert degree to deg-min-sec
	@param
	@返回值
	*/
	static void deg2dms(double deg, double* dms, int ndec);
	//static void eci2ecef(gtime_t tutc, const double* erpv, double* U, double* gmst);
	/**
	@brief  convert degree - minute - second to degree
	@param
	@返回值
	*/
	static double dms2deg(const double* dms);

};



#endif // !CCOORDTRANS_H



