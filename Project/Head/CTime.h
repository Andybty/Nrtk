#ifndef CTIME_H
#define CTIME_H
#include "dpi_types_basic.h"
#include<windows.h>
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define PI          3.1415926535897932  /* pi */


static const double gpst0[] = { 1980,1, 6,0,0,0 }; /* gps time reference */
static const double gst0[] =  { 1999,8,22,0,0,0 }; /* galileo system time reference */
static const double bdt0[] =  { 2006,1, 1,0,0,0 }; /* beidou time reference */
static double leaps[MAXLEAPS + 1][7] = { /* leap seconds (y,m,d,h,m,s,utc-gpst) */
{2017,1,1,0,0,0,-18},
{2015,7,1,0,0,0,-17},
{2012,7,1,0,0,0,-16},
{2009,1,1,0,0,0,-15},
{2006,1,1,0,0,0,-14},
{1999,1,1,0,0,0,-13},
{1997,7,1,0,0,0,-12},
{1996,1,1,0,0,0,-11},
{1994,7,1,0,0,0,-10},
{1993,7,1,0,0,0, -9},
{1992,7,1,0,0,0, -8},
{1991,1,1,0,0,0, -7},
{1990,1,1,0,0,0, -6},
{1988,1,1,0,0,0, -5},
{1985,7,1,0,0,0, -4},
{1983,7,1,0,0,0, -3},
{1982,7,1,0,0,0, -2},
{1981,7,1,0,0,0, -1},
{0}
};
class CTime
{
public:
CTime();
~CTime();
public:
/**
@brief  string to time
@param
@返回值
*/
static int     str2time(const char* s, int i, int n, gtime_t* t);
/**
@brief time to string
@param
@返回值
*/
static void    time2str(gtime_t t, char* str, int n);
/**
@brief convert calendar day/time to time
@param
@返回值
*/
static gtime_t epoch2time(const double* ep);
/**
@brief convert gtime_t struct to calendar day/time
@param
@返回值
*/
static void    time2epoch(gtime_t t, double* ep);
/**
@brief gps time to time 
@param
@返回值
*/
static gtime_t gpst2time(int week, double sec);
/**
@brief time to gps time
@param
@返回值
*/
static double  time2gpst(gtime_t t, int* week);
/**
@brief galileo system time to time
@param
@返回值
*/
static gtime_t gst2time(int week, double sec);
/**
@brief time to galileo system time 
@param
@返回值
*/
static double  time2gst(gtime_t t, int* week);
/**
@brief  beidou time (bdt) to time
@param
@返回值
*/
static gtime_t bdt2time(int week, double sec);
/**
@brief  time to beidouo time (bdt)
@param
@返回值
*/
static double  time2bdt(gtime_t t, int* week);
/**
@brief  get time string 
@param
@返回值
*/
static char* time_str(gtime_t t, int n);
/**
@brief add time
@param
@返回值
*/
static gtime_t timeadd(gtime_t t, double sec);
/**
@brief  time difference 
@param
@返回值
*/
static double  timediff(gtime_t t1, gtime_t t2);
/**
@brief  gpstime to utc
@param
@返回值
*/
static gtime_t gpst2utc(gtime_t t);
/**
@brief get current time in utc
@param
@返回值
*/
static gtime_t timeget(void);
/**
@brief  set current time in utc
@param
@返回值
*/
static void    timeset(gtime_t t);
/**
@brief  * reset current time
@param
@返回值
*/
static void    timereset(void);
/**
@brief   convert utc to gpstime considering leap seconds
@param
@返回值
*/
static gtime_t utc2gpst(gtime_t t);
/**
@brief   convert gpstime to bdt (beidou navigation satellite system time)
@param
@返回值
*/
static gtime_t gpst2bdt(gtime_t t);
/**
@brief    convert bdt (beidou navigation satellite system time) to gpstime
@param
@返回值
*/
static gtime_t bdt2gpst(gtime_t t);
/**
@brief    convert time to day of year
@param
@返回值
*/
static double  time2doy(gtime_t t);
/**
@brief    convert utc to gmst (Greenwich mean sidereal time)
@param
@返回值
*/
static double  utc2gmst(gtime_t t, double ut1_utc);
/**
@brief    adjust gps week number using cpu time
@param
@返回值
*/
static int adjgpsweek(int week);
/**
@brief   get current tick in ms
@param
@返回值
*/
static uint32_t tickget(void);
/**
@brief   sleep ms
@param
@返回值
*/
static void sleepms(int ms);
/*
	
	
	
	
*/
};







#endif // !CTIME_H


