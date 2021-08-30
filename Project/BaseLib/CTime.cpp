#include <direct.h>
#include "CTime.h"
CTime::CTime()
{
}
CTime::~CTime()
{
}
/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
int CTime::str2time(const char* s, int i, int n, gtime_t* t)
{
    double ep[6];
    char str[256], * p = str;

    if (i < 0 || (int)strlen(s) < i || (int)sizeof(str) - 1 < i) return -1;
    for (s += i; *s && --n >= 0;) *p++ = *s++;
    *p = '\0';
    if (sscanf(str, "%lf %lf %lf %lf %lf %lf", ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5) < 6)
        return -1;
    if (ep[0] < 100.0) ep[0] += ep[0] < 80.0 ? 2000.0 : 1900.0;
    *t = epoch2time(ep);
    return 0;
}

/* time to string --------------------------------------------------------------
* convert gtime_t struct to string
* args   : gtime_t t        I   gtime_t struct
*          char   *s        O   string ("yyyy/mm/dd hh:mm:ss.ssss")
*          int    n         I   number of decimals
* return : none
*-----------------------------------------------------------------------------*/
void CTime::time2str(gtime_t t, char* s, int n)
{
    double ep[6];

    if (n < 0) n = 0; else if (n > 12) n = 12;
    if (1.0 - t.dSec < 0.5 / pow(10.0, n)) { t.tiTime++; t.dSec = 0.0; };
    time2epoch(t, ep);
    sprintf(s, "%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f", ep[0], ep[1], ep[2],
        ep[3], ep[4], n <= 0 ? 2 : n + 3, n <= 0 ? 0 : n, ep[5]);
}
/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
 gtime_t CTime::epoch2time(const double* ep)
{
    const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
    gtime_t time = { 0 };
    int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];

    if (year < 1970 || 2099 < year || mon < 1 || 12 < mon) return time;

    /* leap year if year%4==0 in 1901-2099 */
    days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
    sec = (int)floor(ep[5]);
    time.tiTime = (time_t)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
    time.dSec= ep[5] - sec;
    return time;
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
 void CTime::time2epoch(gtime_t t, double* ep)
{
    const int mday[] = { /* # of days in a month */
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days, sec, mon, day;

    /* leap year if year%4==0 in 1901-2099 */
    days = (int)(t.tiTime / 86400);
    sec = (int)(t.tiTime - (time_t)days * 86400);
    for (day = days % 1461, mon = 0; mon < 48; mon++) {
        if (day >= mday[mon]) day -= mday[mon]; else break;
    }
    ep[0] = 1970 + days / 1461 * 4 + mon / 12; ep[1] = mon % 12 + 1; ep[2] = day + 1;
    ep[3] = sec / 3600; ep[4] = sec % 3600 / 60; ep[5] = sec % 60 + t.dSec;
}
/* gps time to time ------------------------------------------------------------
* convert week and tow in gps time to gtime_t struct
* args   : int    week      I   week number in gps time
*          double sec       I   time of week in gps time (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
 gtime_t CTime::gpst2time(int week, double sec)
{
    gtime_t t = CTime::epoch2time(gpst0);
    char tstr[64];
    if (sec < -1E9 || 1E9 < sec) sec = 0.0;
    t.tiTime += (time_t)86400 * 7 * week + (int)sec;
    t.dSec = sec - (int)sec;
    CTime::time2str(t, tstr, 5);
    return t;
}
/* time to gps time ------------------------------------------------------------
* convert gtime_t struct to week and tow in gps time
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gps time (NULL: no output)
* return : time of week in gps time (s)
*-----------------------------------------------------------------------------*/
 double CTime::time2gpst(gtime_t t, int* week)
{
    gtime_t t0 = CTime::epoch2time(gpst0);
    time_t sec = t.tiTime - t0.tiTime;
    int w = (int)(sec / (86400 * 7));

    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7) + t.dSec;
}
/* galileo system time to time -------------------------------------------------
* convert week and tow in galileo system time (gst) to gtime_t struct
* args   : int    week      I   week number in gst
*          double sec       I   time of week in gst (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
 gtime_t CTime::gst2time(int week, double sec)
{
    gtime_t t = CTime::epoch2time(gst0);

    if (sec < -1E9 || 1E9 < sec) sec = 0.0;
    t.tiTime += (time_t)86400 * 7 * week + (int)sec;
    t.dSec = sec - (int)sec;
    return t;
}
/* time to galileo system time -------------------------------------------------
* convert gtime_t struct to week and tow in galileo system time (gst)
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gst (NULL: no output)
* return : time of week in gst (s)
*-----------------------------------------------------------------------------*/
 double CTime::time2gst(gtime_t t, int* week)
{
    gtime_t t0 = CTime::epoch2time(gst0);
    time_t sec = t.tiTime - t0.tiTime;
    int w = (int)(sec / (86400 * 7));

    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7) + t.dSec;
}
/* beidou time (bdt) to time ---------------------------------------------------
* convert week and tow in beidou time (bdt) to gtime_t struct
* args   : int    week      I   week number in bdt
*          double sec       I   time of week in bdt (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
 gtime_t CTime::bdt2time(int week, double sec)
{
    gtime_t t = epoch2time(bdt0);

    if (sec < -1E9 || 1E9 < sec) sec = 0.0;
    t.tiTime += (time_t)86400 * 7 * week + (int)sec;
    t.dSec = sec - (int)sec;
    return t;
}
/* time to beidouo time (bdt) --------------------------------------------------
* convert gtime_t struct to week and tow in beidou time (bdt)
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in bdt (NULL: no output)
* return : time of week in bdt (s)
*-----------------------------------------------------------------------------*/
 double CTime::time2bdt(gtime_t t, int* week)
{
    gtime_t t0 = CTime::epoch2time(bdt0);
    time_t sec = t.tiTime - t0.tiTime;
    int w = (int)(sec / (86400 * 7));

    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7) + t.dSec;
}

 /* get time string -------------------------------------------------------------
* get time string
* args   : gtime_t t        I   gtime_t struct
*          int    n         I   number of decimals
* return : time string
* notes  : not reentrant, do not use multiple in a function
*-----------------------------------------------------------------------------*/
 char* CTime::time_str(gtime_t t, int n)
 {
     static char buff[64];
     time2str(t, buff, n);
     return buff;
 }
/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (s)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
 gtime_t CTime::timeadd(gtime_t t, double sec)
{
    double tt;

    t.dSec += sec; tt = floor(t.dSec); t.tiTime += (int)tt; t.dSec -= tt;
    return t;
}
/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (s)
*-----------------------------------------------------------------------------*/
 double CTime::timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.tiTime, t2.tiTime) + t1.dSec - t2.dSec;
}
/* get current time in utc -----------------------------------------------------
* get current time in utc
* args   : none
* return : current time in utc
*-----------------------------------------------------------------------------*/
 static double timeoffset_ = 0.0;        /* time offset (s) */

 gtime_t CTime::timeget(void)
{
    gtime_t time;
    double ep[6] = { 0 };
#ifdef WIN32
    SYSTEMTIME ts;

    GetSystemTime(&ts); /* utc */
    ep[0] = ts.wYear; ep[1] = ts.wMonth;  ep[2] = ts.wDay;
    ep[3] = ts.wHour; ep[4] = ts.wMinute; ep[5] = ts.wSecond + ts.wMilliseconds * 1E-3;
#else
    struct timeval tv;
    struct tm* tt;

    if (!gettimeofday(&tv, NULL) && (tt = gmtime(&tv.tv_sec))) {
        ep[0] = tt->tm_year + 1900; ep[1] = tt->tm_mon + 1; ep[2] = tt->tm_mday;
        ep[3] = tt->tm_hour; ep[4] = tt->tm_min; ep[5] = tt->tm_sec + tv.tv_usec * 1E-6;
    }
#endif
    time = CTime::epoch2time(ep);

#ifdef CPUTIME_IN_GPST /* cputime operated in gpst */
    time = gpst2utc(time);
#endif
    return CTime::timeadd(time, timeoffset_);
}
/* set current time in utc -----------------------------------------------------
* set current time in utc
* args   : gtime_t          I   current time in utc
* return : none
* notes  : just set time offset between cpu time and current time
*          the time offset is reflected to only timeget()
*          not reentrant
*-----------------------------------------------------------------------------*/
 void CTime::timeset(gtime_t t)
{
    timeoffset_ += CTime::timediff(t, CTime::timeget());
}
/* reset current time ----------------------------------------------------------
* reset current time
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
 void CTime::timereset(void)
{
    timeoffset_ = 0.0;
}

 /* utc to gpstime --------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* args   : gtime_t t        I   time expressed in utc
* return : time expressed in gpstime
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
 gtime_t CTime::utc2gpst(gtime_t t)
 {
     int i;

     for (i = 0; leaps[i][0] > 0; i++) {
         if (CTime::timediff(t, CTime::epoch2time(leaps[i])) >= 0.0) return CTime::timeadd(t, -leaps[i][6]);
     }
     return t;
 }
 /* gpstime to utc --------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in utc
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
  gtime_t CTime::gpst2utc(gtime_t t)
 {
     gtime_t tu;
     int i;

     for (i = 0; leaps[i][0] > 0; i++) {
         tu = CTime::timeadd(t, leaps[i][6]);
         if (CTime::timediff(tu, CTime::epoch2time(leaps[i])) >= 0.0) return tu;
     }
     return t;
 }
 /* gpstime to bdt --------------------------------------------------------------
* convert gpstime to bdt (beidou navigation satellite system time)
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in bdt
* notes  : ref [8] 3.3, 2006/1/1 00:00 BDT = 2006/1/1 00:00 UTC
*          no leap seconds in BDT
*          ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
 gtime_t CTime::gpst2bdt(gtime_t t)
 {
     return CTime::timeadd(t, -14.0);
 }
 /* bdt to gpstime --------------------------------------------------------------
* convert bdt (beidou navigation satellite system time) to gpstime
* args   : gtime_t t        I   time expressed in bdt
* return : time expressed in gpstime
* notes  : see gpst2bdt()
*-----------------------------------------------------------------------------*/
 gtime_t CTime::bdt2gpst(gtime_t t)
 {
     return CTime::timeadd(t, 14.0);
 }

 /* time to day of year ---------------------------------------------------------
* convert time to day of year
* args   : gtime_t t        I   gtime_t struct
* return : day of year (days)
*-----------------------------------------------------------------------------*/
 double CTime::time2doy(gtime_t t)
 {
     double ep[6];

     CTime::time2epoch(t, ep);
     ep[1] = ep[2] = 1.0; ep[3] = ep[4] = ep[5] = 0.0;
     return CTime::timediff(t, CTime::epoch2time(ep)) / 86400.0 + 1.0;
 }
 static double time2sec(gtime_t time, gtime_t* day)
 {
     double ep[6], sec;
     CTime::time2epoch(time, ep);
     sec = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
     ep[3] = ep[4] = ep[5] = 0.0;
     *day = CTime::epoch2time(ep);
     return sec;
 }
 /* utc to gmst -----------------------------------------------------------------
* convert utc to gmst (Greenwich mean sidereal time)
* args   : gtime_t t        I   time expressed in utc
*          double ut1_utc   I   UT1-UTC (s)
* return : gmst (rad)
*-----------------------------------------------------------------------------*/
 double CTime::utc2gmst(gtime_t t, double ut1_utc)
 {
     const double ep2000[] = { 2000,1,1,12,0,0 };
     gtime_t tut, tut0;
     double ut, t1, t2, t3, gmst0, gmst;

     tut = CTime::timeadd(t, ut1_utc);
     ut = time2sec(tut, &tut0);
     t1 = CTime::timediff(tut0, CTime::epoch2time(ep2000)) / 86400.0 / 36525.0;
     t2 = t1 * t1; t3 = t2 * t1;
     gmst0 = 24110.54841 + 8640184.812866 * t1 + 0.093104 * t2 - 6.2E-6 * t3;
     gmst = gmst0 + 1.002737909350795 * ut;

     return fmod(gmst, 86400.0) * PI / 43200.0; /* 0 <= gmst <= 2*PI */
 }

 /* adjust gps week number ------------------------------------------------------
* adjust gps week number using cpu time
* args   : int   week       I   not-adjusted gps week number (0-1023)
* return : adjusted gps week number
*-----------------------------------------------------------------------------*/
int CTime::adjgpsweek(int week)
 {
     int w;
     (void)CTime::time2gpst(CTime::utc2gpst(CTime::timeget()), &w);
     if (w < 1560) w = 1560; /* use 2009/12/1 if time is earlier than 2009/12/1 */
     return week + (w - week + 1) / 1024 * 1024;
 }

/*get tick time-------------------------------------------------------------- -
*get current tick in ms
* args   : none
* return : current tick in ms
* ---------------------------------------------------------------------------- - */
 uint32_t CTime::tickget(void)
{
#ifdef WIN32
    return (uint32_t)timeGetTime();
#else
    struct timespec tp = { 0 };
    struct timeval  tv = { 0 };

#ifdef CLOCK_MONOTONIC_RAW
    /* linux kernel > 2.6.28 */
    if (!clock_gettime(CLOCK_MONOTONIC_RAW, &tp)) {
        return tp.tv_sec * 1000u + tp.tv_nsec / 1000000u;
    }
    else {
        gettimeofday(&tv, NULL);
        return tv.tv_sec * 1000u + tv.tv_usec / 1000u;
    }
#else
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000u + tv.tv_usec / 1000u;
#endif
#endif /* WIN32 */
}

 /*sleep ms--------------------------------------------------------------------
     * sleep ms
     * args   : int   ms         I   miliseconds to sleep(< 0:no sleep)
     * return : none
     * ---------------------------------------------------------------------------- - */
void CTime::sleepms(int ms)
 {
#ifdef WIN32
     if (ms < 5) Sleep(1); else Sleep(ms);
#else
     struct timespec ts;
     if (ms <= 0) return;
     ts.tv_sec = (time_t)(ms / 1000);
     ts.tv_nsec = (long)(ms % 1000 * 1000000);
     nanosleep(&ts, NULL);
#endif
 }