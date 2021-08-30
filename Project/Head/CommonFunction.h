#ifndef COMMONFUNCTION_H
#define COMMONFUNCTION__H
#include"CommonVariables.h"
#include "CRTCM.h"
/**
@brief 通过卫星号计算卫星系统

@param
@返回值
*/
int  SatSys(int sat, int* prn);
/**
@brief 通过卫星号和码通道获取卫星频率

@param
@返回值
*/
double Sat2Freq(int sat, uint8_t code);
/**
@brief 分配观测值结构体内存空间

@param
@返回值
*/
void InitObs(SinStaObsData*& pSinstaobsdata, int n);
/**
@brief 释放观测值的内存空间

@param
@返回值
*/
void FreeObs(SinStaObsData* pSinstaobsdata);
/**
@brief 分配星历观测值内存空间

@param
@返回值
*/
void InitBroadcastEphData(AllBroadcastEphData* &pAllbroadcastephdata);
/**
@brief 释放星历观测值内存空间

@param
@返回值
*/
void FreeBroadcastEphData(AllBroadcastEphData* pAllbroadcastephdata);
/**
@brief 释放站解算信息内存空间

@param
@返回值
*/
void FreeStaSolveInfo(int n,StaSolveInfo*& pStaSolveInfo);
/**
@brief 对基站的解算信息分配内存空间

@param
@返回值
*/
void InitStaSolveInfo(int n,StaSolveInfo*& pStaSolveInfo);
/**
@brief 对基线的解算信息分配内存空间

@param
@返回值
*/
bool InitBLineInfo(int n, BLineInfo*& pBlineInfo);
/**
@brief 释放基线解算信息的内存空间

@param
@返回值
*/
void FreeBLinInfo(BLineInfo* pBlineInfo);
/**
@brief 对基线静态信息分配内存空间

@param
@返回值
*/
bool InitRStaInfo(int n, RStaInfo*& pRStaInfo);
/**
@brief 释放基线静态信息空间

@param
@返回值
*/
void FreeRStaInfo(RStaInfo* pRStaInfo);

/**
@brief 对基线解算信息分配空间

@param
@返回值
*/
bool InitBLineSlon(int n, BLineSlon*& gBlineSlon);
/**
@brief 释放基线解算信息空间

@param
@返回值
*/
void FreeBLineSlon(int n, BLineSlon*& gBlineSlon);

/**
@brief 对三角网信息进行初始化
@param
@返回值
*/
void InitTriNetInfo(int n, TriNetInfo*&  gTrinetInfo);
#endif // !COMMONFUNCTION_H


