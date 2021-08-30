#ifndef CBROADCASTEPH_H
#define CBROADCASTEPH_H
#include "CommonVariables.h"
#include "CommonFunction.h"
#include "dpi_types_basic.h"
#include "CTime.h"
#include "CMath.h"
#define SIN_5        -0.0871557427476582            /* sin(-5.0 deg) */
#define COS_5         0.9961946980917456            /* cos(-5.0 deg) */
#define STD_GAL_NAPA  500.0                         /* error of galileo ephemeris for NAPA (m) */

/*
@@brief 广播星历类，通过卫星星历计算卫星位置和卫星钟差；
*/
class CBroadcastEph
{
public:
	CBroadcastEph();
	~CBroadcastEph();
public:
	/**
	@brief  计算卫星位置
	@param  
	@返回值 
	*/
	int  CalSatPos(gtime_t tiEphTime, SinStaObsData* pSinSatObsData, AllBroadcastEphData* pAllBroadcastEphData,
		           double* rs, double* dts, double* var, int* svh);
	/**
	@brief  计算卫星钟差
	@param  
	@返回值 
	*/
	void CalSatClk(gtime_t tiEphTime, SinStaObsData* pSinSatObsData, AllBroadcastEphData*pALLBroadcastEphData,
		           double* pSatClk, double* pSatClkVar);

	/**
	@brief   选择计算单颗卫星位置模式
	@param
	@返回值
	*/
	int  SatPos(gtime_t time, gtime_t teph, int sat, int ephopt, AllBroadcastEphData* pALLBroadcastEphData,
		        double* rs, double* dts, double* var, int* svh);
private:
	/**
	@brief   选择计算单颗卫星钟差模式
	@param
	@返回值
	*/
	int EphClk(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData, double* dts);
	/**
	@brief   根据星历得到相应的卫星系统
	@param
	@返回值
	*/
	int GetEphSys(int sys);
	/**
	@brief   选择与观测时间最近的星历
	@param
	@返回值
	*/
	BroadcastEphData* SelEph(gtime_t time, int sat, int iode, AllBroadcastEphData* pBroadcastEphData);

	/**
	@brief   计算单颗卫星位置
	@param
	@返回值
	*/
	void Eph2Pos(gtime_t time, BroadcastEphData* eph, double* rs, double* dts,double* var);
	int  EphPos(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData,
		        int iode, double* rs, double* dts, double* var, int* svh);
	/**
	@brief   根据ure参数确定卫星位置的精度
	@param
	@返回值
	*/
	double EphuraVar(int sys, int ura);
private:
	BroadcastEphData* m_pGpsBroadEphData ;    /*GPS广播星历数据结构体指针*/
	BroadcastEphData* m_pGalBroadEphData ;    /*Gal广播星历数据结构体指针*/
	BroadcastEphData* m_pBdsBroadEphData ;    /*BDS广播星历数据结构体指针*/

};

#endif // !CBROADCASTEPH_H


