#ifndef CSTATION_H
#define CSTATION_H
#include "CommonVariables.h"
#include"CommonFunction.h"
#include "CBroadcastEph.h"
#include "CTroposphere.h"
#include "CCoordTrans.h"
#include"CIonosphere.h"
#include "CAntenna.h"
#include "CTide.h"
#include"CMath.h"
/*
@@brief 基站类，处理基站观测信息，将处理好的信息存放在基站处理信息结构体中；
*/
class CStation
{
public:
	CStation();
	~CStation();
public:
	
	/**
	@brief  计算站心距
	@notes    进行了sagnac效应的修正
	*/
	double  CalGeoDist(double* rs, const double* rr, double* e);

	/**
	@brief  计算卫星方位角和高度角
	
	@返回值 卫星高度角 （rad）
	*/
	double CalSatAzEl(double* pos, double* e, double* azel);
	
	/**
	@brief  计算卫星非差相位和伪距残差
	@param
	@返回值
	*/
	void   CalZeroSatRes(double r, ObsData* pObsData, double* azel, double* y_p, double* y_l, double* freq);
	/**
	@brief  处理单站观测数据
	@param
	@返回值
	*/
	bool   Process();
	/**
	@brief  处理n个测站的观测数据
	@param
	@返回值
	*/
	bool   Process(int n);
	/**
	@brief  对测站观测值和星历进行赋值
	@param  SinStaObsData[in]        基站观测数据结构体指针
	@param  pAllbroadcastephdata[in]  观测星历结构体指针
	@返回值
	*/
	void SetObsEphInfo(SinStaObsData* pSinstaobsdata, AllBroadcastEphData* pAllbroadcastephdata);
private:
	/**
	@brief   通过信噪比筛选卫星
	@param
	@返回值
	*/
	bool TestSNRatio();
	/**
	@brief  对观测值进行线性组合
	@param
	@返回值
	*/
	void ObsLineComb();
	
public:
	SinStaObsData       *m_pSinstaobsdata = NULL         ;     /*单站观测数据结构体指针*/
	AllBroadcastEphData *m_pAllbroadcastephdata = NULL   ;     /*星历观测数据结构体指针*/



};

#endif // !CSTATION_H


