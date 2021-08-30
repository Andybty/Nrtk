#ifndef CCYCLESLIP_H
#define CCYCLESLIP_H
#include "CommonVariables.h"
/*
@@brief 周跳类，采用不同的方法探测周跳

*/
class CCycleSlip
{
public:
	CCycleSlip(double dMwThreshold, double dGfThreshold);
	~CCycleSlip();
public:
	/**
	@brief  MW方法探测周跳
	@param
	@返回值
	*/
	void  DetCycleSlipGF(SinStaObsData* pSinSatObsData);
	/**
	@brief  GF方法探测周跳
	@param
	@返回值
	*/
	void  DetCycleSlipMW(SinStaObsData* pSinSatObsData);


private:
    double m_dMWThreshold  ; /*MW探测周跳阈值*/
	double m_dGFthreshold  ; /*GF探测周跳阈值*/
};




#endif // !CCYCLESLIP_H



