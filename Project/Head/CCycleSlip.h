#ifndef CCYCLESLIP_H
#define CCYCLESLIP_H
#include "CommonVariables.h"
/*
@@brief �����࣬���ò�ͬ�ķ���̽������

*/
class CCycleSlip
{
public:
	CCycleSlip(double dMwThreshold, double dGfThreshold);
	~CCycleSlip();
public:
	/**
	@brief  MW����̽������
	@param
	@����ֵ
	*/
	void  DetCycleSlipGF(SinStaObsData* pSinSatObsData);
	/**
	@brief  GF����̽������
	@param
	@����ֵ
	*/
	void  DetCycleSlipMW(SinStaObsData* pSinSatObsData);


private:
    double m_dMWThreshold  ; /*MW̽��������ֵ*/
	double m_dGFthreshold  ; /*GF̽��������ֵ*/
};




#endif // !CCYCLESLIP_H



