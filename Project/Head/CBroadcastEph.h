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
@@brief �㲥�����࣬ͨ������������������λ�ú������Ӳ
*/
class CBroadcastEph
{
public:
	CBroadcastEph();
	~CBroadcastEph();
public:
	/**
	@brief  ��������λ��
	@param  
	@����ֵ 
	*/
	int  CalSatPos(gtime_t tiEphTime, SinStaObsData* pSinSatObsData, AllBroadcastEphData* pAllBroadcastEphData,
		           double* rs, double* dts, double* var, int* svh);
	/**
	@brief  ���������Ӳ�
	@param  
	@����ֵ 
	*/
	void CalSatClk(gtime_t tiEphTime, SinStaObsData* pSinSatObsData, AllBroadcastEphData*pALLBroadcastEphData,
		           double* pSatClk, double* pSatClkVar);

	/**
	@brief   ѡ����㵥������λ��ģʽ
	@param
	@����ֵ
	*/
	int  SatPos(gtime_t time, gtime_t teph, int sat, int ephopt, AllBroadcastEphData* pALLBroadcastEphData,
		        double* rs, double* dts, double* var, int* svh);
private:
	/**
	@brief   ѡ����㵥�������Ӳ�ģʽ
	@param
	@����ֵ
	*/
	int EphClk(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData, double* dts);
	/**
	@brief   ���������õ���Ӧ������ϵͳ
	@param
	@����ֵ
	*/
	int GetEphSys(int sys);
	/**
	@brief   ѡ����۲�ʱ�����������
	@param
	@����ֵ
	*/
	BroadcastEphData* SelEph(gtime_t time, int sat, int iode, AllBroadcastEphData* pBroadcastEphData);

	/**
	@brief   ���㵥������λ��
	@param
	@����ֵ
	*/
	void Eph2Pos(gtime_t time, BroadcastEphData* eph, double* rs, double* dts,double* var);
	int  EphPos(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData,
		        int iode, double* rs, double* dts, double* var, int* svh);
	/**
	@brief   ����ure����ȷ������λ�õľ���
	@param
	@����ֵ
	*/
	double EphuraVar(int sys, int ura);
private:
	BroadcastEphData* m_pGpsBroadEphData ;    /*GPS�㲥�������ݽṹ��ָ��*/
	BroadcastEphData* m_pGalBroadEphData ;    /*Gal�㲥�������ݽṹ��ָ��*/
	BroadcastEphData* m_pBdsBroadEphData ;    /*BDS�㲥�������ݽṹ��ָ��*/

};

#endif // !CBROADCASTEPH_H


