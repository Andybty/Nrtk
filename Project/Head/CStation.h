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
@@brief ��վ�࣬�����վ�۲���Ϣ��������õ���Ϣ����ڻ�վ������Ϣ�ṹ���У�
*/
class CStation
{
public:
	CStation();
	~CStation();
public:
	
	/**
	@brief  ����վ�ľ�
	@notes    ������sagnacЧӦ������
	*/
	double  CalGeoDist(double* rs, const double* rr, double* e);

	/**
	@brief  �������Ƿ�λ�Ǻ͸߶Ƚ�
	
	@����ֵ ���Ǹ߶Ƚ� ��rad��
	*/
	double CalSatAzEl(double* pos, double* e, double* azel);
	
	/**
	@brief  �������Ƿǲ���λ��α��в�
	@param
	@����ֵ
	*/
	void   CalZeroSatRes(double r, ObsData* pObsData, double* azel, double* y_p, double* y_l, double* freq);
	/**
	@brief  ����վ�۲�����
	@param
	@����ֵ
	*/
	bool   Process();
	/**
	@brief  ����n����վ�Ĺ۲�����
	@param
	@����ֵ
	*/
	bool   Process(int n);
	/**
	@brief  �Բ�վ�۲�ֵ���������и�ֵ
	@param  SinStaObsData[in]        ��վ�۲����ݽṹ��ָ��
	@param  pAllbroadcastephdata[in]  �۲������ṹ��ָ��
	@����ֵ
	*/
	void SetObsEphInfo(SinStaObsData* pSinstaobsdata, AllBroadcastEphData* pAllbroadcastephdata);
private:
	/**
	@brief   ͨ�������ɸѡ����
	@param
	@����ֵ
	*/
	bool TestSNRatio();
	/**
	@brief  �Թ۲�ֵ�����������
	@param
	@����ֵ
	*/
	void ObsLineComb();
	
public:
	SinStaObsData       *m_pSinstaobsdata = NULL         ;     /*��վ�۲����ݽṹ��ָ��*/
	AllBroadcastEphData *m_pAllbroadcastephdata = NULL   ;     /*�����۲����ݽṹ��ָ��*/



};

#endif // !CSTATION_H


