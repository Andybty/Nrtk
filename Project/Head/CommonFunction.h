#ifndef COMMONFUNCTION_H
#define COMMONFUNCTION__H
#include"CommonVariables.h"
#include "CRTCM.h"
/**
@brief ͨ�����Ǻż�������ϵͳ

@param
@����ֵ
*/
int  SatSys(int sat, int* prn);
/**
@brief ͨ�����Ǻź���ͨ����ȡ����Ƶ��

@param
@����ֵ
*/
double Sat2Freq(int sat, uint8_t code);
/**
@brief ����۲�ֵ�ṹ���ڴ�ռ�

@param
@����ֵ
*/
void InitObs(SinStaObsData*& pSinstaobsdata, int n);
/**
@brief �ͷŹ۲�ֵ���ڴ�ռ�

@param
@����ֵ
*/
void FreeObs(SinStaObsData* pSinstaobsdata);
/**
@brief ���������۲�ֵ�ڴ�ռ�

@param
@����ֵ
*/
void InitBroadcastEphData(AllBroadcastEphData* &pAllbroadcastephdata);
/**
@brief �ͷ������۲�ֵ�ڴ�ռ�

@param
@����ֵ
*/
void FreeBroadcastEphData(AllBroadcastEphData* pAllbroadcastephdata);
/**
@brief �ͷ�վ������Ϣ�ڴ�ռ�

@param
@����ֵ
*/
void FreeStaSolveInfo(int n,StaSolveInfo*& pStaSolveInfo);
/**
@brief �Ի�վ�Ľ�����Ϣ�����ڴ�ռ�

@param
@����ֵ
*/
void InitStaSolveInfo(int n,StaSolveInfo*& pStaSolveInfo);
/**
@brief �Ի��ߵĽ�����Ϣ�����ڴ�ռ�

@param
@����ֵ
*/
bool InitBLineInfo(int n, BLineInfo*& pBlineInfo);
/**
@brief �ͷŻ��߽�����Ϣ���ڴ�ռ�

@param
@����ֵ
*/
void FreeBLinInfo(BLineInfo* pBlineInfo);
/**
@brief �Ի��߾�̬��Ϣ�����ڴ�ռ�

@param
@����ֵ
*/
bool InitRStaInfo(int n, RStaInfo*& pRStaInfo);
/**
@brief �ͷŻ��߾�̬��Ϣ�ռ�

@param
@����ֵ
*/
void FreeRStaInfo(RStaInfo* pRStaInfo);

/**
@brief �Ի��߽�����Ϣ����ռ�

@param
@����ֵ
*/
bool InitBLineSlon(int n, BLineSlon*& gBlineSlon);
/**
@brief �ͷŻ��߽�����Ϣ�ռ�

@param
@����ֵ
*/
void FreeBLineSlon(int n, BLineSlon*& gBlineSlon);

/**
@brief ����������Ϣ���г�ʼ��
@param
@����ֵ
*/
void InitTriNetInfo(int n, TriNetInfo*&  gTrinetInfo);
#endif // !COMMONFUNCTION_H


