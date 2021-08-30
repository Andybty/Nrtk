#ifndef CVRS_H
#define CVRS_H
#include "dpi_types_basic.h"
#include "CommonVariables.h"
#include "CommonFunction.h"
#include "CCoordTrans.h"
#include "CMath.h"
#include "CBroadcastEph.h"
#include "CStation.h"
/*����۲�ֵ����ֵ*/
typedef struct tagVrsCorrect
{
	double dDTro;/*������*/
	double dDIon;/*�����*/
	double dDPComErr;/*α���ۺ����*/
	double dDCComErr;/*��λ�ۺ����*/
	int CodeType;/*�۲�����*/
	int ref, rov;/*�ο��Ǻ�������*/
	gtime_t  time;/*�ο�ʱ��*/
}VrsCorrect;

/*
@@brief ����վ�࣬��Ҫ������������۲�ֵ��
*/
class CVRS
{
public:
	CVRS();
	~CVRS();
public:
	/**
	@brief  ѡ�����ο�վ
	@param index[in] �������������
	@����ֵ ���ο�վ�ı�����
	*/
	int SelMastReferSta(int index, int& rover1, int& rover2);

	/**
	@brief  ����VRS�۲�ֵ
	@param  Sinstaobsdata[in]        ��վ�۲�����
	@param  master[in] ��վ����
	@param  n[in] �۲�ֵ���͸������Ƿ�Ƶ�ʵģ�
	@����ֵ
	*/
	void GenVRSObs(SinStaObsData Sinstaobsdata, VRSInfo& vrsInfo, int master, int n);

	/**
	@brief  ����VRS�Ĵ�������
	@param  pSinstaobsdata[in]  ��վ�۲����ݽṹ��ָ��
	@param  pTrinetInfo[in]     ��������Ϣ�ṹ��ָ��
	@param  pVrsInfo[in]        ����վ��ṹ��ָ��
	@param  nVrs[in]            ����վ��۲�ֵ����
	@����ֵ
	*/
	void Process(SinStaObsData* pSinstaobsdata, TriNetInfo* pTrinetInfo, VRSInfo* pVrsInfo, int nVrs);
private:
	/**
	@brief  �����վ������վ�ľ���
	@param  sys[in] ����ϵͳ
	@param  dpSta[in] ��վ����
	@param  dpVrs[in] ����վ����
	@����ֵ ��վ������վ����
	*/
	double CalStaVrsDis(double* dpSta, double* dpVrs);
	/**
	@brief  ѡȡ��վ������վ����С����ֵ����
	@param  dpDis[in] ����ϵͳ
	@param  n[in] �������
	@����ֵ ��С������������
	*/
	int SelectMinIndex(double* dpDis, int n);

	/**
	@brief  ȷ����������������
	@param index[in] ����������
	@����ֵ
	*/
	int SelCommonSat(int index);

	/**
	@brief  ���ݻ�վ����ȷ����������
	@param master[in]   ��վ����
	@param bline1[out]  ����1����
	@param bline2[out]  ����2����
	@param flag[out]    ���߷�������
	@����ֵ
	*/
	void SelBlineIndex(int master, int& bline1, int& bline2, int& flag);
	/**
	@brief  ���ݻ�վ����ȷ����������
	@param master[in]   ��վ����
	@param bindex1[in]  ����1����
	@param bindex2[in] ����2����
	@param flag[in]    ���߷�������
	@param rov1[in]    ����վ1����
	@param rov2[in]    ����վ2����
	@����ֵ
	*/
	double* LCM(int master, int rov1, int rov2, VRSInfo& vrsInfo);
	/**
		@brief  ���ݻ�վ����ȷ����������
		@param coeff[in]   LCMϵ������
		@param bindex1[in]  ����1����
		@param bindex2[in] ����2����
		@param flag[in]    ���߷�������
		@����ֵ �۲����Ͷ�Ӧ�Ĵ�������ֵ����
		*/
	int GenVrsAtom(double* coeff, int bindex1, int bindex2, int flag);
	/**
	@brief  ������۲�ֵ����ֵ���г�ʼ��
	@param   n[in]  �ռ�����С
	@����ֵ ����ɹ�����TRUE
	*/
	bool InitVrsCorrect(int n);
	/**
	@brief  �ж��������ߵĹ��������Ƿ�һ��
	@param   ref1[in]  ��һ�����߲ο���
	@param   rov1[in]  ��һ������������
	@param   ref2[in]  �ڶ������߲ο���
	@param   rov2[in]  �ڶ�������������
	@����ֵ ���һ�·���TRUE
	*/
	bool JudgeConsist(int ref1, int rov1, int ref2, int rov2);
private:
	VrsCorrect* m_pVrscorrect=NULL;                        /*����۲�ֵ����ֵ�ṹ��ָ��*/
	CStation      m_cStation;                             /*��վ��*/

};





#endif // !CVRS_H



