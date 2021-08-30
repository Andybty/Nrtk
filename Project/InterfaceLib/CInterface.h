#ifndef CINTERFACE_H
#define CINTERFACE_H
#include "dpi_types_basic.h"
#include "CRTCM.h"
#include "CStation.h"
#include "CBroadcastEph.h"
#include "MacroDefine.h"
#include "CCycleSlip.h"
#include "CBline.h"
#include "CVRS.h"
using namespace std;

/*
@@brief �ӿ��࣬ͨ��RStaDataInput�ӿں��������վ�۲���������VStaOutput�ӿں����������վ�۲�������
*/
class CInterface
{
public:
	CInterface();
	~CInterface();
public:
	/**
	@brief  �����վԭʼ�۲�������
	@param  ui8StaNum  [in] ����Ļ�վ����
	@param  pStaBuffer [in] ����Ļ�վ������ָ��
	@notes  ��������
	*/
	bool RStaDataInput(uint8 ui8StaNum, uint8* pStaBuffer);

	/**
	@brief  �������վ�۲�������
	@param  ui8VRSNum [out] ���������վ����
	@param  pVRSButter [out] ���������վ������ָ��
	@notes  ��������
	*/
	bool VRSDataOutput(uint8 ui8VRSNum, uint8* pVRSButter);

};

#endif