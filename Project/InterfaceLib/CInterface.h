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
@@brief 接口类，通过RStaDataInput接口函数输入基站观测数据流，VStaOutput接口函数输出虚拟站观测数据流
*/
class CInterface
{
public:
	CInterface();
	~CInterface();
public:
	/**
	@brief  输入基站原始观测数据流
	@param  ui8StaNum  [in] 输入的基站个数
	@param  pStaBuffer [in] 输入的基站数据流指针
	@notes  持续调用
	*/
	bool RStaDataInput(uint8 ui8StaNum, uint8* pStaBuffer);

	/**
	@brief  输出虚拟站观测数据流
	@param  ui8VRSNum [out] 输出的虚拟站个数
	@param  pVRSButter [out] 输出的虚拟站数据流指针
	@notes  持续调用
	*/
	bool VRSDataOutput(uint8 ui8VRSNum, uint8* pVRSButter);

};

#endif