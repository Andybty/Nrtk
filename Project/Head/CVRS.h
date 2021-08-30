#ifndef CVRS_H
#define CVRS_H
#include "dpi_types_basic.h"
#include "CommonVariables.h"
#include "CommonFunction.h"
#include "CCoordTrans.h"
#include "CMath.h"
#include "CBroadcastEph.h"
#include "CStation.h"
/*虚拟观测值改正值*/
typedef struct tagVrsCorrect
{
	double dDTro;/*对流层*/
	double dDIon;/*电离层*/
	double dDPComErr;/*伪距综合误差*/
	double dDCComErr;/*相位综合误差*/
	int CodeType;/*观测类型*/
	int ref, rov;/*参考星和流动星*/
	gtime_t  time;/*参考时刻*/
}VrsCorrect;

/*
@@brief 虚拟站类，主要功能生成虚拟观测值；
*/
class CVRS
{
public:
	CVRS();
	~CVRS();
public:
	/**
	@brief  选择主参考站
	@param index[in] 三角网索引编号
	@返回值 主参考站的编译编号
	*/
	int SelMastReferSta(int index, int& rover1, int& rover2);

	/**
	@brief  生成VRS观测值
	@param  Sinstaobsdata[in]        主站观测数据
	@param  master[in] 主站索引
	@param  n[in] 观测值类型个数（是分频率的）
	@返回值
	*/
	void GenVRSObs(SinStaObsData Sinstaobsdata, VRSInfo& vrsInfo, int master, int n);

	/**
	@brief  生成VRS的处理流程
	@param  pSinstaobsdata[in]  单站观测数据结构体指针
	@param  pTrinetInfo[in]     三角网信息结构体指针
	@param  pVrsInfo[in]        虚拟站点结构体指针
	@param  nVrs[in]            虚拟站点观测值个数
	@返回值
	*/
	void Process(SinStaObsData* pSinstaobsdata, TriNetInfo* pTrinetInfo, VRSInfo* pVrsInfo, int nVrs);
private:
	/**
	@brief  计算基站和虚拟站的距离
	@param  sys[in] 卫星系统
	@param  dpSta[in] 基站坐标
	@param  dpVrs[in] 虚拟站坐标
	@返回值 基站与虚拟站距离
	*/
	double CalStaVrsDis(double* dpSta, double* dpVrs);
	/**
	@brief  选取基站和虚拟站的最小距离值索引
	@param  dpDis[in] 卫星系统
	@param  n[in] 距离个数
	@返回值 最小距离的索引编号
	*/
	int SelectMinIndex(double* dpDis, int n);

	/**
	@brief  确定三角网共视卫星
	@param index[in] 三角网索引
	@返回值
	*/
	int SelCommonSat(int index);

	/**
	@brief  根据基站索引确定基线索引
	@param master[in]   主站索引
	@param bline1[out]  基线1索引
	@param bline2[out]  基线2索引
	@param flag[out]    基线方向索引
	@返回值
	*/
	void SelBlineIndex(int master, int& bline1, int& bline2, int& flag);
	/**
	@brief  根据基站索引确定基线索引
	@param master[in]   主站索引
	@param bindex1[in]  基线1索引
	@param bindex2[in] 基线2索引
	@param flag[in]    基线方向索引
	@param rov1[in]    流动站1索引
	@param rov2[in]    流动站2索引
	@返回值
	*/
	double* LCM(int master, int rov1, int rov2, VRSInfo& vrsInfo);
	/**
		@brief  根据基站索引确定基线索引
		@param coeff[in]   LCM系数矩阵
		@param bindex1[in]  基线1索引
		@param bindex2[in] 基线2索引
		@param flag[in]    基线方向索引
		@返回值 观测类型对应的大气改正值个数
		*/
	int GenVrsAtom(double* coeff, int bindex1, int bindex2, int flag);
	/**
	@brief  对虚拟观测值改正值进行初始化
	@param   n[in]  空间分配大小
	@返回值 分配成功返回TRUE
	*/
	bool InitVrsCorrect(int n);
	/**
	@brief  判断两条基线的共视卫星是否一致
	@param   ref1[in]  第一条基线参考星
	@param   rov1[in]  第一条基线流动星
	@param   ref2[in]  第二条基线参考星
	@param   rov2[in]  第二条基线流动星
	@返回值 如果一致返回TRUE
	*/
	bool JudgeConsist(int ref1, int rov1, int ref2, int rov2);
private:
	VrsCorrect* m_pVrscorrect=NULL;                        /*虚拟观测值改正值结构体指针*/
	CStation      m_cStation;                             /*测站类*/

};





#endif // !CVRS_H



