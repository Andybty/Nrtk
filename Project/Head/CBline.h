#ifndef CBLINE_H
#define CBLINE_H
#include "CommonVariables.h"
#include "CommonFunction.h"
#include "CCoordTrans.h"
#include "MacroDefine.h"
#include "CMath.h"
#define ARMODE_OFF       0                                                       /* AR mode: off */
#define ARMODE_CONT      1                                                       /* AR mode: continuous */
#define ARMODE_INST      2                                                       /* AR mode: instantaneous */
#define ARMODE_FIXHOLD   3                                                       /* AR mode: fix and hold */
#define ARMODE_WLNL      4                                                       /* AR mode: wide lane/narrow lane */
#define ARMODE_TCAR      5                                                       /* AR mode: triple carrier ar */
#define VAR_HOLDAMB      0.001                                                   /* constraint to hold ambiguity (cycle^2) */
#define SOLQ_FIX         1                                                       /* solution status: fix */
/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(rtkopt)     ((rtkopt)->ionoopt==IONOOPT_IFLC?1:(rtkopt)->nf)
#define NP(rtkopt)     ((rtkopt)->dynamics==0?0:0)
#define NI(rtkopt)     ((rtkopt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(rtkopt)     ((rtkopt)->tropopt<TROPOPT_EST?0:((rtkopt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(rtkopt)     ((rtkopt)->glomodear!=2?0:NFREQGLO)
#define NB(rtkopt)     ((rtkopt)->mode<=PMODE_DGPS?0:MAXSAT*NF(rtkopt))
#define NR(rtkopt)     (NP(rtkopt)+NI(rtkopt)+NT(rtkopt)+NL(rtkopt))
#define NX(rtkopt)     (NR(rtkopt)+NB(rtkopt))

#define II(s,rtkopt)   (NP(rtkopt)+(s)-1)                                       /* ionos (s:satellite no) */
#define IT(r,rtkopt)   (NP(rtkopt)+NI(rtkopt)+NT(rtkopt)/2*(r))                 /* tropos (r:0=rov,1:ref) */
#define IL(f,rtkopt)   (NP(rtkopt)+NI(rtkopt)+NT(rtkopt)+(f))                   /* receiver h/w bias */
#define IB(s,f,rtkopt) (NR(rtkopt)+MAXSAT*(f)+(s)-1)                            /* phase bias (s:satno,f:freq) */




/*所有基线解算参数结构体*/
typedef struct 
tagStaRrkSolveParam
{
	int n                      ;                                             /*number of baseline*/
	RtkSolveParam* rtksolpara  ;                                             /*baseline solve struct*/
}StaRtkSolPara;


/*
 @@brief CBline类，主要对基线进行解算；
  @brief  基线解算；
  @brief  选择共视卫星；
  @brief  状态参数及其方差更新；
  @brief  卡尔曼滤波；
  @brief  计算双差残差和系数矩阵权矩阵；
  @brief  模糊度固定；
  @brief  计算双差大气（电离层和对流层或综合误差）
*/
class CBline
{
public:
	CBline();
	~CBline();
public:

	/**
	@brief  基线解算
	@param  无
	@返回值 无
	*/
	bool BLineSolution();
public:
	/**
	@brief  选择基线共视卫星
	@param  StaSolveInfo[in]观测的卫星信息
	@param  StaIndex1/StaIndex2[in]基站索引信息
	@param  BLineInfo[in]基线索引信息
	@返回值 基线共视卫星数量
	*/
	int SelCommonSat(StaSolveInfo* pStaSolveInfo, BLineInfo bLineinfo, uint32& StaIndex1, uint32& StaIndex2);

	/**
	@brief  状态参数及其方差更新
	@param  ObsData[in]观测数据
	@param  tt[in]当前历元和前历元的时间差
	@param  ns[in]观测数据有效卫星数
	@param  sat[in]观测数据卫星ID
	@param  index[in]观测数据的基站数
	@返回值 无
	*/
	void UpdateState(ObsData* pSta1ObsData, ObsData* pSta2ObsData, double tt, int ns, int* sat,
		int* iu, int* ir,int index);

	/**
	@brief  卡尔曼滤波
	@param  无
	@返回值 无
	*/
	int Filter();

	/**
	@brief  计算双差残差、系数矩阵和权矩阵
	@param  ObsData[in]观测数据
	@param  x[in/out]前后滤波状态参数
	@param  P[in/out]前后滤波方差参数
	@param  xa[in/out]模糊度固定状态参数
	@param  Pa[in/out]模糊度固定方差参数
	@param  y[in/out]单差残差
	@param  v[in/out]双差残差
	@param  azel[in/out]卫星高度角方位角
	@param  freq[in/out]卫星观测频点
	@param  StaIndex1[in]基站索引
	@返回值
	*/
	int CalDouDiffRes(ObsData* pSta1ObsData, ObsData* pSta2ObsData,double dt, const double* x,
		const double* P, const int* sat, double* y, double* e,
		double* azel, double* freq, const int* iu, const int* ir,
		int ns, double* v, double* H, double* R, int* vflg,int index, uint32 StaIndex1, uint32  StaIndex2);

	/**
	@brief  用lambda方法固定模糊度
	@param  bias[in]相位偏差
	@param  xa[in]浮点解转固定解矩阵
	@返回值 固定解个数
	*/
	int AmbiFixByLambda(double* bias, double* xa,int index);

	/**
	@brief  计算双差大气（电离层和对流层或综合误差）
	@param  ObsData[in]基站观测数据
	@param  sat[in]观测卫星ID
	@param  Index1[in]基站数量
	@param  StaIndex1[in]基站索引
	@返回值 观测卫星数量
	*/
	int CalDouDiffAtmo(ObsData* pSta1ObsData, ObsData* pSta2ObsData, gtime_t  dt, const double* x,
		const int* sat, const int* iu, const int* ir, int ns, int index, uint32 StaIndex1, uint32  StaIndex2);
	/**
	@brief  对单条基线进行处理
	@param  SinStaObsData[in]基站观测数据
	@param  BLineInfo[in]基线信息
	@param  index[in]基站数量
	@返回值 1
	*/
	int Process(SinStaObsData* pSinstaobsdata, BLineInfo bLineinfo,int index);

	/**
	@brief  配置基线解算策略信息
	@param  RtkOpt[in]RTK基线解算策略信息
	@返回值 无
	*/
	void InitRtkopt(RtkOpt rtkopt);
	/**
	@brief  初始化基线解算参数
	@param  n[in]初始化基线数量
	@返回值 true
	*/
	bool InitRtkParam(int n);
	
private:
	/**
	@brief  设置基线解算策略
	@param  nf[in]基线解算频点数量策略
	@param  ionoopt[in]基线电离层模型解算策略
	@param  tropopt[in]基线对流层模型解算策略
	@param  dynamics[in]基线动力模型解算策略
	@param  glomodear[in]基线GLO模糊度解算策略
	@param  navsys[in]基线卫星系统解算策略
	@返回值 无
	*/
	void SetRtkOpt(int nf, int mode, int ionoopt, int tropopt, int dynamics, int glomodear, int navsys);
	/**
	@brief 更新模糊度参数
	@param pSta1ObsData[in]观测数据
	@param  tt[in]当前历元和前历元的时间差
	@param  ns[in]观测数据有效卫星数
	@param  sat[in]观测数据卫星ID
	@param  index[in]观测数据的基站数
	@返回值 无
	*/
	void UpBias(ObsData* pSta1ObsData, ObsData* pSta2ObsData, double tt, int ns, int* sat,
		int* iu, int* ir, int index);
	/**
	@brief  初始化参数
	@param  RtkSolveParam[in]基线解算参数信息
	@param  xi[in]基线解算参数初始状态
	@param  var[in]基线解算参数初始方差
	@param  i[in]基线解算参数数量
	@返回值 无
	*/
	void InitParam(RtkSolveParam* rtksolveparam, double xi, double var, int i);

	/**
	@brief 单差模糊度初值
	@param ObsData[in]观测数据
	@param i[in]观测到第i颗卫星
	@param j[in]观测到的第j个观测值
	@param k[i]观测频点数量
	@返回值 载波单差初始值
	*/
	double SDobs(ObsData* obs1, ObsData* obs2, int i, int j, int k);

	/**
	@brief  计算基线长度
	@param  ru[in]基线端基站坐标
	@param  rb[in]基线另一端基站坐标
	@param  dr[out]基线长度
	@返回值 解算基线长度
	*/
	double Baseline(const double* ru, const double* rb, double* dr);

	/**
	@brief 设置观测值的噪声
	@param sat[in]观测卫星数量
	@param sys[in]观测卫星系统
	@返回值 伪距载波观测噪声
	*/
	double Varerr(int sat, int sys, double el, double bl, double dt, int f,
		RtkOpt* opt);
	/**
	@brief 求双差观测噪声矩阵
	@param Ri/Rj/R[out]双差观测噪声矩阵
	@返回值 无
	*/
	void DDcov(int* nb, int n, double* Ri, double* Rj, int nv, double* R);

	/**
	@brief 单差到双差转换矩阵
	@brief index[in]卫星号索引
	@param ix[out]参考偏差状态矩阵
	@返回值 单差到双差转换矩阵
	*/
	int Ddidx(int* ix,int index);

	/**
	@brief lambda法固定矩阵
	@param QF[out] lambda法模糊度固定矩阵
	@返回值 模糊度固定信息
	*/
	int Lambda(int n, int m, const double* a, const double* Q, double* F,
		double* s);
	/**
	@brief lambda法降相关
	@param LD[int]LD分解因子
	@返回值 无
	*/
	void Reduction(int n, double* L, double* D, double* Z);

	/**
	@brief 排列组合函数
	@param LD[in]LD矩阵排列组合
	@返回值 无
	*/
	void Permutation(int n, double* L, double* D, int j, double del, double* Z);
	/**
	@brief  高斯函数
	@param  LZ[in]高斯变换矩阵
	@返回值 无
	*/

	void Gauss(int n, double* L, double* Z, int i, int j);

	
    /**
	@brief modified lambda(mlambda) search
	@param  修改lambda搜索空间
	@返回值 成功返回0，否则返回-1
	*/
	int Search(int n, int m, const double* L, const double* D,
		const double* zs, double* zn, double* s);

	/**
	@brief 解算矩阵
	@param
	@返回值 
	*/
	int Solve(const char* tr, const double* A, const double* Y, int n,
		int m, double* X);
	/**
	@brief 恢复单差模糊度
	@param  bias[in]相位偏差
	@param  xa[in]双差转单差矩阵
	@param  index[in]卫星索引
	@返回值 无
	*/
	void Restamb(double* bias, int nb, double* xa,int index);

	/**
	@brief 使用hold模糊度方法
	@param xa[in]浮点解转固定解矩阵
	@param  index[in]卫星索引
	@返回值 无
	*/
	void Holdamb(double* xa,int index);

	/**
	@brief 判断系统是否符合要解算的系统
	@param sys[in]卫星系统
	@param m[out]卫星系统对应的标识
	@返回值 返回卫星系统标识
	*/
	int Test_sys(int sys, int m);
	/**
	@brief 保存基线解算参数
	@param stat[in]基线解算状态信息
	@param index[in]基线索引信息
	@返回值 无
	*/
	void SaveSoluStaus(int stat,int index);
	/**
	@brief 简算解算状态的有效性
	@param v[in]解算双差残差
	@返回值 返回1有效，0无效
	*/
	int Valpos(double* v, const double* R,int* vflg,int nv, double thres);
	/**
	@brief 保存每颗卫星的解算状态信息，包括存取历史信息
	@param ObsData[in]观测卫星信息
	@param n1/n2[in]观测卫星数
	@param index[in]观测卫星索引
	@返回值 无
	*/
	void SaveSsat(ObsData* obs1, ObsData* obs2, int n1, int n2, int stat,int index);
public:
	StaSol      m_Stasol         ;             /**所有基线的解算状态信息*/
private:
	int m_naCommonSat[MAXSAT]    ;             /*基线共视卫星编号*/
	int m_naSta1SatID[MAXSAT]    ;             /*起始基线编号对应站的共视卫星索引*/
	int m_naSta2SatID[MAXSAT]    ;             /*终点基线编号对应站的共视卫星索引*/
	int m_numcommonsat           ;             /*单条共视卫星数*/
	RtkOpt m_Rtkopt;
	ambc_t m_ambc[MAXSAT]          ;           /*双差模糊度信息*/
private:
	StaSatSolve m_Stasatsolve      ;           /*多条基线对应的卫星解算状态信息*/
	StaRtkSolPara m_Startksolpara  ;           /*所有基线解算参数结构体*/
	
};

#endif // !CBLINE_H


