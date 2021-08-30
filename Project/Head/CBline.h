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




/*���л��߽�������ṹ��*/
typedef struct 
tagStaRrkSolveParam
{
	int n                      ;                                             /*number of baseline*/
	RtkSolveParam* rtksolpara  ;                                             /*baseline solve struct*/
}StaRtkSolPara;


/*
 @@brief CBline�࣬��Ҫ�Ի��߽��н��㣻
  @brief  ���߽��㣻
  @brief  ѡ�������ǣ�
  @brief  ״̬�������䷽����£�
  @brief  �������˲���
  @brief  ����˫��в��ϵ������Ȩ����
  @brief  ģ���ȹ̶���
  @brief  ����˫������������Ͷ�������ۺ���
*/
class CBline
{
public:
	CBline();
	~CBline();
public:

	/**
	@brief  ���߽���
	@param  ��
	@����ֵ ��
	*/
	bool BLineSolution();
public:
	/**
	@brief  ѡ����߹�������
	@param  StaSolveInfo[in]�۲��������Ϣ
	@param  StaIndex1/StaIndex2[in]��վ������Ϣ
	@param  BLineInfo[in]����������Ϣ
	@����ֵ ���߹�����������
	*/
	int SelCommonSat(StaSolveInfo* pStaSolveInfo, BLineInfo bLineinfo, uint32& StaIndex1, uint32& StaIndex2);

	/**
	@brief  ״̬�������䷽�����
	@param  ObsData[in]�۲�����
	@param  tt[in]��ǰ��Ԫ��ǰ��Ԫ��ʱ���
	@param  ns[in]�۲�������Ч������
	@param  sat[in]�۲���������ID
	@param  index[in]�۲����ݵĻ�վ��
	@����ֵ ��
	*/
	void UpdateState(ObsData* pSta1ObsData, ObsData* pSta2ObsData, double tt, int ns, int* sat,
		int* iu, int* ir,int index);

	/**
	@brief  �������˲�
	@param  ��
	@����ֵ ��
	*/
	int Filter();

	/**
	@brief  ����˫��вϵ�������Ȩ����
	@param  ObsData[in]�۲�����
	@param  x[in/out]ǰ���˲�״̬����
	@param  P[in/out]ǰ���˲��������
	@param  xa[in/out]ģ���ȹ̶�״̬����
	@param  Pa[in/out]ģ���ȹ̶��������
	@param  y[in/out]����в�
	@param  v[in/out]˫��в�
	@param  azel[in/out]���Ǹ߶ȽǷ�λ��
	@param  freq[in/out]���ǹ۲�Ƶ��
	@param  StaIndex1[in]��վ����
	@����ֵ
	*/
	int CalDouDiffRes(ObsData* pSta1ObsData, ObsData* pSta2ObsData,double dt, const double* x,
		const double* P, const int* sat, double* y, double* e,
		double* azel, double* freq, const int* iu, const int* ir,
		int ns, double* v, double* H, double* R, int* vflg,int index, uint32 StaIndex1, uint32  StaIndex2);

	/**
	@brief  ��lambda�����̶�ģ����
	@param  bias[in]��λƫ��
	@param  xa[in]�����ת�̶������
	@����ֵ �̶������
	*/
	int AmbiFixByLambda(double* bias, double* xa,int index);

	/**
	@brief  ����˫������������Ͷ�������ۺ���
	@param  ObsData[in]��վ�۲�����
	@param  sat[in]�۲�����ID
	@param  Index1[in]��վ����
	@param  StaIndex1[in]��վ����
	@����ֵ �۲���������
	*/
	int CalDouDiffAtmo(ObsData* pSta1ObsData, ObsData* pSta2ObsData, gtime_t  dt, const double* x,
		const int* sat, const int* iu, const int* ir, int ns, int index, uint32 StaIndex1, uint32  StaIndex2);
	/**
	@brief  �Ե������߽��д���
	@param  SinStaObsData[in]��վ�۲�����
	@param  BLineInfo[in]������Ϣ
	@param  index[in]��վ����
	@����ֵ 1
	*/
	int Process(SinStaObsData* pSinstaobsdata, BLineInfo bLineinfo,int index);

	/**
	@brief  ���û��߽��������Ϣ
	@param  RtkOpt[in]RTK���߽��������Ϣ
	@����ֵ ��
	*/
	void InitRtkopt(RtkOpt rtkopt);
	/**
	@brief  ��ʼ�����߽������
	@param  n[in]��ʼ����������
	@����ֵ true
	*/
	bool InitRtkParam(int n);
	
private:
	/**
	@brief  ���û��߽������
	@param  nf[in]���߽���Ƶ����������
	@param  ionoopt[in]���ߵ����ģ�ͽ������
	@param  tropopt[in]���߶�����ģ�ͽ������
	@param  dynamics[in]���߶���ģ�ͽ������
	@param  glomodear[in]����GLOģ���Ƚ������
	@param  navsys[in]��������ϵͳ�������
	@����ֵ ��
	*/
	void SetRtkOpt(int nf, int mode, int ionoopt, int tropopt, int dynamics, int glomodear, int navsys);
	/**
	@brief ����ģ���Ȳ���
	@param pSta1ObsData[in]�۲�����
	@param  tt[in]��ǰ��Ԫ��ǰ��Ԫ��ʱ���
	@param  ns[in]�۲�������Ч������
	@param  sat[in]�۲���������ID
	@param  index[in]�۲����ݵĻ�վ��
	@����ֵ ��
	*/
	void UpBias(ObsData* pSta1ObsData, ObsData* pSta2ObsData, double tt, int ns, int* sat,
		int* iu, int* ir, int index);
	/**
	@brief  ��ʼ������
	@param  RtkSolveParam[in]���߽��������Ϣ
	@param  xi[in]���߽��������ʼ״̬
	@param  var[in]���߽��������ʼ����
	@param  i[in]���߽����������
	@����ֵ ��
	*/
	void InitParam(RtkSolveParam* rtksolveparam, double xi, double var, int i);

	/**
	@brief ����ģ���ȳ�ֵ
	@param ObsData[in]�۲�����
	@param i[in]�۲⵽��i������
	@param j[in]�۲⵽�ĵ�j���۲�ֵ
	@param k[i]�۲�Ƶ������
	@����ֵ �ز������ʼֵ
	*/
	double SDobs(ObsData* obs1, ObsData* obs2, int i, int j, int k);

	/**
	@brief  ������߳���
	@param  ru[in]���߶˻�վ����
	@param  rb[in]������һ�˻�վ����
	@param  dr[out]���߳���
	@����ֵ ������߳���
	*/
	double Baseline(const double* ru, const double* rb, double* dr);

	/**
	@brief ���ù۲�ֵ������
	@param sat[in]�۲���������
	@param sys[in]�۲�����ϵͳ
	@����ֵ α���ز��۲�����
	*/
	double Varerr(int sat, int sys, double el, double bl, double dt, int f,
		RtkOpt* opt);
	/**
	@brief ��˫��۲���������
	@param Ri/Rj/R[out]˫��۲���������
	@����ֵ ��
	*/
	void DDcov(int* nb, int n, double* Ri, double* Rj, int nv, double* R);

	/**
	@brief ���˫��ת������
	@brief index[in]���Ǻ�����
	@param ix[out]�ο�ƫ��״̬����
	@����ֵ ���˫��ת������
	*/
	int Ddidx(int* ix,int index);

	/**
	@brief lambda���̶�����
	@param QF[out] lambda��ģ���ȹ̶�����
	@����ֵ ģ���ȹ̶���Ϣ
	*/
	int Lambda(int n, int m, const double* a, const double* Q, double* F,
		double* s);
	/**
	@brief lambda�������
	@param LD[int]LD�ֽ�����
	@����ֵ ��
	*/
	void Reduction(int n, double* L, double* D, double* Z);

	/**
	@brief ������Ϻ���
	@param LD[in]LD�����������
	@����ֵ ��
	*/
	void Permutation(int n, double* L, double* D, int j, double del, double* Z);
	/**
	@brief  ��˹����
	@param  LZ[in]��˹�任����
	@����ֵ ��
	*/

	void Gauss(int n, double* L, double* Z, int i, int j);

	
    /**
	@brief modified lambda(mlambda) search
	@param  �޸�lambda�����ռ�
	@����ֵ �ɹ�����0�����򷵻�-1
	*/
	int Search(int n, int m, const double* L, const double* D,
		const double* zs, double* zn, double* s);

	/**
	@brief �������
	@param
	@����ֵ 
	*/
	int Solve(const char* tr, const double* A, const double* Y, int n,
		int m, double* X);
	/**
	@brief �ָ�����ģ����
	@param  bias[in]��λƫ��
	@param  xa[in]˫��ת�������
	@param  index[in]��������
	@����ֵ ��
	*/
	void Restamb(double* bias, int nb, double* xa,int index);

	/**
	@brief ʹ��holdģ���ȷ���
	@param xa[in]�����ת�̶������
	@param  index[in]��������
	@����ֵ ��
	*/
	void Holdamb(double* xa,int index);

	/**
	@brief �ж�ϵͳ�Ƿ����Ҫ�����ϵͳ
	@param sys[in]����ϵͳ
	@param m[out]����ϵͳ��Ӧ�ı�ʶ
	@����ֵ ��������ϵͳ��ʶ
	*/
	int Test_sys(int sys, int m);
	/**
	@brief ������߽������
	@param stat[in]���߽���״̬��Ϣ
	@param index[in]����������Ϣ
	@����ֵ ��
	*/
	void SaveSoluStaus(int stat,int index);
	/**
	@brief �������״̬����Ч��
	@param v[in]����˫��в�
	@����ֵ ����1��Ч��0��Ч
	*/
	int Valpos(double* v, const double* R,int* vflg,int nv, double thres);
	/**
	@brief ����ÿ�����ǵĽ���״̬��Ϣ��������ȡ��ʷ��Ϣ
	@param ObsData[in]�۲�������Ϣ
	@param n1/n2[in]�۲�������
	@param index[in]�۲���������
	@����ֵ ��
	*/
	void SaveSsat(ObsData* obs1, ObsData* obs2, int n1, int n2, int stat,int index);
public:
	StaSol      m_Stasol         ;             /**���л��ߵĽ���״̬��Ϣ*/
private:
	int m_naCommonSat[MAXSAT]    ;             /*���߹������Ǳ��*/
	int m_naSta1SatID[MAXSAT]    ;             /*��ʼ���߱�Ŷ�Ӧվ�Ĺ�����������*/
	int m_naSta2SatID[MAXSAT]    ;             /*�յ���߱�Ŷ�Ӧվ�Ĺ�����������*/
	int m_numcommonsat           ;             /*��������������*/
	RtkOpt m_Rtkopt;
	ambc_t m_ambc[MAXSAT]          ;           /*˫��ģ������Ϣ*/
private:
	StaSatSolve m_Stasatsolve      ;           /*�������߶�Ӧ�����ǽ���״̬��Ϣ*/
	StaRtkSolPara m_Startksolpara  ;           /*���л��߽�������ṹ��*/
	
};

#endif // !CBLINE_H


