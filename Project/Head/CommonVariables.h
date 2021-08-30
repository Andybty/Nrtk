#ifndef COMMONVARIABLES_H
#define COMMONVARIABLES_H
#include "dpi_types_basic.h"
#include "MacroDefine.h"
using namespace std;

/*卫星解算状态信息(包括上个历元的历史信息)satellite status type*/
typedef struct tagSatSolve {
	uint8_t sys;                                                           /* navigation system */
	uint8_t vs;                                                            /* valid satellite flag single */
	double azel[2];                                                        /* azimuth/elevation angles {az,el} (rad) */
	double resp[NFREQ];                                                    /* residuals of pseudorange (m) */
	double resc[NFREQ];                                                    /* residuals of carrier-phase (m) */
	uint8_t vsat[NFREQ];                                                   /* valid satellite flag */
	uint16_t snr[NFREQ];                                                   /* signal strength (*SNR_UNIT dBHz) */
	uint8_t fix[NFREQ];                                                    /* ambiguity fix flag (1:fix,2:float,3:hold) */
	uint8_t slip[NFREQ];                                                   /* cycle-slip flag */
	uint8_t half[NFREQ];                                                   /* half-cycle valid flag */
	int lock[NFREQ];                                                       /* lock counter of phase */
	uint32_t outc[NFREQ];                                                  /* obs outage counter of phase */
	uint32_t slipc[NFREQ];                                                 /* cycle-slip counter */
	uint32_t rejc[NFREQ];                                                  /* reject counter */
	double gf[NFREQ - 1];                                                  /* geometry-free phase (m) */
	double mw[NFREQ - 1];                                                  /* MW-LC (m) */
	double phw;                                                            /* phase windup (cycle) */
	gtime_t pt[2][NFREQ];                                                  /* previous carrier-phase time */
	double ph[2][NFREQ];                                                   /* previous carrier-phase observable (cycle) */
	double satpos[3];                                                      /* satellite position*/
} ssat_t;

/*多条基线对应的卫星解算状态信息*/
typedef struct tagStaSatSolve
{
	int n;                                                                 /*number of baseline*/
	ssat_t** ssat;

}StaSatSolve;

/*模糊度解算状态信息结构体，*/
typedef struct tagAmb {                                                    /* ambiguity control type */
	gtime_t epoch[4];                                                      /* last epoch */
	int n[4];                                                              /* number of epochs */
	double LC[4];                                                          /* linear combination average */
	double LCv[4];                                                         /* linear combination variance */
	int fixcnt;                                                            /* fix count */
	char flags[MAXSAT];                                                    /* fix flags */
} ambc_t;


/*基线解算结果信息*/
typedef struct tagSol {                                                    /* solution type */
	gtime_t time;                                                          /* time (GPST) */
	double rr[6];                                                          /* position/velocity (m|m/s) */
																		   /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
	float  qr[6];                                                          /* position variance/covariance (m^2) */
																		   /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
																		   /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
	float  qv[6];                                                          /* velocity variance/covariance (m^2/s^2) */
	double dtr[6];                                                         /* receiver clock bias to time systems (s) */
	uint8_t type;                                                          /* type (0:xyz-ecef,1:enu-baseline) */
	uint8_t stat;                                                          /* solution status (SOLQ_???) */
	uint8_t ns;                                                            /* number of valid satellites */
	float age;                                                             /* age of differential (s) */
	float ratio;                                                           /* AR ratio factor for valiation */
	float thres;                                                           /* AR ratio threshold for valiation */
} sol_t;


/*所有基线的解算状态信息*/
typedef struct tagStaSol
{
	int n;                                                                 /*number of baseline*/
	sol_t* sol;

}StaSol;

/*基线解算参数信息*/
typedef struct tagRtkSolveParam {
	double rb[6];                                                          /* base position/velocity (ecef) (m|m/s) */
	double rr[6];                                                          /* rover position/velocity (ecef) (m|m/s) */
	int nx, na;                                                            /* number of float states/fixed states */
	double tt;                                                             /* time difference between current and previous (s) */
	double* x, * P;                                                        /* float states and their covariance */
	double* xa, * Pa;                                                      /* fixed states and their covariance */
	int nfix;                                                              /* number of continuous fixes of ambiguity */
} RtkSolveParam;
/*基线解算配置信息结构体*/
typedef struct tagRtkOpt {
	int nf;                                                               /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
	int mode;                                                             /* positioning mode (PMODE_???) */
	int ionoopt;                                                          /* ionosphere option (IONOOPT_???) */
	int tropopt;                                                          /* troposphere option (TROPOPT_???) */
	int dynamics;                                                         /* dynamics model (0:none,1:velociy,2:accel) */
	int glomodear;                                                        /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
	int bdsmodear;                                                        /* BeiDou AR mode (0:off,1:on) */
	int navsys;                                                           /* navigation system */
	int maxout;                                                           /* obs outage count to reset bias */
	int modear;                                                           /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
	int minlock;                                                          /* min lock count to fix ambiguity */
	int minfix;                                                           /* min fix count to hold ambiguity */
	double prn[6];                                                        /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv [5] pos */
	double std[3];                                                        /* initial-state std [0]bias,[1]iono [2]trop */
	double maxinno;                                                       /* reject threshold of innovation (m) */
	double err[5];                                                        /* measurement error factor */
																		  /* [0]:reserved */
																		  /* [1-3]:error factor a/b/c of phase (m) */
																		  /* [4]:doppler frequency (hz) */
	double sclkstab;                                                      /* satellite clock stability (sec/sec) */
	double eratio[NFREQ];                                                 /* code/phase error ratio */
	double elmaskar;                                                      /* elevation mask of AR for rising satellite (deg) */
	double thresar[8];                                                    /* AR validation threshold */
	double elmaskhold;                                                    /* elevation mask to hold ambiguity (deg) */
}RtkOpt;

/*GNSS时间系统结构体*/
typedef struct tagTimeGnss 
{
	uint16 ui16GpsWeek          ;                                          /*GPSweek*/
	uint32 ui32GpsIntSec        ;                                          /*GPSweek int seconde*/
	double dGpsFracSec          ;                                          /*GPSweek Fractional seconde*/
}TimeGnss;

/*GPS/BDS/GAL广播星历结构体*/
typedef struct  tagBroadcastEphData{      
	int sat                    ;                                          /* satellite number */
	int iode, iodc             ;                                          /* IODE,IODC */
	int sva                    ;                                          /* SV accuracy (URA index) */
	int svh                    ;                                          /* SV health (0:ok) */
	int week                   ;                                          /* GPS/QZS: gps week, GAL: galileo week */
	int code                   ;                                          /* GPS/QZS: code on L2 */
						                                                  /* GAL: data source defined as rinex 3.03 */
						                                                  /* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
	int flag                   ;                                          /* GPS/QZS: L2 P data flag */
						                                                  /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
	gtime_t toe, toc, ttr      ;                                          /* Toe,Toc,T_trans */
  	double A, e, i0, OMG0, omg, M0, deln, OMGd, idot;                     /* SV orbit parameters */
	double crc, crs, cuc, cus, cic, cis;
	double toes                ;                                          /* Toe (s) in week */
	double fit                 ;                                          /* fit interval (h) */
	double f0, f1, f2          ;                                          /* SV clock parameters (af0,af1,af2) */
	double tgd[6]              ;                                          /* group delay parameters */
						                                                  /* GPS/QZS:tgd[0]=TGD */
						                                                  /* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
						                                                  /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
						                                                  /* tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
	double Adot, ndot;                                                    /* Adot,ndot for CNAV */
} BroadcastEphData;
 
typedef struct tagAllBroadcastEphData {
	int nSatNum                         ;                                 /*number of satellite*/
	int MaxNum                          ;                                 /*number of space*/
	BroadcastEphData* pBroadcastEphData ;                                 /*BroadcastEphData struct*/
}AllBroadcastEphData;

typedef struct tagGpsBroadcastEphData {
	int nSatNum;                                                         /*number of satellite*/ 
	int MaxNum;                                                          /*/*number of space*/
	BroadcastEphData* pBroadcastEphData;                                 /*BroadcastEphData struct*/
}GpsBroadcastEphData;
/*观测数据结构*/
typedef struct tagObsData
{
	TimeStamp tTimeStamp                         ;                       /*观测时间*/
	uint8     ui8SatId                           ;                       /*卫星编号索引*/
	uint8     ui8FreqNum                         ;                       /*卫星频率数目*/
	uint8     ui8aCodeType[NFREQ + NEXOBS]       ;                       /*观测码类型*/
	double    daPseRange[NFREQ + NEXOBS]         ;                       /*伪距观测值*/
	double    daCarrPhase[NFREQ + NEXOBS]        ;                       /*载波相位观测值*/
	double    daDoppler[NFREQ + NEXOBS]          ;                       /*多普勒观测值*/
	double    daSigNoiRatio[NFREQ + NEXOBS]      ;                       /*信噪比*/
	double    daPseRanIF[3]                      ;                       /*伪距消电离层组合观测值*/
	double    daPhaseIF[3]                       ;                       /*相位消电离层组合观测值*/
	double    daMW[3]                            ;                       /*MW组合观测值*/
	double    daGF[3]                            ;                       /*GF组合观测值*/
	uint8     ui8aLLI[NFREQ + NEXOBS]            ;                       /* 卫星失锁计数*/
}ObsData;
/*单个测站观测数据结构体*/
typedef struct tagSinStaObsData
{
	uint16    ui16ObsNum         ;                                       /*观测历元个数*/
	uint16    ui16MaxNum         ;                                       /*分配的内存空间*/
	ObsData* pObsData            ;                                       /*观测数据结构体指针*/
}SinStaObsData;
/*卫星解算信息结构体*/
typedef struct  tagSatSlon
{
	uint8  uiObsindex               ;                                   /*在观测结构体中的索引*/
	uint8  ui8SatId                 ;                                   /*卫星号索引*/
	uint8  ui8SatSys                ;                                   /*卫星所属系统*/
	uint8  ui8SatHealSta            ;                                   /*卫星健康状态*/
	uint32 ui32SatTrackNum          ;                                   /*卫星跟踪数目*/
	double dSatClk                  ;                                   /*卫星钟差*/
	double dVar                     ;                                   /*卫星位置和钟差方差*/
	double daAzim[2]                ;                                   /*方位角高度角*/
	double daSigNoiRatio[NFREQ]     ;                                   /*信噪比*/
	double daMultip[6]              ;                                   /*多路径*/
	double daCycleSlipRate[NFREQ]   ;                                   /*周跳比*/
	double daDataInteRate           ;                                   /*数据完整率*/
	double dSta2SatDis              ;                                   /*站星距*/
	double daSatPos[3]              ;                                   /*卫星位置*/
	double daPseuRangeRes[NFREQ]    ;                                   /*伪距残差*/
	double daCarrPhaseRes[NFREQ]    ;                                   /*载波相位残差*/
	double e[3]                     ;                                   /*视线单位向量*/
}SatSoln;

/*基站静态信息结构体*/
typedef struct tagRStaInfo
{
	uint32 ui32StaId		    ;                                      /*基站ID标识*/
	string strName              ;                                      /*基站名称*/
	string strReceAntType	    ;                                      /*接收机天线类型*/
	uint8  ui8ReceAntCorrType   ;                                      /*接收机天线改正类型*/
	double daReveAntCorrVal[3]  ;                                      /*接收机天线改正值*/
	uint8  ui8CoorSysType       ;                                      /*基站坐标系统类型*/
	uint8  ui8CoorType          ;                                      /*基站坐标类型*/
	double daCoorVal[3]         ;                                      /*基站三维坐标值*/
}RStaInfo;


/*基站解算信息结构体*/
typedef struct tagStaSolveInfo
{

	uint32    ui32StaId          ;                                    /*测站ID标识*/
	//TimeGnss  tTimeGnss          ;                                    /*时间*/
	TimeStamp tTimeStamp         ;                                    /*观测时间*/
	uint8     ui8SatNum          ;                                    /*卫星数目*/
	double    dZeniTropDelayDry  ;                                    /*天顶对流层延迟干分量*/
	SatSoln*  pSatSlon           ;                                    /*观测到的卫星信息结构体指针*/
}StaSolveInfo;

/*共视卫星解算信息结构体*/
typedef struct tagComnSatSlon            
{
	uint8  ui8RefSatId          ;                                     /*参考星*/
	uint8  ui8RovSatId          ;                                     /*流动星*/
	uint8  ui8SolStatus         ;                                     /*解算状态*/
	uint32 ui32FixedNum         ;                                     /*累计固定次数*/
	double dTropDelay           ;                                     /*对流层延迟量*/
	double dIonoDelay           ;                                     /*电离层延迟量*/
	double dPCombinErr          ;                                     /*伪距综合误差*/
	double dCCombinErr          ;                                     /*相位综合误差*/
	int    CodeType             ;                                     /*观测值类型*/
}ComnSatSlon;
/*基线静态信息结构体*/
typedef struct tagBLineInfo
{
	uint32 ui32Id                ;                                   /*基线ID*/
	uint32 ui32StartPointId      ;                                   /*基线起点索引*/
	uint32 ui32EndPointId        ;                                   /*基线终点索引*/
}BLineInfo;

/*基线解算信息结构体*/
typedef struct tagBLineSlon               
{
	uint32     ui32Id            ;                                   /*基线ID*/
	TimeStamp  tTimeStamp        ;                                   /*观测时间*/
	double     Ddist             ;                                   /*基线长度*/
	uint8      ui8ComnSatNum     ;                                   /*共视卫星数目*/
	uint8      ui8EquNum         ;                                   /*伪距或相位双差方程个数*/
	uint8      ui8FixedSatNum    ;                                   /*模糊度固定的卫星数目*/
	ComnSatSlon* pComnSatSoln    ;                                   /*共视卫星结构体信息指针*/
}BLineSlon;

/*三角网信息结构体*/
typedef struct tagTriNetInfo
{
	uint32 ui32Id                ;                                   /*三角网ID*/
	uint32 ui32BLine1Id          ;                                   /*基线1索引*/
	uint32 ui32BLine2Id          ;                                   /*基线2索引*/
	uint32 UI32bLine3Id          ;                                   /*基线3索引*/
}TriNetInfo;

typedef struct tagSta {        /* station parameter type */
	char name[MAXANT]; /* marker name */
	char marker[MAXANT]; /* marker number */
	char antdes[MAXANT]; /* antenna descriptor */
	char antsno[MAXANT]; /* antenna serial number */
	char rectype[MAXANT]; /* receiver type descriptor */
	char recver[MAXANT]; /* receiver firmware version */
	char recsno[MAXANT]; /* receiver serial number */
	int antsetup;       /* antenna setup id */
	int itrf;           /* ITRF realization year */
	int deltype;        /* antenna delta type (0:enu,1:xyz) */
	double pos[3];      /* station position (ecef) (m) */
	double del[3];      /* antenna position delta (e/n/u or x/y/z) (m) */
	double hgt;         /* antenna height (m) */
} StaInfo;
/*虚拟站观测信息结构体*/
typedef struct tagVRSInfo               
{
	uint32  un32Id               ;                                   /*虚拟站ID*/
	uint8   un8CoorType          ;                                   /*虚拟站坐标类型*/
	double  daCoorVal[3]         ;                                   /*虚拟站坐标值*/
	uint32  ui32AttTriNerId      ;                                   /*依赖的三角网ID*/
	uint32* ui32AtthBLineId      ;                                   /*生成虚拟站所用的基线索引指针*/
	uint32  ui32MasterStaId      ;                                   /*主基站的索引*/
	int     nNum                 ;                                   /*观测历元个数*/
	int     nMax                 ;                                   /*最大观测历元个数*/
	ObsData* pObsData            ;                                   /*观测数据结构体指针*/
	StaInfo sta                 ;                                   /* station parameters */
}VRSInfo;

extern RStaInfo      *gRStaInfo        ;                             /*基站静态信息*/
extern StaSolveInfo  *gStaSolveInfo    ;                             /*基站解算信息*/
extern BLineInfo     *gBlineInfo       ;                             /*基线静态信息*/
extern BLineSlon     *gBlineSlon       ;                             /*基线解算信息*/
extern TriNetInfo    *gTrinetInfo      ;                             /*三角网信息*/
extern VRSInfo       *gVrsInfo         ;                             /*虚拟站观测信息*/


#endif // !COMMONVARIABLES_H

