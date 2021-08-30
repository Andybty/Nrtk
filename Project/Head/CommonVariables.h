#ifndef COMMONVARIABLES_H
#define COMMONVARIABLES_H
#include "dpi_types_basic.h"
#include "MacroDefine.h"
using namespace std;

/*���ǽ���״̬��Ϣ(�����ϸ���Ԫ����ʷ��Ϣ)satellite status type*/
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

/*�������߶�Ӧ�����ǽ���״̬��Ϣ*/
typedef struct tagStaSatSolve
{
	int n;                                                                 /*number of baseline*/
	ssat_t** ssat;

}StaSatSolve;

/*ģ���Ƚ���״̬��Ϣ�ṹ�壬*/
typedef struct tagAmb {                                                    /* ambiguity control type */
	gtime_t epoch[4];                                                      /* last epoch */
	int n[4];                                                              /* number of epochs */
	double LC[4];                                                          /* linear combination average */
	double LCv[4];                                                         /* linear combination variance */
	int fixcnt;                                                            /* fix count */
	char flags[MAXSAT];                                                    /* fix flags */
} ambc_t;


/*���߽�������Ϣ*/
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


/*���л��ߵĽ���״̬��Ϣ*/
typedef struct tagStaSol
{
	int n;                                                                 /*number of baseline*/
	sol_t* sol;

}StaSol;

/*���߽��������Ϣ*/
typedef struct tagRtkSolveParam {
	double rb[6];                                                          /* base position/velocity (ecef) (m|m/s) */
	double rr[6];                                                          /* rover position/velocity (ecef) (m|m/s) */
	int nx, na;                                                            /* number of float states/fixed states */
	double tt;                                                             /* time difference between current and previous (s) */
	double* x, * P;                                                        /* float states and their covariance */
	double* xa, * Pa;                                                      /* fixed states and their covariance */
	int nfix;                                                              /* number of continuous fixes of ambiguity */
} RtkSolveParam;
/*���߽���������Ϣ�ṹ��*/
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

/*GNSSʱ��ϵͳ�ṹ��*/
typedef struct tagTimeGnss 
{
	uint16 ui16GpsWeek          ;                                          /*GPSweek*/
	uint32 ui32GpsIntSec        ;                                          /*GPSweek int seconde*/
	double dGpsFracSec          ;                                          /*GPSweek Fractional seconde*/
}TimeGnss;

/*GPS/BDS/GAL�㲥�����ṹ��*/
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
/*�۲����ݽṹ*/
typedef struct tagObsData
{
	TimeStamp tTimeStamp                         ;                       /*�۲�ʱ��*/
	uint8     ui8SatId                           ;                       /*���Ǳ������*/
	uint8     ui8FreqNum                         ;                       /*����Ƶ����Ŀ*/
	uint8     ui8aCodeType[NFREQ + NEXOBS]       ;                       /*�۲�������*/
	double    daPseRange[NFREQ + NEXOBS]         ;                       /*α��۲�ֵ*/
	double    daCarrPhase[NFREQ + NEXOBS]        ;                       /*�ز���λ�۲�ֵ*/
	double    daDoppler[NFREQ + NEXOBS]          ;                       /*�����չ۲�ֵ*/
	double    daSigNoiRatio[NFREQ + NEXOBS]      ;                       /*�����*/
	double    daPseRanIF[3]                      ;                       /*α�����������Ϲ۲�ֵ*/
	double    daPhaseIF[3]                       ;                       /*��λ���������Ϲ۲�ֵ*/
	double    daMW[3]                            ;                       /*MW��Ϲ۲�ֵ*/
	double    daGF[3]                            ;                       /*GF��Ϲ۲�ֵ*/
	uint8     ui8aLLI[NFREQ + NEXOBS]            ;                       /* ����ʧ������*/
}ObsData;
/*������վ�۲����ݽṹ��*/
typedef struct tagSinStaObsData
{
	uint16    ui16ObsNum         ;                                       /*�۲���Ԫ����*/
	uint16    ui16MaxNum         ;                                       /*������ڴ�ռ�*/
	ObsData* pObsData            ;                                       /*�۲����ݽṹ��ָ��*/
}SinStaObsData;
/*���ǽ�����Ϣ�ṹ��*/
typedef struct  tagSatSlon
{
	uint8  uiObsindex               ;                                   /*�ڹ۲�ṹ���е�����*/
	uint8  ui8SatId                 ;                                   /*���Ǻ�����*/
	uint8  ui8SatSys                ;                                   /*��������ϵͳ*/
	uint8  ui8SatHealSta            ;                                   /*���ǽ���״̬*/
	uint32 ui32SatTrackNum          ;                                   /*���Ǹ�����Ŀ*/
	double dSatClk                  ;                                   /*�����Ӳ�*/
	double dVar                     ;                                   /*����λ�ú��Ӳ��*/
	double daAzim[2]                ;                                   /*��λ�Ǹ߶Ƚ�*/
	double daSigNoiRatio[NFREQ]     ;                                   /*�����*/
	double daMultip[6]              ;                                   /*��·��*/
	double daCycleSlipRate[NFREQ]   ;                                   /*������*/
	double daDataInteRate           ;                                   /*����������*/
	double dSta2SatDis              ;                                   /*վ�Ǿ�*/
	double daSatPos[3]              ;                                   /*����λ��*/
	double daPseuRangeRes[NFREQ]    ;                                   /*α��в�*/
	double daCarrPhaseRes[NFREQ]    ;                                   /*�ز���λ�в�*/
	double e[3]                     ;                                   /*���ߵ�λ����*/
}SatSoln;

/*��վ��̬��Ϣ�ṹ��*/
typedef struct tagRStaInfo
{
	uint32 ui32StaId		    ;                                      /*��վID��ʶ*/
	string strName              ;                                      /*��վ����*/
	string strReceAntType	    ;                                      /*���ջ���������*/
	uint8  ui8ReceAntCorrType   ;                                      /*���ջ����߸�������*/
	double daReveAntCorrVal[3]  ;                                      /*���ջ����߸���ֵ*/
	uint8  ui8CoorSysType       ;                                      /*��վ����ϵͳ����*/
	uint8  ui8CoorType          ;                                      /*��վ��������*/
	double daCoorVal[3]         ;                                      /*��վ��ά����ֵ*/
}RStaInfo;


/*��վ������Ϣ�ṹ��*/
typedef struct tagStaSolveInfo
{

	uint32    ui32StaId          ;                                    /*��վID��ʶ*/
	//TimeGnss  tTimeGnss          ;                                    /*ʱ��*/
	TimeStamp tTimeStamp         ;                                    /*�۲�ʱ��*/
	uint8     ui8SatNum          ;                                    /*������Ŀ*/
	double    dZeniTropDelayDry  ;                                    /*�춥�������ӳٸɷ���*/
	SatSoln*  pSatSlon           ;                                    /*�۲⵽��������Ϣ�ṹ��ָ��*/
}StaSolveInfo;

/*�������ǽ�����Ϣ�ṹ��*/
typedef struct tagComnSatSlon            
{
	uint8  ui8RefSatId          ;                                     /*�ο���*/
	uint8  ui8RovSatId          ;                                     /*������*/
	uint8  ui8SolStatus         ;                                     /*����״̬*/
	uint32 ui32FixedNum         ;                                     /*�ۼƹ̶�����*/
	double dTropDelay           ;                                     /*�������ӳ���*/
	double dIonoDelay           ;                                     /*������ӳ���*/
	double dPCombinErr          ;                                     /*α���ۺ����*/
	double dCCombinErr          ;                                     /*��λ�ۺ����*/
	int    CodeType             ;                                     /*�۲�ֵ����*/
}ComnSatSlon;
/*���߾�̬��Ϣ�ṹ��*/
typedef struct tagBLineInfo
{
	uint32 ui32Id                ;                                   /*����ID*/
	uint32 ui32StartPointId      ;                                   /*�����������*/
	uint32 ui32EndPointId        ;                                   /*�����յ�����*/
}BLineInfo;

/*���߽�����Ϣ�ṹ��*/
typedef struct tagBLineSlon               
{
	uint32     ui32Id            ;                                   /*����ID*/
	TimeStamp  tTimeStamp        ;                                   /*�۲�ʱ��*/
	double     Ddist             ;                                   /*���߳���*/
	uint8      ui8ComnSatNum     ;                                   /*����������Ŀ*/
	uint8      ui8EquNum         ;                                   /*α�����λ˫��̸���*/
	uint8      ui8FixedSatNum    ;                                   /*ģ���ȹ̶���������Ŀ*/
	ComnSatSlon* pComnSatSoln    ;                                   /*�������ǽṹ����Ϣָ��*/
}BLineSlon;

/*��������Ϣ�ṹ��*/
typedef struct tagTriNetInfo
{
	uint32 ui32Id                ;                                   /*������ID*/
	uint32 ui32BLine1Id          ;                                   /*����1����*/
	uint32 ui32BLine2Id          ;                                   /*����2����*/
	uint32 UI32bLine3Id          ;                                   /*����3����*/
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
/*����վ�۲���Ϣ�ṹ��*/
typedef struct tagVRSInfo               
{
	uint32  un32Id               ;                                   /*����վID*/
	uint8   un8CoorType          ;                                   /*����վ��������*/
	double  daCoorVal[3]         ;                                   /*����վ����ֵ*/
	uint32  ui32AttTriNerId      ;                                   /*������������ID*/
	uint32* ui32AtthBLineId      ;                                   /*��������վ���õĻ�������ָ��*/
	uint32  ui32MasterStaId      ;                                   /*����վ������*/
	int     nNum                 ;                                   /*�۲���Ԫ����*/
	int     nMax                 ;                                   /*���۲���Ԫ����*/
	ObsData* pObsData            ;                                   /*�۲����ݽṹ��ָ��*/
	StaInfo sta                 ;                                   /* station parameters */
}VRSInfo;

extern RStaInfo      *gRStaInfo        ;                             /*��վ��̬��Ϣ*/
extern StaSolveInfo  *gStaSolveInfo    ;                             /*��վ������Ϣ*/
extern BLineInfo     *gBlineInfo       ;                             /*���߾�̬��Ϣ*/
extern BLineSlon     *gBlineSlon       ;                             /*���߽�����Ϣ*/
extern TriNetInfo    *gTrinetInfo      ;                             /*��������Ϣ*/
extern VRSInfo       *gVrsInfo         ;                             /*����վ�۲���Ϣ*/


#endif // !COMMONVARIABLES_H

