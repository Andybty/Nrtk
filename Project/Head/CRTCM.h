#ifndef CRTCM_H
#define CRTCM_H
#include"CommonVariables.h"
#include"CommonFunction.h"
#include"MacroDefine.h"
#include "CTime.h"
#define RTCM3PREAMB 0xD3                                         /* rtcm ver.3 frame preamble */
#define POLYCRC32   0xEDB88320u                                  /* CRC32 polynomial */
/* constants -----------------------------------------------------------------*/

#define PRUNIT_GPS  299792.458                                   /* rtcm ver.3 unit of gps pseudorange (m) */
#define PRUNIT_GLO  599584.916                                   /* rtcm ver.3 unit of glonass pseudorange (m) */
#define RANGE_MS    (CLIGHT*0.001)                               /* range in 1 ms */



#define P2_5        0.03125                                      /* 2^-5 */
#define P2_6        0.015625                                     /* 2^-6 */
#define P2_10       0.0009765625                                 /* 2^-10 */
#define P2_11       4.882812500000000E-04                        /* 2^-11 */
#define P2_15       3.051757812500000E-05                        /* 2^-15 */
#define P2_17       7.629394531250000E-06                        /* 2^-17 */
#define P2_19       1.907348632812500E-06                        /* 2^-19 */

#define P2_24       5.960464477539063E-08                        /* 2^-24 */
#define P2_28       3.725290298461914E-09                        /* 2^-28 */
#define P2_29       1.862645149230957E-09                        /* 2^-29 */
#define P2_30       9.313225746154785E-10                        /* 2^-30 */
#define P2_31       4.656612873077393E-10                        /* 2^-31 */
#define P2_32       2.328306436538696E-10                        /* 2^-32 */
#define P2_33       1.164153218269348E-10                        /* 2^-33 */
#define P2_34       5.820766091346740E-11                        /* 2^-34 */
#define P2_35       2.910383045673370E-11                        /* 2^-35 */

#define P2_38       3.637978807091710E-12                        /* 2^-38 */
#define P2_39       1.818989403545856E-12                        /* 2^-39 */
#define P2_40       9.094947017729280E-13                        /* 2^-40 */
#define P2_41       4.547473508864641E-13                        /* 2^-41 */
#define P2_43       1.136868377216160E-13                        /* 2^-43 */
#define P2_46       1.421085471520200E-14                        /* 2^-46 */
#define P2_48       3.552713678800501E-15                        /* 2^-48 */
#define P2_50       8.881784197001252E-16                        /* 2^-50 */
#define P2_55       2.775557561562891E-17                        /* 2^-55 */
#define P2_59       1.734723475976810E-18                        /* 2^-59 */
#define P2_66       1.355252715606880E-20                        /* 2^-66 */




 /* multi-signal-message header type */
typedef struct tagMesHead {
	uint8_t iod;                                                /* MSM signal ID table */
	uint8_t time_s;                                             /* cumulative session transmitting time */
	uint8_t clk_str;                                            /* clock steering indicator */
	uint8_t clk_ext;                                            /* external clock indicator */
	uint8_t smooth;                                             /* divergence free smoothing indicator */
	uint8_t tint_s;                                             /* soothing interval */
	uint8_t nsat, nsig;                                         /* number of satellites/signals */
	uint8_t sats[64];                                           /* satellites */
	uint8_t sigs[32];                                           /* signals */
	uint8_t cellmask[64];                                       /* cell mask */
} msm_h_t;
/*RTCM数据流结构体*/
typedef struct tagRtcmData {
	int staid;                                                  /* station id */
	int seqno;                                                 /* sequence number for rtcm 2 or iods msm */
	int nbyte;                                                 /* number of bytes in message buffer */
	int nbit;                                                  /* number of bits in word buffer */
	int len;                                                   /* message length (bytes) */
	int obsflag;                                               /* obs data complete flag (1:ok,0:not complete) */
	unsigned char buff[1200];                                  /* message buffer */
	gtime_t time;                                              /* message time */
	gtime_t time_s;                                            /* message start time */
	char msmtype[7][128];                                      /* msm signal types */
	char opt[256];                                             /* RTCM dependent options */
	uint16_t lock[MAXSAT][NFREQ + NEXOBS];                     /* lock time */
	uint16_t loss[MAXSAT][NFREQ + NEXOBS];                     /* loss of lock count */
	gtime_t lltime[MAXSAT][NFREQ + NEXOBS];                    /* last lock time */
	double cp[MAXSAT][NFREQ + NEXOBS]; /* carrier-phase measurement */
	unsigned int nmsg3[400]; /* message count of RTCM 3 (1-299:1001-1299,300-399:2000-2099,0:ohter) */
}RtcmData;

/* stream converter type */
typedef struct tagStrConv
{
	RtcmData out;                                               /* rtcm output data buffer */
	int itype, otype;                                           /* input and output stream type */
	int nmsg;                                                   /* number of output messages */
	int msgs[32];                                               /* output message types */
	double tint[32];                                            /* output message intervals (s) */
	unsigned int tick[32];                                      /* cycle tick of output message */
	int stasel;                                                /* station info selection (0:remote,1:local) */
}StrConv;
/*
@@brief RTCM类，通过RTCM类的方法解析观测数据和星历数据；
*/
class CRTCM
{
public:
	CRTCM();
	~CRTCM();
public:
	/**
	@brief  解码卫星星历数据
	@param  strFileName[in] 输入的数据流文件名
	@param  pui8Stream [in] 输入的数据流指针
	@param  AllBroadcastEphData [out] 输出广播星历数据结构体指针
	@return 函数成功，返回值为2，否则返回-1。
	*/
	int DecodeEphType(AllBroadcastEphData* pAllBroadcastEphData, uint8* pui8Stream, string strFileName);
	/**
	@brief  解码基站观测数据
	@param  pui8Stream[in] 输入数据流指针
	@param  pSinStaObsData [out] 输出基站观测数据结构体指针
	@return 函数成功，返回值为1，否则返回-1。
	*/
	int DecodeObsType(SinStaObsData* pSinStaObsData, uint8* pui8Stream);
	/**
	@brief  初始化解码数据流结构体
	@param	无
	@return 无
	*/
	void InitRtcm();
	/**
	@brief  从文件中读取RTCM数据流
	@param  strFileName1[in] 第一个基站观测数据文件名
	@param  strFileName2[in] 第二个基站观测数据文件名
	@param  strFileName3[in] 第三个基站观测数据文件名
	@return 文件读取成功返回true，否则返回false
	*/
	bool ReadObsFile(string strFileName1, string strFileName2, string strFileName3);
	/**
	@brief  初始化VRS数据流结构体
	@param  无
	@return 无
	*/
	void InitStrCov();
	/**
	@brief  写观测值信息
	@param  time[in]  观测值时间
	@param  vrsInfo[in]  虚拟站观测值结构体
	@param  conv[in]  输出观测值数据流结构体
	@return 无
	*/
	void Write_Obs(gtime_t time,VRSInfo*vrsInfo,StrConv* conv);
	/**
	@brief  写入rtcm3 版本的观测值消息类型
	@param  vrsInfo[in] 虚拟观测值信息
	@param  msg[in]  观测值消息类型
	@param  sync[in]  该条信息是否完成的标识
	@param  conv[out]  输出观测值数据流结构体
	@return 无
	*/
	void Write_rtcm3_msm(RtcmData* out, VRSInfo* vrsInfo, int msg, int sync);
	/**
	@brief  根据消息类型按照RTCM格式生成相应的编码格式
	@param  vrsInfo[in] 虚拟观测值信息
	@param  rtcmdata[out]  rtcm结构体信息
	@param  msg[in]  消息类型
	@param  conv[out]  输出观测值数据流结构体
	@return 存在相应的消息类型则返回1
	*/
	int Gen_Rtcm3(VRSInfo* vrsInfo, RtcmData* rtcmdata, int msg, int sync);
private:
	/**
	@brief  初始化观测值数据结构体
	@param	n[in] 站的数量，默认值为1
	@return 无
	*/
	void InitObs(int n = 1);
	/**
	@brief  初始化RTCM结构体
	@param	rtcmdata[in] rtcm结构体
	@return 无
	*/
	void InitRtcm(RtcmData& rtcmdata);
	/**
	@brief  判断消息类型是不是观测值类型
	@param	msg[in] rtcm结构体
	@return  是观测值类型则返回真
	*/
	bool Is_ObsMsg(int msg);
	/**
	@brief  测试时间间隔
	@param	time[in] 观测时间
	@param	tint[in] 时间间隔
	@return  正确的时间间隔则TRUE
	*/
	bool Is_Tint(gtime_t time, double tint);
	/**
	@brief  把卫星号转换为消息体中的卫星索引
	@param	sys[in] 卫星系统
	@param	sat[in] 卫星编号
	@return  卫星在消息体中的索引
	*/
	int To_SatId(int sys, int sat);
	/**
	@brief  把卫星号转换为消息体中的卫星信号索引
	@param	sys[in] 卫星系统
	@param	code[in] 卫星码类型
	@return  观测值类型在消息体中的信号索引
	*/
	int To_SigId(int sys, unsigned char code);
	/**
	@brief  对无符号的字符串进行编码
	@param	buff[in] 字符串地址
	@param	pos[in]  字节的位置
	@param	len[in]  
	@param	data[in]  要编码的数据
	@return  无
	*/
	void SetBitu(unsigned char* buff, int pos, int len, unsigned int data);
	/**
	@brief  对有符号的字符串进行编码
	@param	buff[in] 字符串地址
	@param	pos[in]  字节的位置
	@param	len[in]
	@param	data[in]  要编码的数据
	@return  无
	*/
	void SetBits(unsigned char* buff, int pos, int len, int data);
	/**
	@brief  对虚拟观测值按照1074电文格式进行编码
	@param	vrsInfo[in] 虚拟观测值结构体指针
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   观测值所属系统
	@param	sysc[in]  消息完成标识
	@return  无
	*/
	int Encode_Msm4(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int sync);
	/**
	@brief  对虚拟观测值按照1074电文格式进行编码
	@param	vrsInfo[in] 虚拟观测值结构体指针
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   观测值所属系统
	@param	sysc[in]  消息完成标识
	@return  无
	*/
	int Encode_Type1005(VRSInfo* vrsInfo, RtcmData* rtcmdata,int sync);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	int Encode_Msm_Head(int type, VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int sync, int* nsat,
		int* ncell, double* rrng, double* rrate,unsigned char* info, double* psrng, double* phrng,
		double* rate, double* lock, unsigned char* half,float* cnr);
	/**
	@brief  对虚拟观测值按照1074电文格式进行编码
	@param	vrsInfo[in] 虚拟观测值结构体指针
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   观测值所属系统
	@param	sysc[in]  消息完成标识
	@return  无
	*/
	void Gen_Msm_Index(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int* nsat, int* nsig,
		int* ncell, unsigned char* sat_ind, unsigned char* sig_ind, unsigned char* cell_ind);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	void Gen_Msm_Sat(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int nsat, const uint8_t* sat_ind,
		double* rrng, double* rrate, uint8_t* info);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	void Gen_Msm_Sig(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int nsat, int nsig, int ncell, const uint8_t* sat_ind,
		const uint8_t* sig_ind, const uint8_t* cell_ind, const double* rrng,
		const double* rrate, double* psrng, double* phrng,
		double* rate, double* lock, uint8_t* half, float* cnr);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	double Locktime_d(gtime_t time, gtime_t* lltime, uint8_t LLI);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	int Encode_Msm_Int_Rrng(RtcmData* rtcm, int i, const double* rrng,
		int nsat);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	int Encode_Msm_Mod_Rrng(RtcmData* rtcm, int i, const double* rrng,
		int nsat);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	int Encode_Msm_Psrng(RtcmData* rtcm, int i, const double* psrng, int ncell);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	int Encode_Msm_Phrng(RtcmData* rtcm, int i, const double* phrng, int ncell);
	/**
	@brief  对虚拟观测值按照1074电文头格式进行编码
	@param	type[in]     电文类型
	@param	vrsInfo[in]  虚拟观测值结构体
	@param	rtcmdata[in]  rtcm结构体
	@param	sys[in]   系统标识
	@param	sysc[in]  消息完成标识
	@param	nsat[in]  卫星数目
	@param	ncell[in]  卫星数目
	@return  无
	*/
	int Encode_Msm_Lock(RtcmData* rtcm, int i, const double* lock, int ncell);
	int To_Msm_Lock(double lock);
	int Encode_Msm_Half_Amb(RtcmData* rtcm, int i, const uint8_t* half,
		int ncell);
	int Encode_Msm_Cnr(RtcmData* rtcm, int i, const float* cnr, int ncell);
protected:
	int Encode_Rtcm3(VRSInfo* vrsInfo, RtcmData* rtcmdata, int type, int sync);
	void set38bits(unsigned char* buff, int pos, double value);
public:
	SinStaObsData* m_pSinstaobsdata;						   //观测数据结构体指针
	AllBroadcastEphData* m_pAllbroadcastephdata;			   //广播星历结构体指针
	

public:
	RtcmData m_RtcmData;									   //RTCM数据流结构体
	StrConv  m_Strconv;                                        //vrs流结构体
	FILE* fp1;												   //第一个基站观测数据文件指针
	FILE* fp2;												   //第二个基站观测数据文件指针
	FILE* fp3;												   //第三个基站观测数据文件指针

};

//校验码对照表
static const unsigned int tbl_CRC24Q[] = {
	0x000000,0x864CFB,0x8AD50D,0x0C99F6,0x93E6E1,0x15AA1A,0x1933EC,0x9F7F17,
	0xA18139,0x27CDC2,0x2B5434,0xAD18CF,0x3267D8,0xB42B23,0xB8B2D5,0x3EFE2E,
	0xC54E89,0x430272,0x4F9B84,0xC9D77F,0x56A868,0xD0E493,0xDC7D65,0x5A319E,
	0x64CFB0,0xE2834B,0xEE1ABD,0x685646,0xF72951,0x7165AA,0x7DFC5C,0xFBB0A7,
	0x0CD1E9,0x8A9D12,0x8604E4,0x00481F,0x9F3708,0x197BF3,0x15E205,0x93AEFE,
	0xAD50D0,0x2B1C2B,0x2785DD,0xA1C926,0x3EB631,0xB8FACA,0xB4633C,0x322FC7,
	0xC99F60,0x4FD39B,0x434A6D,0xC50696,0x5A7981,0xDC357A,0xD0AC8C,0x56E077,
	0x681E59,0xEE52A2,0xE2CB54,0x6487AF,0xFBF8B8,0x7DB443,0x712DB5,0xF7614E,
	0x19A3D2,0x9FEF29,0x9376DF,0x153A24,0x8A4533,0x0C09C8,0x00903E,0x86DCC5,
	0xB822EB,0x3E6E10,0x32F7E6,0xB4BB1D,0x2BC40A,0xAD88F1,0xA11107,0x275DFC,
	0xDCED5B,0x5AA1A0,0x563856,0xD074AD,0x4F0BBA,0xC94741,0xC5DEB7,0x43924C,
	0x7D6C62,0xFB2099,0xF7B96F,0x71F594,0xEE8A83,0x68C678,0x645F8E,0xE21375,
	0x15723B,0x933EC0,0x9FA736,0x19EBCD,0x8694DA,0x00D821,0x0C41D7,0x8A0D2C,
	0xB4F302,0x32BFF9,0x3E260F,0xB86AF4,0x2715E3,0xA15918,0xADC0EE,0x2B8C15,
	0xD03CB2,0x567049,0x5AE9BF,0xDCA544,0x43DA53,0xC596A8,0xC90F5E,0x4F43A5,
	0x71BD8B,0xF7F170,0xFB6886,0x7D247D,0xE25B6A,0x641791,0x688E67,0xEEC29C,
	0x3347A4,0xB50B5F,0xB992A9,0x3FDE52,0xA0A145,0x26EDBE,0x2A7448,0xAC38B3,
	0x92C69D,0x148A66,0x181390,0x9E5F6B,0x01207C,0x876C87,0x8BF571,0x0DB98A,
	0xF6092D,0x7045D6,0x7CDC20,0xFA90DB,0x65EFCC,0xE3A337,0xEF3AC1,0x69763A,
	0x578814,0xD1C4EF,0xDD5D19,0x5B11E2,0xC46EF5,0x42220E,0x4EBBF8,0xC8F703,
	0x3F964D,0xB9DAB6,0xB54340,0x330FBB,0xAC70AC,0x2A3C57,0x26A5A1,0xA0E95A,
	0x9E1774,0x185B8F,0x14C279,0x928E82,0x0DF195,0x8BBD6E,0x872498,0x016863,
	0xFAD8C4,0x7C943F,0x700DC9,0xF64132,0x693E25,0xEF72DE,0xE3EB28,0x65A7D3,
	0x5B59FD,0xDD1506,0xD18CF0,0x57C00B,0xC8BF1C,0x4EF3E7,0x426A11,0xC426EA,
	0x2AE476,0xACA88D,0xA0317B,0x267D80,0xB90297,0x3F4E6C,0x33D79A,0xB59B61,
	0x8B654F,0x0D29B4,0x01B042,0x87FCB9,0x1883AE,0x9ECF55,0x9256A3,0x141A58,
	0xEFAAFF,0x69E604,0x657FF2,0xE33309,0x7C4C1E,0xFA00E5,0xF69913,0x70D5E8,
	0x4E2BC6,0xC8673D,0xC4FECB,0x42B230,0xDDCD27,0x5B81DC,0x57182A,0xD154D1,
	0x26359F,0xA07964,0xACE092,0x2AAC69,0xB5D37E,0x339F85,0x3F0673,0xB94A88,
	0x87B4A6,0x01F85D,0x0D61AB,0x8B2D50,0x145247,0x921EBC,0x9E874A,0x18CBB1,
	0xE37B16,0x6537ED,0x69AE1B,0xEFE2E0,0x709DF7,0xF6D10C,0xFA48FA,0x7C0401,
	0x42FA2F,0xC4B6D4,0xC82F22,0x4E63D9,0xD11CCE,0x575035,0x5BC9C3,0xDD8538
};

/**
@brief  字符串中提取无符号型位数据
@param  buff[in] 字符串指针
@param  pos[in] 位数据读取的起始位置
@param  len[in] 读取的数据长度
@return 提取的无符号型位数据
*/
unsigned int getbitu(const unsigned char* buff, int pos, int len);
/**
@brief  字符串中提取有符号型位数据
@param  buff[in] 字符串指针
@param  pos[in] 位数据读取的起始位置
@param  len[in] 读取的数据长度
@return 提取的有符号型位数据
*/
int getbits(const unsigned char* buff, int pos, int len);
/**
@brief  计算RTCM3格式数据crc24q校验码
@param  buff[in] 字符串指针
@param  len[in] 读取的数据长度
@return crc24q校验码
*/
unsigned int rtk_crc24q(const unsigned char* buff, int len);
/**
@brief  从数据流中输入观测数据流
@param  rtcm[in] RTCM数据流结构体指针
@param  data[in] 输入的数据字符
@param  pSinStaObsData[out] 观测数据结构体指针
@return 函数成功，返回值为1，否则返回-1。
*/
int input_rtcm3obs(SinStaObsData* pSinStaObsData, RtcmData* rtcm, unsigned char data);
/**
@brief  从数据流中输入星历数据流
@param  rtcm[in] RTCM数据流结构体指针
@param  data[in] 输入的数据字符
@param  AllBroadcastEphData[out] 广播星历数据结构体指针
@return 函数成功，返回值为2，否则返回-1。
*/
int input_rtcm3eph(AllBroadcastEphData* pAllBroadcastEphdata, RtcmData* rtcm, unsigned char data);
/**
@brief  从数据流中输入观测数据流
@param  rtcm[in] RTCM数据流结构体指针
@param  pSinStaObsData[out] 观测数据结构体指针
@return 函数成功，返回值为1，否则返回-1。
*/
int decode_rtcm3obs(SinStaObsData* pSinStaObsData, RtcmData* rtcm);
/**
@brief  从数据流中解析星历数据流
@param  rtcm[in] RTCM数据流结构体指针
@param  AllBroadcastEphData[out] 广播星历数据结构体指针
@return 函数成功，返回值为2，否则返回-1。
*/
int decode_rtcm3eph(AllBroadcastEphData* pAllBroadcastEphdata, RtcmData* rtcm);
/**
@brief  从数据流中读取MSM4格式的观测数据
@param  rtcm[in] RTCM数据流结构体指针
@param  sys[in] 卫星所属系统
@param  pSinStaObsData[out] 观测数据结构体指针
@return 函数成功，返回值为1，否则返回-1。
*/
int decode_msm4(SinStaObsData* pSinStaObsData, RtcmData* rtcm, int sys);
/**
@brief  从数据流中读取MSM格式的观测数据头信息
@param  rtcm[in] RTCM数据流结构体指针
@param  sys[in] 卫星所属系统
@param  sync[out] 多系统观测数据标识
@param  iod[out] 数据有效标识
@param  h[out] MSM信息头类型
@param  hsize[out] MSM信息头大小
@return 函数成功，返回值为卫星观测数据个数，否则返回-1。
*/
int decode_msm_head(RtcmData* rtcm, int sys, int* sync, int* iod, msm_h_t* h, int* hsize);
/**
@brief  调整GLONASS观测值时间
@param  rtcm[in,out] RTCM数据流结构体指针
@param  tod[in,out] 周内秒时间
@return 无。
*/
void adjday_glot(RtcmData* rtcm, double tod);
/**
@brief  调整BDS/GPS/GALILEO观测值时间
@param  rtcm[in,out] RTCM数据流结构体指针
@param  tod[in,out] 周内秒时间
@return 无。
*/
void adjweek(RtcmData* rtcm, double tow);
/**
@brief  观测值类型字符转换成整数标识
@param  obs[in] 卫星观测值数据结构体指针
@return 卫星观测值类型标识。
*/
uint8_t obs2code(const char* obs);
/**
@brief  根据卫星系统与观测类型获取频率索引
@param  sys[in] 卫星系统
@param  code[in] 卫星观测数据类型
@return 观测值频率索引。
*/
int code2idx(int sys, uint8_t code);
/**
@brief  观测值类型整数标识转换成字符
@param  obs[in] 卫星观测值类型整数标识
@return 观测值类型字符指针。
*/
char* code2obs(uint8_t code);
/**
@brief  当一个频率有多个类型的观测数据值时，获取每个类型的优先级
@param  sys[in] 卫星系统
@param  code[in] 卫星观测数据类型
@param  opt[in] 观测值类型选项指针
@return 观测值类型优先级，15：最高，-1：最低，0：错误。
*/
int getcodepri(int sys, uint8_t code, const char* opt);
/**
@brief  根据观测值类型获取索引
@param  sys[in] 卫星系统
@param  code[in] 卫星观测数据类型
@param  n[in] 观测值类型总个数
@param  opt[in] 观测值类型选项指针
@param  idx[out] 观测值类型索引
@return 无。
*/
void sigindex(int sys, const uint8_t* code, int n, const char* opt, int* idx);
/**
@brief  存储MSM观测数据到结构体中
@param  rtcm[in] RTCM数据流结构体指针
@param  sys[in] 卫星系统
@param  h[in] MSM信息头类型指针
@param  r[in] 卫星基础观测值
@param  pr[in] 伪距观测值
@param  cp[in] 相位观测值
@param  cnr[in] 观测数据载噪比
@param  idx[out] 观测值类型索引
@return 无。
*/
void save_msm_obs(SinStaObsData* pSinStaObsData, RtcmData* rtcm, int sys, msm_h_t* h,
	const double* r, const double* pr, const double* cp, const double* rr,
	const double* rrf, const double* cnr, const int* lock, const int* ex, const int* half);
/**
@brief  根据卫星系统和prn号计算卫星号
@param  sys[in] 卫星系统
@param  prn[in] 卫星prn号
@return 卫星编号。
*/
int satno(int sys, int prn);
/**
@brief  根据卫星系统、观测值类型和GLONASS FCN号获取卫星频率
@param  sys[in] 卫星系统
@param  code[in] 卫星观测数据类型
@param  fcn[in] GLONASS FCN号
@return 卫星信号频率。
*/
double code2freq(int sys, uint8_t code, int fcn);
/**
@brief  从RTCM数据流中解析GPS星历数据
@param  rtcm[in] RTCM数据流结构体指针
@param  AllBroadcastEphData[out] 广播星历数据结构体指针
@return 函数成功，返回值为2，否则返回-1。
*/
int decode_type1019(RtcmData* rtcm, AllBroadcastEphData* pAllBroadcastEphdata);
#endif // !CRTM_H



