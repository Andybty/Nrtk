#include <direct.h>
#include "CInterface.h"
using namespace std;

char ifile0[1024] = "../../data/BRDC0712.rtcm3";
char ifile1[1024] = "../../data/AQTC0712.rtcm3";
char ifile2[1024] = "../../data/AQYX0712.rtcm3";
char ifile3[1024] = "../../data/CZDZ0712.rtcm3";
int index;
double dt;
RtkOpt rtkopt = { {0} };
SinStaObsData* gpSinstaobsdata;
AllBroadcastEphData* pAllbroadcastephdata;
CRTCM cRtcm;
CStation  cStation;
CBline cBline;
CVRS cVRS;
CInterface::CInterface()
{

}

CInterface::~CInterface()
{

}
bool CInterface::RStaDataInput(uint8 ui8StaNum, uint8* pStaBuffer)
{
	/*初始化测站观测数据*/
	InitObs(gpSinstaobsdata, 3);
	/*初始化星历观测数据*/
	InitBroadcastEphData(pAllbroadcastephdata);
	/*初始化RTCM数据*/
	cRtcm.InitRtcm();
	/*解码卫星星历数据*/
	cRtcm.DecodeEphType(pAllbroadcastephdata, pStaBuffer, ifile0);
	/*初始NRTK基线解算策略*/
	cBline.InitRtkopt(rtkopt);
	/*初始NRTK基线的解算参数*/
	cBline.InitRtkParam(3);
	/*通过读rtcm文件解码观测数据*/
	if (cRtcm.ReadObsFile(ifile1, ifile2, ifile3))
	{
		while (cRtcm.DecodeObsType(gpSinstaobsdata,pStaBuffer) > 0)
		{
			cStation.SetObsEphInfo(gpSinstaobsdata, pAllbroadcastephdata);
			/* 对基站数据进行预处理*/
			cStation.Process(gpSinstaobsdata[0].ui16MaxNum);
			
			for (index = 0; index < 3; index++)
			{
				/*处理基线*/
				cBline.Process(gpSinstaobsdata, gBlineInfo[index],index);
			}
			/*生成虚拟观测值*/
			cVRS.Process(gpSinstaobsdata,gTrinetInfo, gVrsInfo,1);

			for (index = 0; index < 3; index++)
			{
				
				gStaSolveInfo[index].ui8SatNum = 0;
			}
		}
	}
	
	
	





	FreeObs(gpSinstaobsdata);
	FreeBroadcastEphData(pAllbroadcastephdata);
	return false;
}

bool CInterface::VRSDataOutput(uint8 ui8VRSNum, uint8* pStaBuffer)
{
	return false;
}
