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
	/*��ʼ����վ�۲�����*/
	InitObs(gpSinstaobsdata, 3);
	/*��ʼ�������۲�����*/
	InitBroadcastEphData(pAllbroadcastephdata);
	/*��ʼ��RTCM����*/
	cRtcm.InitRtcm();
	/*����������������*/
	cRtcm.DecodeEphType(pAllbroadcastephdata, pStaBuffer, ifile0);
	/*��ʼNRTK���߽������*/
	cBline.InitRtkopt(rtkopt);
	/*��ʼNRTK���ߵĽ������*/
	cBline.InitRtkParam(3);
	/*ͨ����rtcm�ļ�����۲�����*/
	if (cRtcm.ReadObsFile(ifile1, ifile2, ifile3))
	{
		while (cRtcm.DecodeObsType(gpSinstaobsdata,pStaBuffer) > 0)
		{
			cStation.SetObsEphInfo(gpSinstaobsdata, pAllbroadcastephdata);
			/* �Ի�վ���ݽ���Ԥ����*/
			cStation.Process(gpSinstaobsdata[0].ui16MaxNum);
			
			for (index = 0; index < 3; index++)
			{
				/*�������*/
				cBline.Process(gpSinstaobsdata, gBlineInfo[index],index);
			}
			/*��������۲�ֵ*/
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
