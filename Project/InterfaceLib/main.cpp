#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <direct.h>
#include "CInterface.h"
#include "CommonFunction.h"
using namespace std;


RStaInfo* gRStaInfo           ;            /*基站静态信息*/
StaSolveInfo* gStaSolveInfo   ;            /*基站解算信息*/
BLineInfo* gBlineInfo         ;            /*基线静态信息*/
BLineSlon* gBlineSlon         ;            /*基线解算信息*/
TriNetInfo* gTrinetInfo       ;            /*三角网信息*/
VRSInfo*    gVrsInfo          ;            /*虚拟站观测信息*/
 void SetBlineInfo(BLineInfo* pBlineInfo)
{
	 pBlineInfo[0].ui32Id = 1;
	 pBlineInfo[0].ui32StartPointId = 'A';
	 pBlineInfo[0].ui32EndPointId = 'B';
	 pBlineInfo[1].ui32Id = 2;
	 pBlineInfo[1].ui32StartPointId = 'B';
	 pBlineInfo[1].ui32EndPointId = 'C';
	 pBlineInfo[2].ui32Id = 3;
	 pBlineInfo[2].ui32StartPointId = 'C';
	 pBlineInfo[2].ui32EndPointId = 'A';

}

 void SetRStaInfo(int n,RStaInfo* pRStaInfo)
 {
	 for (int i = 0; i < n; i++)
	 {
		 //pRStaInfo[i] = { 0 };
		 memset(pRStaInfo + i, 0, sizeof(RStaInfo));
	 }
	 pRStaInfo[0].ui32StaId = 'A';
	 pRStaInfo[0].daCoorVal[0] = -2537912.4797; pRStaInfo[0].daCoorVal[1] = 4867463.4823;  pRStaInfo[0].daCoorVal[2] = 3236884.0440;
	 pRStaInfo[1].ui32StaId = 'B';
	 pRStaInfo[1].daCoorVal[0] = -2502020.3355; pRStaInfo[1].daCoorVal[1] = 4895812.1999;  pRStaInfo[1].daCoorVal[2] = 3222096.0487;
	 pRStaInfo[2].ui32StaId = 'C';
	 pRStaInfo[2].daCoorVal[0] = -2507959.3203; pRStaInfo[2].daCoorVal[1] = 4919846.2238;  pRStaInfo[2].daCoorVal[2] = 3180879.8902;
 }
 void SetTriNetInfo(TriNetInfo* gTrinetInfo )
 {
	 gTrinetInfo->ui32BLine1Id = 0;
	 gTrinetInfo->ui32BLine2Id = 1;
	 gTrinetInfo->UI32bLine3Id = 2;
	 gTrinetInfo->ui32Id = 0;
 }
 void SetVrsInfo(VRSInfo* &gVrsInfo)
 {
	 int i;
	 gVrsInfo = new VRSInfo;
	 gVrsInfo->daCoorVal[0] = (-2537912.4797 - 2502020.3355 - 2507959.3203) / 3;
	 gVrsInfo->daCoorVal[1] = (4867463.4823 + 4895812.1999 + 4919846.2238) / 3;
	 gVrsInfo->daCoorVal[2] = (3236884.0440 + 3222096.04875 + 3180879.8902) / 3;
	 gVrsInfo->un32Id = 0;
	 gVrsInfo->un8CoorType = 1;
	 gVrsInfo->ui32AttTriNerId = 0;
	 gVrsInfo->ui32MasterStaId = 0;
	 gVrsInfo->nNum = 32;
	 gVrsInfo->pObsData = new ObsData[32];
	
	 gVrsInfo->sta.name[0] = gVrsInfo->sta.marker[0] = '\0';
	 gVrsInfo->sta.antdes[0] = gVrsInfo->sta.antsno[0] = '\0';
	 gVrsInfo->sta.rectype[0] = gVrsInfo->sta.recver[0] = gVrsInfo->sta.recsno[0] = '\0';
	 gVrsInfo->sta.antsetup = gVrsInfo->sta.itrf = gVrsInfo->sta.deltype = 0;
	 for (i = 0; i < 3; i++) {
		 gVrsInfo->sta.pos[i] = gVrsInfo->daCoorVal[i];
	 }
	 gVrsInfo->sta.hgt = 0.0;

 }
int main()
{
	
	uint8* pStaBuffer=NULL;

	CInterface Interface;

	InitBLineInfo(3, gBlineInfo);
	InitStaSolveInfo(3, gStaSolveInfo);
	InitRStaInfo(3, gRStaInfo);
	InitBLineSlon(3, gBlineSlon);
	SetBlineInfo(gBlineInfo);
	SetRStaInfo(3, gRStaInfo);
    InitTriNetInfo(1, gTrinetInfo);
	SetTriNetInfo( gTrinetInfo);
	SetVrsInfo(gVrsInfo);
	if (Interface.RStaDataInput(1, pStaBuffer))
		Interface.VRSDataOutput(1, pStaBuffer);


	
	FreeStaSolveInfo( 3, gStaSolveInfo);
	FreeBLinInfo(gBlineInfo);
	FreeRStaInfo(gRStaInfo);
	return 0;
}