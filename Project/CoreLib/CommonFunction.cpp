#include <direct.h>
#include "CommonFunction.h"
#include "CRTCM.h"
extern int SatSys(int sat, int* prn)
{
    int sys = SYS_NONE;
    if (sat <= 0 || MAXSAT < sat) sat = 0;
    else if (sat <= NSATGPS) {
        sys = SYS_GPS; sat += MINPRNGPS - 1;
    }
    else if ((sat -= NSATGPS) <= NSATGLO) {
        sys = SYS_GLO; sat += MINPRNGLO - 1;
    }
    else if ((sat -= NSATGLO) <= NSATGAL) {
        sys = SYS_GAL; sat += MINPRNGAL - 1;
    }
    else if ((sat -= NSATGAL) <= NSATQZS) {
        sys = SYS_QZS; sat += MINPRNQZS - 1;
    }
    else if ((sat -= NSATQZS) <= NSATCMP) {
        sys = SYS_CMP; sat += MINPRNCMP - 1;
    }
    else if ((sat -= NSATCMP) <= NSATIRN) {
        sys = SYS_IRN; sat += MINPRNIRN - 1;
    }
    else if ((sat -= NSATIRN) <= NSATLEO) {
        sys = SYS_LEO; sat += MINPRNLEO - 1;
    }
    else if ((sat -= NSATLEO) <= NSATSBS) {
        sys = SYS_SBS; sat += MINPRNSBS - 1;
    }
    else sat = 0;
    if (prn) *prn = sat;
    return sys;
}

/* satellite and obs code to frequency -----------------------------------------
* convert satellite and obs code to carrier frequency
* args   : int    sat       I   satellite number
*          uint8_t code     I   obs code (CODE_???)
*          nav_t  *nav_t    I   navigation data for GLONASS (NULL: not used)
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/ 
double Sat2Freq(int sat, uint8_t code)
{
    int i, fcn = 0, sys, prn;

    sys = SatSys(sat, &prn);

    return code2freq(sys, code, fcn);
}

 void FreeObs(SinStaObsData* pSinstaobsdata)
{
    int n_num = pSinstaobsdata[0].ui16MaxNum;
    for (int i = 0; i < n_num; i++)
    {
        delete[] pSinstaobsdata[i].pObsData;
        pSinstaobsdata[i].pObsData = NULL;
    }
    delete[] pSinstaobsdata;
    pSinstaobsdata = NULL;
    

}

 void InitObs(SinStaObsData*& pSinstaobsdata, int n )
{
    if (0 != n)
    {
        pSinstaobsdata = new SinStaObsData[n];

        for (int i = 0; i < n; i++)
        {
            pSinstaobsdata[i].ui16ObsNum = MAXOBS;
            pSinstaobsdata[i].pObsData = new ObsData[MAXOBS];
            memset(pSinstaobsdata[i].pObsData, 0, MAXOBS * sizeof(ObsData));
            pSinstaobsdata[i].ui16MaxNum = n;
        }
    }
}



void InitBroadcastEphData(AllBroadcastEphData* &pAllbroadcastephdata)
{
    pAllbroadcastephdata = new AllBroadcastEphData;
    pAllbroadcastephdata->MaxNum= 2 * MAXSAT;
    pAllbroadcastephdata->nSatNum = 2 * MAXSAT;
    pAllbroadcastephdata->pBroadcastEphData = new BroadcastEphData[2 * MAXSAT];
}

void FreeBroadcastEphData(AllBroadcastEphData* pAllbroadcastephdata)
{
    delete[] pAllbroadcastephdata->pBroadcastEphData;
    pAllbroadcastephdata->pBroadcastEphData = NULL;
    delete[] pAllbroadcastephdata;
    pAllbroadcastephdata = NULL;
}

void FreeStaSolveInfo(int n,StaSolveInfo*& pStaSolveInfo)
{
    
    for (int i = 0; i < n; i++)
    {
        delete[]pStaSolveInfo[i].pSatSlon;
         pStaSolveInfo[i].pSatSlon = NULL;
    }
    if (pStaSolveInfo != NULL)
    {
        delete[]pStaSolveInfo;
    }
    pStaSolveInfo = NULL;
}

void  InitStaSolveInfo(int n, StaSolveInfo*& pStaSolveInfo)
{   
    pStaSolveInfo = new StaSolveInfo[n];
    for (int i = 0; i < n; i++)
    {
        pStaSolveInfo[i].pSatSlon = NULL;
        pStaSolveInfo[i].pSatSlon = new SatSoln[MAXSAT];
        memset(pStaSolveInfo[i].pSatSlon, 0, sizeof(SatSoln)* MAXSAT);
        pStaSolveInfo[i].tTimeStamp = { 0 };
        pStaSolveInfo[i].ui8SatNum = 0;
        pStaSolveInfo[i].dZeniTropDelayDry = 0.0;
    }
   
    
    return ;
}

bool InitBLineInfo(int n, BLineInfo*& pBlineInfo)
{
    if (n < 0)
        return  false;
    else
    {
        pBlineInfo = NULL;
        pBlineInfo = new BLineInfo[3];
        memset(pBlineInfo, 0, 3 * sizeof(BLineInfo));

    }
    return true;
}



void FreeBLinInfo(BLineInfo* pBlineInfo)
{
    delete[]pBlineInfo;
    pBlineInfo = NULL;
}

bool InitRStaInfo(int n, RStaInfo*& pRStaInfo)
{
    if(n < 0)
        return  false;
    else
    {
        pRStaInfo = NULL;
        pRStaInfo = new RStaInfo[3];
        memset(pRStaInfo, 0, 3 * sizeof(RStaInfo));

    }
    return true;
}
void FreeRStaInfo(RStaInfo* pRStaInfo)
{
    delete[]pRStaInfo;
    pRStaInfo = NULL;
}

bool InitBLineSlon(int n, BLineSlon*& gBlineSlon)
{
    int i;
    if (n < 0)
        return  false;
    else
    {
        gBlineSlon = NULL;
        gBlineSlon = new BLineSlon[3];
    }

    for (int i = 0; i < n; i++)
    {
        gBlineSlon[i].pComnSatSoln = NULL;
        gBlineSlon[i].pComnSatSoln = new ComnSatSlon[MAXSAT];
        memset(gBlineSlon[i].pComnSatSoln, 0, sizeof(ComnSatSlon) * MAXSAT);
        gBlineSlon[i].tTimeStamp = { 0 };
        gBlineSlon[i].Ddist = 0.0;
        gBlineSlon[i].ui32Id = 0;
        gBlineSlon[i].ui8ComnSatNum = 0;
        gBlineSlon[i].ui8FixedSatNum = 0;
    }
    return true;
}

void FreeBLineSlon(int n, BLineSlon*& gBlineSlon)
{
    for (int i = 0; i < n; i++)
    {
        delete[]gBlineSlon[i].pComnSatSoln;
        gBlineSlon[i].pComnSatSoln = NULL;
    }
    if (gBlineSlon != NULL)
    {
        delete[]gBlineSlon;
    }
    gBlineSlon = NULL;
}

void InitTriNetInfo(int n, TriNetInfo*& gTrinetInfo)
{
    gTrinetInfo = new TriNetInfo[n];
    for (int i = 0; i < n; i++)
    {
        gTrinetInfo[i] = { 0 };
    }
}