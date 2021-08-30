#include <direct.h>
#include "CVRS.h"
CVRS::CVRS()
{

}

CVRS::~CVRS()
{

}
double CVRS::CalStaVrsDis(double* dpSta, double* dpVrs)
{
	int i;
	double dr[3];
	for (i = 0; i < 3; i++) dr[i] = dpSta[i] - dpVrs[i];
	return CMath::norm(dr, 3);
}

int CVRS::SelectMinIndex(double* dpDis, int n)
{
	double min = dpDis[0];
	int i = 0, index = 0;
	for (i = 1; i < n; i++)
	{
		if (dpDis[i] < min)
		{
			min = dpDis[i];
			index = i;
		}
	}
	return index;
}


int CVRS::SelMastReferSta(int index, int& rover1, int& rover2)
{
	int sta1index = index + index * 2;
	double dis[3];
	int master;
	for (int i = 0; i < 3; i++)
	{
		dis[i] = CalStaVrsDis(gRStaInfo[sta1index + i].daCoorVal, gVrsInfo[index].daCoorVal);
	}
	master = SelectMinIndex(dis, 3);
	switch (master)
	{
	case 0:
		rover1 = sta1index + 1;
		rover2 = sta1index + 2;
		return  sta1index;
	case 1:
		rover1 = sta1index;
		rover2 = sta1index + 2;
		return  sta1index + 1;
	case 2:
		rover1 = sta1index;
		rover2 = sta1index + 1;
		return  sta1index + 2;
	}
}


//int CVRS::SelCommonSat(int index)
//{
//	int sta1index = index + index * 2;
//	int nStaSatNum1 = gStaSolveInfo[sta1index].ui8SatNum;
//	int nStaSatNum2 = gStaSolveInfo[sta1index+1].ui8SatNum;
//	int nStaSatNum3 = gStaSolveInfo[sta1index + 2].ui8SatNum;
//	SatSoln* psatsoln1 = gStaSolveInfo[sta1index].pSatSlon;
//	SatSoln* psatsoln2 = gStaSolveInfo[sta1index+1].pSatSlon;
//	SatSoln* psatsoln3 = gStaSolveInfo[sta1index + 2].pSatSlon;
//	int i, j, k,i_n=0;
//	for (i = 0; i < nStaSatNum1; i++)
//	{
//		for (j = 0; j < nStaSatNum2; j++)
//		{
//			if (psatsoln1[i].ui8SatId != psatsoln2[j].ui8SatId) continue;
//			for (k = 0; k < nStaSatNum3; k++)
//			{
//				if (psatsoln1[i].ui8SatId != psatsoln2[k].ui8SatId) continue;
//				m_naSta1SatID[i_n]=i;            
//				m_naSta2SatID[i_n]=j;             
//				m_naSta3SatID[i_n]=k; 
//				m_naCommonSat[i_n++]= psatsoln1[i].ui8SatId;
//				break;
//			}
//			break;
//		}
//	}
//	return i_n;
//}

void CVRS::GenVRSObs(SinStaObsData Sinstaobsdata, VRSInfo& vrsInfo, int master, int n)
{
	int i, j, k, l,ik=0;
	double* rs, * es, dr, freq;
	ObsData* pObsData = Sinstaobsdata.pObsData;
	StaSolveInfo stasolveinfo = gStaSolveInfo[master];
	int nobs = stasolveinfo.ui8SatNum;
	rs = CMath::mat(1, nobs);
	es = CMath::mat(3, nobs);
	vrsInfo.nNum = 0;
	vrsInfo.ui32MasterStaId = master;
	for (i = 0; i < nobs; i++)
	{
		rs[i] = m_cStation.CalGeoDist(stasolveinfo.pSatSlon[i].daSatPos, vrsInfo.daCoorVal, es + i * 3);
		k = stasolveinfo.pSatSlon[i].uiObsindex;
		memset(&vrsInfo.pObsData[i], 0, sizeof(ObsData));
		for (j = 0; j < n; j++)
		{
			if (pObsData[k].ui8SatId != m_pVrscorrect[j].rov) continue;

			for (l = 0; l < NFREQ + NEXOBS; l++)
			{
				if ((int)pObsData[k].ui8aCodeType[l] != m_pVrscorrect[j].CodeType) continue;
				break;
			}
			if (l == NFREQ + NEXOBS)  break;
			freq = Sat2Freq(pObsData[k].ui8SatId, pObsData[k].ui8aCodeType[l]);
			dr = rs[i] - stasolveinfo.pSatSlon[i].dSta2SatDis- stasolveinfo.pSatSlon[i].dSatClk;
			vrsInfo.pObsData[ik].ui8SatId = pObsData[k].ui8SatId;
			vrsInfo.pObsData[ik].tTimeStamp = pObsData[k].tTimeStamp;
			vrsInfo.pObsData[ik].ui8aCodeType[l] = pObsData[k].ui8aCodeType[l];
			vrsInfo.pObsData[ik].daCarrPhase[l] = pObsData[k].daCarrPhase[l] + dr * freq / CLIGHT -
				m_pVrscorrect[j].dDCComErr * freq / CLIGHT;
			vrsInfo.pObsData[ik].daPseRange[l] = pObsData[k].daPseRange[l] + dr -
				m_pVrscorrect[j].dDPComErr;
			vrsInfo.pObsData[ik].daSigNoiRatio[l] = pObsData[k].daSigNoiRatio[l];
			vrsInfo.pObsData[ik].daDoppler[l] = pObsData[k].daDoppler[l];
			vrsInfo.pObsData[ik].ui8FreqNum++;
			

		}
		if (vrsInfo.pObsData[ik].ui8FreqNum > 0) {
			vrsInfo.nNum++;
			ik++;
		}
	}
	
	free(es); free(rs);
	/*虚拟站与主参考站之差*/



}


void  CVRS::SelBlineIndex(int master, int& bline1, int& bline2, int& flag)
{
	if (master == 3 * master)
	{
		bline1 = master;
		bline2 = master + 2;
		flag = 0;
	}
	else
	{
		bline1 = master - 1;
		bline2 = master;
		flag = 1;
	}
}


double* CVRS::LCM(int master, int rov1, int rov2, VRSInfo& vrsInfo)
{
	int n = 3;
	double* a, * Q, * coef, * W;
	double BLH0[3], BLH1[3], BLH2[3], BLHv[3];
	CCoordTrans::ecef2pos(gRStaInfo[master].daCoorVal, BLH0);
	CCoordTrans::ecef2pos(gRStaInfo[rov1].daCoorVal, BLH1);
	CCoordTrans::ecef2pos(gRStaInfo[rov2].daCoorVal, BLH2);
	CCoordTrans::ecef2pos(vrsInfo.daCoorVal, BLHv);
	a = CMath::zeros(n, 1);      Q = CMath::zeros(n * n, 1);
	coef = CMath::zeros(n * n, 1); W = CMath::zeros(n, 1);
	coef[0] = 1.0; coef[1] = 1.0; coef[2] = 1.0;
	coef[3] = BLH1[0] - BLH0[0];
	coef[4] = BLH2[0] - BLH0[0];
	coef[5] = BLH0[0] - BLH0[0];
	coef[6] = BLH1[1] - BLH0[1];
	coef[7] = BLH2[1] - BLH0[1];
	coef[8] = BLH0[1] - BLH0[1];
	W[0] = 1.0;
	W[1] = BLHv[0] - BLH0[0];
	W[2] = BLHv[1] - BLH0[1];
	CMath::lsq(coef, W, n, n, a, Q);



	//if (a != NULL) { free(a); a = NULL; }
	if (Q != NULL) { free(Q); Q = NULL; }
	if (W != NULL) { free(W); W = NULL; }
	if (coef != NULL)
	{
		free(coef);
		coef = NULL;
	}
	return a;
}

bool CVRS::InitVrsCorrect(int n)
{
	if (n > 0)
	{
		m_pVrscorrect = new VrsCorrect[n];
		memset(m_pVrscorrect, 0, n * sizeof(VrsCorrect));
		return true;

	}
	return false;


}

bool CVRS::JudgeConsist(int ref1, int rov1, int ref2, int rov2)
{
	if (ref1 == ref2 && rov1 == rov2)
		return true;
	return false;
}

int CVRS::GenVrsAtom(double* coeff, int bindex1, int bindex2, int flag)
{
	int n = MIN(gBlineSlon[bindex1].ui8EquNum, gBlineSlon[bindex2].ui8EquNum);
	int k = 0, ik = 0;
	if (InitVrsCorrect(n))
	{
		for (int i = 0; i < gBlineSlon[bindex1].ui8EquNum; i++)
		{
			for (int j = 0; j < gBlineSlon[bindex2].ui8EquNum; j++)
			{
				if ((gBlineSlon[bindex1].pComnSatSoln[i].CodeType == gBlineSlon[bindex2].pComnSatSoln[j].CodeType) &&
					JudgeConsist(gBlineSlon[bindex1].pComnSatSoln[i].ui8RefSatId, gBlineSlon[bindex1].pComnSatSoln[i].ui8RovSatId,
						gBlineSlon[bindex2].pComnSatSoln[j].ui8RefSatId, gBlineSlon[bindex2].pComnSatSoln[j].ui8RovSatId))
				{
					if (flag)
					{
						m_pVrscorrect[k].dDPComErr = -coeff[0] * gBlineSlon[bindex1].pComnSatSoln[i].dPCombinErr +
							coeff[1] * gBlineSlon[bindex2].pComnSatSoln[j].dPCombinErr;
						m_pVrscorrect[k].dDCComErr = -coeff[0] * gBlineSlon[bindex1].pComnSatSoln[i].dCCombinErr +
							coeff[1] * gBlineSlon[bindex2].pComnSatSoln[j].dCCombinErr;
					}
					else
					{
						m_pVrscorrect[k].dDPComErr = coeff[0] * gBlineSlon[bindex1].pComnSatSoln[i].dPCombinErr -
							coeff[1] * gBlineSlon[bindex2].pComnSatSoln[j].dPCombinErr;
						m_pVrscorrect[k].dDCComErr = coeff[0] * gBlineSlon[bindex1].pComnSatSoln[i].dCCombinErr -
							coeff[1] * gBlineSlon[bindex2].pComnSatSoln[j].dCCombinErr;
					}
					m_pVrscorrect[k].CodeType = gBlineSlon[bindex1].pComnSatSoln[i].CodeType;
					m_pVrscorrect[k].ref = gBlineSlon[bindex1].pComnSatSoln[i].ui8RefSatId;
					m_pVrscorrect[k].rov = gBlineSlon[bindex1].pComnSatSoln[i].ui8RovSatId;
					m_pVrscorrect[k++].time = gBlineSlon[bindex1].tTimeStamp;

					break;
				}
			}
		}
	}

	if (coeff)
	{
		free(coeff);
		coeff = NULL;
	}
	return k;
}


void  CVRS::Process(SinStaObsData* pSinstaobsdata,TriNetInfo* pTrinetInfo, VRSInfo* pVrsInfo, int nVrs)
{
	int i, natom;
	int master;
	int rov1, rov2;
	int bline1, bline2, flag;
	double* coeff;
	for (i = 0; i < nVrs; i++)
	{
		//m_nNumCommonSat = SelCommonSat(gTrinetInfo[i].ui32Id);
		master = SelMastReferSta(gTrinetInfo->ui32Id, rov1, rov2);

		SelBlineIndex(master, bline1, bline2, flag);
		coeff = LCM(master, rov1, rov2, gVrsInfo[i]);
		natom = GenVrsAtom(coeff, bline1, bline2, flag);
		GenVRSObs(pSinstaobsdata[master], gVrsInfo[i], master, natom);

		if (m_pVrscorrect != NULL)
		{
			delete[]m_pVrscorrect;
			m_pVrscorrect = NULL;
		}
		

	}
}