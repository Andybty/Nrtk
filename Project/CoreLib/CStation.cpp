#include <direct.h>
#include "CStation.h"
#define elevThreshold 7.0*D2R; 



CStation::CStation()
{

}


CStation::~CStation()
{

}

void CStation::SetObsEphInfo(SinStaObsData* pSinstaobsdata, AllBroadcastEphData* pAllbroadcastephdata)
{
	if ((NULL != pSinstaobsdata) && (NULL != pAllbroadcastephdata) &&
		(pSinstaobsdata->ui16ObsNum > 0) && (pAllbroadcastephdata->nSatNum) > 0)
	{
		m_pAllbroadcastephdata = pAllbroadcastephdata;
		m_pSinstaobsdata = pSinstaobsdata;
	}
	else
		return;
}
/**
 @brief  计算站星距
 @brief  rs[in]卫星位置
 @brief  rr[in]测站位置
 @brief  e[in]测站与观测卫星的矢量矩阵
 返回值  位置与测站的几何距离
 @notes    进行了sagnac效应的修正
*/
double  CStation::CalGeoDist(double* rs, const double* rr, double* e)
{
	double r;
	int i;

	if (CMath::norm(rs, 3) < RE_WGS84) return -1.0;
	for (i = 0; i < 3; i++) e[i] = rs[i] - rr[i];
	r = CMath::norm(e, 3);
	for (i = 0; i < 3; i++) e[i] /= r;
	//返回sagnac效应修正的站星距
	return r + OMGE * (rs[0] * rr[1] - rs[1] * rr[0]) / CLIGHT;
}

/**
@brief  计算卫星方位角和高度角
@brief  pos[in]卫星大地坐标
@brief  e[in]测站与观测卫星的单位向量矩阵
@brief  azel[out]卫星高度角方位角
@返回值 卫星高度角 （rad）
*/
double  CStation::CalSatAzEl( double* pos,  double* e, double* azel)
{
	double az = 0.0, el = PI / 2.0, enu[3];

	if (pos[2] > -RE_WGS84) {
		CCoordTrans::ecef2enu(pos, e, enu);
		//计算卫星方位角
		az = CMath::dot(enu, enu, 2) < 1E-12 ? 0.0 : atan2(enu[0], enu[1]);
		if (az < 0.0) az += 2 * PI;
		//计算卫星高度角
		el = asin(enu[2]);
	}
	if (azel) { azel[0] = az; azel[1] = el; }
	return el;
}

/**
@brief  计算卫星非差相位和伪距残差
@param  pObsData[in]观测卫星数据
@param  r[in]站星几何距离
@brief  azel[in]卫星高度角方位角
@brief  y_p[out]伪距残差
@brief  y_l[out]相位残差
@param  freq[in]观测卫星频点
@返回值
*/
void CStation::CalZeroSatRes(double r, ObsData* pObsData, double* azel,
                             double* y_l, double*y_p,double* freq)
{

	int i;
	for (i = 0; i < NFREQ; i++)
	{
		if ((freq[i] = Sat2Freq(pObsData->ui8SatId, pObsData->ui8aCodeType[i])) == 0.0) continue;

		//伪距残差= 观测值减去改正后的站星几何距离
		if (pObsData->daCarrPhase[i] != 0.0) y_l[i] = pObsData->daCarrPhase[i] * CLIGHT / freq[i] - r ;
		if (pObsData->daPseRange[i] != 0.0) y_p[i] = pObsData->daPseRange[i] - r;
	}
}

bool CStation::TestSNRatio()
{
	return true;
}

void CStation::ObsLineComb()
{

}
/**
@brief  处理单站观测数据
@主要进行非差残差处理
@返回值 处理结束返回true,否则返回false
*/
bool CStation::Process()
{
	int n = m_pSinstaobsdata->ui16ObsNum;
	
	SatSoln* pSatsoln = gStaSolveInfo[0].pSatSlon;
	ObsData* pobsdata = m_pSinstaobsdata->pObsData;

	double rr_[3];
	double pos[3];
	double* rs, * dts, * var,*freq;
	double* e;
	int svh[MAXOBS * 2];
	int i;
	rs = CMath::mat(3, n); 
	dts = CMath::mat(2, n); 
	var = CMath::mat(1, n);
	e  = CMath::mat(3, n); 
	freq = CMath::zeros(NFREQ, n);

	rr_[0] = gRStaInfo->daCoorVal[0];
	rr_[1] = gRStaInfo->daCoorVal[1];
	rr_[2] = gRStaInfo->daCoorVal[2];

	CCoordTrans::ecef2pos(rr_, pos);
	//计算卫星位置
	CBroadcastEph Cbroadcasteph;
	Cbroadcasteph.CalSatPos(m_pSinstaobsdata->pObsData[0].tTimeStamp,&m_pSinstaobsdata[0],
	                            m_pAllbroadcastephdata,rs,dts,var,svh);
	for (i = 0; i < n; i++)
	{   
		if ((pSatsoln[i].dSta2SatDis = CalGeoDist(rs + i * 3, rr_, e + i * 3)) <= 0.0) continue;
		if (CalSatAzEl(pos, e + i * 3, pSatsoln[i].daAzim) < 7.0 * D2R) continue;
		pSatsoln[i].dSta2SatDis+= -CLIGHT * dts[i * 2];
		//计算单站残差,可以加模型进行改正
		CalZeroSatRes(pSatsoln[i].dSta2SatDis, &pobsdata[i], pSatsoln[i].daAzim,  
		              pSatsoln[i].daCarrPhaseRes, pSatsoln[i].daPseuRangeRes, freq);

	}
	free(rs); free(dts); free(var); free(e);free(freq);

	return true;
}
/**
@brief  处理n个观测数据
@主要进行非差残差处理
@返回值 处理结束返回true,否则返回false
*/
bool CStation::Process(int n)
{
	double rr_[3];
	double pos[3];
	int n_sta;
	int j;
	int k;
	CBroadcastEph Cbroadcasteph;
	for (n_sta = 0; n_sta < n; n_sta++)
	{
		int nobs = m_pSinstaobsdata[n_sta].ui16ObsNum;

		SatSoln* pSatsoln = gStaSolveInfo[n_sta].pSatSlon;
		ObsData* pobsdata = m_pSinstaobsdata[n_sta].pObsData;
		double* rs, * dts, * var, * freq;
		double* e;
		int svh[MAXOBS * 2];
		int i;
		k = 0;
		//对卫星位置钟差方差和矩阵初始化
		rs = CMath::mat(3, nobs); 
		dts = CMath::mat(2, nobs); 
		var = CMath::mat(1, nobs);
		e = CMath::mat(3, nobs); 
		freq = CMath::zeros(NFREQ, nobs);
		//测站坐标初始化赋值
		rr_[0] = gRStaInfo[n_sta].daCoorVal[0];
		rr_[1] = gRStaInfo[n_sta].daCoorVal[1];
		rr_[2] = gRStaInfo[n_sta].daCoorVal[2];

		CCoordTrans::ecef2pos(rr_, pos);
		//计算卫星位置
		Cbroadcasteph.CalSatPos(m_pSinstaobsdata[n_sta].pObsData[0].tTimeStamp, 
		                        &m_pSinstaobsdata[n_sta], m_pAllbroadcastephdata, 
								 rs, dts, var, svh);
		for (i = 0; i < nobs; i++)
		{
			if (0.0 == rs[i * 3])
			{
				k++;
				continue;
			}
			if (svh[i])
			{
				k++;
				continue;
			}
			if (((pSatsoln + i-k)->dSta2SatDis = CalGeoDist(rs + i * 3, rr_, e + i * 3)) <= 0.0) continue;
			if (CalSatAzEl(pos, e + i * 3, (pSatsoln + i-k)->daAzim) < 7.0 * D2R) continue;
			(pSatsoln + i-k)->dSta2SatDis += -CLIGHT * dts[i * 2];
			//double dd = CLIGHT * dts[i * 2];
			//计算单站残差,可以加模型进行改正
			CalZeroSatRes((pSatsoln + i-k)->dSta2SatDis, &pobsdata[i], (pSatsoln + i-k)->daAzim,
			              (pSatsoln + i-k)->daCarrPhaseRes, (pSatsoln + i-k)->daPseuRangeRes, freq);



			// 全局变量基站解算信息的存储
			gStaSolveInfo[n_sta].ui8SatNum++;
			(pSatsoln + i-k)->ui8SatId = pobsdata[i].ui8SatId;
			(pSatsoln + i-k)->dSatClk = dts[i * 2]*CLIGHT;
			(pSatsoln + i - k)->uiObsindex = i;
			for (j = 0; j < 3; j++)
			{
				(pSatsoln + i-k)->daSatPos[j] = rs[j + i * 3];
				(pSatsoln + i - k)->e[j] = e[j + i * 3];
			}
		}
		free(rs);free(dts); free(var);free(e); free(freq);
	}

	//通过处理站的个数返回基线的条数
	return true;
}
