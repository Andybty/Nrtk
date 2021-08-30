#include <direct.h>
#include "CBline.h"
CBline::CBline()
{

}

CBline::~CBline()
{

}
/**
 @brief  计算双差残差、系数矩阵和权矩阵
 @param  ObsData[in]观测数据
 @param  x[in/out]前后滤波状态参数
 @param  P[in/out]前后滤波方差参数
 @param  xa[in/out]模糊度固定状态参数
 @param  Pa[in/out]模糊度固定方差参数
 @param  y[in/out]单差残差
 @param  v[in/out]双差残差
 @param  azel[in/out]卫星高度角方位角
 @param  freq[in/out]卫星观测频点
 @param  StaIndex1[in]基站索引
 @返回值
*/
int CBline::CalDouDiffRes(ObsData* pSta1ObsData, ObsData* pSta2ObsData, double dt, const double* x,
	                      const double* P, const int* sat, double* y, double* e, double* azel, 
	                       double* freq, const int* iu, const int* ir, int ns, double* v, double* H,
	                       double* R, int* vflg,int index, uint32 StaIndex1, uint32  StaIndex2)
{
	int b,i,j,m,f,k,nf,nv;
	int sysi,sysj;
	double bl;
	double dr[3], posu[3], posr[3];
	double freqi, freqj;
	double didxi, didxj;
	double* dtdxr, * dtdxu;
	double* Ri, * Rj, * Hi, * im;
	double* tropr, * tropu;
	b = 0;
	nv = 0;
	Hi = NULL;
	didxi = 0.0;
	didxj = 0.0;
	nf = m_Rtkopt.nf;
	int nb[NFREQ * 4 * 2 + 2] = { 0 };
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam= m_Startksolpara.rtksolpara;

	bl = Baseline(gRStaInfo[StaIndex1].daCoorVal, gRStaInfo[StaIndex2].daCoorVal, dr);
	CCoordTrans::ecef2pos(gRStaInfo[StaIndex1].daCoorVal, posu);
	CCoordTrans::ecef2pos(gRStaInfo[StaIndex2].daCoorVal, posr);

	Ri = CMath::mat(ns * nf * 2 + 2, 1);
	Rj = CMath::mat(ns * nf * 2 + 2, 1);
	im = CMath::mat(ns, 1);
	tropu = CMath::mat(ns, 1);
	tropr = CMath::mat(ns, 1); 
	dtdxu = CMath::mat(ns, 3); 
	dtdxr = CMath::mat(ns, 3);

	for (i = 0; i < MAXSAT; i++) 
	{
		for (j = 0; j < NFREQ; j++)
		{
			m_Ssat[i].resp[j] = m_Ssat[i].resc[j] = 0.0;
		}
		
	}
	//计算电离层和对流层延迟
	/* compute factors of ionospheric and tropospheric delay */
	for (i = 0; i < ns; i++) 
	{
		if (m_Rtkopt.ionoopt >= IONOOPT_EST) 
		{
			//im[i] = (ionmapf(posu, azel + iu[i] * 2) + ionmapf(posr, azel + ir[i] * 2)) / 2.0;
			im[i] = 0.0;
		}
		if (m_Rtkopt.tropopt >= TROPOPT_EST) 
		{
			//tropu[i] = prectrop(rtk->sol.time, posu, 0, azel + iu[i] * 2, opt, x, dtdxu + i * 3);
			//tropr[i] = prectrop(rtk->sol.time, posr, 1, azel + ir[i] * 2, opt, x, dtdxr + i * 3);
			tropu[i] = 0.0;
			tropr[i] = 0.0;
		}
	}
	for (m = 0; m < 1; m++)
	{
		for (f = m_Rtkopt.mode > PMODE_DGPS ? 0 : nf; f < nf * 2; f++)
		{

			//根据最大高度角计算参考卫星
			/* search reference satellite with highest elevation */
			for (i = -1, j = 0; j < ns; j++)
			{
				sysi = m_Ssat[sat[j] - 1].sys;
				if (!Test_sys(sysi, m)) continue;
				//if (!validobs(iu[j], ir[j], f, nf, y)) continue;
				//if (i < 0 || azel[1 + iu[j] * 2] >= azel[1 + iu[i] * 2]) i = j;
				if (i<0 || gStaSolveInfo[StaIndex1].pSatSlon[iu[j]].daAzim[1]>gStaSolveInfo[StaIndex1].pSatSlon[iu[i]].daAzim[1]) i = j;

			}
			if (i < 0) continue;

			/* make DD (double difference) */
			for (j = 0; j < ns; j++)
			{
				if (i == j) continue;
				sysi = m_Ssat[sat[i] - 1].sys;
				sysj = m_Ssat[sat[j] - 1].sys;
				if (sysi != sysj) continue;
				if ((freqi = Sat2Freq(pSta2ObsData[iu[i]].ui8SatId, pSta1ObsData[iu[i]].ui8aCodeType[f%2])) == 0.0)continue;
				if ((freqj = Sat2Freq(pSta2ObsData[iu[j]].ui8SatId, pSta1ObsData[iu[j]].ui8aCodeType[f%2])) == 0.0)continue;
				if (H)
				{
					Hi = H + nv * m_RtkSolveparam[index].nx;
					for (k = 0; k < m_RtkSolveparam[index].nx; k++) Hi[k] = 0.0;
				}
				//计算双差残差
				if (f < nf)
				{
					v[nv] = gStaSolveInfo[StaIndex1].pSatSlon[iu[i]].daCarrPhaseRes[f % 2] 
					      - gStaSolveInfo[StaIndex2].pSatSlon[ir[i]].daCarrPhaseRes[f % 2] 
						  - (gStaSolveInfo[StaIndex1].pSatSlon[iu[j]].daCarrPhaseRes[f % 2]
						  - gStaSolveInfo[StaIndex2].pSatSlon[ir[j]].daCarrPhaseRes[f % 2]);
				}
				else
				{
					v[nv] = gStaSolveInfo[StaIndex1].pSatSlon[iu[i]].daPseuRangeRes[f % 2]
						  - gStaSolveInfo[StaIndex2].pSatSlon[ir[i]].daPseuRangeRes[f % 2]
						  - (gStaSolveInfo[StaIndex1].pSatSlon[iu[j]].daPseuRangeRes[f % 2]
						  - gStaSolveInfo[StaIndex2].pSatSlon[ir[j]].daPseuRangeRes[f % 2]);
				}
				/* partial derivatives by rover position */
				/* DD ionospheric delay term */
				if (m_Rtkopt.ionoopt == IONOOPT_EST) 
				{
					didxi = (f < nf ? -1.0 : 1.0) * im[i] * SQR(FREQ1 / freqi);
					didxj = (f < nf ? -1.0 : 1.0) * im[j] * SQR(FREQ1 / freqj);
					v[nv] -= didxi * x[II(sat[i], &m_Rtkopt)] - didxj * x[II(sat[j], &m_Rtkopt)];
					if (H)
					{
						Hi[II(sat[i], &m_Rtkopt)] = didxi;
						Hi[II(sat[j], &m_Rtkopt)] = -didxj;
					}
				}

				/* DD tropospheric delay term */
				if (m_Rtkopt.tropopt == TROPOPT_EST || m_Rtkopt.tropopt == TROPOPT_ESTG) 
				{
					v[nv] -= (tropu[i] - tropu[j]) - (tropr[i] - tropr[j]);
					for (k = 0; k < (m_Rtkopt.tropopt < TROPOPT_ESTG ? 1 : 3); k++) 
					{
						if (!H) continue;
						Hi[IT(0, &m_Rtkopt) + k] = (dtdxu[k + i * 3] - dtdxu[k + j * 3]);
						Hi[IT(1, &m_Rtkopt) + k] = -(dtdxr[k + i * 3] - dtdxr[k + j * 3]);
					}
				}

				//双差相位偏差引起的模糊度更新矩阵H
				/* DD phase-bias term */
				if (f < nf)
				{
					if (m_Rtkopt.ionoopt != IONOOPT_IFLC)
					{
						v[nv] -= CLIGHT / freqi * x[IB(sat[i], f, &m_Rtkopt)] -
							     CLIGHT / freqj * x[IB(sat[j], f, &m_Rtkopt)];
						if (H) 
						{
							Hi[IB(sat[i], f, &m_Rtkopt)] = CLIGHT / freqi;
							Hi[IB(sat[j], f, &m_Rtkopt)] = -CLIGHT / freqj;
						}
					}
					else
					{
						v[nv] -= x[IB(sat[i], f, &m_Rtkopt)] - x[IB(sat[j], f, &m_Rtkopt)];
						if (H) 
						{
							Hi[IB(sat[i], f, &m_Rtkopt)] = 1.0;
							Hi[IB(sat[j], f, &m_Rtkopt)] = -1.0;
						}
					}
				}
				if (f < nf) 
					m_Ssat[sat[j] - 1].resc[f] = v[nv];
				else      
					m_Ssat[sat[j] - 1].resp[f - nf] = v[nv];
				/* test innovation */
				if (m_Rtkopt.maxinno > 0.0 && fabs(v[nv]) > m_Rtkopt.maxinno) {
					if (f < nf) {
						m_Ssat[sat[i] - 1].rejc[f]++;
						m_Ssat[sat[j] - 1].rejc[f]++;
					}
					continue;
				}
				//计算单差测量误差方差
				/* SD (single-differenced) measurement error variances */
				Ri[nv] = Varerr(sat[i], sysi, gStaSolveInfo[1].pSatSlon[iu[i]].daAzim[1], bl, dt, f, &m_Rtkopt);
				Rj[nv] = Varerr(sat[j], sysj, gStaSolveInfo[1].pSatSlon[iu[j]].daAzim[1], bl, dt, f, &m_Rtkopt);
				 /* set valid data flags */
				if (m_Rtkopt.mode > PMODE_DGPS) 
				{
					if (f < nf) 
					m_Ssat[sat[i] - 1].vsat[f] = m_Ssat[sat[j] - 1].vsat[f] = 1;
				}
				else 
				{
					m_Ssat[sat[i] - 1].vsat[f - nf] = m_Ssat[sat[j] - 1].vsat[f - nf] = 1;
				}
				vflg[nv++] = (sat[i] << 16) | (sat[j] << 8) | ((f < nf ? 0 : 1) << 4) | (f % nf);
				nb[b]++;
 			}
			b++;
		}
	}
	/* end of system loop */
	//计算双差测量误差协方差
	  /* DD measurement error covariance */
	DDcov(nb, b, Ri, Rj, nv, R);

	free(Ri); free(Rj); free(im);
	free(tropu); free(tropr); free(dtdxu); free(dtdxr);

	return nv;
}
int CBline::AmbiFixByLambda(double* bias, double* xa,int index)
{
	int* ix; 
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam=m_Startksolpara.rtksolpara;
	sol_t* m_Sol = m_Stasol.sol;
	int nx = m_RtkSolveparam[index].nx, na = m_RtkSolveparam[index].na;
	double* DP, * y, * b, * db, * Qb, * Qab, * QQ, s[2];
	int i,j,nb,info;
	//计算单差模糊度转双差模糊度矩阵
	/* index of SD to DD transformation matrix D */
	ix = CMath::imat(nx, 2);
	nb = Ddidx(ix,index);
	if (nb<= 0)
	{
		free(ix);
		return 0;
	}
	y = CMath::mat(nb, 1);
	DP = CMath::mat(nb, nx - na);
	b = CMath::mat(nb, 2); 
	db = CMath::mat(nb, 1); 
	Qb = CMath::mat(nb, nb);
	Qab = CMath::mat(na, nb); 
	QQ = CMath::mat(na, nb);

	/* y=D*xc, Qb=D*Qc*D', Qab=Qac*D' */
	for (i = 0; i < nb; i++) 
	{
		y[i] = m_RtkSolveparam[index].x[ix[i * 2]] - m_RtkSolveparam[index].x[ix[i * 2 + 1]];
	}

	for (j = 0; j < nx - na; j++)
	{
		for (i = 0; i < nb; i++)
		{
			DP[i + j * nb] =  m_RtkSolveparam[index].P[ix[i*2]+(na+j)*nx]- 
				              m_RtkSolveparam[index].P[ix[i*2+1]+(na+j)*nx];
		}
	}

	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
		{
			Qb[i + j * nb] = DP[i + (ix[j * 2] - na) * nb] - DP[i + (ix[j * 2 + 1] - na) * nb];
		}
	}

	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < na; i++)
		{
			Qab[i + j * na] = m_RtkSolveparam[index].P[i + ix[j * 2] * nx] - 
				              m_RtkSolveparam[index].P[i + ix[j * 2 + 1] * nx];
		}
	}
	//LAMBDA最小二乘估计
	/* LAMBDA/MLAMBDA ILS (integer least-square) estimation */
	if (!(info = Lambda(nb, 2, y, Qb, b, s)))
	{
		m_Sol[index].ratio=s[0]>0? (float)(s[1] / s[0]) : 0.0f;
		if (m_Sol[index].ratio>999.9)
		{
			m_Sol[index].ratio = 999.9f;
		}
		/* validation by popular ratio-test */
		if (s[0] <= 0.0 || s[1] / s[0] >= m_Rtkopt.thresar[0])
		{
			/* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
			for (i = 0; i < na; i++)
			{
				m_RtkSolveparam[index].xa[i] = m_RtkSolveparam[index].x[i];
				for (j = 0; j < na; j++)
				{
					m_RtkSolveparam[index].Pa[i + j * na] = m_RtkSolveparam[index].P[i + j * nx];
				}
			}

			for (i = 0; i < nb; i++)
			{
				bias[i] = b[i];
				y[i] -= b[i];
			}

			if (!CMath::matinv(Qb, nb))
			{
				CMath::matmul("NN", nb, 1, nb, 1.0, Qb, y, 0.0, db);
				CMath::matmul("NN", na, 1, nb, -1.0, Qab, db, 1.0, m_RtkSolveparam[index].xa);

				/* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
				CMath::matmul("NN", na, nb, nb, 1.0, Qab, Qb, 0.0, QQ);
				CMath::matmul("NT", na, na, nb, -1.0, QQ, Qab, 1.0, m_RtkSolveparam[index].Pa);

			  /* restore SD ambiguity */
				Restamb(bias,nb,xa,index);
			}
			else
			{
				nb = 0;
			}
		}
		else
		{
			nb = 0;
		}
	}
	else
	{
		nb = 0;
	}
	free(ix); free(y); free(DP); free(b); 
	free(db); free(Qb); free(Qab); free(QQ);

	return nb;
}
/**
 @brief  计算双差大气（电离层和对流层或综合误差）
 @param  ObsData[in]基站观测数据
 @param  sat[in]观测卫星ID
 @param  Index1[in]基站数量
 @param  StaIndex1[in]基站索引
 @返回值 观测卫星数量
*/
int CBline::CalDouDiffAtmo(ObsData* pSta1ObsData, ObsData* pSta2ObsData, gtime_t dt, const double* x,
	const int* sat, const int* iu, const int* ir, int ns, int index, uint32 StaIndex1, uint32  StaIndex2)
{
	int f, m, i, j, nf;
	int sysi, sysj;
	int nv = 0;
	double bl;
	double dr[3];
	double freqi, freqj;
	double temp, tempn;
	nf = m_Rtkopt.nf;
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam = m_Startksolpara.rtksolpara;

	gBlineSlon[index].ui32Id = index;
	gBlineSlon[index].ui8ComnSatNum = ns;
	gBlineSlon[index].Ddist = Baseline(gRStaInfo[StaIndex1].daCoorVal, gRStaInfo[StaIndex2].daCoorVal, dr);
	gBlineSlon[index].tTimeStamp = dt;
	for (m = 0; m < NSYS; m++)
	{
		for (f = 0; f < nf; f++)
		{
			/* search reference satellite with highest elevation */
			for (i = -1, j = 0; j < ns; j++)
			{
				sysi = m_Ssat[sat[j] - 1].sys;
				if (!Test_sys(sysi, m)) continue;
				//if (!validobs(iu[j], ir[j], f, nf, y)) continue;
				//if (i < 0 || azel[1 + iu[j] * 2] >= azel[1 + iu[i] * 2]) i = j;
				if (i<0 || gStaSolveInfo[StaIndex1].pSatSlon[iu[j]].daAzim[1]>gStaSolveInfo[StaIndex1].pSatSlon[iu[i]].daAzim[1]) i = j;

			}
			if (i < 0) continue;
			for (j = 0; j < ns; j++)
			{
				temp = 0.0;
				//if (i == j) continue;
				sysi = m_Ssat[sat[i] - 1].sys;
				sysj = m_Ssat[sat[j] - 1].sys;
				if (sysi != sysj) continue;
				if ((freqi = Sat2Freq(pSta2ObsData[iu[i]].ui8SatId, pSta1ObsData[iu[i]].ui8aCodeType[f])) == 0.0)continue;
				if ((freqj = Sat2Freq(pSta2ObsData[iu[j]].ui8SatId, pSta1ObsData[iu[j]].ui8aCodeType[f])) == 0.0)continue;
				if (sysi != sysj || fabs(freqi - freqj) > 1e-6) continue;
				gBlineSlon[index].pComnSatSoln[nv].CodeType = pSta1ObsData[iu[j]].ui8aCodeType[f];
				gBlineSlon[index].pComnSatSoln[nv].dCCombinErr = gStaSolveInfo[StaIndex1].pSatSlon[iu[i]].daCarrPhaseRes[f] - gStaSolveInfo[StaIndex2].pSatSlon[ir[i]].daCarrPhaseRes[f] -
					(gStaSolveInfo[StaIndex1].pSatSlon[iu[j]].daCarrPhaseRes[f] - gStaSolveInfo[StaIndex2].pSatSlon[ir[j]].daCarrPhaseRes[f]);

				gBlineSlon[index].pComnSatSoln[nv].dPCombinErr = gStaSolveInfo[StaIndex1].pSatSlon[iu[i]].daPseuRangeRes[f] - gStaSolveInfo[StaIndex2].pSatSlon[ir[i]].daPseuRangeRes[f] -
					(gStaSolveInfo[StaIndex1].pSatSlon[iu[j]].daPseuRangeRes[f] - gStaSolveInfo[StaIndex2].pSatSlon[ir[j]].daPseuRangeRes[f]);


				temp = x[IB(sat[i], f, &m_Rtkopt)] - x[IB(sat[j], f, &m_Rtkopt)];
				tempn = ROUND(temp);
				gBlineSlon[index].pComnSatSoln[nv].dCCombinErr -= CLIGHT / freqi * tempn;
				gBlineSlon[index].pComnSatSoln[nv].ui8RefSatId = sat[i];
				gBlineSlon[index].pComnSatSoln[nv++].ui8RovSatId = sat[j];
			}
		}
	}
	gBlineSlon[index].ui8EquNum = nv;
	return nv;

}
/**
 @brief  选择基线共视卫星
 @param  StaSolveInfo[in]观测的卫星信息
 @param  StaIndex1/StaIndex2[in]基站索引信息
 @param  BLineInfo[in]基线索引信息
 @返回值 基线共视卫星数量
*/
int CBline::SelCommonSat(StaSolveInfo* pStaSolveInfo, BLineInfo bLineinfo, uint32 &StaIndex1, uint32 &StaIndex2)
{
	StaIndex1 = bLineinfo.ui32StartPointId - 65;
	StaIndex2 = bLineinfo.ui32EndPointId - 65;
	int nStaSatNum1 = pStaSolveInfo[StaIndex1].ui8SatNum;
	int nStaSatNum2 = pStaSolveInfo[StaIndex2].ui8SatNum;
	SatSoln* psatsoln1 = pStaSolveInfo[StaIndex1].pSatSlon;
	SatSoln* psatsoln2 = pStaSolveInfo[StaIndex2].pSatSlon;
	
	int i, j, k = 0;
	for (i = 0, j = 0; i < nStaSatNum1 && j < nStaSatNum2; i++, j++)
	{
		if (psatsoln1[i].ui8SatId < psatsoln2[j].ui8SatId) j--;
		else if (psatsoln1[i].ui8SatId > psatsoln2[j].ui8SatId) i--;
		else
		{
			m_naCommonSat[k] = psatsoln1[i].ui8SatId;
			m_naSta1SatID[k] = i;
			m_naSta2SatID[k++] = j;
		}
	}
	m_numcommonsat = k;
	
	return k;
}
/**
 @brief 更新模糊度参数
 @param pSta1ObsData[in]观测数据
 @param  tt[in]当前历元和前历元的时间差
 @param  ns[in]观测数据有效卫星数
 @param  sat[in]观测数据卫星ID
 @param  index[in]观测数据的基站数
 @返回值 无
*/
void CBline::UpBias(ObsData* pSta1ObsData, ObsData* pSta2ObsData, 
	        double tt,int ns, int* sat, int* iu, int* ir,int index)
{
	int nf = m_Rtkopt.nf;
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam;
	m_RtkSolveparam = m_Startksolpara.rtksolpara;
	int k,i,j;
	int slip;
	bool reset;
	double* bias,offset;
	double cp, pr, cp1, cp2, pr1,pr2;
	double freqi, freq1, freq2, C1, C2;
	for (k = 0; k < nf; k++)
	{
		for (i = 1; i <= MAXSAT; i++)
		{
			reset = ++m_Ssat[i - 1].outc[k] > (uint32_t)m_Rtkopt.maxout;

			if (m_Rtkopt.modear == ARMODE_INST && m_RtkSolveparam[index].x[IB(i, k, &m_Rtkopt)] != 0.0)
			{
				InitParam(&m_RtkSolveparam[index], 0.0, 0.0, IB(i, k, &m_Rtkopt));
			}
			else if (reset && m_RtkSolveparam[index].x[IB(i, k, &m_Rtkopt)] != 0.0)
			{
				InitParam(&m_RtkSolveparam[index], 0.0, 0.0, IB(i, k, &m_Rtkopt));
				m_Ssat[i - 1].outc[k] = 0;
			}
			if (m_Rtkopt.modear != ARMODE_INST && reset)
			{
				m_Ssat[i - 1].lock[k] =- m_Rtkopt.minlock;
			}
		
		}
		/* reset phase-bias if detecting cycle slip */
		for (i = 0; i < ns; i++)
		{
			j = IB(sat[i], k, &m_Rtkopt);
			m_RtkSolveparam[index].P[j + j * m_RtkSolveparam[index].nx] += m_Rtkopt.prn[0] * m_Rtkopt.prn[0] * fabs(tt);
			slip = m_Ssat[sat[i] - 1].slip[k];
			if (IONOOPT_IFLC == m_Rtkopt.ionoopt)
				/*默认的是L1 L2的消电离层组合*/
				slip |= m_Ssat[sat[i] - 1].slip[1]; 
			if (ARMODE_INST == m_Rtkopt.modear || !(slip & 1))
				continue;
			m_RtkSolveparam[index].x[j] = 0.0;
			m_Ssat[sat[i] - 1].lock[k] = -m_Rtkopt.minlock;

		}
		bias = CMath::zeros(ns, 1);


		/* estimate approximate phase-bias by phase - code */
		for (i = j = 0, offset = 0.0; i < ns; i++)
		{
			if (m_Rtkopt.ionoopt != IONOOPT_IFLC)
			{
				cp=SDobs(pSta1ObsData, pSta2ObsData, iu[i], ir[i], k); /*cycle*/
				pr = SDobs(pSta1ObsData, pSta2ObsData, iu[i], ir[i], k+ NFREQ);
				freqi = Sat2Freq(sat[i], pSta1ObsData[iu[i]].ui8aCodeType[k]);

				bias[i] = cp - pr * freqi / CLIGHT;
			}
			else
			{
				cp1 = SDobs(pSta1ObsData, pSta2ObsData, iu[i], ir[i], 0);
				cp2 = SDobs(pSta1ObsData, pSta2ObsData, iu[i], ir[i], 1);
				pr1= SDobs(pSta1ObsData, pSta2ObsData, iu[i], ir[i], NFREQ);
				pr2 = SDobs(pSta1ObsData, pSta2ObsData, iu[i], ir[i], NFREQ+1);
				freq1 = Sat2Freq(sat[i], pSta1ObsData[iu[i]].ui8aCodeType[0]);
				freq2 = Sat2Freq(sat[i], pSta1ObsData[iu[i]].ui8aCodeType[1]);
				if (cp1 == 0.0 || cp2 == 0.0 || pr1 == 0.0 || pr2 == 0.0 || freq1 == 0.0 || freq2 <= 0.0) continue;

				C1 = SQR(freq1) / (SQR(freq1) - SQR(freq2));
				C2 = -SQR(freq2) / (SQR(freq1) - SQR(freq2));
				bias[i] = (C1 * cp1 * CLIGHT / freq1 + C2 * cp2 * CLIGHT / freq2) - (C1 * pr1 + C2 * pr2);
			}
			if (m_RtkSolveparam[index].x[IB(sat[i], k, &m_Rtkopt)] != 0.0)
			{
				offset += bias[i] - m_RtkSolveparam[index].x[IB(sat[i], k, &m_Rtkopt)];
				j++;
			}
		}

		/* correct phase-bias offset to enssure phase-code coherency */
		if (j > 0)
		{
			for (i = 1; i <= MAXSAT; i++)
			{
				if (m_RtkSolveparam[index].x[IB(i, k, &m_Rtkopt)] != 0.0)
					m_RtkSolveparam[index].x[IB(i, k, &m_Rtkopt)] += offset / j;
			}
		}
		/* set initial states of phase-bias */
		for (i = 0; i < ns; i++)
		{
			if (bias[i] == 0.0 || m_RtkSolveparam[index].x[IB(sat[i], k, &m_Rtkopt)] != 0.0)
				continue;
			InitParam(&m_RtkSolveparam[index], bias[i], SQR(m_Rtkopt.std[0]), IB(sat[i], k, &m_Rtkopt));
		}
		free(bias);
	}
}

/**
 @brief  状态参数及其方差更新
 @param  ObsData[in]观测数据
 @param  tt[in]当前历元和前历元的时间差
 @param  ns[in]观测数据有效卫星数
 @param  sat[in]观测数据卫星ID
 @param  iu/ir[in]基站卫星ID索引
 @param  index[in]基线的索引
 @返回值 无
*/
void CBline::UpdateState(ObsData* pSta1ObsData, ObsData* pSta2ObsData, double tt, int ns, int* sat,
	int* iu, int* ir,int index)
{
	int i;

	if (m_Rtkopt.ionoopt >= IONOOPT_EST)
	{

	}
	if (m_Rtkopt.tropopt >= TROPOPT_EST)
	{

	}
	if (m_Rtkopt.glomodear == 2 && (m_Rtkopt.navsys & SYS_GLO))
	{

	}
	if (m_Rtkopt.mode > PMODE_DGPS) {
		UpBias(pSta1ObsData, pSta2ObsData,tt,ns,sat,iu,ir,index);
	}
}
int CBline::Filter()
{
	return 0;
}

bool CBline::BLineSolution()
{
	return true;
}
void CBline::SetRtkOpt(int nf, int mode, int ionoopt, int tropopt, int dynamics, int glomodear,int navsys)
{
	
	m_Rtkopt.nf = nf;
	m_Rtkopt.mode = mode;
	m_Rtkopt.ionoopt = ionoopt;
	m_Rtkopt.tropopt = tropopt;
	m_Rtkopt.dynamics = dynamics;
	m_Rtkopt.glomodear = glomodear;
	m_Rtkopt.navsys = navsys;
}

void CBline::InitParam(RtkSolveParam* rtksolveparam, double xi, double var, int i)
{
	int j;
	rtksolveparam->x[i] = xi;
	for (j = 0; j < rtksolveparam->nx; j++)
	{
		rtksolveparam->P[i + j * rtksolveparam->nx] = rtksolveparam->P[j + i * rtksolveparam->nx] = i == j ? var : 0.0;
	}
}
/**
 @brief 单差模糊度初值
 @param ObsData[in]观测数据
 @param i[in]观测到第i颗卫星
 @param j[in]观测到的第j个观测值
 @param k[i]观测频点数量
 @返回值 载波单差初始值
*/
double CBline::SDobs(ObsData* obs1, ObsData* obs2, int i, int j, int k)
{
	double pi = (k < NFREQ) ? obs1[i].daCarrPhase[k] : obs1[i].daPseRange[k - NFREQ];
	double pj = (k < NFREQ) ? obs2[j].daCarrPhase[k] : obs2[j].daPseRange[k - NFREQ];
	return pi == 0.0 || pj == 0.0 ? 0.0 : pi - pj;
}

double CBline::Baseline(const double* ru, const double* rb, double* dr)
{
	int i;
	for (i = 0; i < 3; i++) dr[i] = ru[i] - rb[i];
	return CMath::norm(dr, 3);
}


double CBline::Varerr(int sat, int sys, double el, double bl, double dt, int f,	RtkOpt* opt)
{
	double a, b, c, d, fact, sinel;
	c = opt->err[3] * bl / 1E4;
	d = CLIGHT * opt->sclkstab * dt;
	fact = 1.0;
	sinel = sin(el);
	int nf = opt->nf;
	if (f >= nf) fact = opt->eratio[f - nf];
	if (fact <= 0.0) fact = opt->eratio[0];
	fact *= sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);
	a = fact * opt->err[1];
	b = fact * opt->err[2];
	return 2.0 * (opt->ionoopt == IONOOPT_IFLC ? 3.0 : 1.0) * (a * a + b * b / sinel / sinel + c * c) + d * d;
}
//ddcov(nb, b, Ri, Rj, nv, R);
void CBline:: DDcov(int* nb, int n, double* Ri, double* Rj, int nv, double* R)
{
	int i, j, k, b;
	k = 0;
	for (i = 0; i < nv * nv; i++) R[i] = 0.0;
	for (b = 0; b < n; k += nb[b++])
	{
		for (i = 0; i < nb[b]; i++)
		{
			for (j = 0; j < nb[b]; j++)
			{
				R[k + i + (k + j) * nv] = Ri[k + i] + (i == j ? Rj[k + i] : 0.0);
			}
		}
	}
}
/**
 @brief 单差到双差转换矩阵
 @brief index[in]卫星号索引
 @param ix[out]参考偏差状态矩阵
 @返回值 单差到双差转换矩阵
*/

int CBline::Ddidx( int* ix,int index)
{
	int i, j, k, m, f, nb;
	nb = 0;
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam= m_Startksolpara.rtksolpara;
	int na = m_RtkSolveparam[index].na, nf = m_Rtkopt.nf, nofix;
	for (i = 0; i < MAXSAT; i++)
	{
		for (j = 0; j < NFREQ; j++)
		{
			m_Ssat[i].fix[j] = 0;
		}
	}
	/* m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN */
	for (m = 0; m < NSYS; m++)
	{
		nofix = (m == 1 && m_Rtkopt.glomodear == 0) || (m == 3 && m_Rtkopt.bdsmodear == 0);

		for (f = 0, k = na; f < nf; f++, k += MAXSAT)
		{
			for (i = k; i < k + MAXSAT; i++)
			{
				if (m_RtkSolveparam[index].x[i] == 0 || !m_Ssat[i - k].vsat[f]||!Test_sys(m_Ssat[i - k].sys, m))
				{
					continue;
				}
				if (m_Ssat[i-k].lock[f]> 0 && !(m_Ssat[i - k].slip[f] & 2) &&
					m_Ssat[i - k].azel[1] >= m_Rtkopt.elmaskar && !nofix)
				{
					m_Ssat[i - k].fix[f] = 2;  /*fix*/
					break;
				}
				else
				{
					m_Ssat[i - k].fix[f] = 1;
				}
			}

			for (j = k; j < k + MAXSAT; j++)
			{
				if (i == j || m_RtkSolveparam[index].x[j] == 0.0 || !m_Ssat[j - k].vsat[f])
				{
					continue;
				}
				if (m_Ssat[j - k].lock[f] > 0 && !(m_Ssat[j - k].slip[f] & 2) &&
					m_Ssat[i - k].vsat[f] &&
					m_Ssat[j - k].azel[1] >= m_Rtkopt.elmaskar && !nofix)
				{
					ix[nb * 2] = i; /*state index of ref bias*/
					ix[nb * 2 + 1] = j;
					nb++;
					m_Ssat[j - k].fix[f] = 2;
				}
				else
				{
					m_Ssat[j - k].fix[f] = 1;
				}
			}
		}
	}
	return nb;
}
/**
 @brief lambda法固定矩阵
 @param QF[out] lambda法模糊度固定矩阵
 @返回值 模糊度固定信息
*/
int CBline::Lambda(int n, int m, const double* a, const double* Q, double* F,double* s)
{
	int info;
	double* L, * D, * Z, * z, * E;

	if (n <= 0 || m <= 0) return -1;
	L = CMath::zeros(n, n);
	D = CMath::mat(n, 1);
	Z = CMath::eye(n); 
	z = CMath::mat(n, 1);
	E = CMath::mat(n, m);

	/* LD factorization */
	if (!(info = CMath::LD(n, Q, L, D))) {

		/* lambda reduction */
		Reduction(n, L, D, Z);
		CMath::matmul("TN", n, 1, n, 1.0, Z, a, 0.0, z); /* z=Z'*a */

		/* mlambda search */
		if (!(info = Search(n, m, L, D, z, E, s))) {

			info = Solve("T", Z, E, n, m, F); /* F=Z'\E */
		}
	}
	free(L); free(D); free(Z); free(z); free(E);
	return info;
}

/* permutations --------------------------------------------------------------*/
void CBline::Permutation(int n, double* L, double* D, int j, double del, double* Z)
{
	int k;
	double eta, lam, a0, a1;

	eta = D[j] / del;
	lam = D[j + 1] * L[j + 1 + j * n] / del;
	D[j] = eta * D[j + 1]; D[j + 1] = del;
	for (k = 0; k <= j - 1; k++) {
		a0 = L[j + k * n]; a1 = L[j + 1 + k * n];
		L[j + k * n] = -L[j + 1 + j * n] * a0 + a1;
		L[j + 1 + k * n] = eta * a0 + lam * a1;
	}
	L[j + 1 + j * n] = lam;
	for (k = j + 2; k < n; k++) SWAP(L[k + j * n], L[k + (j + 1) * n]);
	for (k = 0; k < n; k++) SWAP(Z[k + j * n], Z[k + (j + 1) * n]);
}

/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
void CBline::Reduction(int n, double* L, double* D, double* Z)
{
	int i, j, k;
	double del;

	j = n - 2; k = n - 2;
	while (j >= 0) {
		if (j <= k) for (i = j + 1; i < n; i++) Gauss(n, L, Z, i, j);
		del = D[j] + L[j + 1 + j * n] * L[j + 1 + j * n] * D[j + 1];
		if (del + 1E-6 < D[j + 1]) { /* compared considering numerical error */
			Permutation(n, L, D, j, del, Z);
			k = j; j = n - 2;
		}
		else j--;
	}
}

void CBline::Gauss(int n, double* L, double* Z, int i, int j)
{
	int k, mu;

	if ((mu = (int)ROUND(L[i + j * n])) != 0) {
		for (k = i; k < n; k++) L[k + n * j] -= (double)mu * L[k + i * n];
		for (k = 0; k < n; k++) Z[k + n * j] -= (double)mu * Z[k + i * n];
	}
}
/**
 @brief modified lambda(mlambda) search
 @param  修改lambda搜索空间
 @返回值 成功返回0，否则返回-1
*/

int CBline::Search(int n, int m, const double* L, const double* D,
	const double* zs, double* zn, double* s)
{
	int i, j, k, c, nn, imax;
	double newdist, maxdist = 1E99, y;
	double* S = CMath::zeros(n, n);
	double* dist = CMath::mat(n, 1);
	double* zb = CMath::mat(n, 1);
	double* z = CMath::mat(n, 1);
	double* step = CMath::mat(n, 1);
	nn = 0;
	imax = 0;
	k = n - 1;
	dist[k] = 0.0;
	zb[k] = zs[k];
	z[k] = ROUND(zb[k]); y = zb[k] - z[k]; step[k] = SGN(y);
	for (c = 0; c < LOOPMAX; c++) {
		newdist = dist[k] + y * y / D[k];
		if (newdist < maxdist) {
			if (k != 0) {
				dist[--k] = newdist;
				for (i = 0; i <= k; i++)
					S[k + i * n] = S[k + 1 + i * n] + (z[k + 1] - zb[k + 1]) * L[k + 1 + i * n];
				zb[k] = zs[k] + S[k + k * n];
				z[k] = ROUND(zb[k]); y = zb[k] - z[k]; step[k] = SGN(y);
			}
			else {
				if (nn < m) {
					if (nn == 0 || newdist > s[imax]) imax = nn;
					for (i = 0; i < n; i++) zn[i + nn * n] = z[i];
					s[nn++] = newdist;
				}
				else {
					if (newdist < s[imax]) {
						for (i = 0; i < n; i++) zn[i + imax * n] = z[i];
						     s[imax] = newdist;
						for (i = imax = 0; i < m; i++) if (s[imax] < s[i]) 
							 imax = i;
					}
					maxdist = s[imax];
				}
				z[0] += step[0]; y = zb[0] - z[0]; step[0] = -step[0] - SGN(step[0]);
			}
		}
		else {
			if (k == n - 1) break;
			else {
				k++;
				z[k] += step[k]; y = zb[k] - z[k]; step[k] = -step[k] - SGN(step[k]);
			}
		}
	}
	for (i = 0; i < m - 1; i++) { /* sort by s */
		for (j = i + 1; j < m; j++) {
			if (s[i] < s[j]) continue;
			SWAP(s[i], s[j]);
			for (k = 0; k < n; k++) SWAP(zn[k + i * n], zn[k + j * n]);
		}
	}
	free(S); free(dist); free(zb); free(z); free(step);

	if (c >= LOOPMAX) {
		fprintf(stderr, "%s : search loop count overflow\n", __FILE__);
		return -1;
	}
	return 0;
}

int CBline::Solve(const char* tr, const double* A, const double* Y, int n,
	int m, double* X)
{
	double* B = CMath::mat(n, n);
	int info;

	CMath::matcpy(B, A, n, n);
	if (!(info = CMath::matinv(B, n))) CMath::matmul(tr[0] == 'N' ? "NN" : "TN", n, m, n, 1.0, B, Y, 0.0, X);
	free(B);
	return info;
}

/* restore SD (single-differenced) ambiguity ---------------------------------*/
void CBline::Restamb(double* bias, int nb, double* xa,int index)
{
	int i, n, m, f, index1[MAXSAT], nv = 0;
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam = m_Startksolpara.rtksolpara;
	int nf = m_Rtkopt.nf;
	/* init all fixed states to float state values */
	for (i = 0; i < m_RtkSolveparam[index].nx; i++) xa[i] = m_RtkSolveparam[index].x[i]; 
	/* overwrite non phase-bias states with fixed values */
	for (i = 0; i < m_RtkSolveparam[index].na; i++) xa[i] = m_RtkSolveparam[index].xa[i];  

	for (m = 0; m < NSYS; m++)
	{
		for (f = 0; f < nf; f++)
		{
			for (n = i = 0; i < MAXSAT; i++)
			{
				if (m_Ssat[i].fix[f] != 2)
				{
					continue;
				}
				index1[n++] = IB(i + 1, f, &m_Rtkopt);			
			}
			if (n < 2)
				continue;
			xa[index1[0]] = m_RtkSolveparam[index].x[index1[0]];

			for (i = 1; i < n; i++)
			{
				xa[index1[i]] = xa[index1[0]] - bias[nv++];
			}
		}
	}
}

void CBline::Holdamb(double* xa,int index)
{
	double* v, * H, * R;
	int nf = m_Rtkopt.nf;
	int index1[MAXSAT];
	int i,f,m,n,info,nv=0;
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam = m_Startksolpara.rtksolpara;
	sol_t* m_Sol = m_Stasol.sol;
	int nb = m_RtkSolveparam[index].nx - m_RtkSolveparam[index].na;

	v =CMath::mat(nb, 1); H = CMath::zeros(nb, m_RtkSolveparam[index].nx);
	for (m = 0; m < NSYS; m++)
	{
		for (f = 0; f < nf; f++)
		{
			for (n = i = 0; i < MAXSAT; i++)
			{
				if (!Test_sys(m_Ssat[i].sys, m) || m_Ssat[i].fix[f] != 2 ||
					m_Ssat[i].azel[1] < m_Rtkopt.elmaskhold) 
				{
					continue;
				}
				index1[n++] = IB(i + 1, f, &m_Rtkopt);
				m_Ssat[i].fix[f] = 3; /* hold */
			}
			/* constraint to fixed ambiguity */
			for (i = 1; i < n; i++)
			{
				v[nv] = (xa[index1[0]] - xa[index1[i]]) 
					  - (m_RtkSolveparam[index].x[index1[0]]
					  - m_RtkSolveparam[index].x[index1[i]]);

				H[index1[0] + nv * m_RtkSolveparam[index].nx] = 1.0;
				H[index1[i] + nv * m_RtkSolveparam[index].nx] = -1.0;
				nv++;
			}
		}
	}

	if (nv > 0)
	{
		R = CMath::zeros(nv, nv);
		for (i = 0; i < nv; i++) R[i + i * nv] = VAR_HOLDAMB;

		if ((info = CMath::filter(m_RtkSolveparam[index].x, m_RtkSolveparam[index].P, H, v, R, m_RtkSolveparam[index].nx, nv)))
		{
			//errmsg(rtk, "filter error (info=%d)\n", info);
		}
		free(R);
	}
	free(v); free(H);
}

int CBline::Test_sys(int sys, int m)
{
	switch (sys) {
	case SYS_GPS: return m == 0;
	case SYS_SBS: return m == 0;
	case SYS_GLO: return m == 1;
	case SYS_GAL: return m == 2;
	case SYS_CMP: return m == 3;
	case SYS_QZS: return m == 4;
	case SYS_IRN: return m == 5;
	}
	return 0;
}
/**
 @brief 保存基线解算参数
 @param stat[in]基线解算状态信息
 @param index[in]基线索引信息
 @返回值 无
*/
void CBline::SaveSoluStaus(int stat,int index)
{
	int i;
	RtkSolveParam* m_RtkSolveparam = m_Startksolpara.rtksolpara;
	sol_t* m_Sol = m_Stasol.sol;
	if (stat == SOLQ_FIX)
	{
		for (i = 0; i < 3; i++)
		{
			m_Sol[index].rr[i] = m_RtkSolveparam[index].xa[i];
			m_Sol[index].qr[i]= (float)m_RtkSolveparam[index].Pa[i + i * m_RtkSolveparam[index].na];
		}
		m_Sol[index].qr[3] = (float)m_RtkSolveparam[index].Pa[1];
		m_Sol[index].qr[4] = (float)m_RtkSolveparam[index].Pa[1 + 2 * m_RtkSolveparam[index].na];
		m_Sol[index].qr[5] = (float)m_RtkSolveparam[index].Pa[2];

		if (m_Rtkopt.dynamics)
		{
			for (i = 3; i < 6; i++)
			{
				m_Sol[index].rr[i] = m_RtkSolveparam[index].xa[i];
				m_Sol[index].qv[i - 3] = (float)m_RtkSolveparam[index].Pa[i + i * m_RtkSolveparam[index].na];
			}
		}
	}
	else
	{
		for (i = 0; i < 3; i++)
		{
			m_Sol[index].rr[i] = m_RtkSolveparam[index].x[i];
			m_Sol[index].qr[i] = (float)m_RtkSolveparam[index].P[i + i * m_RtkSolveparam[index].nx];
		}
		if (m_Rtkopt.dynamics)
		{
			for (i = 3; i < 6; i++)
			{
				m_Sol[index].rr[i] = m_RtkSolveparam[index].x[i];
				m_Sol[index].qv[i - 3] = (float)m_RtkSolveparam[index].P[i + i * m_RtkSolveparam[index].na];
			}
			m_Sol[index].qv[3] = (float)m_RtkSolveparam[index].P[4 + 3 * m_RtkSolveparam[index].nx];
			m_Sol[index].qv[4] = (float)m_RtkSolveparam[index].P[5 + 4 * m_RtkSolveparam[index].nx];
			m_Sol[index].qv[5] = (float)m_RtkSolveparam[index].P[5 + 3 * m_RtkSolveparam[index].nx];
		}
		m_RtkSolveparam[index].nfix = 0;
	}

}

int CBline::Process(SinStaObsData* pSinstaobsdata, BLineInfo bLineinfo,int index)
{
	int f,i,j,ns,nv;
	uint32 StaIndex1, StaIndex2;
	double tt=1.0;
	double* xp, * Pp, * xa;
	double* H, * v, * R, * bias;
	int vflg[MAXOBS * NFREQ * 2 + 1];
	int ny,info;
	int nf = m_Rtkopt.nf;
	int stat = m_Rtkopt.mode <= PMODE_DGPS ? SOLQ_DGPS : SOLQ_FLOAT;
	double dt = 0.0;
	int niter = 1;         //__________________________________________________________
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	RtkSolveParam* m_RtkSolveparam = m_Startksolpara.rtksolpara;
	sol_t* m_Sol=m_Stasol.sol;
	//_______________________________________________________
	ns = SelCommonSat(gStaSolveInfo, gBlineInfo[index], StaIndex1, StaIndex2);
	ObsData* pSta1ObsData = pSinstaobsdata[StaIndex1].pObsData;
	ObsData* pSta2ObsData = pSinstaobsdata[StaIndex2].pObsData;

	UpdateState( pSta1ObsData,pSta2ObsData,tt, ns, m_naCommonSat, m_naSta1SatID, m_naSta2SatID,index);
	//________
	xp = CMath::mat(m_RtkSolveparam[index].nx, 1); 
	Pp = CMath::zeros(m_RtkSolveparam[index].nx,m_RtkSolveparam[index].nx); 
	xa = CMath::mat(m_RtkSolveparam[index].nx, 1);
	CMath::matcpy(xp, m_RtkSolveparam[index].x, m_RtkSolveparam[index].nx, 1);

	ny= ns * nf * 2 + 2;
	v = CMath::mat(ny, 1); 
	H = CMath::zeros(m_RtkSolveparam[index].nx, ny); 
	R = CMath::mat(ny, ny); 
	bias = CMath::mat(m_RtkSolveparam[index].nx, 1);
	for (i = 0; i < MAXSAT; i++) {
		m_Ssat[i].sys = SatSys(i + 1, NULL);
		for (j = 0; j < NFREQ; j++) m_Ssat[i].vsat[j] = 0;
		for (j = 1; j < NFREQ; j++) m_Ssat[i].snr[j] = 0;
	}
	for (i = 0; i < niter; i++)
	{
		/*CalDouDiffRes(double dt, const double*/
		if ((nv = CalDouDiffRes(pSta1ObsData,pSta2ObsData,dt, xp, Pp, m_naCommonSat, NULL, NULL, NULL, NULL, m_naSta1SatID, m_naSta2SatID, ns, v, H, R,
			vflg,index, StaIndex1, StaIndex2)) < 1)
		{
			//errmsg(rtk, "no double-differenced residual\n");
			stat = SOLQ_NONE;
			break;
		}
		/* Kalman filter measurement update */
		CMath::matcpy(Pp, m_RtkSolveparam[index].P, m_RtkSolveparam[index].nx, m_RtkSolveparam[index].nx);
		if ((info = CMath::filter(xp, Pp, H, v, R, m_RtkSolveparam[index].nx, nv)))
		{
			//errmsg(rtk, "filter error (info=%d)\n", info);
			stat = SOLQ_NONE;
			break;
		}
	}
	if (stat != SOLQ_NONE)
	{
		nv = CalDouDiffRes(pSta1ObsData, pSta2ObsData, dt, xp, Pp, m_naCommonSat, NULL, NULL, NULL, NULL, m_naSta1SatID, m_naSta2SatID, ns, v, NULL, R,
			vflg, index,StaIndex1, StaIndex2);
		/*meaningless*/
		if (Valpos(v, R, vflg, nv, 4.0))
		{
			/* update state and covariance matrix */
			CMath::matcpy(m_RtkSolveparam[index].x, xp, m_RtkSolveparam[index].nx, 1);
			CMath::matcpy(m_RtkSolveparam[index].P, Pp, m_RtkSolveparam[index].nx, m_RtkSolveparam[index].nx);

			/* update ambiguity control struct */
			m_Sol[index].ns = 0;
			for (i = 0; i < ns; i++)
			{
				for (f = 0; f < nf; f++)
				{
					if (!m_Ssat[m_naCommonSat[i] - 1].vsat[f]) continue;
					m_Ssat[m_naCommonSat[i] - 1].lock[f]++;
					m_Ssat[m_naCommonSat[i] - 1].outc[f] = 0;
					if (f == 0) m_Sol[index].ns++; /* valid satellite count by L1 */
				}
			}
			/* lack of valid satellites */
			if (m_Sol[index].ns < 4) 
				stat = SOLQ_NONE;
		}
		else
		{
			stat = SOLQ_NONE;
		}	
	}
	/* resolve integer ambiguity by LAMBDA */
	//if (stat != SOLQ_NONE && AmbiFixByLambda(bias, xa,index) > 1)
	//{
	//	nv = CalDouDiffRes(pSta1ObsData, pSta2ObsData, dt, xp, NULL, m_naCommonSat, NULL, NULL, NULL, NULL, m_naSta1SatID, m_naSta2SatID, ns, v, NULL, R,
	//		vflg, index,StaIndex1, StaIndex2);
	//	/* validation of fixed solution */
	//	if (Valpos(v, R, vflg, nv, 4.0))
	//	{
	//		/* hold integer ambiguity */
	//		if (++m_RtkSolveparam[index].nfix >= m_Rtkopt.minfix &&
	//			m_Rtkopt.modear == ARMODE_FIXHOLD) {
	//			Holdamb(xa,index);
	//		}
	//		stat = SOLQ_FIX;
	//	}
	//}
	
	if (stat != SOLQ_FIX)
	{
		CalDouDiffAtmo(pSta1ObsData, pSta2ObsData, pSta1ObsData[0].tTimeStamp, xp,
			m_naCommonSat, m_naSta1SatID, m_naSta2SatID, ns, index, StaIndex1, StaIndex2);
	}
	else
	{
		CalDouDiffAtmo(pSta1ObsData, pSta2ObsData, pSta1ObsData[0].tTimeStamp, xa,
			m_naCommonSat, m_naSta1SatID, m_naSta2SatID, ns, index, StaIndex1, StaIndex2);
	}
	//SaveSoluStaus(stat,index);
	SaveSsat(pSta1ObsData, pSta2ObsData, pSinstaobsdata[StaIndex1].ui16ObsNum, pSinstaobsdata[StaIndex2].ui16ObsNum, stat,index);
	free(xp); free(Pp);  free(xa);  free(v); free(H); free(R); free(bias);

	if (stat != SOLQ_NONE) m_Sol[index].stat = stat;

	return stat != SOLQ_NONE;
}

void CBline::InitRtkopt(RtkOpt rtkopt)
{
	//m_Rtkopt = rtkopt;
	m_Rtkopt.nf = 2;
	m_Rtkopt.mode = 3;
	m_Rtkopt.ionoopt = 1;
	m_Rtkopt.tropopt = 1;
	m_Rtkopt.dynamics = 0;
	m_Rtkopt.glomodear = 0;
	m_Rtkopt.bdsmodear = 0;
	m_Rtkopt.navsys = 1;
	m_Rtkopt.maxout = 5;
	m_Rtkopt.minlock = 0;
	m_Rtkopt.minfix = 3;
	m_Rtkopt.err[0] = 100;
	m_Rtkopt.err[1] = 0.003;
	m_Rtkopt.err[2] = 0.003;
	m_Rtkopt.err[3] = 0.0;
	m_Rtkopt.err[4] = 1.0;
	m_Rtkopt.prn[0] = 0.0001;
	m_Rtkopt.prn[1] = 0.001;
	m_Rtkopt.prn[2] = 0.0001;
	m_Rtkopt.prn[3] = 10.0;
	m_Rtkopt.prn[4] = 10.0;
	m_Rtkopt.prn[5] = 0.0;
	m_Rtkopt.std[0] = 30.0;
	m_Rtkopt.std[1] = 0.03;
	m_Rtkopt.std[2] = 0.03;
	m_Rtkopt.sclkstab = 5e-12;
	m_Rtkopt.maxinno = 30.0;
	m_Rtkopt.eratio[0] = 100.0;
	m_Rtkopt.eratio[1] = 100.0;
	m_Rtkopt.eratio[2] = 0.0;
	m_Rtkopt.elmaskar = 0.0;
	m_Rtkopt.elmaskhold = 7.0 * D2R;
	m_Rtkopt.thresar[0] = 3.0;
	m_Rtkopt.thresar[1] = 1.0;
	m_Rtkopt.thresar[2] = 0.25;
	m_Rtkopt.thresar[3] = 0.10;
	m_Rtkopt.thresar[4] = 0.05;
	m_Rtkopt.thresar[5] = 0.0;
	m_Rtkopt.thresar[6] = 0.0;
	m_Rtkopt.thresar[7] = 0.0;

}

int CBline::Valpos(double* v, const double* R, int* vflg, int nv, double thres)
{
	double fact = thres * thres;
	int i, stat = 1, sat1, sat2, type, freq;
	char* stype;
	/* post-fit residual test */
	for (i = 0; i < nv; i++) {
		if (v[i] * v[i] <= fact * R[i + i * nv]) continue;
		sat1 = (vflg[i] >> 16) & 0xFF;
		sat2 = (vflg[i] >> 8) & 0xFF;
		type = (vflg[i] >> 4) & 0xF;
		freq = vflg[i] & 0xF;
		stype = type == 0 ? "L" : (type == 1 ? "L" : "C");
		/*errmsg(rtk, "large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
			sat1, sat2, stype, freq + 1, v[i], SQRT(R[i + i * nv]));*/
	}
	return stat;
}

void CBline::SaveSsat(ObsData* obs1, ObsData* obs2, int n1, int n2,int stat,int index)
{
	int i, j,k;
	int nf = m_Rtkopt.nf;
	ssat_t* m_Ssat = m_Stasatsolve.ssat[index];
	
	for (i = 0; i < n1; i++)
	{
		for (j = 0; j < nf; j++)
		{
			if (obs1[i].daCarrPhase[j] == 0.0) continue;
			m_Ssat[obs1[i].ui8SatId - 1].pt[0][j] = obs1[i].tTimeStamp;
			m_Ssat[obs1[i].ui8SatId - 1].ph[0][j] = obs1[i].daCarrPhase[j];
		}
	}
	for (k = 0; k < n2; k++)
	{
		for (j = 0; j < nf; j++)
		{
			if (obs2[k].daCarrPhase[j] == 0.0) continue;
			m_Ssat[obs2[k].ui8SatId - 1].pt[0][j] = obs2[i].tTimeStamp;
			m_Ssat[obs2[k].ui8SatId - 1].ph[0][j] = obs2[i].daCarrPhase[j];
		}
	}
	for (i = 0; i < MAXSAT; i++)
	{
		for (j = 0; j < nf; j++)
		{
			if (m_Ssat[i].fix[j] == 2 && stat != SOLQ_FIX) m_Ssat[i].fix[j] = 1;
			if (m_Ssat[i].slip[j] & 1) m_Ssat[i].slipc[j]++;
		}
	}
}

bool CBline::InitRtkParam(int n)
{
	int i,j;
	ssat_t ssat0 = { 0 };
	if (n < 0) return false;
	m_Stasol.n = m_Stasatsolve.n = m_Startksolpara.n = n;
	m_Stasol.sol = new sol_t[n];
	m_Startksolpara.rtksolpara = new RtkSolveParam[n];
	m_Stasatsolve.ssat = new ssat_t * [n];
	for (i = 0; i < n; i++)
	{
		m_Stasatsolve.ssat[i] = new ssat_t[MAXSAT];
		
		for (j = 0; j < MAXSAT; j++)
		{
			m_Stasatsolve.ssat[i][j] = ssat0;
		}
	}
	for (i = 0; i < n; i++)
	{
		m_Startksolpara.rtksolpara[i].nx= m_Rtkopt.mode <= PMODE_FIXED ? NX(&m_Rtkopt) : NX(&m_Rtkopt);
		m_Startksolpara.rtksolpara[i].na= m_Rtkopt.mode <= PMODE_FIXED ? NR(&m_Rtkopt) : NX(&m_Rtkopt);
		m_Startksolpara.rtksolpara[i].tt = 1.0;
		m_Startksolpara.rtksolpara[i].x=  CMath::zeros(m_Startksolpara.rtksolpara[i].nx, 1);
		m_Startksolpara.rtksolpara[i].P=  CMath::zeros(m_Startksolpara.rtksolpara[i].nx, m_Startksolpara.rtksolpara[i].nx);
		m_Startksolpara.rtksolpara[i].xa= CMath::zeros(m_Startksolpara.rtksolpara[i].na, 1);
		m_Startksolpara.rtksolpara[i].Pa= CMath::zeros(m_Startksolpara.rtksolpara[i].na, m_Startksolpara.rtksolpara[i].na);
		for (j = 0; j < 6; j++)m_Startksolpara.rtksolpara[i].rb[j] = 0.0;
		m_Stasol.sol[i] = { {0} };
	}
	
}




