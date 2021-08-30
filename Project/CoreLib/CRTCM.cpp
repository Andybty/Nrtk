#include <direct.h>
#include "CRTCM.h"

const char* msm_sig_gps[32] = {
	/* GPS: ref [17] table 3.5-91 */
	""  ,"1C","1P","1W",""  ,""  ,""  ,"2C","2P","2W",""  ,""  , /*  1-12 */
	""  ,""  ,"2S","2L","2X",""  ,""  ,""  ,""  ,"5I","5Q","5X", /* 13-24 */
	""  ,""  ,""  ,""  ,""  ,"1S","1L","1X"                      /* 25-32 */
};
const char* msm_sig_glo[32] = {
	/* GLONASS: ref [17] table 3.5-96 */
	""  ,"1C","1P",""  ,""  ,""  ,""  ,"2C","2P",""  ,""  ,""  ,
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char* msm_sig_gal[32] = {
	/* Galileo: ref [17] table 3.5-99 */
	""  ,"1C","1A","1B","1X","1Z",""  ,"6C","6A","6B","6X","6Z",
	""  ,"7I","7Q","7X",""  ,"8I","8Q","8X",""  ,"5I","5Q","5X",
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char* msm_sig_qzs[32] = {
	/* QZSS: ref [17] table 3.5-105 */
	""  ,"1C",""  ,""  ,""  ,""  ,""  ,""  ,"6S","6L","6X",""  ,
	""  ,""  ,"2S","2L","2X",""  ,""  ,""  ,""  ,"5I","5Q","5X",
	""  ,""  ,""  ,""  ,""  ,"1S","1L","1X"
};
const char* msm_sig_sbs[32] = {
	/* SBAS: ref [17] table 3.5-102 */
	""  ,"1C",""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,"5I","5Q","5X",
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char* msm_sig_cmp[32] = {
	/* BeiDou: ref [17] table 3.5-108 */
	""  ,"2I","2Q","2X",""  ,""  ,""  ,"6I","6Q","6X",""  ,""  ,
	""  ,"7I","7Q","7X",""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char* msm_sig_irn[32] = {
	/* NavIC/IRNSS: ref [17] table 3.5-108.3 */
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,"5A",""  ,""  ,
	""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
static char* obscodes[] = {       /* observation code strings */

	""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
	"1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
	"2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
	"6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
	"2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
	"5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
	"6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
};
static char codepris[7][MAXFREQ][16] = {  /* code priority for each freq-index */
   /*    0         1          2          3         4         5     */
	{"CPYWMNSL","PYWCMNDLSX","IQX"     ,""       ,""       ,""      ,""}, /* GPS */
	{"CPABX"   ,"PCABX"     ,"IQX"     ,""       ,""       ,""      ,""}, /* GLO */
	{"CABXZ"   ,"IQX"       ,"IQX"     ,"ABCXZ"  ,"IQX"    ,""      ,""}, /* GAL */
	{"CLSXZ"   ,"LSX"       ,"IQXDPZ"  ,"LSXEZ"  ,""       ,""      ,""}, /* QZS */
	{"C"       ,"IQX"       ,""        ,""       ,""       ,""      ,""}, /* SBS */
	{"IQXDPAN" ,"IQXDPZ"    ,"DPX"     ,"IQXA"   ,"DPX"    ,""      ,""}, /* BDS */
	{"ABCX"    ,"ABCX"      ,""        ,""       ,""       ,""      ,""}  /* IRN */
};
/* GPS obs code to frequency -------------------------------------------------*/
static int code2freq_GPS(uint8_t code, double* freq)
{
	char* obs = code2obs(code);

	switch (obs[0]) {
	case '1': *freq = FREQ1; return 0; /* L1 */
	case '2': *freq = FREQ2; return 1; /* L2 */
	case '5': *freq = FREQ5; return 2; /* L5 */
	}
	return -1;
}
/* GLONASS obs code to frequency ---------------------------------------------*/
static int code2freq_GLO(uint8_t code, int fcn, double* freq)
{
	char* obs = code2obs(code);

	if (fcn < -7 || fcn>6) return -1;

	switch (obs[0]) {
	case '1': *freq = FREQ1_GLO + DFRQ1_GLO * fcn; return 0; /* G1 */
	case '2': *freq = FREQ2_GLO + DFRQ2_GLO * fcn; return 1; /* G2 */
	case '3': *freq = FREQ3_GLO;               return 2; /* G3 */
	case '4': *freq = FREQ1a_GLO;              return 0; /* G1a */
	case '6': *freq = FREQ2a_GLO;              return 1; /* G2a */
	}
	return -1;
}
/* Galileo obs code to frequency ---------------------------------------------*/
static int code2freq_GAL(uint8_t code, double* freq)
{
	char* obs = code2obs(code);

	switch (obs[0]) {
	case '1': *freq = FREQ1; return 0; /* E1 */
	case '7': *freq = FREQ7; return 1; /* E5b */
	case '5': *freq = FREQ5; return 2; /* E5a */
	case '6': *freq = FREQ6; return 3; /* E6 */
	case '8': *freq = FREQ8; return 4; /* E5ab */
	}
	return -1;
}
/* QZSS obs code to frequency ------------------------------------------------*/
static int code2freq_QZS(uint8_t code, double* freq)
{
	char* obs = code2obs(code);

	switch (obs[0]) {
	case '1': *freq = FREQ1; return 0; /* L1 */
	case '2': *freq = FREQ2; return 1; /* L2 */
	case '5': *freq = FREQ5; return 2; /* L5 */
	case '6': *freq = FREQ6; return 3; /* L6 */
	}
	return -1;
}
/* SBAS obs code to frequency ------------------------------------------------*/
static int code2freq_SBS(uint8_t code, double* freq)
{
	char* obs = code2obs(code);

	switch (obs[0]) {
	case '1': *freq = FREQ1; return 0; /* L1 */
	case '5': *freq = FREQ5; return 1; /* L5 */
	}
	return -1;
}
/* BDS obs code to frequency -------------------------------------------------*/
static int code2freq_BDS(uint8_t code, double* freq)
{
	char* obs = code2obs(code);

	switch (obs[0]) {
	case '1': *freq = FREQ1;     return 0; /* B1C */
	case '2': *freq = FREQ1_CMP; return 0; /* B1I */
	case '7': *freq = FREQ2_CMP; return 1; /* B2I/B2b */
	case '5': *freq = FREQ5;     return 2; /* B2a */
	case '6': *freq = FREQ3_CMP; return 3; /* B3 */
	case '8': *freq = FREQ8;     return 4; /* B2ab */
	}
	return -1;
}
/* NavIC obs code to frequency -----------------------------------------------*/
static int code2freq_IRN(uint8_t code, double* freq)
{
	char* obs = code2obs(code);

	switch (obs[0]) {
	case '5': *freq = FREQ5; return 0; /* L5 */
	case '9': *freq = FREQ9; return 1; /* S */
	}
	return -1;
}

static int obsindex(SinStaObsData*obs, gtime_t time, int sat)
{
	int i, j;

	for (i = 0; i < obs->ui16ObsNum; i++) {
		if (obs->pObsData[i].ui8SatId == sat) return i; /* field already exists */
	}
	if (i >= MAXOBS) return -1; /* overflow */

	/* add new field */
	obs->pObsData[i].tTimeStamp = time;
	obs->pObsData[i].ui8SatId = sat;
	/*for (j = 0; j < NFREQ + NEXOBS; j++) {
		obs->data[i].L[j] = obs->data[i].P[j] = 0.0;
		obs->data[i].D[j] = 0.0;
		obs->data[i].SNR[j] = obs->data[i].LLI[j] = obs->data[i].code[j] = 0;
	}
	obs->n++;*/
	for (j = 0; j < NFREQ + NEXOBS; j++)
	{
		obs->pObsData[i].daCarrPhase[j] = obs->pObsData[i].daPseRange[j] = 0.0;
		obs->pObsData[i].daDoppler[j] = 0.0;
		obs->pObsData[i].daSigNoiRatio[j] = obs->pObsData[i].ui8aLLI[j] = obs->pObsData[i].ui8aCodeType[j] = 0;
	}
	obs->ui16ObsNum++;
	return i;
}

/* loss-of-lock indicator ----------------------------------------------------*/
static int lossoflock(RtcmData* rtcm, int sat, int idx, int lock)
{
	int lli = (!lock && !rtcm->lock[sat - 1][idx]) || lock < rtcm->lock[sat - 1][idx];
	rtcm->lock[sat - 1][idx] = (uint16_t)lock;
	return lli;
}

/* satellite system+prn/slot number to satellite number ------------------------
* convert satellite system+prn/slot number to satellite number
* args   : int    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
*          int    prn       I   satellite prn/slot number
* return : satellite number (0:error)
*-----------------------------------------------------------------------------*/
int satno(int sys, int prn)
{
	if (prn <= 0) return 0;
	switch (sys) {
	case SYS_GPS:
		if (prn < MINPRNGPS || MAXPRNGPS < prn) return 0;
		return prn - MINPRNGPS + 1;
	case SYS_GLO:
		if (prn < MINPRNGLO || MAXPRNGLO < prn) return 0;
		return NSATGPS + prn - MINPRNGLO + 1;
	case SYS_GAL:
		if (prn < MINPRNGAL || MAXPRNGAL < prn) return 0;
		return NSATGPS + NSATGLO + prn - MINPRNGAL + 1;
	case SYS_QZS:
		if (prn < MINPRNQZS || MAXPRNQZS < prn) return 0;
		return NSATGPS + NSATGLO + NSATGAL + prn - MINPRNQZS + 1;
	case SYS_CMP:
		if (prn < MINPRNCMP || MAXPRNCMP < prn) return 0;
		return NSATGPS + NSATGLO + NSATGAL + NSATQZS + prn - MINPRNCMP + 1;
	case SYS_IRN:
		if (prn < MINPRNIRN || MAXPRNIRN < prn) return 0;
		return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + prn - MINPRNIRN + 1;
	case SYS_LEO:
		if (prn < MINPRNLEO || MAXPRNLEO < prn) return 0;
		return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN +
			prn - MINPRNLEO + 1;
	case SYS_SBS:
		if (prn < MINPRNSBS || MAXPRNSBS < prn) return 0;
		return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN + NSATLEO +
			prn - MINPRNSBS + 1;
	}
	return 0;
}
/* get signal index ----------------------------------------------------------*/
void sigindex(int sys, const uint8_t* code, int n, const char* opt,
	int* idx)
{
	int i, nex, pri, pri_h[8] = { 0 }, index[8] = { 0 }, ex[32] = { 0 };

	/* test code priority */
	for (i = 0; i < n; i++) {
		if (!code[i]) continue;

		if (idx[i] >= NFREQ) { /* save as extended signal if idx >= NFREQ */
			ex[i] = 1;
			continue;
		}
		/* code priority */
		pri = getcodepri(sys, code[i], opt);

		/* select highest priority signal */
		if (pri > pri_h[idx[i]]) {
			if (index[idx[i]]) ex[index[idx[i]] - 1] = 1;
			pri_h[idx[i]] = pri;
			index[idx[i]] = i + 1;
		}
		else ex[i] = 1;
	}
	/* signal index in obs data */
	for (i = nex = 0; i < n; i++) {
		if (ex[i] == 0);
		else if (nex < NEXOBS) idx[i] = NFREQ + nex++;
		else { /* no space in obs data */
			idx[i] = -1;
		}
#if 0 /* for debug */
#endif
	}
}
/*get code priority---------------------------------------------------------- -
*get code priority for multiple codes in a frequency
* args   : int    sys       I   system(SYS_ ? ? ? )
* uint8_t code     I   obs code(CODE_ ? ? ? )
* char* opt      I   code options(NULL : no option)
* return : priority(15:highest - 1 : lowest, 0 : error)
* ---------------------------------------------------------------------------- - */
int getcodepri(int sys, uint8_t code, const char* opt)
{
	const char* p, * optstr;
	char* obs, str[8] = "";
	int i, j;

	switch (sys) {
	case SYS_GPS: i = 0; optstr = "-GL%2s"; break;
	case SYS_GLO: i = 1; optstr = "-RL%2s"; break;
	case SYS_GAL: i = 2; optstr = "-EL%2s"; break;
	case SYS_QZS: i = 3; optstr = "-JL%2s"; break;
	case SYS_SBS: i = 4; optstr = "-SL%2s"; break;
	case SYS_CMP: i = 5; optstr = "-CL%2s"; break;
	case SYS_IRN: i = 6; optstr = "-IL%2s"; break;
	default: return 0;
	}
	if ((j = code2idx(sys, code)) < 0) return 0;
	obs = code2obs(code);

	/* parse code options */
	for (p = opt; p && (p = strchr(p, '-')); p++) {
		if (sscanf(p, optstr, str) < 1 || str[0] != obs[0]) continue;
		return str[1] == obs[1] ? 15 : 0;
	}
	/* search code priority */
	return (p = strchr(codepris[i][j], obs[1])) ? 14 - (int)(p - codepris[i][j]) : 0;
}
/* obs code to obs code string -------------------------------------------------
* convert obs code to obs code string
* args   : uint8_t code     I   obs code (CODE_???)
* return : obs code string ("1C","1P","1P",...)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
char* code2obs(uint8_t code)
{
	if (code <= CODE_NONE || MAXCODE < code) return "";
	return obscodes[code];
}
/* system and obs code to frequency index --------------------------------------
* convert system and obs code to frequency index
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
* return : frequency index (-1: error)
*                       0     1     2     3     4
*           --------------------------------------
*            GPS       L1    L2    L5     -     -
*            GLONASS   G1    G2    G3     -     -  (G1=G1,G1a,G2=G2,G2a)
*            Galileo   E1    E5b   E5a   E6   E5ab
*            QZSS      L1    L2    L5    L6     -
*            SBAS      L1     -    L5     -     -
*            BDS       B1    B2    B2a   B3   B2ab (B1=B1I,B1C,B2=B2I,B2b)
*            NavIC     L5     S     -     -     -
*-----------------------------------------------------------------------------*/
int code2idx(int sys, uint8_t code)
{
	double freq;

	switch (sys) {
	case SYS_GPS: return code2freq_GPS(code, &freq);
	case SYS_GLO: return code2freq_GLO(code, 0, &freq);
	case SYS_GAL: return code2freq_GAL(code, &freq);
	case SYS_QZS: return code2freq_QZS(code, &freq);
	case SYS_SBS: return code2freq_SBS(code, &freq);
	case SYS_CMP: return code2freq_BDS(code, &freq);
	case SYS_IRN: return code2freq_IRN(code, &freq);
	}
	return -1;
}
/* obs type string to obs code -------------------------------------------------
* convert obs code type string to obs code
* args   : char   *str      I   obs code string ("1C","1P","1Y",...)
* return : obs code (CODE_???)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
uint8_t obs2code(const char* obs)
{
	int i;

	for (i = 1; *obscodes[i]; i++) {
		if (strcmp(obscodes[i], obs)) continue;
		return (uint8_t)i;
	}
	return CODE_NONE;
}
/* save obs data in MSM message----------------------------------------------*/
void save_msm_obs(SinStaObsData* pSinStaObsData,RtcmData* rtcm, int sys, msm_h_t* h, const double* r,
	const double* pr, const double* cp, const double* rr,
	const double* rrf, const double* cnr, const int* lock,
	const int* ex, const int* half)
{
//_____________________________________________________
//______________________________________________________
	const char* sig[32];
	double tt, freq;
	uint8_t code[32];
	char* msm_type = "", * q = NULL;
	int i, j, k, type, prn, sat, fcn, index = 0, idx[32];

	type = getbitu(rtcm->buff, 24, 12);

	switch (sys) {
	case SYS_GPS: msm_type = q = rtcm->msmtype[0]; break;
	case SYS_GLO: msm_type = q = rtcm->msmtype[1]; break;
	case SYS_GAL: msm_type = q = rtcm->msmtype[2]; break;
	case SYS_QZS: msm_type = q = rtcm->msmtype[3]; break;
	case SYS_SBS: msm_type = q = rtcm->msmtype[4]; break;
	case SYS_CMP: msm_type = q = rtcm->msmtype[5]; break;
	case SYS_IRN: msm_type = q = rtcm->msmtype[6]; break;
	}
	/* id to signal */
	for (i = 0; i < h->nsig; i++) {
		switch (sys) {
		case SYS_GPS: sig[i] = msm_sig_gps[h->sigs[i] - 1]; break;
		case SYS_GLO: sig[i] = msm_sig_glo[h->sigs[i] - 1]; break;
		case SYS_GAL: sig[i] = msm_sig_gal[h->sigs[i] - 1]; break;
		case SYS_QZS: sig[i] = msm_sig_qzs[h->sigs[i] - 1]; break;
		case SYS_SBS: sig[i] = msm_sig_sbs[h->sigs[i] - 1]; break;
		case SYS_CMP: sig[i] = msm_sig_cmp[h->sigs[i] - 1]; break;
		case SYS_IRN: sig[i] = msm_sig_irn[h->sigs[i] - 1]; break;
		default: sig[i] = ""; break;
		}
		/* signal to rinex obs type */
		code[i] = obs2code(sig[i]);
		idx[i] = code2idx(sys, code[i]);
	}
	/* get signal index */
	sigindex(sys, code, h->nsig, rtcm->opt, idx);
	for (i = j = 0; i < h->nsat; i++) {

		prn = h->sats[i];
		/*if (sys == SYS_QZS) prn += MINPRNQZS - 1;
		else if (sys == SYS_SBS) prn += MINPRNSBS - 1;*/

		if ((sat = satno(sys, prn))) {
			tt = CTime::timediff(pSinStaObsData->pObsData[0].tTimeStamp, rtcm->time);
			if (rtcm->obsflag || fabs(tt) > 1E-9) {
				pSinStaObsData->ui16ObsNum = rtcm->obsflag = 0;
			}
			index = obsindex(pSinStaObsData, rtcm->time, sat);
		}
		else {
			//trace(2, "rtcm3 %d satellite error: prn=%d\n", type, prn);
		}
		fcn = 0;
		//if (sys == SYS_GLO)
		//{
		//	fcn = -8; /* no glonass fcn info */
		//	if (ex && ex[i] <= 13) {
		//		fcn = ex[i] - 7;
		//		if (!rtcm->nav.glo_fcn[prn - 1]) {
		//			rtcm->nav.glo_fcn[prn - 1] = fcn + 8; /* fcn+8 */
		//		}
		//	}
		//	else if (rtcm->nav.geph[prn - 1].sat == sat) {
		//		fcn = rtcm->nav.geph[prn - 1].frq;
		//	}
		//	else if (rtcm->nav.glo_fcn[prn - 1] > 0) {
		//		fcn = rtcm->nav.glo_fcn[prn - 1] - 8;
		//	}
		//}
		for (k = 0; k < h->nsig; k++) 
		{
			if (!h->cellmask[k + i * h->nsig]) continue;

			if (sat && index >= 0 && idx[k] >= 0) {
				freq = fcn < -7 ? 0.0 : code2freq(sys, code[k], fcn);

				/* pseudorange (m) */
				if (r[i] != 0.0 && pr[j] > -1E12) {
					pSinStaObsData->pObsData[index].daPseRange[idx[k]] = r[i] + pr[j];
				}
				/* carrier-phase (cycle) */
				if (r[i] != 0.0 && cp[j] > -1E12) {
					pSinStaObsData->pObsData[index].daCarrPhase[idx[k]] = (r[i] + cp[j]) * freq / CLIGHT;
				}
				/* doppler (hz) */
				if (rr && rrf && rrf[j] > -1E12) {
					pSinStaObsData->pObsData[index].daDoppler[idx[k]] =
						(float)(-(rr[i] + rrf[j]) * freq / CLIGHT);
				}
				pSinStaObsData->pObsData[index].ui8aLLI[idx[k]] =
					lossoflock(rtcm, sat, idx[k], lock[j]) + (half[j] ? 3 : 0);
				pSinStaObsData->pObsData[index].daSigNoiRatio[idx[k]] = (uint16_t)(cnr[j] / SNR_UNIT + 0.5);
				pSinStaObsData->pObsData[index].ui8aCodeType[idx[k]] = code[k];
			}
			j++;
		}
	}
	//CTime::time2str(rtcm->time, tstr, 2);
}
	
/* system and obs code to frequency --------------------------------------------
* convert system and obs code to carrier frequency
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
*          int    fcn       I   frequency channel number for GLONASS
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
double code2freq(int sys, uint8_t code, int fcn)
{
	double freq = 0.0;

	switch (sys) {
	case SYS_GPS: (void)code2freq_GPS(code, &freq); break;
	case SYS_GLO: (void)code2freq_GLO(code, fcn, &freq); break;
	case SYS_GAL: (void)code2freq_GAL(code, &freq); break;
	case SYS_QZS: (void)code2freq_QZS(code, &freq); break;
	case SYS_SBS: (void)code2freq_SBS(code, &freq); break;
	case SYS_CMP: (void)code2freq_BDS(code, &freq); break;
	case SYS_IRN: (void)code2freq_IRN(code, &freq); break;
	}
	return freq;
}
	
	

	

void adjweek(RtcmData* rtcm, double tow)
{
	double tow_p;
	int week;

	/* if no time, get cpu time */
	if (rtcm->time.tiTime == 0) rtcm->time = CTime::utc2gpst(CTime::timeget());
	tow_p = CTime::time2gpst(rtcm->time, &week);
	if (tow < tow_p - 302400.0) tow += 604800.0;
	else if (tow > tow_p + 302400.0) tow -= 604800.0;
	rtcm->time = CTime::gpst2time(week, tow);
}
/* adjust daily rollover of GLONASS time -------------------------------------*/
void adjday_glot(RtcmData* rtcm, double tod)
{
	gtime_t time,time0;
	double tow, tod_p;
	int week;
	time0 = CTime::gpst2utc(rtcm->time);
	if (rtcm->time.tiTime == 0) rtcm->time = CTime::utc2gpst(CTime::timeget());
	time = CTime::timeadd(time0, 10800.0); /* glonass time */
	tow = CTime::time2gpst(time, &week);
	tod_p = fmod(tow, 86400.0); tow -= tod_p;
	if (tod < tod_p - 43200.0) tod += 86400.0;
	else if (tod > tod_p + 43200.0) tod -= 86400.0;
	time = CTime::gpst2time(week, tow + tod);
	rtcm->time = CTime::utc2gpst(CTime::timeadd(time, -10800.0));
}


int decode_msm_head(RtcmData* rtcm, int sys, int* sync, int* iod,
	msm_h_t* h, int* hsize)
{
	msm_h_t h0 = { 0 };
	double tow, tod;
	char* msg, tstr[64];
	int i = 24, j, dow, mask, staid, type, ncell = 0;

	type = getbitu(rtcm->buff, i, 12); i += 12;

	*h = h0;
	if (i + 157 <= rtcm->len * 8) {
		staid = getbitu(rtcm->buff, i, 12);       i += 12;

		if (sys == SYS_GLO) {
			dow = getbitu(rtcm->buff, i, 3);       i += 3;
			tod = getbitu(rtcm->buff, i, 27) * 0.001; i += 27;
			adjday_glot(rtcm, tod);
		}
		else if (sys == SYS_CMP) {
			tow = getbitu(rtcm->buff, i, 30) * 0.001; i += 30;
			tow += 14.0; /* BDT -> GPST */
			adjweek(rtcm, tow);
		}
		else {
			tow = getbitu(rtcm->buff, i, 30) * 0.001; i += 30;
			adjweek(rtcm, tow);
		}
		*sync = getbitu(rtcm->buff, i, 1);       i += 1;
		*iod = getbitu(rtcm->buff, i, 3);       i += 3;
		h->time_s = getbitu(rtcm->buff, i, 7);       i += 7;
		h->clk_str = getbitu(rtcm->buff, i, 2);       i += 2;
		h->clk_ext = getbitu(rtcm->buff, i, 2);       i += 2;
		h->smooth = getbitu(rtcm->buff, i, 1);       i += 1;
		h->tint_s = getbitu(rtcm->buff, i, 3);       i += 3;
		for (j = 1; j <= 64; j++) {
			mask = getbitu(rtcm->buff, i, 1); i += 1;
			if (mask) h->sats[h->nsat++] = j;
		}
		for (j = 1; j <= 32; j++) {
			mask = getbitu(rtcm->buff, i, 1); i += 1;
			if (mask) h->sigs[h->nsig++] = j;
		}
	}
	else {
		return -1;
	}
	/* test station id */
	//if (!test_staid(rtcm, staid)) return -1;

	if (h->nsat * h->nsig > 64) {
		return -1;
	}
	if (i + h->nsat * h->nsig > rtcm->len * 8) {
		return -1;
	}
	for (j = 0; j < h->nsat * h->nsig; j++) {
		h->cellmask[j] = getbitu(rtcm->buff, i, 1); i += 1;
		if (h->cellmask[j]) ncell++;
	}
	*hsize = i;

	CTime::time2str(rtcm->time, tstr, 2);
	return ncell;
}
/* decode MSM 4: full pseudorange and phaserange plus CNR --------------------*/
 int decode_msm4(SinStaObsData* pSinStaObsData,RtcmData* rtcm, int sys)
{
	 msm_h_t h = { 0 };
	 double r[64], pr[64], cp[64], cnr[64];
	 int i, j, type, sync, iod, ncell, rng, rng_m, prv, cpv, lock[64], half[64];

	 type = getbitu(rtcm->buff, 24, 12);

	 /* decode msm header */
	 if ((ncell = decode_msm_head(rtcm, sys, &sync, &iod, &h, &i)) < 0) return -1;

	 if (i + h.nsat * 18 + ncell * 48 > rtcm->len * 8) {
		 return -1;
	 }
	 for (j = 0; j < h.nsat; j++) r[j] = 0.0;
	 for (j = 0; j < ncell; j++) pr[j] = cp[j] = -1E16;

	 /* decode satellite data */
	 for (j = 0; j < h.nsat; j++) { /* range */
		 rng = getbitu(rtcm->buff, i, 8); i += 8;
		 if (rng != 255) r[j] = rng * RANGE_MS;
	 }
	 for (j = 0; j < h.nsat; j++) {
		 rng_m = getbitu(rtcm->buff, i, 10); i += 10;
		 if (r[j] != 0.0) r[j] += rng_m * P2_10 * RANGE_MS;
	 }
	 /* decode signal data */
	 for (j = 0; j < ncell; j++) { /* pseudorange */
		 prv = getbits(rtcm->buff, i, 15); i += 15;
		 if (prv != -16384) pr[j] = prv * P2_24 * RANGE_MS;
	 }
	 for (j = 0; j < ncell; j++) { /* phaserange */
		 cpv = getbits(rtcm->buff, i, 22); i += 22;
		 if (cpv != -2097152) cp[j] = cpv * P2_29 * RANGE_MS;
	 }
	 for (j = 0; j < ncell; j++) { /* lock time */
		 lock[j] = getbitu(rtcm->buff, i, 4); i += 4;
	 }
	 for (j = 0; j < ncell; j++) { /* half-cycle ambiguity */
		 half[j] = getbitu(rtcm->buff, i, 1); i += 1;
	 }
	 for (j = 0; j < ncell; j++) { /* cnr */
		 cnr[j] = getbitu(rtcm->buff, i, 6) * 1.0; i += 6;
	 }
	 /*save obs data in msm message */
	 
		save_msm_obs(pSinStaObsData,rtcm, sys, &h, r, pr, cp, NULL, NULL, cnr, lock, NULL, half);

	 rtcm->obsflag = !sync;
	 return sync ? 0 : 1;
}
unsigned int getbitu(const unsigned char* buff, int pos, int len)
{
	unsigned int bits = 0;
	int i;
	for (i = pos; i < pos + len; i++) bits = (bits << 1) + ((buff[i / 8] >> (7 - i % 8)) & 1u);
	return bits;
}
int getbits(const unsigned char* buff, int pos, int len)
{
	unsigned int bits = getbitu(buff, pos, len);
	if (len <= 0 || 32 <= len || !(bits & (1u << (len - 1)))) return (int)bits;
	return (int)(bits | (~0u << len)); /* extend sign */
}


unsigned int rtk_crc24q(const unsigned char* buff, int len)
{
	unsigned int crc = 0;
	int i;

	for (i = 0; i < len; i++) crc = ((crc << 8) & 0xFFFFFF) ^ tbl_CRC24Q[(crc >> 16) ^ buff[i]];
	return crc;
}
int decode_rtcm3obs(SinStaObsData* pSinStaObsData,RtcmData* rtcm)
{
	double tow;
	int ret = 0, type = getbitu(rtcm->buff, 24, 12), week;
	switch (type)
	{
	case 1074:
		ret = decode_msm4(pSinStaObsData,rtcm, SYS_GPS);
		break;
	case 1124:
		ret = decode_msm4(pSinStaObsData,rtcm, SYS_CMP); 
		break;
	}
	return ret;
}
int decode_rtcm3eph(AllBroadcastEphData* pAllBroadcastEphdata, RtcmData* rtcm)
{
	double tow;
	int ret = 0, type = getbitu(rtcm->buff, 24, 12), week;
	switch (type)
	{
	case 1019:
		ret = decode_type1019(rtcm, pAllBroadcastEphdata);
		break;
	}
	return ret;
}
int input_rtcm3eph(AllBroadcastEphData* pAllBroadcastEphdata, RtcmData* rtcm, unsigned char data)
{
	
	/* synchronize frame */
	if (rtcm->nbyte == 0) {
		if (data != RTCM3PREAMB) return 0;
		rtcm->buff[rtcm->nbyte++] = data;
		return 0;
	}
	rtcm->buff[rtcm->nbyte++] = data;

	if (rtcm->nbyte == 3) {
		rtcm->len = getbitu(rtcm->buff, 14, 10) + 3; /* length without parity */
	}
	if (rtcm->nbyte < 3 || rtcm->nbyte < rtcm->len + 3) return 0;
	rtcm->nbyte = 0;

	/* check parity */
	if (rtk_crc24q(rtcm->buff, rtcm->len) != getbitu(rtcm->buff, rtcm->len * 8, 24)) {
		return 0;
	}
	/* decode rtcm3 message */
	return decode_rtcm3eph(pAllBroadcastEphdata, rtcm);
}
int input_rtcm3obs(SinStaObsData* pSinStaObsData,RtcmData* rtcm, unsigned char data)
{
	int i;
	/* synchronize frame */
	if (rtcm->nbyte == 0) {
		if (data != RTCM3PREAMB) return 0;
		rtcm->buff[rtcm->nbyte++] = data;
		return 0;
	}
	rtcm->buff[rtcm->nbyte++] = data;

	if (rtcm->nbyte == 3) {
		rtcm->len = getbitu(rtcm->buff, 14, 10) + 3; /* length without parity */
	}
	if (rtcm->nbyte < 3 || rtcm->nbyte < rtcm->len + 3) return 0;
	rtcm->nbyte = 0;

	/* check parity */
	if (rtk_crc24q(rtcm->buff, rtcm->len) != getbitu(rtcm->buff, rtcm->len * 8, 24)) {
		return 0;
	}
	/* decode rtcm3 message */
	   return decode_rtcm3obs(pSinStaObsData,rtcm);
	
}
/*decode type 1019: GPS ephemerides---------------------------------------- - */
int decode_type1019(RtcmData* rtcm, AllBroadcastEphData* pAllBroadcastEphdata)
{
	BroadcastEphData eph = { 0 };
	double toc, sqrtA, tt;
	char* msg;
	int i = 24 + 12, prn, sat, week, sys = SYS_GPS;

	if (i + 476 <= rtcm->len * 8) {
		prn = getbitu(rtcm->buff, i, 6);              i += 6;
		week = getbitu(rtcm->buff, i, 10);              i += 10;
		eph.sva = getbitu(rtcm->buff, i, 4);              i += 4;
		eph.code = getbitu(rtcm->buff, i, 2);              i += 2;
		eph.idot = getbits(rtcm->buff, i, 14) * P2_43 * SC2RAD; i += 14;
		eph.iode = getbitu(rtcm->buff, i, 8);              i += 8;
		toc = getbitu(rtcm->buff, i, 16) * 16.0;         i += 16;
		eph.f2 = getbits(rtcm->buff, i, 8) * P2_55;        i += 8;
		eph.f1 = getbits(rtcm->buff, i, 16) * P2_43;        i += 16;
		eph.f0 = getbits(rtcm->buff, i, 22) * P2_31;        i += 22;
		eph.iodc = getbitu(rtcm->buff, i, 10);              i += 10;
		eph.crs = getbits(rtcm->buff, i, 16) * P2_5;         i += 16;
		eph.deln = getbits(rtcm->buff, i, 16) * P2_43 * SC2RAD; i += 16;
		eph.M0 = getbits(rtcm->buff, i, 32) * P2_31 * SC2RAD; i += 32;
		eph.cuc = getbits(rtcm->buff, i, 16) * P2_29;        i += 16;
		eph.e = getbitu(rtcm->buff, i, 32) * P2_33;        i += 32;
		eph.cus = getbits(rtcm->buff, i, 16) * P2_29;        i += 16;
		sqrtA = getbitu(rtcm->buff, i, 32) * P2_19;        i += 32;
		eph.toes = getbitu(rtcm->buff, i, 16) * 16.0;         i += 16;
		eph.cic = getbits(rtcm->buff, i, 16) * P2_29;        i += 16;
		eph.OMG0 = getbits(rtcm->buff, i, 32) * P2_31 * SC2RAD; i += 32;
		eph.cis = getbits(rtcm->buff, i, 16) * P2_29;        i += 16;
		eph.i0 = getbits(rtcm->buff, i, 32) * P2_31 * SC2RAD; i += 32;
		eph.crc = getbits(rtcm->buff, i, 16) * P2_5;         i += 16;
		eph.omg = getbits(rtcm->buff, i, 32) * P2_31 * SC2RAD; i += 32;
		eph.OMGd = getbits(rtcm->buff, i, 24) * P2_43 * SC2RAD; i += 24;
		eph.tgd[0] = getbits(rtcm->buff, i, 8) * P2_31;        i += 8;
		eph.svh = getbitu(rtcm->buff, i, 6);              i += 6;
		eph.flag = getbitu(rtcm->buff, i, 1);              i += 1;
		eph.fit = getbitu(rtcm->buff, i, 1) ? 0.0 : 4.0; /* 0:4hr,1:>4hr */
	}
	else {
		//trace(2, "rtcm3 1019 length error: len=%d\n", rtcm->len);
		return -1;
	}
	if (prn >= 40) {
		sys = SYS_SBS; prn += 80;
	}
	//trace(4, "decode_type1019: prn=%d iode=%d toe=%.0f\n", prn, eph.iode, eph.toes);

	/*if (rtcm->outtype) {
		msg = rtcm->msgtype + strlen(rtcm->msgtype);
		sprintf(msg, " prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
			prn, eph.iode, eph.iodc, week, eph.toes, toc, eph.svh);
	}*/
	if (!(sat = satno(sys, prn))) {
		//trace(2, "rtcm3 1019 satellite number error: prn=%d\n", prn);
		return -1;
	}
	eph.sat = sat;
	eph.week = CTime::adjgpsweek(week);
	if (rtcm->time.tiTime == 0) 
		rtcm->time = CTime::utc2gpst(CTime::timeget());
	tt = CTime::timediff(CTime::gpst2time(eph.week, eph.toes), rtcm->time);
	if (tt < -302400.0)
		eph.week++;
	else if (tt >= 302400.0) 
		eph.week--;
	eph.toe = CTime::gpst2time(eph.week, eph.toes);
	eph.toc = CTime::gpst2time(eph.week, toc);
	eph.ttr = rtcm->time;
	eph.A = sqrtA * sqrtA;
	if (!strstr(rtcm->opt, "-EPHALL")) {
		if (eph.iode == pAllBroadcastEphdata->pBroadcastEphData[sat - 1].iode) return 0; /* unchanged */
	}
	pAllBroadcastEphdata->pBroadcastEphData[sat - 1] = eph;
	//rtcm->ephsat = sat;
	//rtcm->ephset = 0;
	return 2;
}
CRTCM::CRTCM()
{
	
}

CRTCM::~CRTCM()
{

}

void CRTCM::InitRtcm()
{
	int i, j;
	gtime_t time0 = { 0 };
	m_RtcmData.len = m_RtcmData.nbit = m_RtcmData.nbyte = m_RtcmData.obsflag = 0;
	m_RtcmData.time = m_RtcmData.time_s = time0;
	m_RtcmData.opt[0] = '\0';
	m_RtcmData.staid = 0;
	m_RtcmData.seqno = 0;
	for (i = 0; i < 6; i++) m_RtcmData.msmtype[i][0] = '\0';
	for (i = 0; i < MAXSAT; i++) for (j = 0; j < NFREQ + NEXOBS; j++) {
		m_RtcmData.cp[i][j] = 0.0;
		m_RtcmData.lock[i][j] = m_RtcmData.loss[i][j] = 0;
		m_RtcmData.lltime[i][j] = time0;
	}
	for (i = 0; i < 400; i++) m_RtcmData.nmsg3[i] = 0;
}


int CRTCM::DecodeEphType(AllBroadcastEphData* pAllBroadcastEphData, uint8* str_Stream, string strFileName)
{
	FILE* fp1;
	int info;
	struct stat tagfileinfo;
	char tstr[64];
	fp1 = fopen(strFileName.c_str(), "rb");
	int i, data, ret;
	//const char* c = str_Stream.c_str();
	info = stat(strFileName.c_str(), &tagfileinfo);
	if (info != 0)
	{
		perror("显示文件状态信息出错");
	}
	else
	{
		m_RtcmData.time.tiTime = tagfileinfo.st_ctime;
		m_RtcmData.time.dSec = 0.0;
	}
	m_RtcmData.time_s = m_RtcmData.time;
	CTime::time2str(m_RtcmData.time, tstr, 5);
	while ((data = fgetc(fp1)) != EOF)
	{
		ret = input_rtcm3eph(pAllBroadcastEphData, &m_RtcmData, (unsigned char)data);

	}
	fclose(fp1);
	return 1;
}

int  CRTCM::DecodeObsType(SinStaObsData* pSinStaObsData, uint8* str_Stream)
{
	
	InitObs();
	int i, data, ret;
	double dt;
	//const char* c = str_Stream.c_str();
	for (i = 0; i < 100000; i++)
	{
		if ((data = fgetc(fp1)) == EOF)
		{
			fclose(fp1);
			return -2;
		}
		if ((ret = input_rtcm3obs(&pSinStaObsData[0], &m_RtcmData, (unsigned char)data)))
		{
			break;
		}			
	}
	for (i = 0; i < 10000; i++)
	{
		if ((data = fgetc(fp2)) == EOF)
		{
			fclose(fp2);
			return -2;
		}
		if ((ret = input_rtcm3obs(&pSinStaObsData[1], &m_RtcmData, (unsigned char)data)))
		{
			break;
		}
	}
	dt = CTime::timediff(pSinStaObsData[0].pObsData[0].tTimeStamp, pSinStaObsData[1].pObsData[0].tTimeStamp);
	/*if (dt > 0)
	{
		do {
			for (i = 0; i < 4096; i++)
			{
				if ((data = fgetc(fp2)) == EOF)
				{
					fclose(fp2);
					return -2;
				}
				if ((ret = input_rtcm3obs(&pSinStaObsData[1], &m_RtcmData, (unsigned char)data)))
				{
					break;
				}
			}
			dt = CTime::timediff(pSinStaObsData[0].pObsData[0].tTimeStamp, pSinStaObsData[1].pObsData[0].tTimeStamp);
		} while (dt > 0);
	
	}
	else
	{
		do {
			for (i = 0; i < 4096; i++)
			{
				if ((data = fgetc(fp1)) == EOF)
				{
					fclose(fp1);
					return -2;
				}
				if ((ret = input_rtcm3obs(&pSinStaObsData[0], &m_RtcmData, (unsigned char)data)))
				{
					break;
				}
			}
			dt = CTime::timediff(pSinStaObsData[0].pObsData[0].tTimeStamp, pSinStaObsData[1].pObsData[0].tTimeStamp);
		} while (dt < 0);
	}*/
	



	for (i = 0; i < 10000; i++)
	{
		if ((data = fgetc(fp3)) == EOF)
		{
			fclose(fp3);
			return -2;
		}
		if ((ret = input_rtcm3obs(&pSinStaObsData[2], &m_RtcmData, (unsigned char)data)))
		{
			return ret;
		}
	}
	return 0;
}

bool CRTCM::ReadObsFile(string strFileName1, string strFileName2, string strFileName3)
{
	struct stat tagfileinfo;
	int info;
	if ((fp1 = fopen(strFileName1.c_str(), "rb")) == NULL)
		return false;
	if ((fp2 = fopen(strFileName2.c_str(), "rb")) == NULL)
		return false;
	if ((fp3 = fopen(strFileName3.c_str(), "rb")) == NULL)
		return false;

	info = stat(strFileName1.c_str(), &tagfileinfo);
	if (info != 0)
	{
		perror("显示文件状态信息出错");
	}
	else
	{
		m_RtcmData.time.tiTime = tagfileinfo.st_ctime;
		m_RtcmData.time.dSec = 0.0;
	}
	m_RtcmData.time_s = m_RtcmData.time;
	return true;
}

void CRTCM::InitObs(int n)
{
	if (0 != n)
	{
		m_pSinstaobsdata = new SinStaObsData[n];
	
		for (int i = 0; i < n; i++)
		{
			m_pSinstaobsdata[i].ui16ObsNum = MAXOBS;
			m_pSinstaobsdata[i].pObsData = new ObsData[MAXOBS];
			m_pSinstaobsdata[i].ui16MaxNum = n;
		}
	}
}

void CRTCM::InitRtcm(RtcmData& rtcmdata)
{
	int i, j;
	gtime_t time0 = { 0 };
	rtcmdata.len = rtcmdata.nbit = rtcmdata.nbyte = rtcmdata.obsflag = 0;
	rtcmdata.time = rtcmdata.time_s = time0;
	rtcmdata.opt[0] = '\0';
	rtcmdata.staid = 0;
	rtcmdata.seqno = 0;
	for (i = 0; i < 6; i++) rtcmdata.msmtype[i][0] = '\0';
	for (i = 0; i < MAXSAT; i++) for (j = 0; j < NFREQ + NEXOBS; j++) 
	{
		rtcmdata.cp[i][j] = 0.0;
		rtcmdata.lock[i][j] = rtcmdata.loss[i][j] = 0;
		rtcmdata.lltime[i][j] = time0;
	}
	for (i = 0; i < 400; i++) rtcmdata.nmsg3[i] = 0;
}

void CRTCM::InitStrCov()
{
	m_Strconv.itype = 0; m_Strconv.otype = STRFMT_RTCM3;
	m_Strconv.nmsg =2;/*1074-1084-1094-1124-1005*/
	m_Strconv.msgs[0] = 1074; m_Strconv.msgs[1] = 1005; 

	m_Strconv.tint[0] = 1.0; m_Strconv.tint[1] = 5.0; m_Strconv.tint[2] = 1.0;
	m_Strconv.tint[3] = 1.0; m_Strconv.tint[4] = 1.0; m_Strconv.tint[5] = 5.0;
	InitRtcm(m_Strconv.out);
}

int CRTCM::To_SigId(int sys, unsigned char code)
{
	const char** msm_sig;
	char* sig;
	int i;

	/* signal conversion for undefined signal by rtcm */
	if (sys == SYS_GPS) 
	{
		if (code == CODE_L1Y) code = CODE_L1P;
		else if (code == CODE_L1M) code = CODE_L1P;
		else if (code == CODE_L1N) code = CODE_L1P;
		else if (code == CODE_L2D) code = CODE_L2P;
		else if (code == CODE_L2Y) code = CODE_L2P;
		else if (code == CODE_L2M) code = CODE_L2P;
		else if (code == CODE_L2N) code = CODE_L2P;
	}
	if (!*(sig = code2obs(code))) return 0;
	switch (sys) 
	{
		case SYS_GPS: msm_sig = msm_sig_gps; break;
		case SYS_GLO: msm_sig = msm_sig_glo; break;
		case SYS_GAL: msm_sig = msm_sig_gal; break;
		case SYS_QZS: msm_sig = msm_sig_qzs; break;
		case SYS_SBS: msm_sig = msm_sig_sbs; break;
		case SYS_CMP: msm_sig = msm_sig_cmp; break;
		default: return 0;
	}
	for (i = 0; i < 32; i++) 
	{
		if (!strcmp(sig, msm_sig[i])) return i + 1;
	}
}

int CRTCM::To_SatId(int sys, int sat)
{
	int prn;

	if (SatSys(sat, &prn) != sys) return 0;

	if (sys == SYS_QZS) prn -= MINPRNQZS - 1;
	else if (sys == SYS_SBS) prn -= MINPRNSBS - 1;

	return prn;
}
bool CRTCM::Is_ObsMsg(int msg)
{
	return (1001 <= msg && msg <= 1004) || (1009 <= msg && msg <= 1012) ||
		(1071 <= msg && msg <= 1077) || (1081 <= msg && msg <= 1087) ||
		(1091 <= msg && msg <= 1097) || (1101 <= msg && msg <= 1107) ||
		(1111 <= msg && msg <= 1117) || (1121 <= msg && msg <= 1127);
}

bool CRTCM::Is_Tint(gtime_t time, double tint)
{
	bool a;
	if (tint <= 0.0) return 1;
	double c = fmod(CTime::time2gpst(time, NULL) , tint);
	a=fmod(CTime::time2gpst(time, NULL) + DTTOL, tint) <= 2.0 * DTTOL;
	return a;
}

void CRTCM::SetBitu(unsigned char* buff, int pos, int len, unsigned int data)
{
	unsigned int mask = 1u << (len - 1);
	int i;
	if (len <= 0 || 32 < len) return;
	for (i = pos; i < pos + len; i++, mask >>= 1) {
		if (data & mask) buff[i / 8] |= 1u << (7 - i % 8); else buff[i / 8] &= ~(1u << (7 - i % 8));
	}
}

void CRTCM::SetBits(unsigned char* buff, int pos, int len, int data)
{
	if (data < 0) data |= 1 << (len - 1); else data &= ~(1 << (len - 1)); /* set sign bit */
	SetBitu(buff, pos, len, (unsigned int)data);
}

void CRTCM::Gen_Msm_Index(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int* nsat, int* nsig,
	int* ncell, unsigned char* sat_ind, unsigned char* sig_ind, unsigned char* cell_ind)
{
	int i, j, sat, sig, cell, f;

	*nsat = *nsig = *ncell = 0;
	/* generate satellite and signal index */
	for (i = 0; i < vrsInfo->nNum; i++)
	{
		if (!(sat = To_SatId(sys, vrsInfo->pObsData[i].ui8SatId))) continue;

		for (j = 0; j < NFREQ + NEXOBS; j++) 
		{
			if (!(sig = To_SigId(sys, vrsInfo->pObsData[i].ui8aCodeType[j]))) continue;

			sat_ind[sat - 1] = sig_ind[sig - 1] = 1;
		}
	}
	for (i = 0; i < 64; i++) 
	{
		if (sat_ind[i]) sat_ind[i] = ++(*nsat);
	}
	for (i = 0; i < 32; i++) 
	{
		if (sig_ind[i]) sig_ind[i] = ++(*nsig);
	}

	/* generate cell index */
	for (i = 0; i < vrsInfo->nNum; i++)
	{
		if (!(sat = To_SatId(sys, vrsInfo->pObsData[i].ui8SatId))) continue;

		for (j = 0; j < NFREQ + NEXOBS; j++) 
		{
			if (!(sig = To_SigId(sys, vrsInfo->pObsData[i].ui8aCodeType[j]))) continue;

			cell = sig_ind[sig - 1] - 1 + (sat_ind[sat - 1] - 1) * (*nsig);
			cell_ind[cell] = 1;
		}
	}

	for (i = 0; i < *nsat * (*nsig); i++) 
	{
		if (cell_ind[i] && *ncell < 64) cell_ind[i] = ++(*ncell);
	}
}
/* generate MSM satellite data fields ----------------------------------------*/
void CRTCM::Gen_Msm_Sat(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int nsat, const uint8_t* sat_ind,
	double* rrng, double* rrate, uint8_t* info)
{
	ObsData* data;
	double freq;
	int i, j, k, sat, sig, fcn;

	for (i = 0; i < 64; i++) rrng[i] = rrate[i] = 0.0;

	for (i = 0; i < vrsInfo->nNum; i++)
	{
		data = vrsInfo->pObsData + i;
		//fcn = fcn_glo(data->sat, rtcm); /* fcn+7 */
		fcn = -1;
		if (!(sat = To_SatId(sys, data->ui8SatId))) 
			continue;
		for (j = 0; j < NFREQ + NEXOBS; j++)
		{
			if (!(sig = To_SigId(sys, data->ui8aCodeType[j]))) continue;
			k = sat_ind[sat - 1] - 1;
			freq = code2freq(sys, data->ui8aCodeType[j], fcn - 7);

			/* rough range (ms) and rough phase-range-rate (m/s) */
			if (rrng[k] == 0.0 && data->daPseRange[j] != 0.0) {
				rrng[k] = ROUND(data->daPseRange[j] / RANGE_MS / P2_10) * RANGE_MS * P2_10;
			}
			if (rrate[k] == 0.0 && data->daDoppler[j] != 0.0 && freq > 0.0) {
				rrate[k] = ROUND(-data->daDoppler[j] * CLIGHT / freq) * 1.0;
			}
			/* extended satellite info */
			if (info) info[k] = sys != SYS_GLO ? 0 : (fcn < 0 ? 15 : fcn);
		}
	}
}

double CRTCM::Locktime_d(gtime_t time, gtime_t* lltime, uint8_t LLI)
{
	if (!lltime->tiTime || (LLI & 1)) *lltime = time;
	return CTime::timediff(time, *lltime);
}

int CRTCM::Encode_Msm_Int_Rrng(RtcmData* rtcm, int i, const double* rrng,
	int nsat)
{
	uint32_t int_ms;
	int j;

	for (j = 0; j < nsat; j++)
	{
		if (rrng[j] == 0.0) 
		{
			int_ms = 255;
		}
		else if (rrng[j]<0.0 || rrng[j]>RANGE_MS * 255.0) 
		{
			int_ms = 255;
		}
		else 
		{
			int_ms = ROUND_U(rrng[j] / RANGE_MS / P2_10) >> 10;
		}
		SetBitu(rtcm->buff, i, 8, int_ms); i += 8;
	}
	return i;
}

int CRTCM::Encode_Msm_Mod_Rrng(RtcmData* rtcm, int i, const double* rrng,
	int nsat)
{
	uint32_t mod_ms;
	int j;

	for (j = 0; j < nsat; j++) {
		if (rrng[j] <= 0.0 || rrng[j] > RANGE_MS * 255.0) {
			mod_ms = 0;
		}
		else {
			mod_ms = ROUND_U(rrng[j] / RANGE_MS / P2_10) & 0x3FFu;
		}
		SetBitu(rtcm->buff, i, 10, mod_ms); i += 10;
	}
	return i;
}

int CRTCM::Encode_Msm_Psrng(RtcmData* rtcm, int i, const double* psrng, int ncell)
{
	int j, psrng_val;

	for (j = 0; j < ncell; j++) 
	{
		if (psrng[j] == 0.0) 
		{
			psrng_val = -16384;
		}
		else if (fabs(psrng[j]) > 292.7)
		{
			
			psrng_val = -16384;
		}
		else
		{
			psrng_val = ROUND(psrng[j] / RANGE_MS / P2_24);
		}
		SetBits(rtcm->buff, i, 15, psrng_val); i += 15;
	}
	return i;
}

int CRTCM::Encode_Msm_Phrng(RtcmData* rtcm, int i, const double* phrng, int ncell)
{
	int j, phrng_val;

	for (j = 0; j < ncell; j++) 
	{
		if (phrng[j] == 0.0) 
		{
			phrng_val = -2097152;
		}
		else if (fabs(phrng[j]) > 1171.0)
		{
			
			phrng_val = -2097152;
		}
		else {
			phrng_val = ROUND(phrng[j] / RANGE_MS / P2_29);
		}
		SetBits(rtcm->buff, i, 22, phrng_val); i += 22;
	}
	return i;
}
void CRTCM::set38bits(unsigned char* buff, int pos, double value)
{
	int word_h = (int)floor(value / 64.0);
	unsigned int word_l = (unsigned int)(value - word_h * 64.0);
	SetBits(buff, pos, 32, word_h);
	SetBits(buff, pos + 32, 6, word_l);
}
int CRTCM::To_Msm_Lock(double lock)
{
	if (lock < 0.032) return 0;
	if (lock < 0.064) return 1;
	if (lock < 0.128) return 2;
	if (lock < 0.256) return 3;
	if (lock < 0.512) return 4;
	if (lock < 1.024) return 5;
	if (lock < 2.048) return 6;
	if (lock < 4.096) return 7;
	if (lock < 8.192) return 8;
	if (lock < 16.384) return 9;
	if (lock < 32.768) return 10;
	if (lock < 65.536) return 11;
	if (lock < 131.072) return 12;
	if (lock < 262.144) return 13;
	if (lock < 524.288) return 14;
	return 15;
}

int CRTCM::Encode_Msm_Half_Amb(RtcmData* rtcm, int i, const uint8_t* half,
	int ncell)
{
	int j;

	for (j = 0; j < ncell; j++) 
	{
		SetBitu(rtcm->buff, i, 1, half[j]); i += 1;
	}
	return i;
}

int CRTCM::Encode_Msm_Cnr(RtcmData* rtcm, int i, const float* cnr, int ncell)
{
	int j, cnr_val;

	for (j = 0; j < ncell; j++) {
		cnr_val = ROUND(cnr[j] / 1.0);
		SetBitu(rtcm->buff, i, 6, cnr_val); i += 6;
	}
	return i;
}

int CRTCM::Encode_Msm_Lock(RtcmData* rtcm, int i, const double* lock, int ncell)
{
	int j, lock_val;

	for (j = 0; j < ncell; j++) {
		lock_val = To_Msm_Lock(lock[j]);
		SetBitu(rtcm->buff, i, 4, lock_val); i += 4;
	}
	return i;
}

void CRTCM::Gen_Msm_Sig(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int nsat, int nsig, int ncell, const uint8_t* sat_ind,
	const uint8_t* sig_ind, const uint8_t* cell_ind, const double* rrng,
	const double* rrate, double* psrng, double* phrng,
	double* rate, double* lock, uint8_t* half, float* cnr)
{
	ObsData* data;
	double freq, lambda, psrng_s, phrng_s, rate_s, lt;
	int i, j, k, sat, sig, fcn, cell, LLI;
	for (i = 0; i < ncell; i++) 
	{
		if (psrng) psrng[i] = 0.0;
		if (phrng) phrng[i] = 0.0;
		if (rate) rate[i] = 0.0;
	}

	for (i = 0; i < vrsInfo->nNum; i++)
	{
		data = vrsInfo->pObsData + i;
		//fcn = fcn_glo(data->sat, rtcm); /* fcn+7 */
		fcn = -1;

		if (!(sat = To_SatId(sys, data->ui8SatId))) continue;

		for (j = 0; j < NFREQ + NEXOBS; j++)
		{
			if (!(sig = To_SigId(sys, data->ui8aCodeType[j]))) continue;

			k = sat_ind[sat - 1] - 1;
			if ((cell = cell_ind[sig_ind[sig - 1] - 1 + k * nsig]) >= 64) continue;
	
			freq = code2freq(sys, data->ui8aCodeType[j], fcn - 7);
			lambda = freq == 0.0 ? 0.0 : CLIGHT / freq;
			psrng_s = data->daPseRange[j] == 0.0 ? 0.0 : data->daPseRange[j] - rrng[k];
			phrng_s = data->daCarrPhase[j] == 0.0 || lambda <= 0.0 ? 0.0 : data->daCarrPhase[j] * lambda - rrng[k];
			rate_s = data->daDoppler[j] == 0.0 || lambda <= 0.0 ? 0.0 : -data->daDoppler[j] * lambda - rrate[k];
			/* subtract phase - psudorange integer cycle offset */
			LLI = data->ui8aLLI[j];
			if ((LLI & 1) || fabs(phrng_s - rtcmdata->cp[data->ui8SatId - 1][j]) > 1171.0) 
			{
				rtcmdata->cp[data->ui8SatId - 1][j] = ROUND(phrng_s / lambda) * lambda;
				LLI |= 1;
			}

			phrng_s -= rtcmdata->cp[data->ui8SatId - 1][j];

			lt = Locktime_d(data->tTimeStamp, rtcmdata->lltime[data->ui8SatId - 1] + j, LLI);

			if(psrng && psrng_s != 0.0) psrng[cell - 1] = psrng_s;
			if (phrng && phrng_s != 0.0) phrng[cell - 1] = phrng_s;
			if (rate && rate_s != 0.0) rate[cell - 1] = rate_s;
			if (lock) lock[cell - 1] = lt;
			if (half) half[cell - 1] = (data->ui8aLLI[j] & 2) ? 1 : 0;
			if (cnr) cnr[cell - 1] = (float)(data->daSigNoiRatio[j] * SNR_UNIT);
		}
	}
}

int CRTCM::Encode_Msm_Head(int type, VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int sync, int* nsat,
	int* ncell, double* rrng, double* rrate, unsigned char* info, double* psrng, double* phrng,
	double* rate, double* lock, unsigned char* half, float* cnr)
{
	double tow;
	unsigned char sat_ind[64] = { 0 }, sig_ind[32] = { 0 }, cell_ind[32 * 64] = { 0 };
	unsigned int dow, epoch;
	int i = 24, j, nsig = 0;
	char tstr[64];
	
	switch (sys) 
	{
		case SYS_GPS: type += 1070; break;
		case SYS_GLO: type += 1080; break;
		case SYS_GAL: type += 1090; break;
		case SYS_QZS: type += 1110; break;
		case SYS_SBS: type += 1100; break;
		case SYS_CMP: type += 1120; break;
		default: return 0;
	}
	/* generate msm satellite, signal and cell index */
	Gen_Msm_Index(vrsInfo, rtcmdata,sys, nsat, &nsig, ncell, sat_ind, sig_ind, cell_ind);
	CTime::time2str(rtcmdata->time, tstr, 5);
	if (sys == SYS_GLO)
	{
		/* GLONASS time (dow + tod-ms) */
		tow = CTime::time2gpst(CTime::timeadd(CTime::gpst2utc(rtcmdata->time), 10800.0), NULL);
		dow = (uint32_t)(tow / 86400.0);
		epoch = (dow << 27) + ROUND_U(fmod(tow, 86400.0) * 1E3);
	}
	else if (sys == SYS_CMP) 
	{
		/* BDS time (tow-ms) */
		epoch = ROUND_U(CTime::time2gpst(CTime::gpst2bdt(rtcmdata->time), NULL) * 1E3);
	}
	else 
	{
		/* GPS, QZSS, Galileo and IRNSS time (tow-ms) */
		epoch = ROUND_U(CTime::time2gpst(rtcmdata->time, NULL) * 1E3);
	}
	/* encode msm header (ref [15] table 3.5-78) */
	SetBitu(rtcmdata->buff, i, 12, type); i += 12; /* message number */
	SetBitu(rtcmdata->buff, i, 12, rtcmdata->staid); i += 12; /* reference station id */
	SetBitu(rtcmdata->buff, i, 30, epoch); i += 30; /* epoch time */
	SetBitu(rtcmdata->buff, i, 1, sync); i += 1; /* multiple message bit */
	SetBitu(rtcmdata->buff, i, 3, rtcmdata->seqno); i += 3; /* issue of data station */
	SetBitu(rtcmdata->buff, i, 7, 0); i += 7; /* reserved */
	SetBitu(rtcmdata->buff, i, 2, 0); i += 2; /* clock streering indicator */
	SetBitu(rtcmdata->buff, i, 2, 0); i += 2; /* external clock indicator */
	SetBitu(rtcmdata->buff, i, 1, 0); i += 1; /* smoothing indicator */
	SetBitu(rtcmdata->buff, i, 3, 0); i += 3; /* smoothing interval */

	 /* satellite mask */
	for (j = 0; j < 64; j++) 
	{
		SetBitu(rtcmdata->buff, i, 1, sat_ind[j] ? 1 : 0); i += 1;
	}
	/* signal mask */
	for (j = 0; j < 32; j++) 
	{
		SetBitu(rtcmdata->buff, i, 1, sig_ind[j] ? 1 : 0); i += 1;
	}
	/* cell mask */
	for (j = 0; j < *nsat * nsig && j < 64; j++) 
	{
		SetBitu(rtcmdata->buff, i, 1, cell_ind[j] ? 1 : 0); i += 1;
	}

	/* generate msm satellite data fields */
	Gen_Msm_Sat(vrsInfo, rtcmdata, sys, *nsat, sat_ind, rrng, rrate, info);

	/*generate msm signal data fields*/
	Gen_Msm_Sig(vrsInfo, rtcmdata, sys, *nsat, nsig, *ncell, sat_ind, sig_ind, cell_ind, rrng, rrate,
		psrng, phrng, rate, lock, half, cnr);

	return i;
}


int CRTCM::Encode_Type1005(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sync)
{
	double* p = vrsInfo->daCoorVal;
	int i = 24;
	SetBitu(rtcmdata->buff, i, 12, 1005); i += 12; /* message no */
	SetBitu(rtcmdata->buff, i, 12, rtcmdata->staid); i += 12; /* ref station id */
	SetBitu(rtcmdata->buff, i, 6, 0); i += 6; /* itrf realization year */
	SetBitu(rtcmdata->buff, i, 1, 1); i += 1; /* gps indicator */
	SetBitu(rtcmdata->buff, i, 1, 1); i += 1; /* glonass indicator */
	SetBitu(rtcmdata->buff, i, 1, 1); i += 1; /* galileo indicator */
	SetBitu(rtcmdata->buff, i, 1, 0); i += 1; /* ref station indicator */
	set38bits(rtcmdata->buff, i, p[0] / 0.0001); i += 38; /* antenna ref point ecef-x */

	return 0;
}



int CRTCM::Encode_Msm4(VRSInfo* vrsInfo, RtcmData* rtcmdata, int sys, int sync)
{
	double rrng[64], rrate[64], psrng[64], phrng[64], lock[64];
	float cnr[64];
	unsigned char half[64];
	int i, nsat, ncell;
	/* encode msm header */
	if (!(i = Encode_Msm_Head(4, vrsInfo,rtcmdata, sys, sync, &nsat, &ncell, rrng, rrate, NULL, psrng,
		phrng, NULL, lock, half, cnr))) 
	{
		return 0;
	}

	/* encode msm satellite data */
	i = Encode_Msm_Int_Rrng(rtcmdata, i, rrng, nsat); /* rough range integer ms */
	i = Encode_Msm_Mod_Rrng(rtcmdata, i, rrng, nsat); /* rough range modulo 1 ms */

	 /* encode msm signal data */
	i = Encode_Msm_Psrng(rtcmdata, i, psrng, ncell); /* fine pseudorange */
	i = Encode_Msm_Phrng(rtcmdata, i, phrng, ncell); /* fine phase-range */
	i = Encode_Msm_Lock(rtcmdata, i, lock, ncell); /* lock-time indicator */
	i = Encode_Msm_Half_Amb(rtcmdata, i, half, ncell); /* half-cycle-amb indicator */
	i = Encode_Msm_Cnr(rtcmdata, i, cnr, ncell); /* signal cnr */
	rtcmdata->nbit = i;
	return 1;

}
int CRTCM::Encode_Rtcm3(VRSInfo* vrsInfo, RtcmData* rtcmdata, int type, int sync)
{
	int ret = 0;
	switch (type)
	{
		case 1005: ret = Encode_Type1005(vrsInfo,rtcmdata, sync);     break;
		//case 1006: ret = encode_type1006(rtcm, sync);     break;
		case 1074: ret = Encode_Msm4(vrsInfo,rtcmdata, SYS_GPS, sync); break;

	}
	if (ret > 0) 
	{
		type -= 1000;
		if (1 <= type && type <= 299) rtcmdata->nmsg3[type]++; /* 1001-1299 */
		else if (1000 <= type && type <= 1099) rtcmdata->nmsg3[type - 700]++; /* 2000-2099 */
		else rtcmdata->nmsg3[0]++;
	}
	return ret;
}


int CRTCM::Gen_Rtcm3(VRSInfo* vrsInfo,RtcmData *rtcm, int type, int sync)
{
	unsigned int crc;
	int i = 0;
	rtcm->nbit = rtcm->len = rtcm->nbyte = 0;
	rtcm->time = vrsInfo->pObsData[0].tTimeStamp;//____0723____
	/* set preamble and reserved */
	SetBitu(rtcm->buff, i, 8, RTCM3PREAMB); i += 8;
	SetBitu(rtcm->buff, i, 6, 0); i += 6;
	SetBitu(rtcm->buff, i, 10, 0); i += 10;
	/* encode rtcm 3 message body */
	if (!Encode_Rtcm3(vrsInfo,rtcm, type, sync)) return 0;
	/* padding to align 8 bit boundary */
	for (i = rtcm->nbit; i % 8; i++) 
	{
		SetBitu(rtcm->buff, i, 1, 0);
	}
	/* message length (header+data) (bytes) */
	if ((rtcm->len = i / 8) >= 3 + 1024) 
	{
		rtcm->nbit = rtcm->len = 0;
		return 0;
	}
	/* message length without header and parity */
	SetBitu(rtcm->buff, 14, 10, rtcm->len - 3);

	/* crc-24q */
	crc = rtk_crc24q(rtcm->buff, rtcm->len);
	SetBitu(rtcm->buff, i, 24, crc);

	/* length total (bytes) */
	rtcm->nbyte = rtcm->len + 3;

	return 1;
}

void CRTCM::Write_rtcm3_msm(RtcmData* out,VRSInfo* vrsInfo, int msg, int sync)
{
	ObsData* data, buff[MAXOBS];
	int i, j, n, ns, sys, nobs, code, nsat = 0, nsig = 0, nmsg, mask[MAXCODE] = { 0 };
	if (1071 <= msg && msg <= 1077) sys = SYS_GPS;
	else if (1081 <= msg && msg <= 1087) sys = SYS_GLO;
	else if (1091 <= msg && msg <= 1097) sys = SYS_GAL;
	else if (1101 <= msg && msg <= 1107) sys = SYS_SBS;
	else if (1111 <= msg && msg <= 1117) sys = SYS_QZS;
	else if (1121 <= msg && msg <= 1127) sys = SYS_CMP;
	else return;

	data = vrsInfo->pObsData;
	nobs = vrsInfo->nNum;
	/* count number of satellites and signals */
	for (i = 0; i < nobs && i < MAXOBS; i++) {
		if (SatSys(data[i].ui8SatId, NULL) != sys) continue;
		nsat++;
		for (j = 0; j < NFREQ + NEXOBS; j++) {
			if (!(code = data[i].ui8aCodeType[j]) || mask[code - 1]) continue;
			mask[code - 1] = 1;
			nsig++;
		}
	}
	if (nsig <= 0 || nsig > 64) return;

	/* pack data to multiple messages if nsat x nsig > 64 */
	ns = 64 / nsig;         /* max number of sats in a message */
	nmsg = (nsat - 1) / ns + 1; /* number of messages */

	vrsInfo->pObsData = buff;
	for (i = j = 0; i < nmsg; i++) 
	{
		for (n = 0; n < ns && j < nobs && j < MAXOBS; j++) 
		{
			if (SatSys(data[j].ui8SatId, NULL) != sys) continue;
			vrsInfo->pObsData[n++] = data[j];
		}
		vrsInfo->nNum = n;

		if (Gen_Rtcm3(vrsInfo,out,msg, i < nmsg - 1 ? 1 : sync))
		{
			//strwrite(str, out->buff, out->nbyte);
		}
	}
	vrsInfo->pObsData = data;
	vrsInfo->nNum = nobs;
}
void CRTCM::Write_Obs(gtime_t time, VRSInfo* vrsInfo, StrConv* conv)
{
	int i, j = 0;
	int nbyte = 0;
	char tstr[64];
	CTime::time2str(time, tstr, 5);
	for (i = 0; i < conv->nmsg; i++) {
		if (!Is_ObsMsg(conv->msgs[i]) ||! Is_Tint(time, conv->tint[i])) 
			continue;

		j = i; /* index of last message */
	}
	for (i = 0; i < conv->nmsg; i++)
	{
		if ( !Is_Tint(time, conv->tint[i]))
			continue;
		if (conv->otype == STRFMT_RTCM3)
		{
			if (conv->msgs[i] <= 1012)
			{
				if(!Gen_Rtcm3(vrsInfo,&conv->out,conv->msgs[i],i!=j))//VRSInfo* vrsInfo, RtcmData* rtcm, int type, int sync
					continue;
			}
			else
			{
				Write_rtcm3_msm(&conv->out, vrsInfo, conv->msgs[i], i != j);
			}
		}
	}
}

