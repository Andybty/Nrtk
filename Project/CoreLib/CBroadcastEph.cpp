
#include "CBroadcastEph.h"
#include "CommonFunction.h"

/* ephemeris selections */
static int eph_sel[] = {
   // GPS,GLO,GAL,QZS,BDS,IRN,SBS 
    0,0,0,0,0,0,0
};
/** 
  @brief 广播星历类的构造函数
  */

CBroadcastEph::CBroadcastEph()
{

}
/**
  @brief 广播星历类的析构函数
  */

CBroadcastEph::~CBroadcastEph()
{

}

/**
 @brief   根据星历得到相应的卫星系
 @param   sys[in]卫星系统数量
 @返回值  卫星系统标识
*/

int CBroadcastEph::GetEphSys(int sys)
{
    //判断卫星系统
    switch (sys) {
    case SYS_GPS: return eph_sel[0];
    case SYS_GLO: return eph_sel[1];
    case SYS_GAL: return eph_sel[2];
    case SYS_QZS: return eph_sel[3];
    case SYS_CMP: return eph_sel[4];
    case SYS_IRN: return eph_sel[5];
    case SYS_SBS: return eph_sel[6];
    }
    return 0;
}
 /**
  @brief   根据ure参数确定卫星位置的精度方差
  @param   sys[in]卫星系统
  @param   ura[in]卫星位置精度参数
  @返回值  卫星位置的精度方差
  */
double CBroadcastEph::EphuraVar(int sys, int ura)
{
    //初始化URA值
    const double ura_value[] = {
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,
        96.0,192.0,384.0,768.0,1536.0,3072.0,6144.0
    };
    if (sys == SYS_GAL) { 
        //根据GLO系统ura值计算GLO卫星位置精度方差
        if (ura <= 49) return SQR(ura * 0.01);
        if (ura <= 74) return SQR(0.5 + (ura - 50) * 0.02);
        if (ura <= 99) return SQR(1.0 + (ura - 75) * 0.04);
        if (ura <= 125) return SQR(2.0 + (ura - 100) * 0.16);
        return SQR(STD_GAL_NAPA);
    }
       // 否则返回GPS系统ura值的卫星位置精度方差
    else { 
        return ura < 0 || 14 < ura ? SQR(6144.0) : SQR(ura_value[ura]);
    }
}

	
void CBroadcastEph::CalSatClk(gtime_t tiEphTime, SinStaObsData* pSinSatObsData, 
     AllBroadcastEphData* pBroadcastEphData, double* pSatClk, double* pSatClkVar)
{

}
 /**
  @brief   计算单颗卫星位置，并给出健康状态
  @param   time[in]GPS时间
  @param   eph[in]广播星历时间
  @param   sat[in]卫星ID
  @param   pBroadcastEphData[in]广播星历
  @param   iode[in]广播星历数据龄期
  @param   rs[in]卫星位置
  @param   dts[in]卫星钟差
  @param   var[in]卫星位置和钟方差
  @param   svh[in]卫星健康状态
  @返回值  解算完成返回1，否则返回0
  */
int CBroadcastEph::EphPos(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData,
                           int iode, double* rs, double* dts, double* var, int* svh)
{
    BroadcastEphData* eph;
    double rst[3], dtst[1], tt = 1E-3;
    int i, sys;
    if (!(eph = SelEph(teph, sat, -1, pBroadcastEphData))) return 0;
    Eph2Pos(time, eph, rs, dts, var);
    time = CTime::timeadd(time, tt);
    Eph2Pos(time, eph, rst, dtst, var);
    *svh = eph->svh;
    return 1;
}
 /**
  @brief  选择计算单颗卫星位置模式
  @param  time[in]GPS时间
  @param  teph[in]GPS星历时间
  @param  sat[in]卫星数
  @param  ephopt[in]星历选项
  @param  AllBroadcastEphData[in]星历数据
  @param  rs[out]卫星位置
  @param  dts[out]卫星钟差
  @param  var[out]卫星位置钟差方差
  @param  svh[out]卫星健康状态
  @返回值 解算完成返回1，否则返回0
  */

int  CBroadcastEph::SatPos(gtime_t time, gtime_t teph, int sat, int ephopt, AllBroadcastEphData* pALLBroadcastEphData,
                           double* rs, double* dts, double* var, int* svh)
{
    *svh = 0;
    //判断是广播星历还是精密星历
    switch (ephopt)
    {
    case EPHOPT_BRDC: 
    //利用广播星历计算卫星位置
      return EphPos(time, teph, sat, pALLBroadcastEphData, -1, rs, dts, var, svh);
    }
     *svh = -1;
     return 0;
}
 /**
  @brief   选择与观测时间最近的星历
  @param  time[in]GPS时间
  @param  sat[in]卫星数
  @param  iode[out]卫星数据龄期
  @param  AllBroadcastEphData[in]星历数据
  @返回值  卫星广播星历结构体
 */
BroadcastEphData* CBroadcastEph::SelEph(gtime_t time, int sat, int iode, AllBroadcastEphData* pAllBroadcastEphData)
{
    double t, tmax, tmin;
    int i, j , sys, sel;
    j = -1;
    sys = SatSys(sat, NULL);
    //判断不同卫星系统观测时间历元间隔
    switch (sys) {
    case SYS_GPS: tmax = MAXDTOE + 1.0; sel = eph_sel[0]; break;
    case SYS_GAL: tmax = MAXDTOE_GAL; sel = eph_sel[2]; break;
    case SYS_QZS: tmax = MAXDTOE_QZS + 1.0; sel = eph_sel[3]; break;
    case SYS_CMP: tmax = MAXDTOE_CMP + 1.0; sel = eph_sel[4]; break;
    case SYS_IRN: tmax = MAXDTOE_IRN + 1.0; sel = eph_sel[5]; break;
    default: tmax = MAXDTOE + 1.0; break;
    }
    tmin = tmax + 1.0;
    for (i = 0; i < pAllBroadcastEphData->nSatNum; i++) {
        if (pAllBroadcastEphData->pBroadcastEphData[i].sat != sat) continue;
        if (iode >= 0 && pAllBroadcastEphData->pBroadcastEphData[i].iode != iode) continue;
        if (sys == SYS_GAL) {
            sel = GetEphSys(SYS_GAL);
            if (sel == 0 && !(pAllBroadcastEphData->pBroadcastEphData[i].code & (1 << 9))) continue; /* I/NAV */
            if (sel == 1 && !(pAllBroadcastEphData->pBroadcastEphData[i].code & (1 << 8))) continue; /* F/NAV */
            if (CTime::timediff(pAllBroadcastEphData->pBroadcastEphData[i].toe, time) >= 0.0) continue; /* AOD<=0 */
        }
        if ((t = fabs(CTime::timediff(pAllBroadcastEphData->pBroadcastEphData[i].toe, time))) > tmax) continue;
        if (iode >= 0) return pAllBroadcastEphData->pBroadcastEphData + i;
        //判断与toe最近的时间
        if (t <= tmin) { j = i; tmin = t; } 
    }
    if (iode >= 0 || j < 0) {
        return NULL;
    }
    //返回最近的观测历元星历
    return pAllBroadcastEphData->pBroadcastEphData + j;;
}
 /**
  @brief   选择计算单颗卫星钟差模式
  @param  time[in]GPS时间
  @param  teph[in]GPS星历时间
  @param  sat[in]卫星数
  @param  AllBroadcastEphData[in]星历数据
  @param  dts[out]卫星钟差
  @返回值 解算完成返回1，否则返回0
 */
int  CBroadcastEph::EphClk(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData, double* dts)
{
    BroadcastEphData * eph;
    double t, ts;
    int i;
    if (!(eph = SelEph(teph, sat, -1, pBroadcastEphData))) return 0;
    t = ts = CTime::timediff(time, eph->toc);
    //根据f0/f1/f2星历参数计算卫星钟差
    for (i = 0; i < 2; i++) 
    {
        t = ts - (eph->f0 + eph->f1 * t + eph->f2 * t * t);
    }
    *dts= eph->f0 + eph->f1 * t + eph->f2 * t * t;
    return 1;
}
 /**
  @brief   计算单颗卫星位置
  @param   time[in]GPS时间
  @param   eph[in]广播星历时间
  @param   rs[in]卫星位置
  @param   dts[in]卫星钟差
  @param   var[in]卫星位置和钟方差
  @返回值  无
 */
void CBroadcastEph::Eph2Pos(gtime_t time, BroadcastEphData* eph, double* rs, double* dts, double* var)
{
    double tk, M, E, Ek, sinE, cosE, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, cosi, mu, omge;
    double xg, yg, zg, sino, coso;
    int n, sys, prn;

    if (eph->A <= 0.0) {
        rs[0] = rs[1] = rs[2] = *dts = *var = 0.0;
        return;
    }
    //归化时间
    tk = CTime::timediff(time, eph->toe);
    switch ((sys = SatSys(eph->sat, &prn))) 
    {
        case SYS_GAL: mu = MU_GAL; omge = OMGE_GAL; break;
        case SYS_CMP: mu = MU_CMP; omge = OMGE_CMP; break;
        default:      mu = MU_GPS; omge = OMGE;     break;
    }
    //计算卫星平近点角
    M = eph->M0 + (sqrt(mu / (eph->A * eph->A * eph->A)) + eph->deln) * tk;
    //迭代计算偏近点角
    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
        Ek = E; E -= (E - eph->e * sin(E) - M) / (1.0 - eph->e * cos(E));
    }
    if (n >= MAX_ITER_KEPLER) {
        return;
    }
    sinE = sin(E); cosE = cos(E);
    //计算真近点角和升叫距角
    u = atan2(sqrt(1.0 - eph->e * eph->e) * sinE, cosE - eph->e) + eph->omg;
    //摄动改正计算升交距角/卫星矢量/轨道倾角
    r = eph->A * (1.0 - eph->e * cosE);
    i = eph->i0 + eph->idot * tk;
    sin2u = sin(2.0 * u); cos2u = cos(2.0 * u);
    u += eph->cus * sin2u + eph->cuc * cos2u;
    r += eph->crs * sin2u + eph->crc * cos2u;
    i += eph->cis * sin2u + eph->cic * cos2u;
    // 卫星在轨道面坐标系中的坐标
    x = r * cos(u); y = r * sin(u); 
    //计算北斗GEO卫星坐标
    cosi = cos(i);
    if (sys == SYS_CMP && (prn <= 5 || prn >= 59)) { /* ref [9] table 4-1 */
        O = eph->OMG0 + eph->OMGd * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        xg = x * cosO - y * cosi * sinO;
        yg = x * sinO + y * cosi * cosO;
        zg = y * sin(i);
        sino = sin(omge * tk); coso = cos(omge * tk);
        rs[0] = xg * coso + yg * sino * COS_5 + zg * sino * SIN_5;
        rs[1] = -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5;
        rs[2] = -yg * SIN_5 + zg * COS_5;
    }
    else {
        // 观测时刻t的升交点经度
        O = eph->OMG0 + (eph->OMGd - omge) * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        //卫星在WGS-84坐标系中的坐标
        rs[0] = x * cosO - y * cosi * sinO;
        rs[1] = x * sinO + y * cosi * cosO;
        rs[2] = y * sin(i);
    }
    tk = CTime::timediff(time, eph->toc);
    *dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;
    //卫星钟差
    *dts -= 2.0 * sqrt(mu * eph->A) * eph->e * sinE / SQR(CLIGHT);

   //卫星卫星和钟差方差
    *var = EphuraVar(sys, eph->sva);
}
 /**
  @brief  计算卫星位置钟差
  @param  gtime_t[in]GPS星历时间
  @param  SinStaObsData[in]观测数据时刻
  @param  AllBroadcastEphData[in]星历数据
  @param  rs[out]卫星位置
  @param  dts[out]卫星钟差
  @param  var[out]卫星位置钟差方差
  @param  svh[out]卫星健康状态
  @返回值 共视卫星数量
 */
int CBroadcastEph::CalSatPos(gtime_t tiEphTime, SinStaObsData* pSinSatObsData,AllBroadcastEphData* pBroadcastEphData,
                             double* rs, double* dts, double* var, int* svh)
{
   
    int i, stat;
    gtime_t time[2 * MAXOBS] = { 0 };
    double dt, pr;
    int j,k=0;
    int n = pSinSatObsData->ui16ObsNum;

    for (i = 0; i < n && i < 2 * MAXOBS; i++) 
    {
        for (j = 0; j < 3; j++) 
            rs[j + i * 3] = 0.0;
        for (j = 0; j < 2; j++) 
            dts[j + i * 2] = 0.0;

        var[i] = 0.0; svh[i] = 0;

        //判断观测数据伪距不为空
        for (j = 0, pr = 0.0; j < NFREQ; j++)
          if ((pr = pSinSatObsData->pObsData[i].daPseRange[j]) != 0.0) break;

        if (j >= NFREQ) {
            continue;
        }
        //根据卫星钟计算发射时间
        time[i] = CTime::timeadd(pSinSatObsData->pObsData[i].tTimeStamp, -pr / CLIGHT);

        //根据广播星历计算卫星钟差
        if (!EphClk(time[i], tiEphTime, pSinSatObsData->pObsData[i].ui8SatId, pBroadcastEphData, &dt)) {
            continue;
        }
        time[i] = CTime::timeadd(time[i], -dt);

        //计算卫星发射时刻的位置和钟
        if (!SatPos(time[i], tiEphTime, pSinSatObsData->pObsData[i].ui8SatId, EPHOPT_BRDC, 
            pBroadcastEphData, rs + i * 3, dts + i * 2, var + i, svh + i)) {
            continue;
        }
        //共视卫星个数
        k++;
    }
    return k;
}