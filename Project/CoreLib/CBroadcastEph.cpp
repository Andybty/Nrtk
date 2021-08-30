
#include "CBroadcastEph.h"
#include "CommonFunction.h"

/* ephemeris selections */
static int eph_sel[] = {
   // GPS,GLO,GAL,QZS,BDS,IRN,SBS 
    0,0,0,0,0,0,0
};
/** 
  @brief �㲥������Ĺ��캯��
  */

CBroadcastEph::CBroadcastEph()
{

}
/**
  @brief �㲥���������������
  */

CBroadcastEph::~CBroadcastEph()
{

}

/**
 @brief   ���������õ���Ӧ������ϵ
 @param   sys[in]����ϵͳ����
 @����ֵ  ����ϵͳ��ʶ
*/

int CBroadcastEph::GetEphSys(int sys)
{
    //�ж�����ϵͳ
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
  @brief   ����ure����ȷ������λ�õľ��ȷ���
  @param   sys[in]����ϵͳ
  @param   ura[in]����λ�þ��Ȳ���
  @����ֵ  ����λ�õľ��ȷ���
  */
double CBroadcastEph::EphuraVar(int sys, int ura)
{
    //��ʼ��URAֵ
    const double ura_value[] = {
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,
        96.0,192.0,384.0,768.0,1536.0,3072.0,6144.0
    };
    if (sys == SYS_GAL) { 
        //����GLOϵͳuraֵ����GLO����λ�þ��ȷ���
        if (ura <= 49) return SQR(ura * 0.01);
        if (ura <= 74) return SQR(0.5 + (ura - 50) * 0.02);
        if (ura <= 99) return SQR(1.0 + (ura - 75) * 0.04);
        if (ura <= 125) return SQR(2.0 + (ura - 100) * 0.16);
        return SQR(STD_GAL_NAPA);
    }
       // ���򷵻�GPSϵͳuraֵ������λ�þ��ȷ���
    else { 
        return ura < 0 || 14 < ura ? SQR(6144.0) : SQR(ura_value[ura]);
    }
}

	
void CBroadcastEph::CalSatClk(gtime_t tiEphTime, SinStaObsData* pSinSatObsData, 
     AllBroadcastEphData* pBroadcastEphData, double* pSatClk, double* pSatClkVar)
{

}
 /**
  @brief   ���㵥������λ�ã�����������״̬
  @param   time[in]GPSʱ��
  @param   eph[in]�㲥����ʱ��
  @param   sat[in]����ID
  @param   pBroadcastEphData[in]�㲥����
  @param   iode[in]�㲥������������
  @param   rs[in]����λ��
  @param   dts[in]�����Ӳ�
  @param   var[in]����λ�ú��ӷ���
  @param   svh[in]���ǽ���״̬
  @����ֵ  ������ɷ���1�����򷵻�0
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
  @brief  ѡ����㵥������λ��ģʽ
  @param  time[in]GPSʱ��
  @param  teph[in]GPS����ʱ��
  @param  sat[in]������
  @param  ephopt[in]����ѡ��
  @param  AllBroadcastEphData[in]��������
  @param  rs[out]����λ��
  @param  dts[out]�����Ӳ�
  @param  var[out]����λ���Ӳ��
  @param  svh[out]���ǽ���״̬
  @����ֵ ������ɷ���1�����򷵻�0
  */

int  CBroadcastEph::SatPos(gtime_t time, gtime_t teph, int sat, int ephopt, AllBroadcastEphData* pALLBroadcastEphData,
                           double* rs, double* dts, double* var, int* svh)
{
    *svh = 0;
    //�ж��ǹ㲥�������Ǿ�������
    switch (ephopt)
    {
    case EPHOPT_BRDC: 
    //���ù㲥������������λ��
      return EphPos(time, teph, sat, pALLBroadcastEphData, -1, rs, dts, var, svh);
    }
     *svh = -1;
     return 0;
}
 /**
  @brief   ѡ����۲�ʱ�����������
  @param  time[in]GPSʱ��
  @param  sat[in]������
  @param  iode[out]������������
  @param  AllBroadcastEphData[in]��������
  @����ֵ  ���ǹ㲥�����ṹ��
 */
BroadcastEphData* CBroadcastEph::SelEph(gtime_t time, int sat, int iode, AllBroadcastEphData* pAllBroadcastEphData)
{
    double t, tmax, tmin;
    int i, j , sys, sel;
    j = -1;
    sys = SatSys(sat, NULL);
    //�жϲ�ͬ����ϵͳ�۲�ʱ����Ԫ���
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
        //�ж���toe�����ʱ��
        if (t <= tmin) { j = i; tmin = t; } 
    }
    if (iode >= 0 || j < 0) {
        return NULL;
    }
    //��������Ĺ۲���Ԫ����
    return pAllBroadcastEphData->pBroadcastEphData + j;;
}
 /**
  @brief   ѡ����㵥�������Ӳ�ģʽ
  @param  time[in]GPSʱ��
  @param  teph[in]GPS����ʱ��
  @param  sat[in]������
  @param  AllBroadcastEphData[in]��������
  @param  dts[out]�����Ӳ�
  @����ֵ ������ɷ���1�����򷵻�0
 */
int  CBroadcastEph::EphClk(gtime_t time, gtime_t teph, int sat, AllBroadcastEphData* pBroadcastEphData, double* dts)
{
    BroadcastEphData * eph;
    double t, ts;
    int i;
    if (!(eph = SelEph(teph, sat, -1, pBroadcastEphData))) return 0;
    t = ts = CTime::timediff(time, eph->toc);
    //����f0/f1/f2�����������������Ӳ�
    for (i = 0; i < 2; i++) 
    {
        t = ts - (eph->f0 + eph->f1 * t + eph->f2 * t * t);
    }
    *dts= eph->f0 + eph->f1 * t + eph->f2 * t * t;
    return 1;
}
 /**
  @brief   ���㵥������λ��
  @param   time[in]GPSʱ��
  @param   eph[in]�㲥����ʱ��
  @param   rs[in]����λ��
  @param   dts[in]�����Ӳ�
  @param   var[in]����λ�ú��ӷ���
  @����ֵ  ��
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
    //�黯ʱ��
    tk = CTime::timediff(time, eph->toe);
    switch ((sys = SatSys(eph->sat, &prn))) 
    {
        case SYS_GAL: mu = MU_GAL; omge = OMGE_GAL; break;
        case SYS_CMP: mu = MU_CMP; omge = OMGE_CMP; break;
        default:      mu = MU_GPS; omge = OMGE;     break;
    }
    //��������ƽ�����
    M = eph->M0 + (sqrt(mu / (eph->A * eph->A * eph->A)) + eph->deln) * tk;
    //��������ƫ�����
    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
        Ek = E; E -= (E - eph->e * sin(E) - M) / (1.0 - eph->e * cos(E));
    }
    if (n >= MAX_ITER_KEPLER) {
        return;
    }
    sinE = sin(E); cosE = cos(E);
    //���������Ǻ����о��
    u = atan2(sqrt(1.0 - eph->e * eph->e) * sinE, cosE - eph->e) + eph->omg;
    //�㶯���������������/����ʸ��/������
    r = eph->A * (1.0 - eph->e * cosE);
    i = eph->i0 + eph->idot * tk;
    sin2u = sin(2.0 * u); cos2u = cos(2.0 * u);
    u += eph->cus * sin2u + eph->cuc * cos2u;
    r += eph->crs * sin2u + eph->crc * cos2u;
    i += eph->cis * sin2u + eph->cic * cos2u;
    // �����ڹ��������ϵ�е�����
    x = r * cos(u); y = r * sin(u); 
    //���㱱��GEO��������
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
        // �۲�ʱ��t�������㾭��
        O = eph->OMG0 + (eph->OMGd - omge) * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        //������WGS-84����ϵ�е�����
        rs[0] = x * cosO - y * cosi * sinO;
        rs[1] = x * sinO + y * cosi * cosO;
        rs[2] = y * sin(i);
    }
    tk = CTime::timediff(time, eph->toc);
    *dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;
    //�����Ӳ�
    *dts -= 2.0 * sqrt(mu * eph->A) * eph->e * sinE / SQR(CLIGHT);

   //�������Ǻ��Ӳ��
    *var = EphuraVar(sys, eph->sva);
}
 /**
  @brief  ��������λ���Ӳ�
  @param  gtime_t[in]GPS����ʱ��
  @param  SinStaObsData[in]�۲�����ʱ��
  @param  AllBroadcastEphData[in]��������
  @param  rs[out]����λ��
  @param  dts[out]�����Ӳ�
  @param  var[out]����λ���Ӳ��
  @param  svh[out]���ǽ���״̬
  @����ֵ ������������
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

        //�жϹ۲�����α�಻Ϊ��
        for (j = 0, pr = 0.0; j < NFREQ; j++)
          if ((pr = pSinSatObsData->pObsData[i].daPseRange[j]) != 0.0) break;

        if (j >= NFREQ) {
            continue;
        }
        //���������Ӽ��㷢��ʱ��
        time[i] = CTime::timeadd(pSinSatObsData->pObsData[i].tTimeStamp, -pr / CLIGHT);

        //���ݹ㲥�������������Ӳ�
        if (!EphClk(time[i], tiEphTime, pSinSatObsData->pObsData[i].ui8SatId, pBroadcastEphData, &dt)) {
            continue;
        }
        time[i] = CTime::timeadd(time[i], -dt);

        //�������Ƿ���ʱ�̵�λ�ú���
        if (!SatPos(time[i], tiEphTime, pSinSatObsData->pObsData[i].ui8SatId, EPHOPT_BRDC, 
            pBroadcastEphData, rs + i * 3, dts + i * 2, var + i, svh + i)) {
            continue;
        }
        //�������Ǹ���
        k++;
    }
    return k;
}