#ifndef MACRODEFINE_H
#define MACRODEFINE_H
#define STRFMT_RTCM3 1                  /* stream format: RTCM 3 */
#define DTTOL       0.025               /* tolerance of time difference (s) */
#define SQR(x)   ((x)* (x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define MAX(x,y)    ((x)>=(y)?(x):(y))
#define LOOPMAX     10000                                                /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define ROUND(x)    ((int)floor((x)+0.5))
#define ROUND_U(x)  ((uint32_t)floor((x)+0.5))
#define RE_GLO   6378136.0                                              /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14                                           /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14                                           /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14                                         /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14                                         /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3                                           /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5                                            /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5                                        /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5                                            /* earth angular velocity (rad/s) ref [9] */
#define OMGE     7.2921151467E-5                                        /* earth angular velocity (IS-GPS) (rad/s) */
#define RTOL_KEPLER 1E-13                                               /* relative tolerance for Kepler equation */
#define MAX_ITER_KEPLER 30                                              /* max number of iteration of Kelpler */

#define NFREQGLO    2                                                   /* number of carrier frequencies of GLONASS */

#define MINPRNGPS   1                                                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                                                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1)                             /* number of GPS satellites */
#define NSYSGPS     1

#define ENAGLO
#ifdef ENAGLO
#define MINPRNGLO   1                                                  /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                                                 /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1)                            /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif 
#define ENAGAL
#ifdef ENAGAL
#define MINPRNGAL   1                                                 /* min satellite PRN number of Galileo */
#define MAXPRNGAL   36                                                /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1)                            /* number of Galileo satellites */
#define NSYSGAL     1
#else
#define MINPRNGAL   0
#define MAXPRNGAL   0
#define NSATGAL     0
#define NSYSGAL     0
#endif
#ifdef ENAQZS
#define MINPRNQZS   193                                              /* min satellite PRN number of QZSS */
#define MAXPRNQZS   202                                              /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                                              /* min satellite PRN number of QZSS L1S */
#define MAXPRNQZS_S 191                                              /* max satellite PRN number of QZSS L1S */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1)                          /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif
#define ENACMP
#ifdef ENACMP
#define MINPRNCMP   1                                               /* min satellite sat number of BeiDou */
#define MAXPRNCMP   63                                              /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1)                         /* number of BeiDou satellites */
#define NSYSCMP     1
#else
#define MINPRNCMP   0
#define MAXPRNCMP   0
#define NSATCMP     0
#define NSYSCMP     0
#endif
#ifdef ENAIRN
#define MINPRNIRN   1                                              /* min satellite sat number of IRNSS */
#define MAXPRNIRN   14                                             /* max satellite sat number of IRNSS */
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1)                        /* number of IRNSS satellites */
#define NSYSIRN     1
#else
#define MINPRNIRN   0
#define MAXPRNIRN   0
#define NSATIRN     0
#define NSYSIRN     0
#endif
#ifdef ENALEO
#define MINPRNLEO   1                                             /* min satellite sat number of LEO */
#define MAXPRNLEO   10                                            /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1)                       /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif
#define NSYS        (NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP+NSYSIRN+NSYSLEO) /* number of systems */

#define MINPRNSBS   120                                          /* min satellite PRN number of SBAS */
#define MAXPRNSBS   158                                          /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1)                      /* number of SBAS satellites */

#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
#define FREQ1       1.57542E9                                   /* L1/E1/B1C  frequency (Hz) */
#define FREQ2       1.22760E9                                   /* L2         frequency (Hz) */
#define FREQ5       1.17645E9                                   /* L5/E5a/B2a frequency (Hz) */
#define FREQ6       1.27875E9                                   /* E6/L6  frequency (Hz) */
#define FREQ7       1.20714E9                                   /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9                                  /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9                                  /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9                                   /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6                                   /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9                                   /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6                                   /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9                                  /* GLONASS G3 frequency (Hz) */
#define FREQ1a_GLO  1.600995E9                                  /* GLONASS G1a frequency (Hz) */
#define FREQ2a_GLO  1.248060E9                                  /* GLONASS G2a frequency (Hz) */
#define FREQ1_CMP   1.561098E9                                  /* BDS B1I     frequency (Hz) */
#define FREQ2_CMP   1.20714E9                                   /* BDS B2I/B2b frequency (Hz) */
#define FREQ3_CMP   1.26852E9                                   /* BDS B3      frequency (Hz) */

#define NFREQ       3                                           /* number of carrier frequencies */
#define NEXOBS      2                                           /* number of extended obs codes */
#define MAXFREQ     7                                           /* max NFREQ */

#define EFACT_GPS   1.0                                         /* error factor: GPS */
#define EFACT_GLO   1.5                                         /* error factor: GLONASS */
#define EFACT_GAL   1.0                                         /* error factor: Galileo */
#define EFACT_QZS   1.0                                         /* error factor: QZSS */
#define EFACT_CMP   1.0                                         /* error factor: BeiDou */
#define EFACT_SBS   3.0                                         /* error factor: SBAS */

#define MAXSTA      255                                         /* max satellite number (1 to MAXSAT) */

#ifndef MAXOBS
#define MAXOBS      96                                         /* max number of obs in an epoch */
#endif
#define MAXRCV      64                                         /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                                         /* max number of obs type in RINEX */

#define SNR_UNIT    0.001                                     /* SNR unit (dBHz) */

#define SC2RAD      3.1415926535898                           /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0                            /* 1 AU (m) */
#define AS2R        (D2R/3600.0)                              /* arc sec to radian */

#define D2R         (PI/180.0)                                /* deg to rad */
#define R2D         (180.0/PI)                                /* rad to deg */
#define CLIGHT      299792458.0                               /* speed of light (m/s) */
#define SYS_NONE    0x00                                      /* navigation system: none */
#define SYS_GPS     0x01                                      /* navigation system: GPS */
#define SYS_SBS     0x02                                      /* navigation system: SBAS */
#define SYS_GLO     0x04                                      /* navigation system: GLONASS */
#define SYS_GAL     0x08                                      /* navigation system: Galileo */
#define SYS_QZS     0x10                                      /* navigation system: QZSS */
#define SYS_CMP     0x20                                      /* navigation system: BeiDou */
#define SYS_IRN     0x40                                      /* navigation system: IRNS */
#define SYS_LEO     0x80                                      /* navigation system: LEO */
#define SYS_ALL     0xFF                                      /* navigation system: all */

#define CODE_NONE   0                                         /* obs code: none or unknown */
#define CODE_L1C    1                                         /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                                         /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
#define CODE_L1W    3                                         /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                                         /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                                         /* obs code: L1M        (GPS) */
#define CODE_L1N    6                                         /* obs code: L1codeless,B1codeless (GPS,BDS) */
#define CODE_L1S    7                                         /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                                         /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                                         /* (not used) */
#define CODE_L1A    10                                        /* obs code: E1A,B1A    (GAL,BDS) */
#define CODE_L1B    11                                        /* obs code: E1B        (GAL) */
#define CODE_L1X    12                                        /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
#define CODE_L1Z    13                                        /* obs code: E1A+B+C,L1S (GAL,QZS) */
#define CODE_L2C    14                                        /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                                        /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                                        /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                                        /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                                        /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
#define CODE_L2P    19                                        /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                                        /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                                        /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                                        /* obs code: L2M        (GPS) */
#define CODE_L2N    23                                        /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                                        /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                                        /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                                        /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
#define CODE_L7I    27                                        /* obs code: E5bI,B2bI  (GAL,BDS) */
#define CODE_L7Q    28                                        /* obs code: E5bQ,B2bQ  (GAL,BDS) */
#define CODE_L7X    29                                        /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
#define CODE_L6A    30                                        /* obs code: E6A,B3A    (GAL,BDS) */
#define CODE_L6B    31                                        /* obs code: E6B        (GAL) */
#define CODE_L6C    32                                        /* obs code: E6C        (GAL) */
#define CODE_L6X    33                                        /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
#define CODE_L6Z    34                                        /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
#define CODE_L6S    35                                        /* obs code: L6S        (QZS) */
#define CODE_L6L    36                                        /* obs code: L6L        (QZS) */
#define CODE_L8I    37                                        /* obs code: E5abI      (GAL) */
#define CODE_L8Q    38                                        /* obs code: E5abQ      (GAL) */
#define CODE_L8X    39                                        /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
#define CODE_L2I    40                                        /* obs code: B1_2I      (BDS) */
#define CODE_L2Q    41                                        /* obs code: B1_2Q      (BDS) */
#define CODE_L6I    42                                        /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                                        /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                                        /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                                        /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                                        /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                                        /* obs code: B1I        (BDS) (obsolute) */
#define CODE_L1Q    48                                        /* obs code: B1Q        (BDS) (obsolute) */
#define CODE_L5A    49                                        /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                                        /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                                        /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                                        /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                                        /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                                        /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                                        /* obs code: SB+C       (IRN) */
#define CODE_L1D    56                                        /* obs code: B1D        (BDS) */
#define CODE_L5D    57                                        /* obs code: L5D(L5S),B2aD (QZS,BDS) */
#define CODE_L5P    58                                        /* obs code: L5P(L5S),B2aP (QZS,BDS) */
#define CODE_L5Z    59                                        /* obs code: L5D+P(L5S) (QZS) */
#define CODE_L6E    60                                        /* obs code: L6E        (QZS) */
#define CODE_L7D    61                                        /* obs code: B2bD       (BDS) */
#define CODE_L7P    62                                        /* obs code: B2bP       (BDS) */
#define CODE_L7Z    63                                        /* obs code: B2bD+P     (BDS) */
#define CODE_L8D    64                                        /* obs code: B2abD      (BDS) */
#define CODE_L8P    65                                        /* obs code: B2abP      (BDS) */
#define CODE_L4A    66                                        /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4B    67                                        /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4X    68                                        /* obs code: G1al1OCd+p (GLO) */
#define MAXCODE     68                                        /* max number of obs code */

#define PMODE_SINGLE 0                                        /* positioning mode: single */
#define PMODE_DGPS   1                                        /* positioning mode: DGPS/DGNSS */
#define PMODE_KINEMA 2                                        /* positioning mode: kinematic */
#define PMODE_STATIC 3                                        /* positioning mode: static */
#define PMODE_MOVEB  4                                        /* positioning mode: moving-base */
#define PMODE_FIXED  5                                        /* positioning mode: fixed */

#define SOLQ_NONE   0                                         /* solution status: no solution */
#define SOLQ_FIX    1                                         /* solution status: fix */
#define SOLQ_FLOAT  2                                         /* solution status: float */
#define SOLQ_SBAS   3                                         /* solution status: SBAS */
#define SOLQ_DGPS   4                                         /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 5                                         /* solution status: single */
#define SOLQ_PPP    6                                         /* solution status: PPP */
#define SOLQ_DR     7                                         /* solution status: dead reconing */
#define MAXSOLQ     7                                         /* max number of solution status */

#define MAXDTOE     7200.0                                    /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0                                    /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0                                   /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0                                   /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0                                    /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_IRN 7200.0                                    /* max time difference to IRNSS Toe (s) */

#define EPHOPT_BRDC 0                                         /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC 1                                         /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS 2                                         /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC 3                                       /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM 4                                       /* ephemeris option: broadcast + SSR_COM */

#define RE_WGS84    6378137.0                                 /* earth semimajor axis (WGS84) (m) */

#define IONOOPT_OFF 0                                         /* ionosphere option: correction off */
#define IONOOPT_BRDC 1                                        /* ionosphere option: broadcast model */
#define IONOOPT_SBAS 2                                        /* ionosphere option: SBAS model */
#define IONOOPT_IFLC 3                                        /* ionosphere option: L1/L2 iono-free LC */
#define IONOOPT_EST 4                                         /* ionosphere option: estimation */
#define IONOOPT_TEC 5                                         /* ionosphere option: IONEX TEC model */
#define IONOOPT_QZS 6                                         /* ionosphere option: QZSS broadcast model */
#define IONOOPT_STEC 8                                        /* ionosphere option: SLANT TEC model */

#define TROPOPT_OFF 0                                         /* troposphere option: correction off */
#define TROPOPT_SAAS 1                                        /* troposphere option: Saastamoinen model */
#define TROPOPT_SBAS 2                                        /* troposphere option: SBAS model */
#define TROPOPT_EST 3                                         /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG 4                                        /* troposphere option: ZTD+grad estimation */
#define TROPOPT_ZTD 5                                         /* troposphere option: ZTD correction */
#define MAXANT      64                  /* max length of station name/antenna type */
#endif // !MACRODEFINE_H

