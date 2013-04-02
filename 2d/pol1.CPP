#include "Poloz1.h"
#include "mf.h"
/*================strip_waveguide=====================*/
/*-----------------Planar------------------------*/
double
nref_max = 2.5,
nref_substrat = 2.0,
al_nref,
h_waveguide = 2,
w_waveguide = 2,

/* ========Koshiba==
      nref_max=3.44,
       nref_substrat=3.34,
       al_nref,
       h_waveguide=5,
       w_waveguide=5,
===========*/
lam = 1.55,
klam;
/*Test Koshiba iee proceeding-j. vol 139, No.2, April, 1992*/
// rib structure
double  W_test = 2,
        t_test = 0.2,
        h_test = 1.1,
        Xs_test = 3.0,
        Ys_test = 5.0,
        Yc_test = 1.0;
int     Mode = 0;
double  *RCyl, *FCyl;
double   HRCyl;
int NCyl;
int FlagCyl = 0;
double cl =  137.0388, //скорость света
       zMin, zMax, // пространство интегрирования в переменной z
       zh,   //  шаг интегрирования
       zhMap, chMap; /*шаг по изображению контуров*/
double *r, *U, *w, *w1; //массивы расстояния от ядра и потенциала
double *HE[MaxM],//матрица энергии
       *S[MaxM], // матрица перекрывания
       *D,       // массив для значений энергии
       *DLbuf, *Ebuf, *BufK; // вспомогательные массивы для диагонализации
unsigned  *ind;   /*индексы для разложения функций*/
int *Map[N_MAP_STR];
int DebugMain = 0, MfeDebug = 0;
int LL, //длина базиса
    nf2,
    Nh = 2, // число эрмитовых компонентов на узле
    iiz, // число точек в массивах пространства интегрирования
    NNz = 10, // число узлов для сплайнов
    Nf = 36, //число точек между узлами
    NcMap = 324, NzMap = 324; /*размер области для контура*/
long timt, tim0;
Datime  Dat;
/*===============масивы для работы с финитными функциями===========*/
double *f,     //массив для интегрирования матричных элементов
       *fi[2], //массив для сплайнов
       *d_fi[2], // массив для производных сплайнов
       *ff[2][2][3],// массив для произведений сплайнов
       *fdf[2][2][3],// массив для произведения сплайнов и производных
       *dfdf[2][2][3], // произвкдения производной
       *v, *f_v; // вектор коэффициентов разложения и функция
double zzN[2][2],
       zdzN[2][2];
/*================для двумерного случая==================*/
double *v_ground, *fz, *Uzc[2000], *buf_cz[2000],
       *zband;//массив для блокирования отражения от границ (по супергауссу)

unsigned long *ind3p;   /*индексы для разложения функций*/
double dzdzN[2][2],
       zcN[2][2],
       ccN[2][2], cdcN[2][2], dcdcN[2][2];

double *z, *c, cMin, cMax,
       M1, M2, me, //массы ядер, электрона
       MfeZ1 = 1, MfeZ2 = 1, // заряды ядер
       *Uz, *Uc, /*потенциал по z*/
       zav0,
       MfeDab0,
       ch;/*шаг по межъядерному расстоянию*/
double MfeDab = 2, // межъядерноe расстоянию
       E_Space = 30; // пространство для электрона

int    NNc = 20, // число узлов для сплайнов
       iic;// число точек в массивах пространства интегрирования
int MAX_v,//длина вспомогательных массивов
    iX, iY; // положение контурных графиков плотности
double delta_Z = 1; /*устранение сингулярности в кулоне*/

time_t  timcur(void)
{
    return (time_t)(time(&timt) - tim0);
}
void dump(double** Arr, int maxX, int maxY, char* filename) {
    //задампим HE
    FILE* Fwr = fopen(filename, "wt");
    int i, j;
    for (i = 0; i < maxY; i++) {
        if(i > 0) {
            fprintf(Fwr, "\n");
        }
        for (j = 0; j < maxX; j++) {
            if(j > 0) {
                fprintf(Fwr, ",");
            }
            fprintf(Fwr, "%f", Arr[j][i]);
        }
    }
    fclose(Fwr);
}
int Contur(int  *Map[], int ramxr , int ramyr, int NLevel,
           int XoX, int YoY)
{
    int err,
        ramxl, ramyl,
        MapMax, MapMin,
        LtGtLevgen,
        levgen0, levgen, HLevel ;
    TColor ColorForConture;
    int *buf3[3];
    int i, kk, jj, j, i0, i1, i2, ip, im,
        k00, k01, k02,
        k10, k11, k12,
        k20, k21, k22,
        k33, l,
        stcon = 0.0,
        sdef = 0.0;
    LtGtLevgen = 1;
    //XoX=200;YoY=300;
    err = 0;
    ramxl = 1;           /*Adding 13.10.91*/
    ramyl = 1;
    MapMax = MapMin = Map[1][1];
    for (j = 0; j < ramyr; j++)for (i = 0; i < ramxr; i++)
        {
            if (Map[j][i] > MapMax)MapMax = Map[j][i];
            if (Map[j][i] < MapMin)MapMin = Map[j][i];
        }
    setcolor(clBlack);
    setfillstyle(0, clWhite);
    rect(ramxl + XoX, ramyl + YoY, ramxr + XoX, ramyr + YoY);
    /* ...........  USER CODE  ............. */

    HLevel = (MapMax - MapMin) / NLevel;
    /*
    */
    levgen0 = HLevel;
    for (levgen = levgen0; levgen < MapMax - MapMin; levgen += HLevel)
    {
        if (levgen < 0)break;
        im = ip = 0;
        for (jj = 0; jj < 3; jj++)    buf3[jj] = Map[jj];
        jj = 2;
        i0 = 0; i1 = 1; i2 = 2; stcon = sdef = 0.0;
        for (l = 1; l < ramyr; l++)
        {

            if (l > ramyl && l < ramyr )
            {
                for (i = ramxl; i < ramxr - 1; i++)
                {
                    stcon = stcon + 1.0;
                    ip = i + 1; im = i - 1;
                    kk  = (int)buf3[i1][i] ;
                    if ((kk - levgen)*LtGtLevgen <= 0 )
                    {
                        k00 =  (int)buf3[i0][im];
                        k01 =  (int)buf3[i0][i] ;
                        k02 =  (int)buf3[i0][ip];
                        k10 =  (int)buf3[i1][im];
                        k12 =  (int)buf3[i1][ip];
                        k20 =  (int)buf3[i2][im];
                        k21 =  (int)buf3[i2][i] ;
                        k22 =  (int)buf3[i2][ip];
                        if ((k00 - levgen)*LtGtLevgen > 0) k00 = 0; else k00 = 1;
                        if ((k01 - levgen)*LtGtLevgen > 0) k01 = 0; else k01 = 1;
                        if ((k02 - levgen)*LtGtLevgen > 0) k02 = 0; else k02 = 1;
                        if ((k10 - levgen)*LtGtLevgen > 0) k10 = 0; else k10 = 1;
                        if ((k12 - levgen)*LtGtLevgen > 0) k12 = 0; else k12 = 1;
                        if ((k20 - levgen)*LtGtLevgen > 0) k20 = 0; else k20 = 1;
                        if ((k21 - levgen)*LtGtLevgen > 0) k21 = 0; else k21 = 1;
                        if ((k22 - levgen)*LtGtLevgen > 0) k22 = 0; else k22 = 1;

                        k33 = k00 * k01 * k02 * k10 * k12 * k20 * k21 * k22;
                        k11 = k00 + k01 + k02 + k10 + k12 + k20 + k21 + k22;
                        if (k11 == 1 || k11 == 2 || k11 == 3 )
                        {
                            buf3[i1][i] = /*(unsigned char)*/ (levgen + LtGtLevgen);
                            continue;
                        }
                        if ((k11 > 1 && k11 < 9 && k33 == 0) && (
                                    k00 + k01 + k10 == 3 ||
                                    k01 + k02 + k12 == 3 ||
                                    k10 + k20 + k21 == 3 ||
                                    k21 + k22 + k12 == 3
                                ))
                        {
                            if (kk >/*=-MapMin*/0) ColorForConture = clBlack; //clWhite;
                            else
                            {
                                return -1;
                                ColorForConture = clBlue;
                            }
                            putpixel(XoX + i, YoY + l, ColorForConture);
                            //                       putpixel(XoX+i,YoY-l,ColorForConture);
                        }
                    }
                }/* i */
            }/* if(l) */

            jj++; if (jj == 3) jj = 0;
            buf3[jj] = Map[l + 2];

            i0++; if (i0 == 3) i0 = 0;
            i1++; if (i1 == 3) i1 = 0;
            i2++; if (i2 == 3) i2 = 0;

        }/* l */
    }/*levgen*/

Mexit:
    return err;
} /* contur() */


/* аналитические выражения для финитных функций начало */
double Fb(double x)// эрмитова компонента 0
{
    double a;
    if (fabs(x) > 1.0 ) return 0.0;
    if (x >= 0.0)
    {
        a = (1.0 - x);
        return a * a * (1.0 + 2 * x);
    }
    if (x < 0.0)
    {
        a = (1.0 + x);
        return a * a * (1.0 - 2 * x);
    }
    return 0.0;
}/*F(x)*/
double FbH(double x) // эрмитова компонента 1
{
    static double a;
    if (fabs(x) > 1.0 ) return 0.0;
    if (x >= 0.0)
    {
        a = (1.0 - x);
        return a * a * x;
    }
    if (x < 0.0)
    {
        a = (1.0 + x);
        return a * a * x;
    }
    return 0.0;
}/*FH(x)*/

double DFb(double x)//  производная эрмитовой компонента 0
{
    static double a;
    if (fabs(x) > 1.0 ) return 0;
    if (x >= 0.0)
    {
        a = (1.0 - x);
        return -2 * a * (1 + 2 * x) + a * a * 2;
    }
    if (x < 0.0)
    {
        a = (1.0 + x);
        return  2 * a * (1 - 2 * x) - a * a * 2;
    }
    return 0.0;
}/*DF(x)*/

double DFbH(double x)  //  производная эрмитовой компонента 1
{
    static double a;
    if (fabs(x) > 1.0 ) return 0.0;
    if (x >= 0.0)
    {
        a = (1.0 - x);
        return (-2 * a * x + a * a) ;
    }
    if (x < 0.0)
    {
        a = (1.0 + x);
        return  (2 * a * x + a * a) ;
    }
    return 0.0;
}/*DFH(x)*/
//====================================================
Boolean Region(int Max, int r, int r1, int *ir0, int *ir1)
{
    int pr, pr1, k, Maxr;
    /* r,r1 - номера узлов , где центрируются базисные функции*/
    /* ir0,ir1 - полученая   область интегрирования           */
    /* если Maxr < 0  вычисляется область для переменной  ro */
    Maxr = abs(Max);
    k  = abs(r - r1);
    if (k > 1)
     return False;
    pr = r * Nf;
    pr1 = (r + 2) * Nf;
    if (k == 0) {
        *ir0 = pr;
        *ir1 = pr1;
    }/*k=0*/
    if (r < r1)  {
        *ir0 = pr + Nf;
        *ir1 = pr1;
    }/*r<r1*/
    if (r > r1)
    {
        *ir0 = pr;
        *ir1 = pr1 - Nf;
    }/*r>r1*/
    return True;
}/*Region*/
double Simps(int nn, double hh, double f[])
{
    static double a, b, h;
    static int i, j, L, N, n;
    h = hh; n = nn;
    L = 0; N = n;
    a = b = 0.0;
    /*Ньютона-Котеса n=6
      (Дьяконов В.П. Справочник по алгоритмам и программам
      на языке бейсик для персональных ЭВМ. М:"Наука",1987, стр.106)*/
    b = 1.0 / 140;

    for (j = 0; j < nn; j += 6)
    {
        a += (41 * (f[j] + f[j + 6]) + 216 * (f[j + 1] + f[j + 5]) + 27 * (f[j + 2] + f[j + 4]) + 272 * f[j + 3]);
    }/*j*/
    a = a * h * b;

    return a;
}/*Simps*/
//======================================

//======================================
int  TREQLAB(int N, double *A[], double *B[], double *D,
             double *DL, double *E)
{
    /*================C-like Arrays[0<=i<MM]              ==================   *
    *       GENERAL MATRIX PROBLEM AX=DBX                                      *
    *       COMPLETELY DOUBLE PRECISION                                        *
    *                                                                          *
    *       n - THE MATRIX ORDER                                               *
    *       NMAX - IF THE PROGRAM IS TO WORK WITH DIFFERENT N, NMAX            *
    *               SHOULD BE THE LARGEST N. IT SHOULD COINSIDE WITH           *
    *               DIMENSION OR COMMON DECLARATION                            *
    *       A[NMAX,NMAX], B[NMAX,NMAX], D[NMAX], DL[NMAX], E[NMAX]             *
    *               DOUBLE PRECISION ARRAYS                                    *
    *               AT THE BEGINNING                                           *
    *                  A AND B CONTAIN INITIAL MATRIX. ONLY ELEMENTS OF UPPER  *
    *                       RIGHT TRIANGLE ARE USED                            *
    *               IN THE END                                                 *
    *                  D[NMAX] CONTAINS EIGENVALUES                            *
    *                  A[NMAX,NMAX] CONTAINS EIGENVECTORS ARRANGED IN COLUMNS  *
    *                    SO A[I,K] IS I-TH ELEMENT OF K-TH VECTOR              *
    *       return - ERROR INDEX : return  0  - SUCCESS                        *
    *                              return -1  - MATRIX B IS NOT POSITIVE       *
    *                              return  1  - MORE THAN 30 ITERATIONS        *
    *=====================================================================*/
    static double ETA, MACHEP, TOL;
    static double BB, C, F, G, H, HH, P, R, S, T, TSTA, TSTB;
    static int I, J, K, L, M, NERR, ST = 0, j;
    /*     DebugMain = 0;*/
    /*       PARAMETERS ETA AND MACHEPS FOR UKNC DOUBLE PRECISION */
    ETA = pow(1.00 / 2.00 , 1023);
    MACHEP = pow(1.00 / 2.00, 53);

    TOL = ETA / MACHEP;
    ST = 0;
    if (N == 1)
    {
        D[1 - 1] = A[1 - 1][1 - 1] / B[1 - 1][1 - 1];
        A[1 - 1][1 - 1] = 1.0;
        return  0;
    }/*endif*/

    /* --------------------      reduc1   -------------------------*/

    for (I = 1; I <= N; I++)
    {
        if ((I / 10) * 10 == I)
        {
            sprintf (str, " (TREQLAB) reduc1= %d %5ld s.", I, timcur());
            FastWriteA (str, 1, 50, clBlack, clWhite);
        }
        //     if(kbhit()){if ((ch = getch()) == ESC) return -777;}
        //        if(StopMain){ return -777;}

        for (J = 1; J <= N; J++)
        {
            G = B[I - 1][J - 1];
            if (I == 1) goto met102;
            for (L = 1; L <= I - 1; L++)
            {
                K = I - L;
                G = G - B[I - 1][K - 1] * B[J - 1][K - 1];
            }/*L*/
met102:     if (I != J) goto met106;
            if (G >  0.00 ) goto met104;
            sprintf (str, " (TREQLAB) reduc= B[%d,%d] %5ld s.[G<0]", I - 1, J - 1, timcur());
            return -1;
met104:     H = sqrt(G);
            DL[I - 1] = H;
            continue;
met106:     B[J - 1][I - 1] = G / H;
        }/*J 108*/
    }/*I 110*/


    for (I = 1; I <= N; I++)
    {
        if ((I / 10) * 10 == I)
        {
            sprintf (str, " (TREQLAB) reduc2= %d %5ld s. ", I, timcur());
            FastWriteA (str, 1, 50, clBlack, clWhite);
        }

        //     if(kbhit()){if ((ch = getch()) == ESC) return -777;}
        //        if(StopMain){ return -777;}

        H = DL[I - 1];
        for (J = 1; J <= N; J++)
        {
            G = A[I - 1][J - 1];
            if (I == 1) goto met122;
            for (L = 1; L <= I - 1; L++)
            {
                K = I - L;
                G = G - B[I - 1][K - 1] * A[J - 1][K - 1];
            }/*L*/
met122:      A[J - 1][I - 1] = G / H;
        }/*J*/
    }/* I*/
    //   if(ST){getch();
    //           TSTA=TSTB=0.0;for(j=0;j<N;j++){TSTA+=A[10][j];TSTB+=B[10][j];}
    //         }

    for (J = 1; J <= N; J++)
    {
        if ((J / 10) * 10 == J)
        {
            sprintf (str, " (TREQLAB) reduc3= %d %5ld s. ", J, timcur());
            FastWriteA (str, 1, 50, clBlack, clWhite);
        }
        //     if(kbhit()){if ((ch = getch()) == ESC) return -777;}
        //        if(StopMain){ return -777;}

        for (I = J; I <= N; I++)
        {
            G = A[I - 1][J - 1];
            if (I <= J) goto met132;
            for (L = J; L <= I - 1; L++)
            {
                K = I + J - L - 1;
                G = G - A[K - 1][J - 1] * B[I - 1][K - 1];
            }/*L*/
met132:        if (J == 1) goto met136;
            for (L = 1; L <= J - 1; L++)
            {
                K = J - L;
                G = G - A[J - 1][K - 1] * B[I - 1][K - 1];
            }/*L*/
met136:        A[I - 1][J - 1] = G / DL[I - 1];
        }/*I*/
    }/*J*/

    /*-----------------       tred2 ----------------*/
    //         if(ST){getch();
    //           TSTA=TSTB=0.0;for(j=0;j<N;j++){TSTA+=A[10][j];TSTB+=B[10][j];}
    //         }

    I = N + 1;
met200:     I = I - 1;
    if ((I / 10) * 10 == I)
    {
        sprintf (str, " (TREQLAB) tred2.0=%d %5ld s. ", I, timcur());
        FastWriteA (str, 1, 50, clBlack, clWhite);
    }
    //     if(kbhit()){if ((ch = getch()) == ESC) return -777;}
    //        if(StopMain){ return -777;}

    if (I <  2)goto met600;
    L = I - 2;
    F = A[I - 1][I - 1 - 1];
    G = 0.0;
    K = 0;
met205:     K = K + 1;
    if (K >  L) goto met210;
    P = A[I - 1][K - 1];
    G = G + P * P;
    goto met205;
met210:     H = G + F * F;
    if (G >  TOL) goto met220;
    E[I - 1] = F;
    H = 0.00;
    goto met500;
met220:     L = L + 1;
    P = sqrt(H);
    G = P;
    if (F >= ETA) G = -P;
    E[I - 1] = G;
    H = H - F * G;
    A[I - 1][I - 1 - 1] = F - G;
    F = 0.00;


    for (J = 1; J <= L; J++)
    {
        A[J - 1][I - 1] = A[I - 1][J - 1] / H;
        G = 0.00;
        for (K = 1; K <= J; K++)
        {
            G = G + A[J - 1][K - 1] * A[I - 1][K - 1];
        }/*K*/
        K = J + 1;
met320:   if (K >  L) goto met330;
        G = G + A[K - 1][J - 1] * A[I - 1][K - 1];
        K = K + 1;
        goto met320;
met330:   E[J - 1] = G / H;
        F = F + G * A[J - 1][I - 1];
    }/*J*/


    HH = F / (H + H);
    for (J = 1; J <= L; J++)
    {
        F = A[I - 1][J - 1];
        G = E[J - 1] - HH * F;
        E[J - 1] = G;
        for (K = 1; K <= J; K++)
        {
            P = F * E[K - 1] + G * A[I - 1][K - 1];
            A[J - 1][K - 1] = A[J - 1][K - 1] - P;
        }/*K*/
    }/*J*/
met500: D[I - 1] = H;
    goto met200;

met600: D[1 - 1] = 0.00;
    E[1 - 1] = 0.00;
    //         if(ST){getch();
    //           TSTA=TSTB=0.0;for(j=0;j<N;j++){TSTA+=A[10][j];TSTB+=B[10][j];}
    //         }
    for (I = 1; I <= N; I++)
    {
        if ((I / 10) * 10 == I)
        {
            sprintf (str, " (TREQLAB) tred2.1=%d %5ld s. ", I, timcur());
            FastWriteA (str, 1, 50, clBlack, clWhite);
        }
        //     if(kbhit()){if ((ch = getch()) == ESC) return -777;}
        //        if(StopMain){ return -777;}
        if (I == 1) goto met790;
        L = I - 1;
        if (fabs(D[I - 1]) <  ETA) goto met790;
        for (J = 1; J <= L; J++)
        {
            G = 0.00;
            for (K = 1; K <= L; K++)
            {
                G = G + A[I - 1][K - 1] * A[K - 1][J - 1];
            }/*K*/
            for (K = 1; K <= L; K++)
            {
                A[K - 1][J - 1] = A[K - 1][J - 1] - G * A[K - 1][I - 1];
            }/*K*/
        }/*J*/
met790:   D[I - 1] = A[I - 1][I - 1];
        A[I - 1][I - 1] = 1.00;
        J = 1;
met800:   if (J >  L) continue;;
        A[I - 1][J - 1] = 0.00;
        A[J - 1][I - 1] = 0.00;
        J = J + 1;
        goto met800;
    }/*I*/

    /*---------------------       tql2  ---------------------*/

    NERR = 0;
    for (I = 2; I <= N; I++)
    {
        E[I - 1 - 1] = E[I - 1];
    }/*I*/
    E[N - 1] = 0.00;
    BB = 0.00;
    F = 0.00;

    //         if(ST){getch();
    //           TSTA=TSTB=0.0;for(j=0;j<N;j++){TSTA+=A[10][j];TSTB+=B[10][j];}
    //         }

    for (L = 1; L <= N; L++)
    {
        J = 0;
        H = MACHEP * (fabs(D[L - 1]) + fabs(E[L - 1]));
        if (BB <  H) BB = H;

        for (M = L; M <= N; M++)
        {
            if (fabs(E[M - 1]) <= BB) break;
        }/*M*/
        /*--------------         CONT - 1110 ---------------*/

met1110:   if (M == L) goto met1300;
        /*----------------         NEXTIT - 1200-------------*/
met1200:  if (J <  30) goto met1210;
        return 1;
met1210:  J = J + 1;

        G = D[L - 1];
        P = (D[L + 1 - 1] - G) / (2.0 * E[L - 1]);
        R = sqrt(P * P + 1.00 );
        T = P + R;
        if (P <  0.0) T = P - R;
        D[L - 1] = E[L - 1] / T;
        H = G - D[L - 1];
        I = L + 1;
met1220:  if (I >  N) goto met1230;
        D[I - 1] = D[I - 1] - H;
        I = I + 1;
        goto met1220;
met1230:  F = F + H;
        P = D[M - 1];
        C = 1.00;
        S = 0.00;
        I = M - 1;
met1240:  if (I <  L) goto met1280;
        T = E[I - 1];
        G = C * T;
        H = C * P;
        if (fabs(P) <  fabs(T)) goto met1250;
        C = T / P;
        R = sqrt(C * C + 1.00 );
        E[I + 1 - 1] = S * P * R;
        S = C / R;
        C = 1.00 / R;
        goto met1260;
met1250:  C = P / T;
        R = sqrt(C * C + 1.00 );
        E[I + 1 - 1] = S * T * R;
        S = 1.00 / R;
        C = C / R;
met1260:  P = C * D[I - 1] - S * G;
        D[I + 1 - 1] = H + S * (C * G + S * D[I - 1]);
        for (K = 1; K <= N; K++)
        {
            T = A[K - 1][I - 1];
            H = A[K - 1][I + 1 - 1];
            A[K - 1][I + 1 - 1] = S * T + C * H;
            A[K - 1][I - 1] = C * T - S * H;
        }/*K*/
        I = I - 1;
        goto met1240;
met1280:
        T = S * P;
        E[L - 1] = T;
        D[L - 1] = C * P;
        if (fabs(T) >  BB) goto met1200;
        /*------------------------         ROOT - 1300----------------*/
met1300:
        D[L - 1] = D[L - 1] + F;
        if ((L / 10) * 10 == L)
        {
            sprintf (str, " (TREQLAB) D[%3d]=%e [%5ld s]. ", L, D[L - 1], timcur());
            FastWriteA (str, 1, 50, clBlack, clWhite);
        }
        //     if(kbhit()){if ((ch = getch()) == ESC) return -777;}
        //        if(StopMain){ return -777;}
    }/*L*/

    for (I = 1; I <= N; I++)
    {
        K = I;
        P = D[I - 1];
        J = I + 1;
met1510:  if (J >  N) goto met1530;
        if (D[J - 1] >= P) goto met1520;
        K = J;
        P = D[J - 1];
met1520:  J = J + 1;
        goto met1510;
met1530:  if (K == I) continue;
        D[K - 1] = D[I - 1];
        D[I - 1] = P;
        for (J = 1; J <= N; J++)
        {
            P = A[J - 1][I - 1];
            A[J - 1][I - 1] = A[J - 1][K - 1];
            A[J - 1][K - 1] = P;
        }/*J*/
    }/*I*/

    /*--------------------       rebaka ------------------*/

    for (J = 1; J <= N; J++)
    {
        if ((J / 10) * 10 == J)
        {
            sprintf (str, " (TREQLAB) rebak =%d [%5ld s]. ", J, timcur());
            FastWriteA (str, 1, 50, clBlack, clWhite);
        }
        for (L = 1; L <= N; L++)
        {
            I = N - L + 1;
            G = A[I - 1][J - 1];
            if (I == N) goto met1610;
            for (K = I + 1; K <= N; K++)
            {
                G = G - B[K - 1][I - 1] * A[K - 1][J - 1];
            }/*K*/
met1610:     A[I - 1][J - 1] = G / DL[I - 1];

        }/*L*/
    }/*J*/
    return 0;

}/*TREQLAB*/
//=======================================
//===============================================
/*===================h2p1d();===========================*/
int DosGetDateTime(Datime *DT)
{
    struct time t;
    struct date d;
    Datime Date;
    gettime(&t);
    getdate(&d);
    Date.day = d.da_day;
    Date.month = d.da_mon;
    Date.year = d.da_year;
    Date.hours = t.ti_hour;
    Date.minutes = t.ti_min;
    Date.seconds = t.ti_sec;
    *DT = Date;
}/*DosGetDateTime*/
//================================================
double Nref_Planar(double y1)
{
    double y;
    y = y1;   if (Cyl == 1)y = -y1;
    if (y > 0)return 1;
    else
    {
        //экспонента
        if (Mode == 0) return  nref_substrat + (nref_max - nref_substrat) * exp(al_nref * y);
        // ступенька
        if (Mode == 1)
        {
            if (y <= 0 && fabs(y) < h_waveguide)return  nref_max;
            else  return  nref_substrat;
        }
        // типа функции Гаусса
        if (Mode == 2) return  nref_substrat + (nref_max - nref_substrat) * exp(-al_nref * y * y);
        // парабола
        if (Mode == 3)
        {
            if (y <= 0 && fabs(y) < h_waveguide)
                return  y * y * ((nref_substrat - nref_max) / sqr(h_waveguide)) + nref_max;
            else     return  nref_substrat;
        }
    }
}/*Nref*/



double FunPot_cz(double c, double z)
{
    double y, a;
    y = c;
    //   a= nref_substrat+(nref_max-nref_substrat)*exp(-sqr(y/w_waveguide));
    a = (nref_max - nref_substrat) * exp(-sqr(y / w_waveguide));
    if (z > 0) return  1;
    else  return  (Nref_Planar(z) - nref_substrat) * a  +  nref_substrat;
}/*FunPot_cz*/


/*
double FunPot_cz(double c,double z) // test from Koshiba
{
    double nc, ns, nf, a;
    double t, h, W;
    int var;
    // Mode=1;
    nc = 1;
    ns = nref_substrat;
    nf = nref_max;
    t = t_test;
    h = h_test;
    W = W_test;
    a = nc;

    if (z <= -t)a = ns;
    if (z <= 0 && z > -t)a = nf;
    if ((z > 0 && z < h) && (c > -W / 2 && c < W / 2)) a = nf;
    return a;

}//FunPot_cz test from Koshiba

*/
//----------------------------
double MatElS_1d_3p(int ig, int jg)
{
    /*интеграл перекрывания*/
    static int i, j, z0, z1, ix, iy,
           c0, c1 ;
    static double a, rez, za, ca;
    static int iz, k, iz1, hz, hz1, d, dl,
           ic,   ic1, hc, hc1, d1, dl1;
    UUU3p u;
    rez = 0.0;
    /*-----------*/
    u.i = ind3p[ig]; iz  = u.zrc.iz; hz  = u.zrc.hz;
    ic  = u.zrc.ic; hc  = u.zrc.hc;
    u.i = ind3p[jg]; iz1 = u.zrc.iz; hz1 = u.zrc.hz;
    ic1 = u.zrc.ic; hc1 = u.zrc.hc;
    rez = a = 0.0;
    d = iz1 - iz + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d - 1) > 1) return 0.0;
    dl = 0; if (d == 2 ) dl = Nf;

    d1 = ic1 - ic + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d1 - 1) > 1) return 0.0;
    dl1 = 0; if (d1 == 2 ) dl1 = Nf;

    if ( Region(-(NNz - 1), iz, iz1, &z0, &z1) && Region(-(NNc - 1), ic, ic1, &c0, &c1))
    {
        /*
             sprintf(str,"S[%d][%d]",ig,jg);
            FastWriteA (str, 5,5 , clWhite, clRed);
        */
        for (i = z0; i <= z1; i++)
        {
            f[i - z0] = ff[hz][hz1][d][i - z0 + dl];
        }
        za = Simps(z1 - z0, zh, f);

        for (i = c0; i <= c1; i++)
        {
            f[i - c0] = ff[hc][hc1][d1][i - c0 + dl1];
        }
        ca = Simps(c1 - c0, ch, f);

        a = zzN[hz][hz1] * ccN[hc][hc1];
        rez = za * ca * a; goto Mexit;
    }/*Region*/
    /* -------------*/
Mexit:
    if (DebugMain == 1)
    {
        sprintf(str,
                "S {iz,iz1[%d-%d]z0,z1[%d-%d]}{(ic,ic1[%d,%d],c0,c1[%d,%d]}za=%e,ca=%e,a=%e ch=%f ",
                iz, iz1, z0, z1,
                ic, ic1, c0, c1, za, ca, a, ch);
        FastWriteA (str, 5, 5 , clWhite, clRed);
        sprintf(str, "{hz=%d,hz1=%d,d=%d,dl=%d} {hc=%d,hc1=%d,d1=%d,dl1=%d}rez=%e",
                hz, hz1, d, dl,
                hc,
                hc1,
                d1,
                dl1, rez);
        FastWriteA (str, 5, 5 , clWhite, clRed);
    }
    return rez;
}/*MatElS_1d_3p*/

double MatElV_1d_3p(int ig, int jg)
{
    /*матричные элементы потенциала*/
    static int i, j, z0, z1, ix, iy,
           c0, c1 ;
    static double a, bz, bc, rez, za, ca, cza;
    static int iz, k, iz1, hz, hz1, d, dl,
           ic,   ic1, hc, hc1, d1, dl1;
    UUU3p u;
    rez = 0.0;
    /*-----------*/
    u.i = ind3p[ig]; iz  = u.zrc.iz; hz  = u.zrc.hz;
    ic  = u.zrc.ic; hc  = u.zrc.hc;
    u.i = ind3p[jg]; iz1 = u.zrc.iz; hz1 = u.zrc.hz;
    ic1 = u.zrc.ic; hc1 = u.zrc.hc;
    rez = a = 0.0;
    d = iz1 - iz + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d - 1) > 1) return 0.0;
    dl = 0; if (d == 2 ) dl = Nf;

    d1 = ic1 - ic + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d1 - 1) > 1) return 0.0;
    dl1 = 0; if (d1 == 2 ) dl1 = Nf;

    if ( Region(-(NNz - 1), iz, iz1, &z0, &z1) && Region(-(NNc - 1), ic, ic1, &c0, &c1))
    {
        sprintf(str, "MatElV_1d_3p[%d][%d]", ig, jg);
        FastWriteA (str, 5, 5 , clWhite, clRed);
        /*
               for(i=z0;i<=z1;i++){
               sprintf(str,"i=%d ,z0=%d,z1=%d,",i ,z0,z1);
               f[i-z0] =ff[hz][hz1][d][i-z0+dl]*Uz[i];
                     //    ix=i*msx*zh+XoX;
                     //    iy=bottom-(f[i-z0]*100+YoY);
                     //    putpixel(ix,iy,clLime);
               }
               za=Simps(z1-z0,zh,f);
               for(i=z0;i<=z1;i++){
                f[i-z0] =ff[hz][hz1][d][i-z0+dl];
               }
               bz=Simps(z1-z0,zh,f);

        for (i = c0; i <= c1; i++)
        {
            f[i - c0] = ff[hc][hc1][d1][i - c0 + dl1] * Uc[i];
        }
        ca = Simps(c1 - c0, ch, f);
        for (i = c0; i <= c1; i++)
        {
            f[i - c0] = ff[hc][hc1][d1][i - c0 + dl1];
        }
        bc = Simps(c1 - c0, ch, f);

        */
        for (i = z0; i <= z1; i++)
        {
            for (j = c0; j <= c1; j++)
            {
                f[j - c0] = ff[hz][hz1][d][i - z0 + dl] *
                            ff[hc][hc1][d1][j - c0 + dl1] * Uzc[i][j];
            }
            fz[i - z0] = Simps(c1 - c0, ch, f);
        }/*i*/
        cza = Simps(z1 - z0, zh, fz);
        /*-------------------*/
        a = zzN[hz][hz1] * ccN[hc][hc1];
        //rez= (za*bc +ca*bz +  cza)*a; goto Mexit;
        rez =  cza * a; goto Mexit;
    }/*Region*/
    /* -------------*/
Mexit:
    if (DebugMain == 1)
    {
        sprintf(str, "V iz[%d-%d]z[%d-%d](%d[%d,%d],%d[%d,%d])=%e ",
                iz, iz1, z0, z1, ig, iz, hz, jg, iz1, hz1, rez);
        FastWriteA (str, 5, 5 , clWhite, clRed);
    }
    return rez;

}/*MatElV_1d_3p*/

double MatElT_1d_3p(int ig, int jg)
{
    /*матричные элементы кинетической энергии*/
    static int i, j, z0, z1, ix, iy,
           c0, c1 ;
    static double a, ba, rez, za, ca, bz, bc, dza, dca;
    static int iz, k, iz1, hz, hz1, d, dl,
           ic,   ic1, hc, hc1, d1, dl1;
    UUU3p u;
    rez = 0.0;
    /*-----------*/
    u.i = ind3p[ig];   iz  = u.zrc.iz; hz  = u.zrc.hz;
    ic  = u.zrc.ic; hc  = u.zrc.hc;
    u.i = ind3p[jg]; iz1 = u.zrc.iz; hz1 = u.zrc.hz;
    ic1 = u.zrc.ic; hc1 = u.zrc.hc;
    rez = a = 0.0;
    d = iz1 - iz + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d - 1) > 1) return 0.0;
    dl = 0; if (d == 2 ) dl = Nf;

    d1 = ic1 - ic + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d1 - 1) > 1) return 0.0;
    dl1 = 0; if (d1 == 2 ) dl1 = Nf;

    if ( Region(-(NNz - 1), iz, iz1, &z0, &z1) && Region(-(NNc - 1), ic, ic1, &c0, &c1))
    {

        for (i = z0; i <= z1; i++)
        {
            f[i - z0] = ff[hz][hz1][d][i - z0 + dl];
        }
        za = Simps(z1 - z0, zh, f);

        for (i = z0; i <= z1; i++)
        {
            f[i - z0] = -dfdf[hz][hz1][d][i - z0 + dl];
        }
        bz = Simps(z1 - z0, zh, f);

        /*------------------------*/
        for (i = c0; i <= c1; i++)
        {
            f[i - c0] = ff[hc][hc1][d1][i - c0 + dl1];
        }
        ca = Simps(c1 - c0, ch, f);

        for (i = c0; i <= c1; i++)
        {
            f[i - c0] = -dfdf[hc][hc1][d1][i - c0 + dl1];
        }
        bc = Simps(c1 - c0, ch, f);

        /*--------------------------*/
        for (i = z0; i <= z1; i++)
        {
            f[i - z0] = fdf[hz][hz1][d][i - z0 + dl];
        }
        dza = Simps(z1 - z0, zh, f);

        for (i = c0; i <= c1; i++)
        {
            f[i - c0] = fdf[hc][hc1][d1][i - c0 + dl1];
        }
        dca = Simps(c1 - c0, ch, f);


        /*--------------------------*/
        a = zzN[hz][hz1] * dcdcN[hc][hc1];
        rez = (1 / M1 + 1 / M2) * za * bc * a;

        a = dzdzN[hz][hz1] * ccN[hc][hc1];
        rez += (1 / M2 + 1 / me) * bz * ca * a;

        //  a=zdzN[hz][hz1]*cdcN[hc][hc1];
        // rez+=-(2/M2)*dza*dca*a;

        goto Mexit;
    }/*Region*/
    /* -------------*/
Mexit:
    if (DebugMain == 1)
    {
        sprintf(str, "T iz[%d-%d]z[%d-%d](%d[%d,%d],%d[%d,%d])=%e ",
                iz, iz1, z0, z1, ig, iz, hz,
                jg, iz1, hz1, rez);
        FastWriteA (str, 5, 5 , clWhite, clRed);
    }
    return 0.5 * rez;
}/*MatElT_1d_3p_3p*/
double MatElX_1d_3p(int ig, int jg, char key)
{
    /*ch=(z,c,r)- матричные элементы (z,c+MfeDab,r) */

    static int i, j, z0, z1, ix, iy,
           c0, c1 ;
    static double a, aa, rez, za, ca;
    static int iz, k, iz1, hz, hz1, d, dl,
           ic,   ic1, hc, hc1, d1, dl1;
    UUU3p u;
    rez = 0.0;
    /*-----------*/
    u.i = ind3p[ig]; iz  = u.zrc.iz; hz  = u.zrc.hz;
    ic  = u.zrc.ic; hc  = u.zrc.hc;
    u.i = ind3p[jg]; iz1 = u.zrc.iz; hz1 = u.zrc.hz;
    ic1 = u.zrc.ic; hc1 = u.zrc.hc;
    rez = a = 0.0;
    d = iz1 - iz + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d - 1) > 1) return 0.0;
    dl = 0; if (d == 2 ) dl = Nf;

    d1 = ic1 - ic + 1; /*сочетание узлов -1(d==0),0(d==1),+1(d==2)*/
    if (abs(d1 - 1) > 1) return 0.0;
    dl1 = 0; if (d1 == 2 ) dl1 = Nf;

    if ( Region(-(NNz - 1), iz, iz1, &z0, &z1) && Region(-(NNc - 1), ic, ic1, &c0, &c1))
    {
        //     sprintf(str,"[%d][%d]",ig,jg);
        //   FastWriteA (str, 5,5 , clWhite, clRed);
        for (i = z0; i <= z1; i++)
        {
            aa = 1; if (key == 'z')aa = z[i];
            f[i - z0] = ff[hz][hz1][d][i - z0 + dl] * aa;
        }
        za = Simps(z1 - z0, zh, f);
        for (i = c0; i <= c1; i++)
        {
            aa = 1; if (key == 'c')aa = c[i];
            f[i - c0] = ff[hc][hc1][d1][i - c0 + dl1] * aa;
        }
        ca = Simps(c1 - c0, ch, f);
        a = zzN[hz][hz1] * ccN[hc][hc1];
        rez = za * ca * a; goto Mexit;
    }/*Region*/
    /* -------------*/
Mexit:
    if (DebugMain == 1)
    {
        sprintf(str,
                "X {iz,iz1[%d-%d]z0,z1[%d-%d]}{(ic,ic1[%d,%d],c0,c1[%d,%d]}za=%e,ca=%e,a=%e ch=%f ",
                iz, iz1, z0, z1,
                ic, ic1, c0, c1, za, ca, a, ch);
        FastWriteA (str, 5, 5 , clWhite, clRed);
        sprintf(str, "{hz=%d,hz1=%d,d=%d,dl=%d} {hc=%d,hc1=%d,d1=%d,dl1=%d}rez=%e",
                hz, hz1, d, dl,
                hc,
                hc1,
                d1,
                dl1, rez);
        FastWriteA (str, 5, 5 , clWhite, clRed);
    }
    return rez;
}/*MatElX_1d_3p*/
//функция моды от номера узла
double Mfe_fun_1d_3p(int iz, int  ic, double *v)
{
    /*значение функции в номерах точек интегрирования   iz,ic*/
    static  double a, b1, b2, b3, b4, tau, aa;
    static  int i, j, m, n, zh, ch, l, k;
    /*зачение рел. компон. функции как сумма от ближайших соседей*/
    a = 0.0;
    l = -77;

    i = iz / Nf - 1;          j = ic / Nf - 1; /*номера узлов*/

    if (iz < 0    || ic < 0   ) goto Mexit;
    if (iz >= iiz || ic >= iic) goto Mexit;
    if (i >= NNz  || j >= NNc ) goto Mexit;
    /* zh,ch - эрмитовы компоненты ; i + m, j + n  - , ближайшие узлы */

    for (zh = 0; zh <= 1; zh++)
        for (ch = 0; ch <= 1; ch++)
            for (m = 0; m <= 1; m++)
                for (n = 0; n <= 1; n++)
                {
                    aa = 1;
                    if (i + m >= NNz || j + n >= NNc ) continue;
                    if (i + m < 0 || j + n < 0 ) continue;
                    //       l= rtkh((m+i),(n+j),k,rh,th)  ;
                    l = ch+(n + j)*2+zh*NNc*2+(m + i)*NNc*4;
                    //l = zhch((m + i), zh, (n + j), ch);
                    if (l >= LL || l < 0) goto Mexit;
                    b1 =   fi[zh][iz - (i + m) * Nf];
                    b2 =   fi[ch][ic - (j + n) * Nf];
                    b3 =   v[l];
                    b4 =   zcN[zh][ch];
                    a += b1 * b2 * b3 * b4 * aa;
                }/*rh,kh,m,n*/
Mexit:
    return a;
}/*Mfe_fun_1d_3p*/

//функция моды от координат
double Mfe_fun_1d_3p_zc(double z, double c, double *v) {
    /*значение функции координатах z,c */
    if (z <= zMin || z >= zMax || c <= cMin || c >= cMax) return 0.0;
    int ic, iz;
    iz = (z - zMin) / zh; if (iz >= iiz)iz = iiz - 1;
    ic = (c - cMin) / ch; if (ic >= iic)ic = iic - 1;
    return Mfe_fun_1d_3p(iz, ic, v);
}/*Mfe_fun_1d_3p_zc*/


int Init_3p(void)
{
    int err, l, i, j, k, d, l0, l1,
        iz, hz, ic, hc;
    double x, a, az, b, del_ab;
    UUU3p u;
    err = 0;
    // Nf = Mfe_Nf;
    /**/
    nf2 = 2 * Nf;
    iiz = (NNz + 1) * Nf; /*число точек интегрирования*/
    iic = (NNc + 1) * Nf; /*число точек интегрирования*/
    del_ab = 0.9;
    /*
    cMin=MfeDab*(1-del_ab); //минимальное расст. между ядрами
    cMax= MfeDab*(1+2*del_ab); // максимальное расст. между ядрами
         zMin=-E_Space/2-MfeDab;
         zMax=E_Space/2;
    */
    zh = (zMax - zMin) / iiz; /*шаг интегрирования*/
    ch = (cMax - cMin) / iic; /*шаг интегрирования*/
    LL = (NNz) * Nh * (NNc) * Nh;
    if ((ind3p = (unsigned long *)calloc( sizeof(unsigned long), LL + 1 )) == NULL)
    {
        sprintf(str, "No memory for ind3p[] ");
        err = 1; return err ;
    }
    MAX_v = 0;
    if (LL > MAX_v)MAX_v = LL;
    if (iiz > MAX_v)MAX_v = iiz;
    if (iic > MAX_v)MAX_v = iic;
    /*индексы для разложения функций*/
    l = 0;
    for (iz = 0; iz < NNz; iz++) {
        for (hz = 0; hz < 2; hz++) {
            for (ic = 0; ic < NNc; ic++) {
                for (hc = 0; hc < 2; hc++) {
                    u.zrc.hz = hz;
                    u.zrc.iz = iz;
                    u.zrc.hc = hc;
                    u.zrc.ic = ic;
                    ind3p[l] = u.i;
                    sprintf(str, "ind3p[%3d,%2d,%d,%2d,%d]",
                            l, u.zrc.iz, u.zrc.hz, u.zrc.ic, u.zrc.hc);
                    FastWriteA (str, 5, 5 , clWhite, clRed);
                    l++;
                }/*hc*/
            }/*ic*/
        }/*hz*/
    }/*iz*/


    /*calloc begin*/
    if ((    z    = (double *)calloc( sizeof(double), iiz + 2)) == NULL)
    {
        sprintf(str, "No memory for z... ");         err = 1; return err ;
    }
    if ((    c    = (double *)calloc( sizeof(double), iic + 2)) == NULL)
    {
        sprintf(str, "No memory for c... ");         err = 1; return err ;
    }

    if ((    Uz    = (double *)calloc( sizeof(double), iiz + 2)) == NULL)
    {
        sprintf(str, "No memory for z... ");         err = 1; return err ;
    }

    if ((    zband    = (double *)calloc( sizeof(double), iiz + 2)) == NULL)
    {
        sprintf(str, "No memory for zband... ");      err = 1; return err ;
    }

    if ((    Uc    = (double *)calloc( sizeof(double), iic + 2)) == NULL)
    {
        sprintf(str, "No memory for Uc... ");         err = 1; return err ;
    }
    for (i = 0; i <= iiz; i++)
    {
        if ((    Uzc[i]    = (double *)calloc( sizeof(double), iic + 2)) == NULL)
        {
            sprintf(str, "No memory for Uzc[%d]... ", i);         err = 1; return err ;
        }
        if ((    buf_cz[i]    = (double *)calloc( sizeof(double), iic + 2)) == NULL)
        {
            sprintf(str, "No memory for buf_zc[%d]... ", i);         err = 1; return err ;
        }

    }/*i*/


    if ((    v    = (double *)calloc( sizeof(double), MAX_v + 2)) == NULL)
    {
        sprintf(str, "No memory for v... ");         err = 1; return err ;
    }
    if ((    f_v    = (double *)calloc( sizeof(double), MAX_v + 2)) == NULL)
    {
        sprintf(str, "No memory for f_v... ");         err = 1; return err ;
    }
    if ((    v_ground    = (double *)calloc( sizeof(double), MAX_v + 2)) == NULL)
    {
        sprintf(str, "No memory for f_v... ");         err = 1; return err ;
    }

    if ((    f    = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
    {
        sprintf(str, "No memory for f... ");         err = 1; return err ;
    }

    if ((    fz    = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
    {
        sprintf(str, "No memory for fr... ");         err = 1; return err ;
    }

    for (i = 0; i < 2; i++)
    {
        if ((fi[i] = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
        {
            sprintf(str, "No memory for ff ");         err = 1; return err ;
        }
    }

    for (i = 0; i < 2; i++)
    {
        if ((d_fi[i] = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
        {
            sprintf(str, "No memory for ff ");         err = 1; return err ;
        }
    }

    for (i = 0; i < 2; i++)  for (j = 0; j < 2; j++)  for (k = 0; k < 3; k++)
                if ((ff[i][j][k] = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
                {
                    sprintf(str, "No memory for ff ");         err = 1; return err ;
                }
    for (i = 0; i < 2; i++)  for (j = 0; j < 2; j++)  for (k = 0; k < 3; k++)
                if ((fdf[i][j][k] = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
                {
                    sprintf(str, "No memory for fdf ");         err = 1; return err ;
                }
    for (i = 0; i < 2; i++)  for (j = 0; j < 2; j++)  for (k = 0; k < 3; k++)
                if ((dfdf[i][j][k] = (double *)calloc( sizeof(double), 2 * Nf + 2)) == NULL)
                {
                    sprintf(str, "No memory for dfdf ");         err = 1; return err ;
                }


    /*calloc end*/
    /*заполнение массива координат*/
    for (i = 0; i <= iiz; i++)
    {
        z[i] = zMin + i * zh;
        for (j = 0; j <= iic; j++)
        {
            c[j] = cMin + j * ch;
            Uzc[i][j] = sqr(klam * FunPot_cz(c[j], z[i]));
        }/*j*/
    }/*i*/
    for (i = 0; i <= iic; i++)
    {
        c[i] = cMin + i * ch;
    }/*i*/

    /*----------заполнение массива финитных функций начало---------*/
    for (i = 0; i <= nf2; i++) {
        x = i * 1.0 / Nf;
        fi[0][i] =      Fb(x - 1.0);
        fi[1][i] =     FbH(x - 1.0);
        d_fi[0][i] =   DFb(x - 1.0);
        d_fi[1][i] =  DFbH(x - 1.0);
    }/*i*/
    /*----------заполнение массива финитных функций конец---------*/

    /*----------иллюстрация вида ff и  fdf dfdf --------*/
    for (i = 0; i < 2; i++)
      for (j = 0; j < 2; j++)
        for (k = 0; k < 3; k++) {
                l0 = 0; l1 = nf2; d = 0;
                if (k == 0) {
                    l1 = Nf;
                    d =  Nf;
                }
                if (k == 2) {
                    l0 = Nf;
                    d = -Nf;
                }
                for (l = l0; l <= l1; l++) {
                    a =  ff[i][j][k][l] = fi[i][l] * fi[j][l + d];
                    b = fdf[i][j][k][l] = fi[i][l] * d_fi[j][l + d];
                    b = dfdf[i][j][k][l] = d_fi[i][l] * d_fi[j][l + d];
                }/*l*/
            }/*i,j,k*/
    /*
      y=x/b b-расстояние между узлами
      [f(y)],[fh(y)*b]- финитные функции
      d/dx[f(y)]=1/b*d/dy[f(y)],
      d/dx[fh(y)]=d/dy[f(y)],
      нормировка производится из требования равенства производной
      fh(0)=1 в узловых точках (необходимо для быстрого разложения
      функций по известным значениям функций и производных в узлах )
    */
    /*-------------z---------------*/
    az = 1 / (zh * Nf);
    a = az;
    zzN[0][0] = 1;
    zzN[0][1] = 1 * (1 / a);
    zzN[1][0] = (1 / a) * 1;
    zzN[1][1] = (1 / a) * (1 / a);

    zdzN[0][0] = 1 * a;
    zdzN[0][1] = 1 * 1;
    zdzN[1][0] = (1 / a) * a;
    zdzN[1][1] = (1 / a) * 1;

    dzdzN[0][0] = a * a;
    dzdzN[0][1] = a * 1;
    dzdzN[1][0] = 1 * a;
    dzdzN[1][1] = 1 * 1;
    /*------------zc---------------*/
    b = 1 / (ch * Nf);
    zcN[0][0] = 1;
    zcN[0][1] = 1 * (1 / b);
    zcN[1][0] = (1 / a) * 1;
    zcN[1][1] = (1 / a) * (1 / b);

    /*-----------c----------------*/
    a = 1 / (ch * Nf);
    ccN[0][0] = 1;
    ccN[0][1] = 1 * (1 / a);
    ccN[1][0] = (1 / a) * 1;
    ccN[1][1] = (1 / a) * (1 / a);

    cdcN[0][0] = 1 * a;
    cdcN[0][1] = 1 * 1;
    cdcN[1][0] = (1 / a) * a;
    cdcN[1][1] = (1 / a) * 1;

    dcdcN[0][0] = a * a;
    dcdcN[0][1] = a * 1;
    dcdcN[1][0] = 1 * a;
    dcdcN[1][1] = 1 * 1;
    /*------------r---------------*/
    return err;
}/* Init_3p*/
void Mfe_free()
{
    int i, j, k;
    free(ind);
    free(z);
    free(c);
    free(Uz);
    free(Uc);
    free(zband);
    for (i = 0; i <= iiz; i++)
    {
        free(Uzc[i]);
        free(buf_cz[i]);
    }
    free(v);
    free(f_v);
    free(v_ground);
    free(f);
    free(fz);
    for (i = 0; i < 2; i++)
    {
        free(fi[i]);
        free(d_fi[i]);
    }/*i*/

    for (i = 0; i < 2; i++)  for (j = 0; j < 2; j++)  for (k = 0; k < 3; k++)
            {
                free(ff[i][j][k]);
                free(fdf[i][j][k]);
                free(dfdf[i][j][k]);
            }/*i,j,k,*/
}/*Mfe_free*/

void strip_waveguide(void) {
    FILE *Inf;
    int ix, iy, ix1, iy1;
    int i, k, j, l, n, N, NN, M, err, ig, jg;
    double eps, a, b, x, nm, beta, beta_max, nm_koshiba, beta_koshiba;

    double *HE[MaxM], *A[MaxM], *B[MaxM],
           *Astore[MaxM],//запом.матрица дип. мом.на первом шагу
           *ST[MaxM],//хранится треугольное разложение матрицы перекрывания
           *BasTime[MaxM],// список функций для временного развития
           *X, *Y, *D, *DLbuf, *Ebuf;

    double L, Lw, Pol, A00, A01,
           W, U, Nr,
           MulMass = 1,
           norm,
           a1, aa, aa0, aa1, bb,
           cc, zz,
           d, hx, h_tim,
           zav, cav,
           msy, mst,
           pozcl,/*положение свободной классической частицы в поле*/
           pozcl0, pozcl1, pozcl2, pozcl3,
           h_mark_tim,
           fcl, fcl1, fcl2 ;
    TColor ColAr[16];

    int kk, k1, il, is,
        s_in, s_fin,// начальный и конечный номер баз. функций
        only_Stationary,
        shift_graph,
        kc, kz,
        iz, izh, ic, ich,
        count,
        z0, z1,
        c0, c1,
        ix0, iy0,
        ixa0, iya0,
        Width, Height, ixcl, iycl;

    err = 0;

    tim0 =  time(&timt);
    err = 0;
    top = left = 0; PanelSize (&right , &bottom); //уст. размер граф.панели
    setfillstyle(1, clWhite);
    bar(left, top, right, bottom);
    XoX = (left + right) / 11;
    YoY = (top + bottom) / 2;
    XoX = (right) / 2;
    xmin = -4; xmax = 4;
    ymin = -2 ; ymax = 3;
    Mode = Form1->Edit1->Text.ToInt();
    nref_substrat = Form1->Substrat->Text.ToDouble();
    nref_max =     Form1->Nref_max->Text.ToDouble();
    h_waveguide =     Form1->H_waveguide->Text.ToDouble();
    w_waveguide =     Form1->W_waveguide->Text.ToDouble();
    lam =     Form1->Lam->Text.ToDouble();

    if (Mode == 1) {
        W_test = 2,
        t_test = 0.2,
        h_test = 1.1,
        Xs_test = 3.0,
        Ys_test = 5.0,
        Yc_test = 1.0;
        beta_koshiba = 0.4981;
        nm_koshiba =  sqrt(beta_koshiba * (sqr(nref_max) - sqr(nref_substrat)) + sqr(nref_substrat));
    }
    if (Mode == 2) {
        W_test = 3,
        t_test = 0.9,
        h_test = 0.1,
        Xs_test = 3.0,
        Ys_test = 5.0,
        Yc_test = 1.0;
        beta_koshiba = 0.4380;
        nm_koshiba =  sqrt(beta_koshiba * (sqr(nref_max) - sqr(nref_substrat)) + sqr(nref_substrat));
    }
    if (Mode == 3) {
        W_test = 4,
        t_test = 3.5,
        h_test = 2.5,
        Xs_test = 4.5,
        Ys_test = 7.5,
        Yc_test = 1.0;
        beta_koshiba = 0.3974;
        nm_koshiba =  sqrt(beta_koshiba * (sqr(nref_max) - sqr(nref_substrat)) + sqr(nref_substrat));
    }

    if (Mode == 4) {
        W_test = 3,
        t_test = 0.7,
        h_test = 0.3,
        Xs_test = 3.0,
        Ys_test = 5,
        Yc_test = 1.0;
        beta_koshiba = 0.3974;
        //эта штука 4 раза одна и та же
        nm_koshiba =  sqrt(beta_koshiba * (sqr(nref_max) - sqr(nref_substrat)) + sqr(nref_substrat));
    }
    /*--------
    //===============Koshiba
    xmax=Xs_test+W_test/2; xmin=-xmax;
    ymax=Yc_test+h_test;
    ymin=-(t_test+Ys_test/3);
    //===========Koshiba
    ------------*/

    xmin = -20;
    xmax = 20;
    ymin = -20;
    ymax = 20;



    cMin = xmin; cMax = xmax;
    zMin = ymin; zMax = ymax;
    hrx = 1.5;
    hry = 2;
    //    Mode=1;

    //    nref_substrat=1.2;
    //    nref_max=1.9;

    al_nref = 1 / h_waveguide;
    klam = 2 * M_PI / lam;


    /*стационарная задача*/
    only_Stationary = 1;


    //  M1=1836.154; //масса первого ядра
    //  M2=1836.154; //масса второго ядра
    M1 = M2 = 1;
    me = 1;      //масса электрона
    L = E_Space;

    NNc = 11;
    NNz = 11;
    s_in = 0; s_fin = 100;
    DebugMain = 0;
    tim0 = time(&tim0);
    //clearviewport();
    /*===========*/
    if ((Inf = fopen("mfe.rez", "at")) == NULL) {
        sprintf(str, "error write File");  goto Mexit;
    }
    //===============================
    Init_3p();
    //=============================
    NN = LL;
    /*===========*/
    DosGetDateTime(&Dat);
    sprintf (str,
             "**********Polozok: %2d-%2d-%4d   [%2d:%2d:%2d] **********",
             Dat.day, Dat.month,
             Dat.year,//-(Dat.year/100)*100,
             Dat.hours, Dat.minutes, Dat.seconds );
    /*------------*/
    fprintf(Inf, "\n\n%s", str);
    sprintf(str,
            "Test: W=%6.3f t=%6.3f h=%6.3f Xs=%6.3f Ys=%f Yc=%f NN=%d(%dx%d) ",
            W_test, t_test, h_test,
            Xs_test, Ys_test, Yc_test,
            NN, NNz, NNc);
    FastWriteA (str, 2, 1, clBlack, clWhite);
    fprintf(Inf, "\n n_max=%f nref_substrat=%f",
            nref_max, nref_substrat);
    fprintf(Inf, "\n%s", str);
    fprintf(Inf, "\n xmin=%f xmax=%f", xmin, xmax);
    fprintf(Inf, "\n ymin=%f ymax=%f", ymin, ymax);
    fprintf(Inf, "\n Mode=%d nm_koshiba=%f", Mode, nm_koshiba);
    setcolor(clBlack);
    line(left, YoY, right, YoY);
    line(XoX, top, XoX, bottom);

    if (NN >= MaxM)  {
        sprintf(str, "%d >= MaxM=%d (mfe_supp.h)", NN, MaxM);
        err = -1; goto Mexit;
    }
    //=========== Стационарная задача============
    for (i = 0; i <= NN; i++)
    {
        if ((HE[i] = (double *)calloc( sizeof(double), NN + 1 )) == NULL)
        {
            sprintf(str, "No memory for HE[] ");
            err = 1; goto Mexit ;
        }

        if ((B[i] = (double *)calloc( sizeof(double), NN + 1 )) == NULL)
        {
            sprintf(str, "No memory for B[] ");
            err = 1; goto Mexit ;
        }
    }/*i*/


    if ((DLbuf = (double *)calloc( sizeof(double), NN + 1 )) == NULL)
    {
        sprintf(str, "No memory for DLbuf[] ");
        err = 1; goto Mexit ;
    }
    if ((D = (double *)calloc( sizeof(double), NN + 1 )) == NULL)
    {
        sprintf(str, "No memory for DLbuf[] ");
        err = 1; goto Mexit ;
    }
    if ((Ebuf = (double *)calloc( sizeof(double), NN + 1 )) == NULL)
    {
        sprintf(str, "No memory for Ebuf[] ");
        err = 1; goto Mexit ;
    }


    /*
    */
    /*--------------------------------------------*/
    msx = 600 / L;
    msy = 10;
    msx = (right - left) / (xmax - xmin);
    msy = (bottom - top) / (ymax - ymin);
    //Axes();
    //mst=right*0.75/(Mul_tim_fin*tim_fin-tim_in);
    //Hfont=(bottom-12)/40;

    hx = L / 500;
    setcolor(clGray);
    //разметка сплайнов
    for (i = 0; i < NNz; i++) {
        Region(-(NNz - 1), i, i, &z0, &z1);
        if (z1 > iiz) {
            sprintf(str, "Nf=%d,iiz=%d,z0=%d,z1=%d", Nf, iiz, z0, z1);
            FastWriteA (str, 10, 10, clBlack, clWhite);
            goto Mexit;
        }
        ix = (zMin + zh * ((z0 + z1) / 2)) * msx + XoX;
        iy = YoY;
        line(ix, bottom - iy, ix, bottom - (iy + 10));
        Mark(ix, bottom - iy, 2, clBlue);
    }/**/
    ix = zMin * msx + XoX; iy = YoY; line(ix, bottom - iy, ix, bottom - (iy + 100));
    ix = zMax * msx + XoX; iy = YoY; line(ix, bottom - iy, ix, bottom - (iy + 100));
    setcolor(clGray);
    for (x = zMin; x < zMax; x += 1)
    {
        ix = x * msx + XoX;
        iy = YoY;
        line(ix, bottom - iy, ix, bottom - iy - 150);
    }/*x*/
    i = 0;
    for (x = zMin; x < zMax; x += hx) {
        ix = x * msx + XoX;
        iy = (FunPot_cz(0, x) - nref_substrat) * msy * 10 + YoY;
        putpixel(ix, bottom - iy, clBlue);
        iy = (FunPot_cz(cMax, x) - nref_substrat) * msy * 10 + YoY;
        putpixel(ix, bottom - iy, clBlue);
        i++;
    }/*x*/
    /*-----прорисовка потенциала конец---------*/
    /*-----заполнение матрицы энергии начало--*/
    DebugMain = 0;
    k = 0;
    for (i = 0; i < NN; i++) {
        if (GetAsyncKeyState(VK_ESCAPE)) {
            sprintf(str, "Operator  Break !!!");
            err = -777; goto Mexit;
        }/*if(ESC)*/
        sprintf(str, "HE[%d][]", i);
        FastWriteA (str, 4, 30, clBlack, clWhite);
        for (j = i; j < NN; j++) {
            HE[i][j] = MatElT_1d_3p( i, j) + MatElV_1d_3p( i, j);
            B[i][j] = B[j][i] = MatElS_1d_3p( i, j);
            HE[j][i] = HE[i][j];
        }/*j*/
    }/*i*/
    /*-----заполнение матрицы энергии конец--*/

    err = TREQLAB(LL, HE, B, D, DLbuf, Ebuf); if (err != 0)  goto Mexit;

    /*восстановление B[i][j]=S*/
    for (i = 0; i < NN; i++)
    {
        for (j = i; j < NN; j++)
        {
            B[i][j] = B[j][i] = MatElS_1d_3p( i, j);
        }/*j*/
    }/*i*/
    only_Stationary = 1;
    clearviewport();
    sprintf(str, "c=(%f, %f)  z=(%f  ,%f) ",
            cMin, cMax, zMin, zMax);
    FastWriteA (str, 1, 1, clWhite, clBlack);
    /*----------прорисовка контуров----начало------*/
    if (only_Stationary == 1) {
        for (i = 0; i <= NzMap; i++) {
            if ((Map[i] = (int *)calloc( sizeof(int), NcMap + 1 )) == NULL) {
                sprintf(str, "No memory for Map[%d] ", i);
                err = 1; goto Mexit ;
            }
        }/*i*/
        clearviewport();
        dump(HE, LL, LL, "dump2d.csv");
        iY = 25;
        iX = 20;
        //for(is=s_in;is<=s_fin;is++){
        for (is = NN - 1; is >= 0; is--) {
            if (is < NN - 1)break;
            //  is=17;

            iX = 20;
            if ((is / 2) * 2 != is)iX = 400;
            sprintf(str, "F[%3d]", NN - 1 - is);
            ix = iX / Hfont; iy = 1;
            FastWriteA (str, iy, ix, clBlack, clWhite);
            // clearviewport();
            //setfillstyle(0, clBlack);
            //rect(iX,iY,iic+iX,iiz+iY);
            for (i = 0; i < LL; i++) {
              v[i] = (HE[i][is]);
            }
            sprintf(str, "--%3f--", Mfe_fun_1d_3p_zc(-0.5, 0, v));
            FastWriteA (str, iy, ix, clBlack, clWhite);
            //нахождение максимума функции
            aa = 0;
            for (i = 0; i < iiz; i++) {
                for (j = 0; j < iic; j++)
                {
                    a = Mfe_fun_1d_3p(i, j, v); a = a * a;
                    if (a > aa)aa = a;
                }/*j*/
            }/*i*/
            if (fabs(a) < 1.0e-10)aa = 1;

            zhMap = (zMax - zMin) / NzMap;
            chMap = (cMax - cMin) / NcMap;
            for (i = 0; i < NzMap; i++) {
                zz = i * zhMap + zMin;
                for (j = 0; j < NcMap; j++) {
                    cc = j * chMap + cMin;
                    a = Mfe_fun_1d_3p_zc(zz, cc, v); a = a * a;
                    b = (a / aa) * 32000;
                    //                   if(b<-32000)b=-32000;
                    Map[(NzMap - 1) - i][j] = b;
                }/*j*/
            }/*i*/
            Contur(  Map, NcMap, NzMap, 10, iX, iY);
            setcolor(clGray);
            for (x = -20; x < zMax; x++)if (x > zMin && x < zMax)
                {
                    i = (zMax - zMin) / zhMap - (x - zMin) / zhMap;
                    if (i > NzMap)break;
                    line(iX, iY + i, iX + NcMap, iY + i);
                    sprintf(str, "%4.1f", x);
                    i += iY - Hfont / 2;
                    ix = iX + NcMap;
                    gprintfColor(&ix, &i , str, clBlack);
                }/*x*/
            for (x = -20; x < cMax; x++)if (x > cMin && x < cMax)
                {
                    i = (x - cMin) / chMap;
                    if (i > NcMap)break;
                    line(iX + i, iY, iX + i, iY + NzMap);
                    sprintf(str, "%3.1f", x);
                    i += iX;
                    iy = iY + NzMap;
                    gprintfColor(&i, &iy , str, clBlack);
                }/*x*/
        }/*is*/
        for (i = 0; i <= NzMap; i++)free(Map[i]);
    }/*only_Stationary==1*/
    /*----------вывод спектра начало----------------*/
    j = 0;
    for (k = NN - 1; k >= 0; k--) {
        nm = sqrt(fabs(D[k])) / klam;
        beta = (sqr(nm) - sqr(nref_substrat)) /
               (sqr(nref_max) - sqr(nref_substrat));
        if (nm < nref_substrat )break;

        //     if(k<NN-2)break;
        sprintf(str, "[%2d]    Nm=%f beta=%f ",  j, nm, beta);
        FastWriteA (str, NzMap / (Hfont + 4) + 4 + j, 1, clBlack, clWhite);
        fprintf(Inf, "\n%s", str);
        j++;
    }/*k*/
    sprintf(str, "max_Mode=%2d",  j);
    FastWriteA (str, NzMap / (Hfont + 4) + 4, 30, clBlack, clWhite);
    //=====================================
    goto Mexit;
    int fl;
    double r, hy, fi, hr, pc, pz, Sb, Sa, Sab;
    //*RCyl,*FCyl;
    is = NN - 1;

    // HRCyl=RCyl[NCyl-1]/NCyl;
    for (i = 0; i < 256; i++)
    {
        rgbcol.RGB[3] = (unsigned char)0;
        rgbcol.RGB[0] = (unsigned char)(i);
        rgbcol.RGB[1] = (unsigned char)(i);
        rgbcol.RGB[2] = (unsigned char)(i);
        ColSet[i] = rgbcol.c;
    }

    for (i = 0; i < LL; i++)v[i] = (HE[i][is]);

    fl = 0;
    for (pc = 0; pc < 2; pc += 0.1)
        for (pz = 0; pz > -1; pz -= 0.1)
        {
            if (GetAsyncKeyState(VK_ESCAPE))
            {
                sprintf(str, "Operator  Break !!!");
                err = -777; goto Mexit;
            }/*if(ESC)*/

            Sa = Sb = Sab = 0;
            if (fl == 0)aa = 0;
            for (i = 0; i < iiz; i++)
            {
                for (j = 0; j < iic; j++)
                {
                    zz = i * zh + zMin;
                    cc = j * ch + cMin;
                    r = sqrt(sqr(zz - pz) + sqr(cc - pc));
                    k = r / HRCyl;
                    if (k > NCyl - 1)k = NCyl - 1;
                    b = fabs(FCyl[k]); Sb += b * b * zh * ch;
                    a = Mfe_fun_1d_3p(i, j, v); Sa += a * a * zh * ch;
                    //a=1;
                    a = a * b;          Sab += a * zh * ch;
                    if (a > aa && fl == 0)aa = a;
                }/*j*/
            }/*i*/
            fl = 1;
            sprintf(str, "aa=%f",  aa);
            FastWriteA (str, 1, 1, clBlack, clWhite);
            for (i = 0; i < iiz; i++)
            {
                for (j = 0; j < iic; j++)
                {
                    zz = i * zh + zMin;
                    cc = j * ch + cMin;
                    r = sqrt(sqr(zz - pz) + sqr(cc - pc));
                    k = r / HRCyl;
                    if (k > NCyl - 1)k = NCyl - 1;
                    b = fabs(FCyl[k]);
                    a = Mfe_fun_1d_3p(i, j, v);
                    //a=1;
                    a = a * b ;

                    k = 0.7 * ((a) / aa) * 255;
                    putpixel(j, bottom - i, ColSet[k]);
                    // if(fabs((zz-pz)*(cc-pc))<0.01)putpixel(j,bottom- i,clBlue);
                }/*j*/
            }/*i*/
            setcolor(clBlue);
            j = (pz - zMin) / zh;
            i = (pc - cMin) / ch;
            line(i - 5, bottom - j, i + 5, bottom - j);
            line(i, bottom - (j - 5), i, bottom - (j + 5));
            Sab = sqr(Sab) / (Sa * Sb);
            sprintf(str, "Sab=%f Sa=%f Sb=%f", Sab, Sa, Sb / (2 * M_PI));
            FastWriteA (str, 1, 20, clBlack, clWhite);


        }/*pc,pz*/


    /*----------вывод спектра конец----------------*/
    /*----------прорисовка контуров----конец------*/



    /*


        for(x=-h_waveguide;x<h_waveguide;x+=hx)
       for(y=-h_waveguide;y<h_waveguide;y+=hy){
          r=sqrt(x*x+y*y);
          if(fabs(r)<1.0e-10)fi=0;
          else{
             if(x>=0)fi=asin(y/r);
             else fi=M_PI - asin(y/r);
          }
           ix=x*msx*0.5+XoX+k*2.1*h_waveguide*msx*0.5;
           iy=bottom-(y*msx*0.5+YoY-h_waveguide*msx);
           a=ModeFun_cyl(m, M-k, M,HE,r)*cos(m*fi);
           a=a*a;
           j=(a/aa)*255;
           putpixel(ix,iy,ColSet[j]);
       }

    */


    /*--------------------------------------------*/
Mexit:;
    for (i = 0; i <= NN; i++)
    {
        free(HE[i]);
        free(B[i]);
    }/*i*/
    free(DLbuf);
    free(D);
    free(Ebuf);
    Mfe_free();

    if (err != 0)FastWriteA(str, 2, 2, clBlack, clRed);
    fclose(Inf);
}/*strip_waveguide()*/



double Ftest(int k, double x)
{
    double a;
    a = 1;
    /*
    if(k==0) return sin(x*a);
    if(k==1)  return a*cos(x*a);
    if(k==2) return -a*a*sin(x*a);
    */
    if (k == 0) return exp(-x * x);
    if (k == 1)  return -2 * x * exp(-x * x);
    if (k == 2) return -2 * exp(-x * x) + sqr(2 * x) * exp(-x * x);

    return    0;

}/*Ftest*/
void Test_splines(void)
{
    double C0[100], C1[100];

    int N, ix1, iy1, k, i, j, n;

    double hx, a, y1, df, d0, d1;
    xmin = -4;
    xmax = 4;
    ymin = -1.5;
    ymax = 1.5;
    top = left = 0; PanelSize (&right , &bottom); //уст. размер граф.панели
    //устанавливаем цвет заливки и заливаем панель
    setfillstyle(1, clWhite);
    bar(left, top, right, bottom);

    XoX = (left + right) / 2; // начало координат в
    YoY = (top + bottom) / 2; // в центре экрана
    // вычисляем масштабные коэффициенты
    msx = (right - left) / (xmax - xmin);
    msy = (bottom - top) / (ymax - ymin);
    // устанавливаем шаг разметки осей
    hrx = 1.0;
    hry = 0.5;
    // размечаем ось x с шагом hrx
    Axes();
    hx = hrx / 100;
    for ( k = 0; k < 3; k++)
    {
        N = 0;
        for (x = xmin ; x < xmax; x += hx)
        {
            y = Ftest(k, x);
            ix = x * msx + XoX;
            iy = bottom - (y * msy + YoY);
            if (x == xmin)
            {
                ix1 = ix;
                iy1 = iy;
            }
            line(ix, iy, ix1, iy1);
            ix1 = ix; iy1 = iy;
            if ((N / 10) * 10 == N)Mark(ix, iy, 1, clBlack);
            N++;
        }/*y*/
    }/*k*/
    n = N;
    N = 12;
    for (j = 0; j <= N; j++)
    {
        C0[j] = Ftest(0, j - N / 2); C1[j] = Ftest(1, j - N / 2);
    }/*j*/
    d0 = 0;
    setcolor(clBlue);
    for (x = xmin ; x < xmax; x += hx)
    {
        a = 0;
        for (j = 0; j <= N; j++)
        {
            a += C0[j] * Fb(x - (j - N / 2)) + C1[j] * FbH(x - (j - N / 2));
        }/*j*/
        y = a; d0 += fabs(Ftest(0, x) - a);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        if (x == xmin)
        {
            ix1 = ix;
            iy1 = iy;
        }
        line(ix, iy, ix1, iy1);
        ix1 = ix; iy1 = iy;
    }/*y*/
    d0 /= n;
    d1 = 0;
    for (x = xmin ; x < xmax; x += hx)
    {
        a = 0;
        for (j = 0; j <= N; j++)
        {
            a += C0[j] * DFb(x - (j - N / 2)) + C1[j] * DFbH(x - (j - N / 2));
        }/*j*/
        y = a;  d1 += fabs(Ftest(1, x) - a);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        if (x == xmin)
        {
            ix1 = ix;
            iy1 = iy;
            y1 = y;
        }
        line(ix, iy, ix1, iy1);
        df = (y - y1) / hx;
        iy1 = bottom - (df * msy + YoY);
        putpixel(ix, iy1, clBlue);
        ix1 = ix; iy1 = iy; y1 = y;
    }/*y*/
    d1 /= n;
    sprintf(str, "d0=%e d1=%e",  d0, d1);
    FastWriteA (str, 1 + 1, 30, clBlack, clWhite);

}/*Test_splines*/
//=================================
//=================================
