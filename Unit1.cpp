//---------------------------------------------------------------------------
#include <vcl\vcl.h>
#pragma hdrstop

#include "Unit1.h"
#include <math.h>
#include <stdio.h>
//---------------------------------------------------------------------------
#pragma resource "*.dfm"

#define  sqr(x)  ((x)*(x))
#define  MAT    100
#define  MaxM   200
#define  N_MAP_STR  1000
//------------------------------
typedef struct
{
    double p1, p2;
} PP2;


typedef struct
{
    char obj;
    double Poz;
    double diam;
    double foc;
    double app;
} OPT;

typedef struct
{
    double x, y, z;
} XYZ;

typedef   union
{
    TColor c;
    unsigned char RGB[4];
} RGB_Col;

typedef struct
{
    int i, j;
} IJ;




//----------------------------
TForm1 *Form1;
int Hfont = 8;
int kl = 1;
int left, right, top, bottom, XoX, YoY,
    ix, iy, ix1, iy1;
TCanvas *Canva = 0;
double hrx, hry, msx, msy, x, y,
       xmin, ymin, xmax, ymax, hx, hy;
char str[100];
//-------Volume1-----------------
double t, t0, p0 = 5, p10 = 4, s = 0.1,
              s10 = 0.001, hreg = 0.5, P10;
//---------Volume2----------------------
PP2 P, P0, D, B;
double s12 = 0.1, s20 = 0.001, p20 = 3.5;

/*===============Planar=================*/
FILE *Fwr, Frd;
double hw = 7,
       nsubstrat = 1.2,
       ncover = 1,
       nmax = 1.5;
double Lbas;
int Mode = 0,
    Mtot = 20,
    N_tot,
    M_bas;
double *Cbuf;
int Plan1 = 0, NxMap, NyMap = 300;
double xMax, xMin, yMax, yMin;
double *Fi_bas[N_MAP_STR], *n_xy[N_MAP_STR];
double *yy1;
double h_x, h_y;
IJ *ind;

//============Cylynder==================================
int M_fi = 0, M_r = 0;
TColor ColSet[256];
RGB_Col rgbcol;
//================Potential Well=======================
double max_pot = -5;

//-----------------------------
void putpixel(int x, int y, TColor cl)
{
    if (Canva)
    {
        Canva->Brush->Color = cl;
        Canva->Brush->Style = bsSolid;
        Canva->FillRect(Rect(x, y, x + 1, y + 1));
    }
};

//Установка цвета линии
void setcolor(TColor x)
{
    if (Canva)    Canva->Pen->Color = x;
};
//Прорисовка линии из точки (x,y) в точку (x1,y1)
void line(int x, int y, int x1, int y1)
{
    if (Canva)
    {
        Canva->MoveTo(x, y);
        Canva->LineTo(x1, y1);
    }
};
//Вывод текста начиная с графических координат (x,y)
void gprintf(int *x, int *y, char *s)
{
    int x1, y1;
    x1 = *x; y1 = *y;
    if (Canva)
    {
        Canva->Font->Color = clBlack;
        SetBkMode(Canva->Handle, TRANSPARENT);
        Canva->TextOut(x1, y1, s);
    }
};

void FastWriteA(char *s, int y, int x,
                TColor clLet, TColor clFon)
{
    if (Canva)
    {
        Canva->Font->Color = clLet;
        SetBkColor(Canva->Handle, clFon);
        SetBkMode(Canva->Handle, OPAQUE);
        Canva->TextOut((x)*Hfont, (y)*Hfont, s);
    }
};

void Mark(int i, int j, int h, TColor col)
{
    int k, l;
    for (k = -h; k <= h; k++)
        for (l = -h; l <= h; l++)
            putpixel(i + k, j + l, col);
}/*Mark*/

void setfillstyle(int x, TColor cl)
{
    if (Canva)    Canva->Brush->Color = cl;
};
void bar(int x, int y, int x1, int y1)
{
    if (Canva) Canva->FillRect(Rect(x, y, x1, y1));
};

void rect(int x, int y, int x1, int y1)
{
    if (Canva)    Canva->Rectangle(x, y, x1, y1);
};

//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent *Owner)
    : TForm(Owner)
{
    DecimalSeparator = '.';
    sprintf(str, "%d", Mode); Edit1->Text = str;

    sprintf(str, "%d", M_fi); Edit2->Text = str;
    sprintf(str, "%d", M_r); Edit3->Text = str;
    sprintf(str, "%g", hw); Edit4->Text = str;
    sprintf(str, "%g", nsubstrat); Edit5->Text = str;

}
//---------------------------------------------------------------------------
double FF(double x)
{
    return sin(x * x);

}/*FF(x)*/
/*Button1*/
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
            fprintf(Fwr, "%d", Arr[j][i]);
        }
    }
    fclose(Fwr);
}
//---------------------------------------------------------------------------
double PP(double t)
{
    double C, al, b;
    al = -(s10 + s * kl);
    b = -p0 * s10 / al;
    C = P10 - b;
    return  C * exp(al * t) + b;
}/*PP*/

double Freg(double t, double hreg)
{
    return -(p10 / 500) * t + p10 + hreg;
} /*Freg*/
//---------------------------------------------------------------------------/*Volume1*/
//---------------------------------------------------------------------------
void root2(double A, double B, double C,
           double *x1, double *x2)
{
    double d;
    d = B * B - 4 * A * C;
    *x1 = (-B + sqrt(d)) / (2 * A);
    *x2 = (-B - sqrt(d)) / (2 * A);
}/*root2*/

PP2   P2(double t)
{
    PP2 p;
    double al1, al2, C1, C2, k1, k2,
           m11, m12, m21, m22, d, d1, d2;
    m11 = -(s10 + s12); m12 = s12;
    m21 = s12; m22 = -(s * kl + s12 + s20);
    d = m11 * m22 - m12 * m21;
    D.p1 = p0 * s10; D.p2 = p0 * s20;
    d1 = -D.p1 * m22 + D.p2 * m12;
    d2 = -m11 * D.p2 + m21 * D.p1;
    B.p1 = d1 / d; B.p2 = d2 / d;
    // B.p1=0.1; B.p2=0.5;
    root2(1, -(m11 + m22), m11 * m22 - m12 * m21, &al1, &al2);
    // al1=-0.01;al2=-0.02;
    k1 = -(m11 - al1) / m12;
    k2 = -(m11 - al2) / m12;

    //   k1=1;k2=1;
    d = k2 - k1;
    d1 = (P0.p1 - B.p1) * k2 - (P0.p2 - B.p2);
    d2 = (P0.p2 - B.p2) - k1 * (P0.p1 - B.p1);
    C1 = d1 / d; C2 = d2 / d;

    // C1=1;C2=1;
    p.p1 = C1 * exp(al1 * t) + C2 * exp(al2 * t) + B.p1;
    p.p2 = C1 * k1 * exp(al1 * t) + C2 * k2 * exp(al2 * t) + B.p2;
    return p;
}/*P2*/

//-------------------------------------/*Volume2*/
//---------------------------------------------------------------------------

/**/
int ancip (char buf[], double digit[])
{
    //        extern double atof();
    int count, i;
    static char dubl[20];
    char ch, *pntb, *pntd;
    int flag, fldig;
    count = flag = fldig = 0;
    pntb = buf; pntd = dubl;
    do
    {
        ch = *(pntb++);
        if ((ch == '9') || (ch == '8') ||
                (ch == '7') || (ch == '6') || (ch == '5') || (ch == '4') ||
                (ch == '3') || (ch == '2') || (ch == '1') || (ch == '0') ||
                (ch == '-') || (ch == '+') || ( ch == '.'))
        {
            flag = 1; fldig = 1;
            *(pntd++) = ch;
        }
        else
        {
            if (((ch == 'e') || (ch == 'E')) && (fldig != 0))
            {
                *(pntd++) = ch;
                fldig = 0;
                continue;
            }
            if (flag != 0)
            {
                *pntd = '\0';
                count++;
                digit [count] = atof(dubl);
                pntd = dubl;
                flag = fldig = 0;
            }
        }
    }
    while (ch != '\0');
    return (count);
}/* end ancip */

void Lens(double F, XYZ V0, XYZ *v1, XYZ p)
/*
 F-фокусное расст., V0 - исходный вeктор,
    который попал
  в точку линзы p.x,p.y оптическая ось линзы
  совпадает с осью x,
  V1 - вектор, выходящий из точки (p.x,p.y)
*/
{
    XYZ  V1;
    V1.x = 2 * F;
    V1.y = (2 * F / V0.x * V0.y - p.y) - p.y;
    *v1 = V1;
}/*Lens*/

//==================================/*OptDesk*/
//---------------------------------------------------------------------------/*Heat*/
//---------------------------------------------------------------------------
double syleqmax(int N, double A[MAT][MAT], double *X, double *B)
{
    int i, j, k, l, m, indi[MAT];
    double a, b, det = 1;
    //  DebugMain=1;
    for (i = 0; i < N; i++)indi[i] = -1;
    for (i = 0; i < N; i++)
    {
        X[i] = 0;
        a = A[i][0]; k = 0; indi[i] = k;
        for (j = 1; j < N; j++)if (fabs(A[i][j]) > fabs(a))
            {
                a = A[i][j];
                k = j;
            }
        if (k >= 0)
        {
            indi[i] = k;

            if (fabs(A[i][k]) < 1.0e-200)return 0.0;
            b = 1 / A[i][k];
            if (fabs(det) < 1.0e300)det = det * A[i][k];

            for (l = i + 1; l < N; l++)
            {
                a = A[l][k];
                for (m = 0; m < N; m++)A[l][m] = A[l][m] - b * A[i][m] * a;
                B[l] = B[l] - b * B[i] * a;
            }/*l*/
        }/*k>=0*/
        if (fabs(det) < 1.0e-200) goto Mexit;
    }/*i*/
    /*
    i=N;
    X[indi[i-1]]=B[i-1]/A[i-1][indi[i-1]];
    */
    for (i = N - 1; i >= 0; i--)
    {
        b = B[i]; for (j = 0; j < N; j++)if (j != indi[i])b -= A[i][j] * X[j];
        X[indi[i]] = b / A[i][indi[i]];
    }/*i*/
    for (i = 0; i < N; i++)
    {
        for (j = i; j < N; j++)if (i == indi[j] && indi[i] != i)
            {
                indi[j] = indi[i]; det = -det; break;
            }/*j*/
    }
Mexit:
    return det;
}/*syleqj*/


double F(double x)
{
    double s;
    s = 5;
    return 7 * exp(-sqr(x / s));
}/*F(x)*/

double fi(int k, double x)
{
    if (k == 0)return 1;
    else return pow(x, k);
}/*fi*/

double APR(int M, double *C, double x)
{
    int k;
    double a;
    a = 0; for (k = 0; k < M; k++)a += C[k] * fi(k, x);
    return a;
}/*apr*/


/*MNK*/
//--------------------------------
double fun_init(double x)
{
    if (x < 0) return 0;
    else
    {
        if (x < 3)return 5.0 / 3 * x;
        else   return  (x - 3) * (-5 / (xmax - 3)) + 5;
    }
}/*fun_init*/
double fis(int k, double x)
{
    return sin((k + 1) * M_PI * x / xmax);
}/*fi*/
double APRS(int M, double *C, double x, double t)
{
    int k;
    double a;
    a = 0; for (k = 0; k < M; k++)a += C[k] * fis(k, x) *
                                           cos((k + 1) * M_PI * t / xmax) * exp(-k / 10.0 * t);
    return a;
}/*apr*/


//---------------------------------------------------------------------------/*String*/
//Редукция Хаусхолдера действительной симметричной матрицы a[1...n][1...n]
// на выходе получается трехдиагональная матрица
// d[1...n] возвращает диагональ трехдиагональной матрицы. 
// e[1...n] возвращает внедиагональные элементы, причем e[1]=0. 
void tred2(int n, double tol, double *a[], double d[], double e[]) {
    static int i, j, k, l, nm, inm, jnm, knm;
    static double f, g, h, h1;
    nm = n;
    for (i = n; i >= 2; i += -1)
    {
        inm = nm * i;
        l = i - 2; f = a[i][ i - 1]; g = 0.0;
        for (k = 1; k <= l; k += 1)g = g + a[i][ k] * a[i][ k]; h = g + f * f;
        if (g <= tol)
        {
            e[i] = f; h = 0.0;
            goto skip;
        };
        l = l + 1;
        if (f >= 0)g = e[i] = -sqrt(h) ; else g = e[i] = sqrt(h);
        h = h - f * g; a[i][ i - 1] = f - g; f = 0.0;
        for (j = 1; j <= l; j += 1)
        {
            jnm = j * nm;
            a[j][ i] = a[i][ j] / h; g = 0.0;
            for (k = 1; k <= j; k += 1)g = g + a[j][ k] * a[i][ k];
            for (k = j + 1; k <= l; k += 1)g = g + a[k ][ j] * a[i][ k];
            e[j] = g / h; f = f + g * a[j][ i];
        };
        h1 = f / (h + h);
        for (j = 1; j <= l; j += 1)
        {
            jnm = j * nm;
            f = a[i][ j]; g = e[j] = e[j] - h1 * f;
            for (k = 1; k <= j; k += 1)
                a[j][ k] = a[j][ k] - f * e[k] - g * a[i][ k];
        };
skip:           d[i] = h;
    };
    d[1] = e[1] = 0.0;
    for (i = 1; i <= n; i += 1)
    {
        inm = i * nm;
        l = i - 1; if (d[i] != 0.0)
            for (j = 1; j <= l; j += 1)
            {
                g = 0.0;
                for (k = 1; k <= l; k += 1)g = g + a[i][ k] * a[k ][ j];
                for (k = 1; k <= l; k += 1)
                {
                    knm = k * nm; a[k][ j] = a[k][ j] - g * a[k][ i];
                }
            };
        d[i] = a[i][ i];
        a[i][ i] = 1.0;
        for (j = 1; j <= l; j += 1)a[i][ j] = a[j ][ i] = 0.0;
    }
}/*'eop'*/

//Вычиcление coбcтвенных знaчений и coбcтвенных вектopoв матрицы a[1...n][1...n]
// d[1...n] Нa выхoде coдеpжит coбcтвенные знaчения этoй мaтpицы в вoзpacтaющем  пopядке. 
// e[1...n] мaccив paзмеpнocти N, пocледние N-1 элементoв кoтopoгo нa  вхoде coдеpжaт внедиaгoнaльные элементы  cимметpичеcкoй тpехдиaгoнaльнoй мaтpицы, величинa  Е(1) пpoизвoльнa.
void imtql2(int n, double macheps, double *a[], double d[], double e[])
{
    static int i, i1, j, k, l, m, nm, i1nm, jnm;
    static double b, c, f, g, p, r, s;
    nm = n;
    for (i = 2; i <= n; i += 1)e[i - 1] = e[i];
    e[n] = 0.0; k = n - 1;
    for (l = 1; l <= n; l += 1)
    {
        j = 0;
tst:            for (m = l; m <= k; m += 1)
            if (fabs(e[m]) <= macheps * (fabs(d[m]) + fabs(d[m + 1])))
                goto cont1;
        m = n;
cont1:          p = d[l];
        if (m == l)goto root;
        if (j == 30)goto fail;
        j = j + 1; g = (d[l + 1] - p) / (2 * e[l]); r = sqrt(1 + g * g);
        if (g < 0.0)      g = d[m] - p - e[l] / (g - r) ;
        else           g = d[m] - p - e[l] / (g + r) ;
        s = c = 1.0; p = d[m];
        for (i = m - 1; i >= l; i += -1)
        {
            f = s * e[i]; b = c * e[i];
            if (fabs(f) >= fabs(g))
            {
                c = g / f; r = sqrt(c * c + 1.0);
                e[i + 1] = f * r; s = 1.0 / r; c = c / r;
            }
            else
            {
                c = f / g; r = sqrt(c * c + 1.0);
                e[i + 1] = g * r; s = c / r; c = 1.0 / r;
            };
            f = c * d[i] - s * b; g = c * b - s * p; r = d[i] + p; p = c * f - s * g;
            g = s * f + c * g; d[i + 1] = r - p;
            for (i1 = 1; i1 <= n; i1 += 1)
            {
                i1nm = i1 * nm;
                f = a[i1][ i + 1]; a[i1][ i + 1] = s * a[i1][ i] + c * f;
                a[i1][ i] = c * a[i1][ i] - s * f;
            }
        };
        d[l] = p; e[l] = g; e[m] = 0.0;
        goto tst;
root :;
    }
    for (i = 1; i <= n; i += 1)
    {
        k = i; p = d[i];
        for (j = i + 1; j <= n; j += 1)
            if (d[j] < p)
            {
                k = j; p = d[j];
            };
        if (k != i)
        {
            d[k] = d[i]; d[i] = p;
            for (j = 1; j <= n; j += 1)
            {
                jnm = nm * j;
                p = a[j][ i]; a[j][ i] = a[j][ k]; a[j][ k] = p;
            }
        }
    };
fail :;
}/*'eop'*/



double nRef(double y)
{
    if (y > 0)return ncover;
    else
    {
        if (Mode == 0)
        {
            if (y < -hw)return nsubstrat;
            else return nmax;
        }
        if (Mode == 1)
            return (nmax - nsubstrat) * exp(y / hw)
                   + nsubstrat;
        if (Mode == 2)
            return (nmax - nsubstrat) * exp(-sqr(y / hw))
                   + nsubstrat;
        if (Mode == 3)
        {
            if (y < -hw)return nsubstrat;
            else return nmax - y * y * (nmax - nsubstrat) / (hw * hw);
        }
    }
}/*nRef*/

//элементарная гармоника
double fi_p(int k, double y) {
    double nr;
    nr = 1 / sqrt((ymax - ymin) / 2);
    return nr * sin((k + 1) * M_PI * (y - ymin) / (ymax - ymin));
}/*fi*/

//нормированная гармоника
double fid2_p(int k, double y) {
    double nr;
    nr = 1 / sqrt((ymax - ymin) / 2);
    return -nr * sqr(   (k + 1) * M_PI / (ymax - ymin)  ) *
           sin((k + 1) * M_PI * (y - ymin) / (ymax - ymin));
}/*fi*/

double ModeFun (int m, int M, double *HE[], double x)
{
    int k; double a;
    a = 0; for (k = 0; k < M; k++) a += HE[k][m] * fi_p(k, x);
    return a;
} /*ModeFun*/



//---------------------------------------------------------------------------
void __fastcall TForm1::PlanarClick(TObject *Sender)
{
    int N, M, i, j, k, l, err;

    double *HE[MaxM], *D, *DLbuf, *n_y, *f, *f1, *yy1;

    double a, lam, klam, nm;

    Plan1 = 0;
    err = 0;
    Canva = PaintBox1->Canvas;
    Mode = Edit1->Text.ToInt();
    hw = Edit4->Text.ToDouble();
    nsubstrat = Edit5->Text.ToDouble();
    xmin = -20; xmax = 20;
    ymin = -20; ymax = 20;
    top = left = 0;
    right = PaintBox1->Width;
    bottom = PaintBox1->Height;

    setfillstyle(1, clWhite);
    bar(left, top, right, bottom);

    XoX = (right - left) / 5;
    YoY = (bottom - top) / 2;
    setcolor(clGray);
    line(left, bottom - YoY, right, bottom - YoY);
    line(XoX, top, XoX, bottom);

    msx = (right - left) / (xmax - xmin);
    msy = (bottom - top) / (ymax - ymin);
    hrx = 1;
    hry = 1;

    for (x = xmin; x < xmax; x += hrx)
    {
        y = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix, iy + 5);
        sprintf(str, "%g", x);
        gprintf(&ix, &iy, str);

    }/*x*/

    for (y = ymin; y < ymax; y += hry)
    {
        x = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix + 5, iy);
        sprintf(str, "%g", y);
        gprintf(&ix, &iy, str);
    }/*y*/
    hx = hrx / 100;
    hy = hry / 100;
    N = 0;
    for (y = ymin; y < ymax; y += hy) {
        x = nRef(y);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        putpixel(ix, iy, clBlack);
        N++;
    } /*y*/

    /*
     for(k=0;k<4;k++)
      for(y=ymin;y<ymax;y+=hy)
     {
      x=fi_p(k,y);;
      ix=x*msx+XoX;
      iy=bottom-(y*msy+YoY);
      putpixel(ix,iy,clBlack);
     }/*y*/
    */

    /*
    i=0;
    for(l=0;l<5;l++)
       for(k=0;k<5;k++){

     a=0;for(y=ymin;y<ymax;y+=hy)a+=fi_p(l,y)*fi_p(k,y);
     a=a*hy;
     sprintf(str,"%d %d %e",l,k,a);
     FastWriteA(str,1+2*i,30,clBlack, clWhite);
     i++;
    }/*k,l*/
    */

    if ((n_y = (double *)calloc( sizeof(double), N + 1 )) == NULL)
    {
        sprintf(str, "No memory for n_y[] ");
        err = 1; goto Mexit ;
    }

    if ((f = (double *)calloc( sizeof(double), N + 1 )) == NULL)
    {
        sprintf(str, "No memory for f[] ");
        err = 1; goto Mexit ;
    }
    if ((yy1 = (double *)calloc( sizeof(double), N + 1 )) == NULL)
    {
        sprintf(str, "No memory for yy1[] ");
        err = 1; goto Mexit ;
    }

    i = 0;
    for (y = ymin; y < ymax; y += hy)
    {
        n_y[i] = nRef(y); yy1[i] = y;
        i++;
    }

    M = Mtot;   lam = 1.0;   klam = 2 * M_PI / lam;
    /*===============================*/
    for (i = 0; i <= M; i++)
    {
        if ((HE[i] = (double *)calloc( sizeof(double), M + 1 )) == NULL)
        {
            sprintf(str, "No memory for HE[] ");
            err = 1; goto Mexit ;
        }
    }/*i*/

    if ((DLbuf = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for DLbuf[] ");
        err = 1; goto Mexit ;
    }

    if ((D = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for D ");
        err = 1; goto Mexit ;
    }

    if ((Cbuf = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for Cbuf ");
        err = 1; goto Mexit ;
    }


    /*==============================================================*/

    for (i = 1; i <= M; i++)
    {
        for (k = 0; k < N; k++)
            f[k] = fi_p(i - 1, yy1[k]);
        for (j = 1; j <= i; j++)
        {
            a = 0; 
            for (k = 0; k < N; k++)  {
                a += f[k] * (fid2_p(j - 1, yy1[k]) +
                             sqr(klam * n_y[k]) * fi_p(j - 1, yy1[k]));
            }
            HE[j][i] = HE[i][j] = a * hy;
        }/*j*/
    }/*i*/
    /*-------*/

    tred2(M, 1.0e-10, HE, D, DLbuf);
    imtql2(M, 1.0e-8, HE, D, DLbuf);

    for (i = 1; i <= M; i++)  {
        D[i - 1] = D[i];
        for (j = 1; j <= M; j++)
        {
            HE[i - 1][j - 1] = HE[i][j];
        }/*j*/
    }/*i*/

    k = 0;
    for (i = M - 1; i >= 0; i--)
    {
        nm = sqrt(fabs(D[i])) / klam;
        if (nm > nmax  || nm < nsubstrat )break;
        sprintf(str, "[%2d]    nm=%f ",  k, nm);
        FastWriteA (str, 4 + (M - i) * 2, 20, clBlack, clWhite);

        ix = nm * msx + XoX;
        iy = bottom - YoY;
        line(ix, iy + 10, ix, iy - 10);

        k++;
    }
    for (i = 0; i < M; i++)Cbuf[i] = HE[i][M - 1];

    dump(HE, M, M, "dump1d.csv");
    for (k = 0; k < 2; k++)
        for (y = ymin; y < ymax; y += hy)
        {
            x = ModeFun (M - 1 - k, M, HE, y);
            ix = x * msx + XoX;
            iy = bottom - (y * msy + YoY);
            if (k == 0)Mark(ix, iy, 1, clRed); else  putpixel(ix, iy, clBlack);

        }/*y*/


    Plan1 = 1;
Mexit:;
    free(n_y);
    free(f);
    free(yy1);

    for (i = 0; i <= M; i++)free(HE[i]);
    free(DLbuf);
    free(D);

}/*Planar*/

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
            sprintf (str, " (TREQLAB) reduc1= %d %5ld s.", I);
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
            sprintf (str, " (TREQLAB) reduc= B[%d,%d] %5ld s.[G<0]", I - 1, J - 1);
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
            sprintf (str, " (TREQLAB) reduc2= %d %5ld s. ", I);
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
            sprintf (str, " (TREQLAB) reduc3= %d %5ld s. ", J);
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
        sprintf (str, " (TREQLAB) tred2.0=%d %5ld s. ", I);
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
            sprintf (str, " (TREQLAB) tred2.1=%d %5ld s. ", I);
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
            sprintf (str, " (TREQLAB) D[%3d]=%e [%5ld s]. ", L, D[L - 1]);
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
            sprintf (str, " (TREQLAB) rebak =%d [%5ld s]. ", J);
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





double nRef_cyl(double  r)
{
    if (r < 0)return 0;
    else return nRef(-r);
}/*nRef_cyl*/

double Bas_fun_cyl(int m, int k, double x)
{
    double r;
    if (m == 0) return cos((k + 0.5) * M_PI * x / Lbas);
    else return       sin((k +  1) * M_PI * x / Lbas);
}/*Bas_fun_cyl*/

double dBas_fun_cyl(int m, int k, double x)
{
    double r;
    if (m == 0) return -(k + 0.5) * M_PI / Lbas * sin((k + 0.5) * M_PI * x / Lbas);
    else return          (k + 1) * M_PI / Lbas * cos((k +  1) * M_PI * x / Lbas);
}/*dBas_fun_cyl*/

double d2Bas_fun_cyl(int m, int k, double x)
{
    double r;
    if (m == 0) return -sqr((k + 0.5) * M_PI / Lbas) * cos((k + 0.5) * M_PI * x / Lbas);
    else       return -sqr((k +  1) * M_PI / Lbas) * sin((k +  1) * M_PI * x / Lbas);
}/*d2Bas_fun_cyl*/

double ModeFun_cyl(int mfi, int m, int M, double *HE[], double x)
{
    int k; double a;

    a = 0; for (k = 0; k < M; k++)a += HE[k][m] * Bas_fun_cyl(mfi, k, x);
    return a;
}/*ModeFun*/



//---------------------------------------------------------------------------
void __fastcall TForm1::CYLINDERClick(TObject *Sender)
{
    int N, M, i, j, k, l, err, m;

    double *HE[MaxM], *S[MaxM],
           *D, *DLbuf, *Ebuf, *n_r, *f, *f1, *rr1;
    double a, aa, lam, klam, nm, ro, fi;

    TColor ColSet[256];
    RGB_Col rgbcol;

    for (i = 0; i < 256; i++)
    {
        rgbcol.RGB[3] = (unsigned char)0;
        rgbcol.RGB[0] = (unsigned char)(i);
        rgbcol.RGB[1] = (unsigned char)(i);
        rgbcol.RGB[2] = (unsigned char)(i);
        ColSet[i] = rgbcol.c;
    }



    err = 0;
    lam = 1;
    klam = 2 * M_PI / lam;
    M = 50;



    err = 0;
    Canva = PaintBox1->Canvas;
    Mode = Edit1->Text.ToInt();
    M_fi = Edit2->Text.ToInt();
    M_r = Edit3->Text.ToInt();

    xmin = -2; xmax = 10;
    ymin = -1; ymax = 3;
    top = left = 0;
    right = PaintBox1->Width;
    bottom = PaintBox1->Height;

    setfillstyle(1, clWhite);
    bar(left, top, right, bottom);

    XoX = (right - left) / 6;
    YoY = (bottom - top) / 4;
    setcolor(clGray);
    line(left, bottom - YoY, right, bottom - YoY);
    line(XoX, top, XoX, bottom);

    msx = (right - left) / (xmax - xmin);
    msy = (bottom - top) / (ymax - ymin);
    hrx = 1;
    hry = 0.5;

    for (x = xmin; x < xmax; x += hrx)
    {
        y = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix, iy + 5);
        sprintf(str, "%g", x);
        gprintf(&ix, &iy, str);

    }/*x*/

    for (y = ymin; y < ymax; y += hry)
    {
        x = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix + 5, iy);
        sprintf(str, "%g", y);
        gprintf(&ix, &iy, str);
    }/*y*/
    hx = hrx / 100;
    hy = hry / 100;
    N = 0;
    for (x = 0; x <= xmax; x += hx)
    {
        y = nRef_cyl(x);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        putpixel(ix, iy, clBlue);
        N++;
    }/*x*/
    Lbas = xmax;
    m = M_fi;
    /*
    for(k=0;k<3;k++){
      for(x=0;x<=xmax;x+=hx)
          {
          y=Bas_fun_cyl(m,k,x);
          ix=x*msx+XoX;
          iy=bottom-(y*msy+YoY);
           putpixel(ix,iy,clBlue);
      }/*x*/
}/*k*/
*/
/*==============================================================*/

if ((n_r = (double *)calloc( sizeof(double), N + 1 )) == NULL)
{
    sprintf(str, "No memory for n_r[] ");
    err = 1; goto Mexit ;
}
if ((f = (double *)calloc( sizeof(double), N + 1 )) == NULL)
{
    sprintf(str, "No memory for f[] ");
    err = 1; goto Mexit ;
}
if ((rr1 = (double *)calloc( sizeof(double), N + 1 )) == NULL)
{
    sprintf(str, "No memory for rr1[] ");
    err = 1; goto Mexit ;
}
N = 0;
for (x = 0; x <= xmax; x += hx)
{
    y = nRef_cyl(x);
    n_r[N] = y;
    rr1[N] = x;
    N++;
}/*x*/

/*
 i=0;
 for(l=0;l<5;l++)
    for(k=0;k<5;k++){

  a=0;for(x=0;x<xmax;x+=hx)a+=d2Bas_fun_cyl(m,l,x)
                                 *d2Bas_fun_cyl(m,k,x)*x;
  a=a*hx;
  sprintf(str,"%d %d %e",l,k,a);
  FastWriteA(str,1+2*i,30,clBlack, clWhite);
  i++;
 }/*k,l*/

*/
/**/
for (i = 0; i <= M; i++)
{
    if ((HE[i] = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for HE[] ");
        err = 1; goto Mexit ;
    }
    if ((S[i] = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for S[] ");
        err = 1; goto Mexit ;
    }
}/*i*/

if ((DLbuf = (double *)calloc( sizeof(double), M + 1 )) == NULL)
{
    sprintf(str, "No memory for DLbuf[] ");
    err = 1; goto Mexit ;
}
if ((Ebuf = (double *)calloc( sizeof(double), M + 1 )) == NULL)
{
    sprintf(str, "No memory for Ebuf[] ");
    err = 1; goto Mexit ;
}
if ((D = (double *)calloc( sizeof(double), M + 1 )) == NULL)
{
    sprintf(str, "No memory for D[] ");
    err = 1; goto Mexit ;
}
/*==============================================================*/

for (i = 0; i < M; i++)
{
    for (k = 1; k < N; k++)f[k] = Bas_fun_cyl(m, i, rr1[k]);

    for (j = 0; j <= i; j++)
    {
        HE[i][j] = 0;
        aa = 0; for (k = 1; k < N; k++)
        {
            aa += f[k] *
                  (rr1[k] * d2Bas_fun_cyl(m, j, rr1[k]) +
                   dBas_fun_cyl(m, j, rr1[k]) -
                   m * m / rr1[k] * Bas_fun_cyl(m, j, rr1[k]) +
                   rr1[k] * sqr(klam * n_r[k]) * Bas_fun_cyl(m, j, rr1[k]));

            a = 0;
        }
        HE[j][i] = HE[i][j] = aa * hx;

        aa = 0; for (k = 1; k < N; k++)aa += f[k] *
                                                 Bas_fun_cyl(m, j, rr1[k]) * rr1[k];

        S[j][i] = S[i][j] = aa * hx;
    }/*j*/
}/*i*/


err = TREQLAB(M, HE, S, D, DLbuf, Ebuf);

k = 0;
for (i = M - 1; i >= 0; i--)
{
    nm = sqrt(fabs(D[i])) / klam;
    if (nm > nmax  || nm < nsubstrat )break;
    sprintf(str, "[%2d]    nm=%f ",  k, nm);
    FastWriteA (str, 4 + (M - i) * 2, 20, clBlack, clWhite);

    iy = bottom - (nm * msy + YoY);
    ix = 0 + XoX;
    line(ix, iy, ix + 15, iy);
    k++;
}

for (k = 0; k < 2; k++)
    for (x = 0; x <= xmax; x += hx)
    {
        y = 0.5 * ModeFun_cyl( m, M - 1 - k, M, HE, x);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        if (k == 0)Mark(ix, iy, 1, clBlue); else putpixel(ix, iy, clBlack);

    }/*x*/
k = M_r; m = M_fi; aa = 0;
for (x = 0; x <= xmax; x += hx)
{
    y = ModeFun_cyl( m, M - 1 - k, M, HE, x); y = y * y;
    if (y > aa)aa = y;
}/*x*/

for (x = -hw; x < hw; x += hx)
{
    for (y = -hw; y < hw; y += hx)
    {
        ix = msx * x + XoX;
        iy = bottom - (msx * y + YoY);
        ro = sqrt(x * x + y * y);

        if (x >= 0)fi = asin(y / ro);
        else fi = M_PI - asin(y / ro);
        a = ModeFun_cyl(M_fi, M - 1 - k, M, HE, ro) * cos(M_fi * fi); a = a * a;

        i = a / aa * 256; if (i < 0)i = 0; if (i > 255)i = 255;

        putpixel(ix, iy, ColSet[i]);

    }/*y*/
}/*x*/




Mexit:;

}/*CYLINDER*/

double Pot(double x)
{
    if (Mode == 0)
    {
        if (fabs(x) > hw)return  0;
        else return max_pot;
    }
    if (Mode == 1) return  max_pot * exp(-fabs(x / hw)) ;

    if (Mode == 2) return  max_pot * exp(-sqr(x / hw)) ;

    if (Mode == 3)
        if (fabs(x) > hw)return  0;
        else  return  -max_pot / sqr(hw) * x * x + max_pot;

}/*Pot*/

double fi_w(int k, double x)
{
    double nr;
    nr = 1 / sqrt((xmax - xmin) / 2);
    return nr * sin((k + 1) * M_PI * (x - xmin) / (xmax - xmin));
}/*fi*/

double d2fi_w(int k, double x)
{
    double nr;
    nr = 1 / sqrt((xmax - xmin) / 2);
    return -nr * sin((k + 1) * M_PI * (x - xmin) / (xmax - xmin)) *
           sqr((k + 1) * M_PI / (xmax - xmin));
}/*fi*/

double StateFun (int m, int M, double *HE[], double x)
{
    int k; double a;
    a = 0; for (k = 0; k < M; k++) a += HE[k][m] * fi_w(k, x);
    return a;
} /*ModeFun*/

//---------------------------------------------------------------------------/*Potential_Well*/
//---------------------------------------------------------------------------
void __fastcall TForm1::Planar1Click(TObject *Sender)
{
    int N, M, i, j, k, l, err;

    double *HE[MaxM], *D, *DLbuf, *n_y, *f, *f1, *yy1;

    double a, lam, klam, nm, pp, P;


    err = 0;
    Canva = PaintBox1->Canvas;
    Mode = Edit1->Text.ToInt();
    hw = Edit4->Text.ToDouble();
    nsubstrat = Edit5->Text.ToDouble();
    if (Plan1 == 0)
    {
        sprintf(str, "  Planar !!! ");
        err = -777;  goto Mexit;
    }

    xmin = -1; xmax = 4;
    ymin = -10; ymax = 10;
    top = left = 0;
    right = PaintBox1->Width;
    bottom = PaintBox1->Height;

    setfillstyle(1, clWhite);
    bar(left, top, right, bottom);

    XoX = (right - left) / 5;
    YoY = (bottom - top) / 2;
    setcolor(clGray);
    line(left, bottom - YoY, right, bottom - YoY);
    line(XoX, top, XoX, bottom);

    msx = (right - left) / (xmax - xmin);
    msy = (bottom - top) / (ymax - ymin);
    hrx = 1;
    hry = 1;

    for (x = xmin; x < xmax; x += hrx)
    {
        y = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix, iy + 5);
        sprintf(str, "%g", x);
        gprintf(&ix, &iy, str);

    }/*x*/

    for (y = ymin; y < ymax; y += hry)
    {
        x = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix + 5, iy);
        sprintf(str, "%g", y);
        gprintf(&ix, &iy, str);
    }/*y*/
    hx = hrx / 100;
    hy = hry / 100;
    N = 0;
    for (y = ymin; y < ymax; y += hy)
    {
        x = nRef(y);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        putpixel(ix, iy, clBlack);
        N++;
    }/*y*/

    /*
     for(k=0;k<4;k++)
      for(y=ymin;y<ymax;y+=hy)
     {
      x=fi_p(k,y);;
      ix=x*msx+XoX;
      iy=bottom-(y*msy+YoY);
      putpixel(ix,iy,clBlack);
     }/*y*/
    */

    /*
    i=0;
    for(l=0;l<5;l++)
       for(k=0;k<5;k++){

     a=0;for(y=ymin;y<ymax;y+=hy)a+=fi_p(l,y)*fi_p(k,y);
     a=a*hy;
     sprintf(str,"%d %d %e",l,k,a);
     FastWriteA(str,1+2*i,30,clBlack, clWhite);
     i++;
    }/*k,l*/
    */

    if ((n_y = (double *)calloc( sizeof(double), N + 1 )) == NULL)
    {
        sprintf(str, "No memory for n_y[] ");
        err = 1; goto Mexit ;
    }

    if ((f = (double *)calloc( sizeof(double), N + 1 )) == NULL)
    {
        sprintf(str, "No memory for f[] ");
        err = 1; goto Mexit ;
    }
    if ((yy1 = (double *)calloc( sizeof(double), N + 1 )) == NULL)
    {
        sprintf(str, "No memory for yy1[] ");
        err = 1; goto Mexit ;
    }

    i = 0;
    for (y = ymin; y < ymax; y += hy)
    {
        n_y[i] = nRef(y); yy1[i] = y;
        i++;
    }

    M = Mtot;   lam = 1.0;   klam = 2 * M_PI / lam;
    /*===============================*/
    for (i = 0; i <= M; i++)
    {
        if ((HE[i] = (double *)calloc( sizeof(double), M + 1 )) == NULL)
        {
            sprintf(str, "No memory for HE[] ");
            err = 1; goto Mexit ;
        }
    }/*i*/

    if ((DLbuf = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for DLbuf[] ");
        err = 1; goto Mexit ;
    }

    if ((D = (double *)calloc( sizeof(double), M + 1 )) == NULL)
    {
        sprintf(str, "No memory for D ");
        err = 1; goto Mexit ;
    }



    /*==============================================================*/

    for (i = 1; i <= M; i++)
    {
        for (k = 0; k < N; k++)f[k] = fi_p(i - 1, yy1[k]);
        for (j = 1; j <= i; j++)
        {
            a = 0; for (k = 0; k < N; k++)
            {
                a += f[k] * (fid2_p(j - 1, yy1[k]) +
                             sqr(klam * n_y[k]) * fi_p(j - 1, yy1[k]));
            }
            HE[j][i] = HE[i][j] = a * hy;
        }/*j*/
    }/*i*/
    /*-------*/

    tred2(M, 1.0e-10, HE, D, DLbuf);
    imtql2(M, 1.0e-8, HE, D, DLbuf);

    for (i = 1; i <= M; i++)
    {
        D[i - 1] = D[i];
        for (j = 1; j <= M; j++)
        {
            HE[i - 1][j - 1] = HE[i][j];
        }/*j*/
    }/*i*/

    k = 0; P = 0;
    for (i = M - 1; i >= 0; i--)
    {
        nm = sqrt(fabs(D[i])) / klam;
        if (nm > nmax  || nm < nsubstrat )break;
        sprintf(str, "[%2d]    nm=%f ",  k, nm);
        pp = 0;
        for (j = 0; j < M; j++)
        {
            pp += (Cbuf[j] * HE[j][i]);
        }

        sprintf(str, "[%3d] pp=%g", i, sqr(pp));
        FastWriteA(str, 4 + 2 * (M - i), 1, clBlack, clYellow);

        P += sqr(pp);
        sprintf(str, "[%2d]    nm=%f ",  k, nm);
        FastWriteA (str, 4 + (M - i) * 2, 20, clBlack, clWhite);
        ix = nm * msx + XoX;
        iy = bottom - YoY;
        line(ix, iy + 10, ix, iy - 10);

        k++;
    }
    sprintf(str, "P=%g", P);
    FastWriteA(str, 1, 1, clBlack, clYellow);

    for (k = 0; k < 2; k++)
        for (y = ymin; y < ymax; y += hy)
        {
            x = ModeFun (M - 1 - k, M, HE, y);
            ix = x * msx + XoX;
            iy = bottom - (y * msy + YoY);
            if (k == 0)Mark(ix, iy, 1, clRed); else  putpixel(ix, iy, clBlack);

        }/*y*/



Mexit:;
    if (err != 0 )FastWriteA(str, 10, 10, clBlack, clRed);
}/*Planar1*/

//функция отрисовки линий уровня для заданного массива
int Contur(int  *Map[N_MAP_STR], int ramxr , int ramyr, int NLevel,
           int XoX, int YoY) {
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
    /*
     for(i=0;i<3;i++){
       if ((buf3[i] = (int *)farcalloc( sizeof(int),ramxr+1 )) == NULL) {
          sprintf(str,"No memory for buf3[%d] ",i);
          err=1;goto Mexit ;
       }
    }
    */
    //XoX=200;YoY=300;
    err = 0;
    ramxl = 1;           /*Adding 13.10.91*/
    ramyl = 1;
    MapMax = MapMin = Map[1][1];
    //поиск максимума показателя преломления и минимума
    //не очень понятно, потому что максимум всегда 32000
    for (j = 0; j < ramyr; j++) {
        for (i = 0; i < ramxr; i++) {
            if (Map[j][i] > MapMax) MapMax = Map[j][i];
            if (Map[j][i] < MapMin) MapMin = Map[j][i];
        }
    }

    //нарисуем черный квадрат
    setcolor(clWhite);
    setfillstyle(0, clBlack);
    rect(ramxl + XoX, ramyl + YoY, ramxr + XoX, ramyr + YoY);

    /* ...........  USER CODE  ............. */
    //шаг изменения показателя преломления
    HLevel = (MapMax - MapMin) / NLevel;
    /*
        DebugMain =1;
            sprintf(str,"MapMax=%d MapMin=%d NLevel=%d ",
            MapMax,MapMin,NLevel);
        FastWriteA (str, 6, 40, clWhite, clRed);
        DebugMain =0;
    */
    levgen0 = HLevel;
    //
    for (levgen = levgen0 + MapMin; levgen < MapMax; levgen += HLevel)  {
        //в ходе увеличения показатель преломления не может стать больше нуля
        //скорее всего это костыль
        if (levgen < 0) break;
        im = ip = 0;
        for (jj = 0; jj < 3; jj++)    buf3[jj] = Map[jj];
        jj = 2;
        i0 = 0; i1 = 1; i2 = 2; stcon = sdef = 0.0;
        for (l = 1; l < ramyr; l++)
        {
            //             if(kbhit())  if ((ch = getch()) == ESC) break;
            if (l > ramyl && l < ramyr )
            {
                for (i = ramxl; i < ramxr - 1; i++)
                {
                    stcon = stcon + 1.0;
                    ip = i + 1; im = i - 1;
                    kk  = (int)buf3[i1][i];
                    //если текущий показатель преломления меньше рассматриваемого уровня
                    if ((kk - levgen)*LtGtLevgen <= 0 ) {
                        //запишем каждую из соседних клеток в отдельную переменную
                        k00 =  (int)buf3[i0][im];
                        k01 =  (int)buf3[i0][i] ;
                        k02 =  (int)buf3[i0][ip];
                        k10 =  (int)buf3[i1][im];
                        k12 =  (int)buf3[i1][ip];
                        k20 =  (int)buf3[i2][im];
                        k21 =  (int)buf3[i2][i] ;
                        k22 =  (int)buf3[i2][ip];
                        //пронормируем все соседние клетки:
                        //0 - там больше уровня
                        //1 - там меньше уровня
                        if ((k00 - levgen)*LtGtLevgen > 0) k00 = 0; else k00 = 1;
                        if ((k01 - levgen)*LtGtLevgen > 0) k01 = 0; else k01 = 1;
                        if ((k02 - levgen)*LtGtLevgen > 0) k02 = 0; else k02 = 1;
                        if ((k10 - levgen)*LtGtLevgen > 0) k10 = 0; else k10 = 1;
                        if ((k12 - levgen)*LtGtLevgen > 0) k12 = 0; else k12 = 1;
                        if ((k20 - levgen)*LtGtLevgen > 0) k20 = 0; else k20 = 1;
                        if ((k21 - levgen)*LtGtLevgen > 0) k21 = 0; else k21 = 1;
                        if ((k22 - levgen)*LtGtLevgen > 0) k22 = 0; else k22 = 1;

                        // равно 1, когда хотя бы одна соседняя клетка меньше уровня
                        k33 = k00 * k01 * k02 * k10 * k12 * k20 * k21 * k22;
                        // общее число клеток меньше уровня
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
                            if (kk >= MapMin) ColorForConture = clWhite;
                            else    //return -1;
                            {
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
    // for(i=0;i<ramyr;i++)farfree(Map[i]);
    // for(i=0;i<3;i++)free(buf3[i]);
    return err;
} /* contur() */



//функция показателя преломления от двухмерных координат (x,y)
double nRef_2D(double x, double y) {
    double r;
    r = sqrt(sqr(x) + sqr(y));
    if (y > 0)
        return ncover;
    else {
        if (Mode == 0) {
            if (r > hw) return nsubstrat;
            else return nmax;
        }
        if (Mode == 1)
            return (nmax - nsubstrat) * exp(-r / hw)
                   + nsubstrat;
        if (Mode == 2)
            return (nmax - nsubstrat) * exp(-sqr(r / hw))
                   + nsubstrat;
        if (Mode == 3) {
            if (r > hw) 
                return nsubstrat;
            else 
                return nmax - r * r * (nmax - nsubstrat) / (hw * hw);
        }
    }
}/*nRef*/

// возвращает 0 при k != k1 и какое-то значение, когда они равны
double Mat_elT(int k, int k1) {
    double a, b; int i;
    if (k != k1) {
        return 0;
    }
    else {
        a = b = 0; 
        for (i = 0; i < N_tot; i++) {
            a += Fi_bas[ind[k].i][i] * fid2_p( ind[k].i, yy1[i]);
            b += Fi_bas[ind[k].j][i] * fid2_p( ind[k].j, yy1[i]);
        }
        return (a + b) * h_x * h_y;
    }
}/*Mat_elT*/
double Mat_elN(int k, int k1)
{
    double a; 
    int i, i1, j, j1, m, m1;
    i = ind[k].i; i1 = ind[k1].i;
    j = ind[k].j; j1 = ind[k1].j;
    a = 0; 
    for (m = 0; m < N_tot; m++) {
        for (m1 = 0; m1 < N_tot; m1++) {
            a += Fi_bas[i][m] * Fi_bas[j][m1] * n_xy[m][m1] *
                 Fi_bas[i1][m] * Fi_bas[j1][m1];
        }
    }
    return a * h_x * h_y;
}/*Mat_elN*/
double ModeFun_2d(int m, double *HE[], int i, int j)
{
    double a; int k, i1, j1;
    a = 0; for (k = 0; k < M_bas; k++)
    {
        i1 = ind[k].i; j1 = ind[k].j;
        a += Fi_bas[i1][i] * Fi_bas[j1][j] * HE[k][m];
    }/*k*/
    return a;

}/*ModeFun_2d*/
//---------------------------------------------------------------------------
void __fastcall TForm1::Waveguide_2DClick(TObject *Sender)
{
    int N, M, i, i1, j, j1, k, l, err, mm;
    int iX, iY;
    int *Map[N_MAP_STR];
    double maxDens, yhMap, xhMap;

    double *HE[N_MAP_STR], *D, *DLbuf, *n_y, *f, *f1; //*yy1;

    double a, b, lam, klam, nm;

    Plan1 = 0;
    err = 0;
    Canva = PaintBox1->Canvas;
    Mode = Edit1->Text.ToInt();
    hw = Edit4->Text.ToDouble();
    nsubstrat - Edit5->Text.ToDouble();
    xmin = -20; xmax = 20;
    ymin = -20; ymax = 20;

    xMax = xmax; xMin = xmin;
    yMax = ymax; yMin = ymin;
    top = left = 0;
    right = PaintBox1->Width;
    bottom = PaintBox1->Height;

    setfillstyle(1, clWhite);
    bar(left, top, right, bottom);

    XoX = (right - left) / 2;
    YoY = (bottom - top) / 2;
    setcolor(clGray);
    line(left, bottom - YoY, right, bottom - YoY);
    line(XoX, top, XoX, bottom);

    msx = (right - left) / (xmax - xmin);
    msy = (bottom - top) / (ymax - ymin);
    hrx = 1;
    hry = 1;

    //ось x
    for (x = xmin; x < xmax; x += hrx) {
        y = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix, iy + 5);
        sprintf(str, "%g", x);
        gprintf(&ix, &iy, str);

    }/*x*/

    //ось y
    for (y = ymin; y < ymax; y += hry) {
        x = 0;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        line(ix, iy, ix + 5, iy);
        sprintf(str, "%g", y);
        gprintf(&ix, &iy, str);
    }/*y*/

    h_x = hx = hrx / 10;
    h_y = hy = hry / 10;
    N = 0;

    //показатель преломления на поверхности
    for (y = ymin; y < ymax; y += hy) {
        x = nRef_2D(0, y);
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        putpixel(ix, iy, clBlack);
        N++;
    }/*y*/
    
    //запишем число проходов в счетчик
    N_tot = N;

    /*прорисовка контуров*/
    //   if(((count_tim/nt)/1)*1==(count_tim/nt)){
    NyMap = 600; 
    NxMap = NyMap * (xMax - xMin) / (yMax - yMin); /*размер области для контура*/
    //выделим память под показатели преломления в каждой точке
    for (i = 0; i <= NyMap + 2; i++) {
        if ((Map[i] = (int *)calloc( sizeof(int), NxMap + 2 )) == NULL)
        {
            sprintf(str, "No memory for Map[%d] ", i);
            err = 1; goto Mexit ;
        }
    }/*i*/

    yhMap = (yMax - yMin) / NyMap;
    xhMap = (xMax - xMin) / NxMap;
    k = 0;

    //максимальный показатель преломления
    maxDens = 0;
    for (i = 0; i < NyMap; i++) {
        //zz=i*zhMap+zMin;
        y = i * yhMap + yMin;
        for (j = 0; j < NxMap; j++) {
            x = j * xhMap + xMin;
            b = nRef_2D(x, y);
            if (b > maxDens) maxDens = b;
        }/*j*/
    }/*i*/

    //заполним Map нормированными по максимуму показателями преломления. Макс значение = 32000
    for (i = 0; i < NyMap; i++) {
        //zz=i*zhMap+zMin;
        y = i * yhMap + yMin;
        for (j = 0; j < NxMap; j++) {
            x = j * xhMap + xMin;
            b = nRef_2D(x, y);
            b = (b / maxDens) * 32000;
            Map[(NyMap - 1) - i][j] = (int)b;
        }/*j*/
    }/*i*/

    //функция отрисовки контура показателя преломления
    iX = 0; iY = 0;
    Contur(Map, NxMap, NyMap, 10, iX, iY); // цикл по времени

    //напишем число итераций по оси y
    sprintf(str, "[N=%2d] ",  N);
    FastWriteA (str, 4, 20, clBlack, clWhite);


    //Fi_bas[MaxM];

    // return;
    /*
     for(k=0;k<4;k++)
      for(y=ymin;y<ymax;y+=hy)
     {
      x=fi_p(k,y);;
      ix=x*msx+XoX;
      iy=bottom-(y*msy+YoY);
      putpixel(ix,iy,clBlack);
     }/*y*/
    */

    /*
    i=0;
    for(l=0;l<5;l++)
       for(k=0;k<5;k++){

     a=0;for(y=ymin;y<ymax;y+=hy)a+=fi_p(l,y)*fi_p(k,y);
     a=a*hy;
     sprintf(str,"%d %d %e",l,k,a);
     FastWriteA(str,1+2*i,30,clBlack, clWhite);
     i++;
    }/*k,l*/
    */

    //раздача памяти - массив из N+1 эл-тов'
    if ((n_y = (double *)calloc( sizeof(double), N + 1 )) == NULL) {
        sprintf(str, "No memory for n_y[] ");
        err = 1; goto Mexit ;
    }

    if ((f = (double *)calloc( sizeof(double), N + 1 )) == NULL) {
        sprintf(str, "No memory for f[] ");
        err = 1; goto Mexit ;
    }
    if ((yy1 = (double *)calloc( sizeof(double), N + 1 )) == NULL) {
        sprintf(str, "No memory for yy1[] ");
        err = 1; goto Mexit ;
    }
    // массив (N+1)x(N+1)
    // n_xy - матрица квадратных корней из показателя преломления в этой точке
    for (i = 0; i <= N; i++) {
        if ((n_xy[i] = (double *)calloc( sizeof(double), N + 1 )) == NULL) {
            sprintf(str, "No memory for n_y[] ");
            err = 1; goto Mexit ;
        }
    }/*i*/

    //заполним yy1 и n_xy
    i = 0;
    for (y = ymin; y < ymax; y += hy) {
        yy1[i] = y; j = 0;
        for (x = xmin; x < xmax; x += hx) {
            n_xy[j][i] = sqr(nRef_2D(x, y));
            j++;
        }
        i++;
    }
    M = Mtot;   lam = 1.0;   klam = 2 * M_PI / lam;
    M_bas = M * M;
    /*===============================*/
    //раздача памяти под HE и Fi_bas - двухмерные массивы (M * M)x(M * M)
    for (i = 0; i <= M_bas; i++) {
        if ((HE[i] = (double *)calloc( sizeof(double), M_bas + 1 )) == NULL) {
            sprintf(str, "No memory for HE[] ");
            err = 1; goto Mexit ;
        }
        if ((Fi_bas[i] = (double *)calloc( sizeof(double), M_bas + 1 )) == NULL) {
            sprintf(str, "No memory for Fi_bas[] ");
            err = 1; goto Mexit ;
        }
    }/*i*/

    // DLbuf - одномерный массив длиной (M * M)
    if ((DLbuf = (double *)calloc( sizeof(double), M_bas + 1 )) == NULL) {
        sprintf(str, "No memory for DLbuf[] ");
        err = 1; goto Mexit ;
    }
    // D - одномерный массив длиной (M * M)
    if ((D = (double *)calloc( sizeof(double), M_bas + 1 )) == NULL)
    {
        sprintf(str, "No memory for D ");
        err = 1; goto Mexit ;
    }
    // Cbuf - одномерный массив длиной (M * M)
    if ((Cbuf = (double *)calloc( sizeof(double), M_bas + 1 )) == NULL)
    {
        sprintf(str, "No memory for Cbuf ");
        err = 1; goto Mexit ;
    }
    // ind - одномерный массив длиной (M * M) из струкуры с полями i,j
    if (( ind = (IJ *)calloc( sizeof(IJ), M_bas + 1 )) == NULL)
    {
        sprintf(str, "No memory for D ");
        err = 1; goto Mexit ;
    }

    /*==============================================================*/
    //Fi_bas - кэш значений fi_p
    for (k = 0; k < M; k++) {
        for (i = 0; i < N; i++) {
            Fi_bas[k][i] = fi_p(k - 1, yy1[i]);
        }
    }

    k = 0;
    //ind = [{i:0, j:0}, {i:0, j:1}, .. , {i:0, j:M-1}, {i:1, j:0}, .. , {i:M-1, j:M-1}]
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            ind[k].i = i;
            ind[k].j = j;
            k++;
        }
    }
    //выведем размерность самой большой матрицы 
    sprintf(str, "M_bas= %d  k= %d ", M_bas, k);
    FastWriteA(str, 1, 30, clBlack, clWhite);

    /*==============================================================*/
    //  return;

    for (i = 1; i <= M_bas; i++)
    {
        sprintf(str, "i= %d ", i);
        FastWriteA(str, 3, 30, clBlack, clWhite);

        //заполним матрицу HE - симметрично относительно диагонали
        for (j = 1; j <= i; j++) {
            HE[j][i] = HE[i][j] = Mat_elT(i - 1, j - 1) + sqr(klam) * Mat_elN(i - 1, j - 1);
        }/*j*/
    }/*i*/
    /*-------*/

    tred2(M_bas, 1.0e-10, HE, D, DLbuf);
    imtql2(M_bas, 1.0e-8, HE, D, DLbuf);
    //D - массив собственных чисел HE

    //отбросим первую строку и колонку
    for (i = 1; i <= M_bas; i++) {
        D[i - 1] = D[i];
        for (j = 1; j <= M_bas; j++) {
            HE[i - 1][j - 1] = HE[i][j];
        }/*j*/
    }/*i*/

    k = 0;
    setcolor(clBlack);
    for (i = M_bas - 1; i >= 0; i--) {
        nm = sqrt(fabs(D[i])) / klam;
        if (nm > nmax  || nm < nsubstrat) break;
        sprintf(str, "[%2d]    nm=%f ",  k, nm);
        FastWriteA (str, 4 + (M_bas - i) * 2, 20, clBlack, clWhite);

        ix = nm * msx + XoX;
        iy = bottom - YoY;
        line(ix, iy + 10, ix, iy - 10);
        
        k++;
    }

    dump(HE, M_bas, M_bas, "dump2d.csv");

    //найдем максимальное значение модовой функции
    maxDens = 0; mm = 0;
    for (i = 0; i < NyMap; i++) {
        i1 = i * yhMap / hy;
        for (j = 0; j < NxMap; j++) {
            j1 = j * xhMap / hx;
            b = ModeFun_2d(mm, HE, j1, i1); 
            b = b * b;
            if (b > maxDens) maxDens = b;
        }/*j*/
    }/*i*/

    //построим карту уровней значений функции
    k = 0;
    for (i = 0; i < NyMap; i++) {
        i1 = (i) * yhMap / hy;
        for (j = 0; j < NxMap; j++) {
            j1 = (j) * xhMap / hx;
            b = ModeFun_2d(mm, HE, j1, i1); b = b * b;
            b = (b / maxDens) * 32000;
            Map[(NyMap - 1) - i][j] = (int)b;
        }/*j*/
    }/*i*/
    iX = 0; iY = 0;
    //======================
    // нарисуем этот контур
    Contur(Map, NxMap, NyMap, 10, iX, iY); // цикл по времени

    //нарисуем сечение функции в середине
    for (i = 0; i < N_tot; i++) {
        x = xmin + hx * i;
        y = ModeFun_2d(mm, HE, N_tot / 2, i) * 100;
        ix = x * msx + XoX;
        iy = bottom - (y * msy + YoY);
        Mark(ix, iy, 1, clBlue);
    }/*i*/
    /*==============================================*/
    sprintf(str, "maxDens=%e ", maxDens);
    FastWriteA (str, 6, 40, clBlack, clWhite);

    return;
    //никогда не исполняемый код

    for (i = 0; i < M; i++)Cbuf[i] = HE[i][M - 1];

    for (k = 0; k < 2; k++)
        for (y = ymin; y < ymax; y += hy)
        {
            x = ModeFun (M - 1 - k, M, HE, y);
            ix = x * msx + XoX;
            iy = bottom - (y * msy + YoY);
            if (k == 0)Mark(ix, iy, 1, clRed); else  putpixel(ix, iy, clBlack);

        }/*y*/


    Plan1 = 1;

    //освобождение памяти
    Mexit:;
    free(n_y);
    free(f);
    free(yy1);

    for (i = 0; i <= M; i++) {
        free(HE[i]);    
    }
    free(DLbuf);
    free(D);

}/*Waveguide_2D*/
//---------------------------------------------------------------------------


