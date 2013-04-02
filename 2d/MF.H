
// описание необходимых стандартных библиотек
#include <vcl.h>
#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <time.h>
#include <dos.h>
//определение возведения в квадрат
#define  sqr(x)    ((x)*(x))
// Трехкомпонентные переменные
typedef struct{
  double x,y,z;
}XYZ;

//двухкомпонентная структура давлений в объемах
typedef struct{
   double p1,p2;
}PP;
typedef struct{
     char def;
     double posx,
           diam,
           foc,
           apert;
}OPT;
typedef struct{
 double re,im;
}complex;
typedef union {/*общий нерелятивистский трехмерный случай*/
       unsigned long i;
       char c[4];
       struct bit_field{
                unsigned hz   :1;//эрмитова компонента по z
                unsigned hr   :1;//эрмитова компонента по r
                unsigned hc   :1;//эрмитова компонента движения ядер
                unsigned iz   :10;//номер узла по z
                unsigned ir   :8;//номер узла по r
                unsigned ic   :8;//номер узла по c
        } zrc;

}UUU3p;


typedef union {
       unsigned i;
       char c[2];
       struct bit_field{
                unsigned k   :1;
                unsigned h   :1;
                unsigned iz  :14;
        } khz;
}UUU;

 typedef   union{
          TColor c;
          unsigned char RGB[4];
          }RGB_Col;
typedef struct{
   int     day, month, year,  hours,   minutes,    seconds;

}Datime;
typedef union{
  unsigned long i;
     struct bit_field{
          unsigned i   :16;//
          unsigned j   :16;//
     }ij;
}NONZERO;      // расшифровка ненулевых элементов матрицы
          

#define MaxN  1000 //макс. число значений в таблице
#define MaxM  1000  //макс. число базисных функций
#define N_MAP_STR 1000
#define zhch(iz,izh,ic,ich)  (ich+ic*2+izh*NNc*2+iz*2*NNc*2)
extern char str[];
extern int ix,iy, XoX,YoY,left,right,top,bottom;
extern double xmin,xmax,ymin,ymax,msx,msy,x,y,hrx,hry;
extern    long timt,tim0;
extern Hfont;
void Dirac_Coulomb(void);
void strip_waveguide(void);
double FunPot_cz(double c,double z);
//=======================
extern double xmin,xmax;
void putpixel(int x, int y, TColor cl);
//Установка цвета линии
void setcolor(TColor x);
//Прорисовка линии из точки (x,y) в точку (x1,y1)
void line(int x, int y, int x1, int y1);
//Вывод текста начиная с графических координат (x,y)
void gprintf(int *x, int *y, char *s);
//Установка цвета заливки
void setfillstyle(int x, TColor cl);
//Прорисовка и заливка прямоугольника с  диагональю (x,y,x1,y1)
void bar(int x, int y, int x1, int y1);
void rect(int x, int y, int x1, int y1);
void gprintfColor(int *x, int *y, char *s,TColor cl);
void clearviewport(void);
//Вывод текста из строки s (цвет clLet) с предварительной очисткой фона
// (цвет clFon) в строчные (не графические) позиции (x,y). При этом должна
//быть описана глобальная переменная Hfont и ей присвоено необходимое
//значение высоты букв в пикселях (удобно положить Hfont = 8)
void FastWriteA(char *s, int y, int x, TColor clLet, TColor clFon);
void FastWriteG(char *s, int y, int x, TColor clLet, TColor clFon);
//вывод маркерных квадратов размером 2*h+1 пикселей
void Mark(int i,int j,int h,TColor col);
// функция определения размеров графической панели
void PanelSize (int *right , int *bottom);
//------------------------------
// описание функции для разметки осей
void Axes(void);
 int  TREQLAB(int N,double *A[], double *B[],double *D,
                          double *DL,double *E);

extern   double  nref_max,nref_substrat,h_waveguide,w_waveguide,al_nref,lam;
extern  int Mode;
double Nref_Planar(double y1);
void Test_splines(void);
void justir(void);
extern double
       nref_max,
       nref_substrat,
       al_nref,
       h_waveguide,
       w_waveguide,

/* ========Koshiba==
      nref_max=3.44,
       nref_substrat=3.34,
       al_nref,
       h_waveguide=5,
       w_waveguide=5,
===========*/
       lam,
       klam;
       /*Test Koshiba iee proceeding-j. vol 139, No.2, April, 1992*/
       // rib structure
extern  int Cyl,FlagCyl;
extern  double  *RCyl,*FCyl;
extern  double   HRCyl;
extern  int NCyl;
extern     TColor ColSet[256];
extern     RGB_Col rgbcol;


