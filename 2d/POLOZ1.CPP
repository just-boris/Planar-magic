//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Poloz1.h"
#include "mf.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------
//блок глобальных переменных
TCanvas *Canva=0;
//TPaintBox *PaintBox1;
char str[200];
int Hfont = 8;
int ix,iy, XoX,YoY,left,right,top,bottom;
double xmin,xmax,ymin,ymax,msx,msy,x,y,hrx,hry;
double sq2pi;

double al,msz,zmin,zmax,z;


double  AxonXx,AxonXy,AxonXz,	AxonYx,AxonYy,AxonYz;
double eal,ebt,egm;
double EuM[3][3];
double EulerAl=80,EulerBt=70,EulerGm=-10;/*углы Эйлера в градусах*/
double delAl=0,delBt=0,delGm=0;/* добавка к углам Эйлера в градусах*/
//double EulerAl=0,EulerBt=0,EulerGm=0;/*углы Эйлера в градусах*/
//=================Planar=================
int Cyl=0; //Признак цилиндрического волновода
   TColor ColSet[256];
   RGB_Col rgbcol;
double Lbas;
//Вывод точки цвета cl
void putpixel(int x, int y, TColor cl)
{
  if (Canva) {
    Canva->Brush->Color = cl;
    Canva->Brush->Style = bsSolid;
    Canva->FillRect(Rect(x,y,x+1,y+1));
  }
};

//Установка цвета линии
void setcolor(TColor x)
{
  if (Canva) 	Canva->Pen->Color = x;
};

//Прорисовка линии из точки (x,y) в точку (x1,y1)
void line(int x, int y, int x1, int y1)
{
  if (Canva) {
    Canva->MoveTo(x,y);
    Canva->LineTo(x1,y1);
  }
};

//Вывод текста начиная с графических координат (x,y)
void gprintf(int *x, int *y, char *s)
{
   int x1,y1;
   x1=*x; y1=*y;
   if (Canva) {
     Canva->Font->Color = clBlack;
     SetBkMode(Canva->Handle, TRANSPARENT);
     Canva->TextOut(x1,y1,s);
   }

};
void gprintfColor(int *x, int *y, char *s,TColor cl)
{
   int x1,y1;
   x1=*x; y1=*y;
   if (Canva) {
     Canva->Font->Color = cl;
     SetBkMode(Canva->Handle, TRANSPARENT);
     Canva->TextOut(x1,y1,s);
   }
};

//Установка цвета заливки
void setfillstyle(int x, TColor cl)
{
  if (Canva)    Canva->Brush->Color = cl;
};

//Прорисовка и заливка прямоугольника с  диагональю (x,y,x1,y1)
void bar(int x, int y, int x1, int y1)
{
  if (Canva)    Canva->FillRect(Rect(x,y,x1,y1));
};
void clearviewport(void)
{
 setfillstyle(1,clWhite);
 bar(left,top,right,bottom);
};
void rect(int x, int y, int x1, int y1)
{
  if (Canva)    Canva->Rectangle(x,y,x1,y1);
};

//Вывод текста из строки s (цвет clLet) с предварительной очисткой фона
// (цвет clFon) в строчные (не графические) позиции (x,y). При этом должна
//быть описана глобальная переменная Hfont и ей присвоено необходимое
//значение высоты букв в пикселях (удобно положить Hfont = 8)
void FastWriteA(char *s, int y, int x, TColor clLet, TColor clFon)
{
  if (Canva) {
    Canva->Font->Color = clLet;
    SetBkColor(Canva->Handle, clFon);
    SetBkMode(Canva->Handle, OPAQUE);
    Canva->TextOut((x)*(Hfont),(y)*(Hfont+4),s);
  }
};
void FastWriteG(char *s, int y, int x, TColor clLet, TColor clFon)
{
  if (Canva) {
    Canva->Font->Color = clLet;
    SetBkColor(Canva->Handle, clFon);
    SetBkMode(Canva->Handle, OPAQUE);
    Canva->TextOut((x),(y),s);
  }
};

//вывод маркерных квадратов размером 2*h+1 пикселей
void Mark(int i,int j,int h,TColor col)
{  int k,l;
      for(k=-h;k<=h;k++)
      for(l=-h;l<=h;l++)
      putpixel(i+k,j+l,col);
}/*Mark*/
// функция определения размеров графической панели
void PanelSize (int *right , int *bottom)
{
    Canva=Form1->PaintBox1->Canvas;
    *right = Form1->PaintBox1->Width; // из атрибутов панели
    *bottom =Form1->PaintBox1->Height;// из атрибутов панели
}/*PanelSize*/


//------------------------------
// описание функции для разметки осей
void Axes(void)
{
//проводим оси координат
    line(left,bottom-YoY,right,bottom-YoY);
    line(XoX,top,XoX,bottom);
		// размечаем ось x с шагом hrx
    for(x=xmin;x<=xmax;x+=hrx){
        y=0;
        //блок перевода физических координат в экранные
        ix=msx*x+XoX;
        iy=bottom-(msy*y+YoY);
        // отмечаем позиции на оси x
        line(ix,iy,ix,iy+5);
        //выводим значения отмеченных позиций
        sprintf(str,"%4.1f",x);
        gprintf(&ix,&iy,str);
    }/*x*/

	// размечаем ось y с шагом hry
    for(y=ymin;y<=ymax;y+=hry){
        x=0;
        ix=msx*x+XoX;
        iy=bottom-(msy*y+YoY);
        line(ix,iy,ix+5,iy);
        sprintf(str,"%4.1f",y);
        gprintf(&ix,&iy,str);
    }/*x*/

}/*Axes*/
//==============================================
double FF(double x, double y)
{

return FunPot_cz(x,y)-nref_substrat;

}

double F(double x, double y)
{ double sx,sy;
    sx = 4; sy = 3;
   return  6*exp(-(sqr(x/sx)+sqr(y/sy)));

}/*F(x,y)*/

void Axon (XYZ p, XYZ *p1) /*перевод в плоскость экрана*/
{   XYZ a;
   a.x=AxonXx*p.x+AxonXy*p.y+AxonXz*p.z;
   a.y=AxonYx*p.x+AxonYy*p.y+AxonYz*p.z;
   a.z=0.0;
   *p1=a;
}/*Axon*/
//заполнение матрицы Эйлера
void EulerMat(double al,double bt, double gm,double A[][3])
{   static double sa,sb,sg,ca,cb,cg;
    sa=sin(al);sb=sin(bt);sg=sin(gm);
    ca=cos(al);cb=cos(bt);cg=cos(gm);
    A[0][0] =(ca*cb*cg-sa*sg);A[0][1]=-(sa*cb*cg+ca*sg);A[0][2]=sb*cg;
    A[1][0] =(ca*cb*sg+sa*cg);A[1][1]=(-sa*cb*sg+ca*cg);A[1][2]=sb*sg;
    A[2][0] =    -ca*sb      ;A[2][1]=     sa*sb       ;A[2][2]=  cb;

}/*EulerMat*/;


void Ax(double A[][3],XYZ p, XYZ *p1)/*поворот в сотв. с углами Эйлера*/
{ static  XYZ a;
     a.x= A[0][0]*p.x+A[0][1]*p.y+A[0][2]*p.z;
     a.y= A[1][0]*p.x+A[1][1]*p.y+A[1][2]*p.z;
     a.z= A[2][0]*p.x+A[2][1]*p.y+A[2][2]*p.z;
   *p1=a;
}/*Ax*/
// Вычисление экранных координат(ix,iy) по заданным физическим(x,y,z)
// в аксонометрической проекции
void xyz_to_ixiyA(double x,double y,double z,int *ix, int *iy)
{ XYZ p,p1;
  p.x=x*msx;p.y=y*msy;p.z=z*msz;
  Axon (p, &p1);
  *ix=p1.x+XoX;
  *iy=bottom-(p1.y+YoY);
}/*xyz_to_ixiy*/


//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
 DecimalSeparator = '.';
 sprintf(str,"%d", Mode);    Edit1->Text=str;
 sprintf(str,"%g", nref_substrat);    Substrat->Text=str;
 sprintf(str,"%g", nref_max);Nref_max->Text=str;
 sprintf(str,"%g", h_waveguide);H_waveguide->Text=str;
 sprintf(str,"%g", w_waveguide);W_waveguide->Text=str;
 sprintf(str,"%g", lam);Lam->Text=str;
}
//---------------------------------------------------------------------------
void Axon_3d(void)
{
     // описание переменных
   int ix1,iy1;
   int Nx,Ny;
   double hx,hy;
   XYZ p,p1;
   Mode=Form1->Edit1->Text.ToInt();
       xmin=-5;
       xmax=5;
       ymin=-5;
       ymax=5;
       zmin=-1;
       zmax=5;
       al=30.0/180*M_PI; //
     top=left=0; PanelSize (&right , &bottom);//уст. размер граф.панели
    msx=(right-left)/(xmax-xmin)*0.6;
    msy=(bottom-top)/(ymax-ymin);
    msz=(bottom-top)/(zmax-zmin);;
    //================
/*
    nref_max=1.9;
    nref_substrat=1.2;
    h_waveguide=5;
    w_waveguide=3;
*/
    Mode=Form1->Edit1->Text.ToInt();
    nref_substrat=Form1->Substrat->Text.ToDouble();
    nref_max=     Form1->Nref_max->Text.ToDouble();
    h_waveguide=     Form1->H_waveguide->Text.ToDouble();
    w_waveguide=     Form1->W_waveguide->Text.ToDouble();



    al_nref=1/h_waveguide;
     //====================

    //устанавливаем цвет заливки и заливаем панель
    setfillstyle(1,clWhite);
	bar(left,top,right,bottom);

    XoX=(left+right)/2;       // начало координат в
    YoY=(top+bottom)/2;       // в центре экрана

     eal =  M_PI*  (EulerAl+delAl)/180;
     ebt =  M_PI*  (EulerBt+delBt)/180;
     egm =  M_PI*  (EulerGm+delGm)/180;
     EulerMat(eal,ebt, egm, EuM);/*заполнение матрицы Эйлера для зад. угл.*/
          /*определение проекций на экран*/
     p.x=1.0;p.y=0.0;p.z=0.0;
     Ax(EuM, p, &p1);
       AxonXx=p1.x;AxonXy=p1.y;AxonXz=p1.z;
     p.x=0.0;p.y=1.0;p.z=0.0;
     Ax(EuM, p, &p1);
       AxonYx=p1.x;AxonYy=p1.y;AxonYz=p1.z;





        setcolor(clGray);
        line(left,bottom-YoY,right,bottom-YoY);
        //пририсовка 3D осей
        setcolor(clBlack);
        xyz_to_ixiyA(0,0,0,&ix, &iy);
        xyz_to_ixiyA(0,0,8,&ix1, &iy1);
        line(ix,iy,ix1,iy1); gprintf(&ix1,&iy1," F(x,y)");
        xyz_to_ixiyA(0,8,0,&ix1, &iy1);
        line(ix,iy,ix1,iy1);  iy1-=20; gprintf(&ix1,&iy1,"     Y");
        xyz_to_ixiyA(8,0,0,&ix1, &iy1);
        line(ix,iy,ix1,iy1);  ix1-=30; gprintf(&ix1,&iy1," X");

        //Построение поверхности F(x,y)
         Nx=30;
         Ny=30;
         hx=(xmax-xmin)/Nx;
         hy=(ymax-ymin)/Ny;
      for(x=xmin;x<xmax;x+=hx){
         for(y=ymin;y<ymax;y+=hy){
            z=FF(x,y);   xyz_to_ixiyA(x,y,z,&ix, &iy);
            z=FF(x+hx,y);xyz_to_ixiyA(x+hx,y,z,&ix1, &iy1);
            line(ix,iy,ix1,iy1);
            z=FF(x,y+hy);xyz_to_ixiyA(x,y+hy,z,&ix1, &iy1);
            line(ix,iy,ix1,iy1);
         }/*y*/
      }/*x*/

        //=======================================
} /* Axon_3d*/
double erf(double x)
{  double xx,a,hx,t;
     xx=fabs(x);
     hx=xx/500;
 if(xx< 1.0e-10)a=0;
 else {
   a=0;for(t=0;t<=xx;t+=hx)a+=exp(-sqr(t));a*=hx;
   if(x<0)a=-a;
 }
 return a*sq2pi;
}/*erf*/

double Fun(double x)
{  double W,dz,z;
 W=5;//mkm
 dz=7.22;//mkm
 z=x;
 return 0.5*(erf((1+2*z/W)*W/(2*dz))+erf((1-2*z/W)*W/(2*dz)));


  //return sin(x);



}/**/
void graph_1d(void)
{
     double hx;
       xmin=-15;
       xmax=15;
       ymin=-2;
       ymax=2;
     top=left=0; PanelSize (&right , &bottom);//уст. размер граф.панели
    //устанавливаем цвет заливки и заливаем панель
    setfillstyle(1,clWhite);
	bar(left,top,right,bottom);

    XoX=(left+right)/2;       // начало координат в
    YoY=(top+bottom)/2;       // в центре экрана
	// вычисляем масштабные коэффициенты
    msx=(right-left)/(xmax-xmin);
    msy=(bottom-top)/(ymax-ymin);
	// устанавливаем шаг разметки осей
    hrx=1.0;
    hry=0.5;
	// размечаем ось x с шагом hrx
   Axes();
//==============блок вывода графиков==================
   sq2pi=2/(sqrt(M_PI));
	// выводим график функции точками с шагом hx
	hx=hrx/100;
    for(x=xmin;x<=xmax;x+=hx){
        y=Fun(x);
        // y=erf(x);
        //блок перевода физических координат в экранные
        ix=msx*x+XoX;
        iy=bottom-(msy*y+YoY);
        //выводим значения функции по точкам
        putpixel(ix,iy,clBlue);
    }/*x*/
//=======================================

}/*graph_1d()*/
//================================
double Bas_fun_cyl(int m,int k, double x)
{ double r,ll;
if(m == 0) return cos((k+0.5)*M_PI*x/Lbas);
else return sin((k+1)*M_PI*x/Lbas);
}/*Bas_fun_cyl*/
double dBas_fun_cyl(int m,int k, double x)
{ double r,ll;
if(m == 0) return -(k+0.5)*M_PI/Lbas*sin((k+0.5)*M_PI*x/Lbas);
else return (k+0.5)*M_PI/Lbas*cos((k+1)*M_PI*x/Lbas);
}/*Bas_fun_cyl*/

double d2Bas_fun_cyl(int m,int k, double x)
{ double r,ll;
if(m == 0) return -sqr((k+0.5)*M_PI/Lbas)*cos((k+0.5)*M_PI*x/Lbas);
else       return -sqr((k+1)*M_PI/Lbas)*sin((k+1)*M_PI*x/Lbas);
}/*Bas_fun_cyl*/

double ModeFun_cyl(int mfi,int m,int M, double *HE[],double x)
{   int k; double a;

    a=0;for(k=0;k<M;k++)a+=HE[k][m]*Bas_fun_cyl(mfi,k, x);
    return a;
}/*ModeFun*/




void cylindr(void)
{

 int ix,iy,ix1,iy1;
 int i,k,j,l,N,M,m,err;
 double hx,L,aa,a,t,ht,hy,lam,klam,nm,r,fi;
 double *HE[MaxM],*S[MaxM],*D,*DLbuf,*Ebuf,*n_y,*f,*f1,*yy1;
  err=0;
       for(i=0;i<256;i++){
       rgbcol.RGB[3]=(unsigned char)0;
       rgbcol.RGB[0]=(unsigned char)(i);
       rgbcol.RGB[1]=(unsigned char)(i);
       rgbcol.RGB[2]=(unsigned char)(i);
       ColSet[i]=rgbcol.c;
    }
 tim0 =  time(&timt);
 top=left=0; PanelSize (&right , &bottom);//уст. размер граф.панели
 setfillstyle(1,clWhite);
 bar(left,top,right,bottom);
 XoX=(left+right)/11;
 YoY=(top+bottom)/3;
  setcolor(clGray);
 line(left,bottom-YoY,right, bottom-YoY);
 line(XoX, top, XoX,bottom);
 ymin=-1; ymax=3;
 xmin=-1; xmax=10;
 hrx=1;
 hry=0.1;
 msx=(right-left)/(xmax-xmin);
 msy=(bottom-top)/(ymax-ymin);
    Axes();//проводим и размечаем оси
    Mode=3;
        Mode=Form1->Edit1->Text.ToInt();
    m=0;
//=================================
    Cyl=1;FlagCyl=0;
    nref_max=1.9;
    nref_substrat=1.8;
    h_waveguide=1;
    Lbas=xmax;
    al_nref=1/h_waveguide;
//================================    
    //lam=2;
    lam=     Form1->Lam->Text.ToDouble();
    klam=2*M_PI/lam;
   sprintf(str,"Cylindr: nref_max=%g nref_substrat=%g h_waveguide=%g lam=%g ",
      nref_max,nref_substrat, h_waveguide, lam);
   FastWriteA (str, 1,1, clBlack, clWhite);



   // выводим профиль показателя преломления
  setcolor (clBlack);
    hx=hrx/100;
    N=0;
    for(x=0 ;x< xmax;x+=hx){
     y=Nref_Planar(x);
     ix=x*msx+XoX;
     iy=bottom-(y*msy+YoY);
     if(N==0){ix1=ix;iy1=iy;}
     line(ix,iy,ix1,iy1);
     ix1=ix;iy1=iy;
     N++;
    }/*y*/
    NCyl=N;
    HRCyl=hx;
/*
     m=4;k=0;
    for(k=0;k<10;k++)
    for(x=0 ;x< xmax;x+=hx){
     y=Bas_fun_cyl(m,k,  x);
     ix=x*msx+XoX;
     iy=bottom-(y*msy+YoY);
     if(x==0){ix1=ix;iy1=iy;}
     line(ix,iy,ix1,iy1);
     ix1=ix;iy1=iy;
    }/*y*/
*/
    M =100;
    if(M>MaxM)M=MaxM;
  /*==============================================================*/
   if ((n_y = (double *)calloc( sizeof(double),N+1 )) == NULL) {
      sprintf(str,"No memory for n_z[] ");
         err=1;goto Mexit ;
   }
   if ((f = (double *)calloc( sizeof(double),N+1 )) == NULL) {
      sprintf(str,"No memory for f[] ");
         err=1;goto Mexit ;
   }
/*   */
   if ((yy1 = (double *)calloc( sizeof(double),N+1 )) == NULL) {
      sprintf(str,"No memory for yy1[] ");
         err=1;goto Mexit ;
   }
   if ((FCyl = (double *)calloc( sizeof(double),N+1 )) == NULL) {
      sprintf(str,"No memory for FCyl[] ");
         err=1;goto Mexit ;
   }
   if ((RCyl = (double *)calloc( sizeof(double),N+1 )) == NULL) {
      sprintf(str,"No memory for RCyl[] ");
         err=1;goto Mexit ;
   }

/**/
   for(i=0;i<=M;i++){
     if ((HE[i] = (double *)calloc( sizeof(double),M+1 )) == NULL) {
      sprintf(str,"No memory for HE[] ");
         err=1;goto Mexit ;
     }
     if ((S[i] = (double *)calloc( sizeof(double),M+1 )) == NULL) {
      sprintf(str,"No memory for S[] ");
         err=1;goto Mexit ;
     }
   }/*i*/
   if ((DLbuf = (double *)calloc( sizeof(double),M+1 )) == NULL) {
      sprintf(str,"No memory for DLbuf[] ");
         err=1;goto Mexit ;
   }
      if ((Ebuf = (double *)calloc( sizeof(double),M+1 )) == NULL) {
      sprintf(str,"No memory for Ebuf[] ");
         err=1;goto Mexit ;
   }

   if ((D = (double *)calloc( sizeof(double),M+1 )) == NULL) {
      sprintf(str,"No memory for DLbuf[] ");
         err=1;goto Mexit ;

   }
  /*==============================================================*/
/*==========заполнение профиля n_y ======*/
  i=0;
  for(x=0;x<xmax;x+=hx)
  {
    n_y[i]=Nref_Planar(x);yy1[i]=x;i++;
  }/*y*/
/*==========заполнение матрицы и вычисление собств. знач. начало======*/
//m=0;
  for(i=0;i<M;i++){
    sprintf(str,"N=%4d Mbas=%3d m_fi=%d",N,i,m);
    FastWriteA (str, 2,20, clBlack, clWhite);
      for(k=1;k<N;k++)f[k]=Bas_fun_cyl(m,i,yy1[k]);
     for(j=0;j<=i;j++){
        HE[i][j] =0;
       aa=0; for(k=1;k<N;k++)aa+=f[k]*(yy1[k]*d2Bas_fun_cyl(m,j,yy1[k])+
                                              dBas_fun_cyl(m,j,yy1[k])-
                                            m*m/yy1[k]*Bas_fun_cyl(m,j,yy1[k])+
                             yy1[k]*sqr(klam*n_y[k])*Bas_fun_cyl(m,j,yy1[k]));
        HE[j][i]=HE[i][j]=aa*hx;
       aa=0; for(k=1;k<N;k++)aa+=f[k]*Bas_fun_cyl(m,j,yy1[k])*yy1[k];
               S[j][i]=S[i][j]=aa*hx;

     }/*j*/
  }/*i*/
 /*-----заполнение матрицы энергии конец--*/



    err=TREQLAB(M,HE,S,D,DLbuf,Ebuf);

  k=0;
  for(i=M-1;i>=0;i--){
   nm=sqrt(fabs(D[i]))/klam;
   if(nm<Nref_Planar(ymin) || nm<Nref_Planar(ymax) )break;
   sprintf(str,"[%2d]    Nm=%f ",  k,nm);
   FastWriteA (str, 4+(M-i)*2,20, clBlack, clWhite);
   ix=0+XoX;
   iy=bottom-(msy*nm+YoY);
    setcolor(clBlue);
    line (ix+20,iy,ix-20,iy);
    ix-=30;
    sprintf(str,"%d",k);
    gprintf(&ix,&iy,str);
    k++;
  }
   i=k;
/*-------------------------*/
double Sc;
int NCyl ;
hy=hx;
for(k=1;k<=i;k++ ){
    aa=0;j=0;
    for(x=0 ;x< xmax;x+=hx){
     y=ModeFun_cyl(m, M-k, M,HE,x); //x=sqr(x);
     if(m==0 && k==1){
      FCyl[j]=y; Sc+=x*y*y*hx;
      RCyl[j]=x;
     }
     y=y*y;
     if(fabs(y)>aa)aa=fabs(y);
     ix=x*msx+XoX;
     iy=bottom-(y*msy*0.5+YoY);
     if(k==1)Mark(ix,iy,1,clGray);
      else   putpixel(ix,iy,clBlack);
      j++;
    }/*y*/
    NCyl=j;HRCyl=hx;
    ix=2*h_waveguide*msx;
    hx=2*h_waveguide/ix;
   if(k<5)
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
   }/**/
    if(k>1)break;
 }/*k*/

   sprintf(str,"Sc=%f ",Sc  );
   FastWriteA (str, 10,20, clBlack, clWhite);

/*-----------------------------*/
    Mexit:;
if(err != 0)FastWriteA(str,2,2,clBlack,clRed);

   for(i=0;i<=M;i++){
     free(HE[i]);
     free(S[i]);
   }/*i*/
   free(DLbuf);
   free(Ebuf);
   free(D);
   free(n_y);
   free(f);
   free(yy1);
   Cyl=0;FlagCyl=1;

}/*planar_mfe*/

void __fastcall TForm1::Axon_3DClick(TObject *Sender)
{
Axon_3d();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::ScrollBar1Change(TObject *Sender)
{
delAl=ScrollBar1->Position;
Axon_3d();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::ScrollBar2Change(TObject *Sender)
{
delBt=ScrollBar2->Position;
Axon_3d();

}
//---------------------------------------------------------------------------

void __fastcall TForm1::ScrollBar3Change(TObject *Sender)
{
delGm=ScrollBar3->Position;
Axon_3d();
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Strip_waveguideClick(TObject *Sender)
{
 strip_waveguide();
}
//---------------------------------------------------------------------------
/*
*/
//---------------------------------------------------------------------------

void __fastcall TForm1::Graph_1dClick(TObject *Sender)
{
graph_1d();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::CylindrClick(TObject *Sender)
{
cylindr();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::TestClick(TObject *Sender)
{
Test_splines();        
}
//---------------------------------------------------------------------------


