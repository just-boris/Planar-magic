#include <math.h>

typedef struct{
 double p1, p2;
}PP;

typedef struct{
 double p1, p2, p3;
}PPP;

#define sqr(x) ((x)*(x))

typedef struct {
               char def;/*S-sourse, L- lens, D-diaphragm, R-screen*/
               double foc,posx,diam,apert;
                  /* foc  - ������. ����������, ��
                     posx - ������� �� ���. ���, ��
                     diam - �������, ��
                     apert- ��������, �� */
               } OPT_OBJ;

typedef struct {
double x,y,z;
} XYZ;