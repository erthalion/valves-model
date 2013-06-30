#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <csignal>
#include <fstream>
#include <string.h>
#include <iostream>
#include <time.h>

#include "constants.hpp"
#include "utils.hpp"
#include "groups.h"
#include "output.h"

#include "ndarray.h"

using namespace std;

/*
 * Variable type for U array
 */
const int VELOCITY_U = 0;
const int VELOCITY_V = 1;
const int VELOCITY_W = 2;
const int PRESSURE = 3;

// Глобальные переменные
long double
    *Hx, *Hy, *Hz, // шаги
    *Cx, *Cy, *Cz, // координаты
    Rn, R0,
    Rn_1, Rn_2,
    ***U, // Вектор неизвестных
    ***U_1, ***U_2, // для ускорения
    ***R, // Невязка
    ***Z; // для метода н/а

int
    iters,
    ***G; // Маска узлов

// Работа с группами и матрицы влияния
matrix_ind ***arg, ***func;

// Массив содержит количество элементов, которые изменяться (и будут пересчитаны), если изменить i,j,k элемент
int ***global_carg;

Utils *utils;
GroupsGenerator *groupGenerator;
Output *output;

// Методы проекта
long double norm(long double*** v);
long double norm(long double*** v1, long double*** v2);
void residual();
void residual(int i, int j, int k);
long double A1(long double ***U1, long double ***U2, int i, int j, int k);
long double A2(long double ***U, int i, int j, int k);
extern "C" {
    int main();
}
long double A(long double ***U1, long double ***U2, int i, int j, int k);
void eval_scalars(long double ***u, long double ***R, int i1, int j1, int k1, long double &Rn_F1, long double &Rn_F2, long double &F1_F1, long double &F2_F2, long double &F1_F2);
long double calc_alpha(long double Rn_F1, long double Rn_F2, long double F1_F1, long double F2_F2, long double F1_F2);
int cubic(long double *x, long double a, long double b, long double c);
void alpha_iter();


// Норма вектора
long double norm(long double*** v)
{
    long double s = 0;
    for(int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                s += v[i][j][k]*v[i][j][k]*Hx[i]*Hy[j]*Hz[k];
    return sqrt(s);
}

long double norm(long double*** v1, long double*** v2)
{
    long double s = 0, a;
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                if(G[i][j][k])
                {
                    a = v1[i][j][k]-v2[i][j][k];
                    s += a*a*Hx[i]*Hy[j]*Hz[k];
                }
            }
        }
    }
    return sqrt(s);
}

/*
 * Полный пересчет невязки
 */
void residual()
{
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                if(G[i][j][k])
                {
                    R[i][j][k] = A(U,U,i,j,k);
                }
            }
        }
    }
}

/*
 * Пересчет невязки для конкретного элемента.
 * Сначала выбирается количество элементов, которые изменятся,
 * затем индексы для них - по ним происходит пересчет.
 */
void residual(int i, int j, int k)
{
    int count_changed = global_carg[i][j][k];
    for(int oc=0; oc<count_changed; ++oc)
    {
        indexes &changed=arg[i][j][k][oc];
        int ii=changed.i;
        int jj=changed.j;
        int kk=changed.k;
        if(G[ii][jj][kk])
            R[ii][jj][kk] = A(U,U,ii,jj,kk);
    }
}


/////////////////////////////////////////////////////////////////////////////////
//  нелинейная часть оператора
long double A1(long double ***U1, long double ***U2, int i, int j, int k)
{
    int g = G[i][j][k];

    // check mask
    if( g == 0 || k>=dNz*3 ) return 0;

    int uk = k%dNz;
    int vk = uk + dNz;
    int wk = vk + dNz;

    long double ddx, ddy, ddz;
    if(k==uk)
    {
        long double appv, appw;
        switch(g)
        {
            case 2:
                ddx = (U2[i+1][j][uk]-U2[i][j][uk])*2/(Hx[i+1]+Hx[i])*U1[i][j][uk];
                break;
            case 3:
                ddx = (U2[i][j][uk]-U2[i-1][j][uk])*2/(Hx[i]+Hx[i-1])*U1[i][j][uk];
                break;
            default:
                ddx = (U2[i+1][j][uk]-U2[i-1][j][uk])*2/(Hx[i+1]+2*Hx[i]+Hx[i-1])*U1[i][j][uk];
        }

        appv = 0.5*( Hy[j]/(Hy[j]+Hy[j+1])*(U1[i][j+1][vk]+U1[i-1][j+1][vk]) +
                Hy[j+1]/(Hy[j]+Hy[j+1])*(U1[i][j][vk]+U1[i-1][j][vk]) );
        ddy = (U2[i][j+1][uk]-U2[i][j-1][uk])/(Hy[j+1]+Hy[j])*appv;

        appw = 0.5*( Hz[k]/(Hz[k]+Hz[k+1])*(U1[i][j][wk+1]+U1[i-1][j][wk+1]) +
                Hz[k+1]/(Hz[k]+Hz[k+1])*(U1[i][j][wk]+U1[i-1][j][wk]) );
        ddz = (U2[i][j][uk+1]-U2[i][j][uk-1])/(Hz[k+1]+Hz[k])*appw;
    }
    else if(k==vk)
    {
        long double appu, appw;
        switch(g)
        {
            case 2:
                appu = 0.5*( U1[i+1][j][uk] + U1[i+1][j-1][uk] );
                ddx = (U2[i+1][j][vk]-U2[i][j][vk])/Hx[i+1]*appu;
                break;
            case 3:
                appu = 0.5*( U1[i][j][uk] + U1[i][j-1][uk] );
                ddx = (U2[i][j][vk]-U2[i-1][j][vk])/Hx[i]*appu;
                break;
            default:
                appu = 0.5*( Hx[i+1]/(Hx[i]+Hx[i+1])*(U1[i][j][uk]+U1[i][j-1][uk]) +
                        Hx[i]/(Hx[i]+Hx[i+1])*(U1[i+1][j][uk]+U1[i+1][j-1][uk]) );
                ddx = (U2[i+1][j][vk]-U2[i-1][j][vk])/(Hx[i+1]+Hx[i])*appu;
        }

        ddy = (U2[i][j+1][vk]-U2[i][j-1][vk])*2/(Hy[j+1]+2*Hy[j]+Hy[j-1])*U1[i][j][vk];

        appw = 0.5*( Hz[k]/(Hz[k]+Hz[k+1])*(U1[i][j][wk+1]+U1[i][j-1][wk+1]) +
                Hz[k+1]/(Hz[k]+Hz[k+1])*(U1[i][j][wk]+U1[i][j-1][wk]) );
        ddz = (U2[i][j][vk+1]-U2[i][j][vk-1])/(Hz[k+1]+Hz[k])*appw;
    }
    else if(k==wk)
    {
        long double appu, appv;
        switch(g)
        {
            case 2:
                appu = 0.5*( U1[i+1][j][uk] + U1[i+1][j][uk-1] );
                ddx = (U2[i+1][j][wk]-U2[i][j][wk])/Hx[i+1]*appu;
                break;
            case 3:
                appu = 0.5*( U1[i][j][uk] + U1[i][j][uk-1] );
                ddx = (U2[i][j][wk]-U2[i-1][j][wk])/Hx[i]*appu;
                break;
            default:
                appu = 0.5*( Hx[i+1]/(Hx[i]+Hx[i+1])*(U1[i][j][uk]+U1[i][j][uk-1]) +
                        Hx[i]/(Hx[i]+Hx[i+1])*(U1[i+1][j][uk]+U1[i+1][j][uk-1]) );
                ddx = (U2[i+1][j][wk]-U2[i-1][j][wk])/(Hx[i+1]+Hx[i])*appu;
        }

        appv = 0.5*( Hy[j]/(Hy[j]+Hy[j+1])*(U1[i][j+1][vk]+U1[i][j+1][vk-1]) +
                Hy[j+1]/(Hy[j]+Hy[j+1])*(U1[i][j][vk]+U1[i][j][vk-1]) );
        ddy = (U2[i][j+1][wk]-U2[i][j-1][wk])/(Hy[j+1]+Hy[j])*appv;

        ddz = (U2[i][j][wk+1]-U2[i][j][wk-1])*2/(Hz[k+1]+2*Hz[k]+Hz[k-1])*U1[i][j][wk];
    }

    return (ddx+ddy+ddz);
}
// линейная часть оператора
long double A2(long double ***U, int i, int j, int k)
{
    int g=G[i][j][k];

    // check mask
    if(!g) return 0;

    int uk = k%dNz;
    int vk = uk + dNz;
    int wk = vk + dNz;
    int pk = wk + dNz;

    if(k>=3*dNz) //давление
    {
        switch(g)
        {
            case 5:
                return (U[i][j][pk+1]-U[i][j][pk])/Hz[k+1];
            case 6:
                return (U[i][j][pk]-U[i][j][pk-1])/Hz[k];
            case 7:
                return (U[i][j+1][pk]-U[i][j][pk])/Hy[j+1];
            case 8:
                return (U[i][j][pk]-U[i][j-1][pk])/Hy[j];
            default:
                return (U[i+1][j][uk]-U[i][j][uk])*2/(Hx[i+1]+Hx[i]) +
                    (U[i][j+1][vk]-U[i][j][vk])*2/(Hy[j+1]+Hy[j]) +
                    (U[i][j][wk+1]-U[i][j][wk])*2/(Hz[k+1]+Hz[k]);
        }
    }

    long double lap_u;
    int i1=i;
    switch(g)
    {
        case 2:
            i1+=1;
            break;
        case 3:
            i1-=1;
    }
    if(k==uk)
    {
        lap_u=nu*(
                ((U[i1+1][j][k]-U[i1][j][k])*2/(Hx[i1+1]+Hx[i1])-
                 (U[i1][j][k]-U[i1-1][j][k])*2/(Hx[i1]+Hx[i1-1]))*2/(0.5*Hx[i1+1]+Hx[i1]+0.5*Hx[i1-1]) +
                ((U[i][j+1][k]-U[i][j][k])/Hy[j+1]-
                 (U[i][j][k]-U[i][j-1][k])/Hy[j])/(Hy[j+1]+Hy[j])*2.0 +
                ((U[i][j][k+1]-U[i][j][k])/Hz[k+1]-
                 (U[i][j][k]-U[i][j][k-1])/Hz[k])/(Hz[k+1]+Hz[k])*2.0);
    }
    else if(k==vk)
    {
        lap_u=nu*(
                ((U[i1+1][j][k]-U[i1][j][k])/Hx[i1+1]-
                 (U[i1][j][k]-U[i1-1][j][k])/Hx[i1])/(Hx[i1+1]+Hx[i1])*2.0+
                ((U[i][j+1][k]-U[i][j][k])*2/(Hy[j+1]+Hy[j])-
                 (U[i][j][k]-U[i][j-1][k])*2/(Hy[j]+Hy[j-1]))*2/(0.5*Hy[j+1]+Hy[j]+0.5*Hy[j-1])+
                ((U[i][j][k+1]-U[i][j][k])/Hz[k+1]-
                 (U[i][j][k]-U[i][j][k-1])/Hz[k])/(Hz[k+1]+Hz[k])*2.0);
    }
    else if(k==wk)
    {
        lap_u=nu*(
                ((U[i1+1][j][k]-U[i1][j][k])/Hx[i1+1]-
                 (U[i1][j][k]-U[i1-1][j][k])/Hx[i1])/(Hx[i1+1]+Hx[i1])*2.0+
                ((U[i][j+1][k]-U[i][j][k])/Hy[j+1]-
                 (U[i][j][k]-U[i][j-1][k])/Hy[j])/(Hy[j+1]+Hy[j])*2.0 +
                ((U[i][j][k+1]-U[i][j][k])*2/(Hz[k+1]+Hz[k])-
                 (U[i][j][k]-U[i][j][k-1])*2/(Hz[k]+Hz[k-1]))*2/(0.5*Hz[k+1]+Hz[k]+0.5*Hz[k-1]));
    }

    long double grad_p = 0;
    if(k==uk)
        grad_p=(U[i][j][pk]-U[i-1][j][pk])/Hx[i];
    else if(k==vk)
        grad_p=(U[i][j][pk]-U[i][j-1][pk])/Hy[j];
    else if(k==wk)
        grad_p=(U[i][j][pk]-U[i][j][pk-1])/Hz[k];

    return grad_p/rho - lap_u;
}


long double A(long double ***U1, long double ***U2, int i, int j, int k)
{
    return A1(U1,U2,i,j,k)+A2(U2,i,j,k);
}


/////////////////////////////////////////////////////////////////////////////////
void eval_scalars(long double ***u, long double ***R, int i1, int j1,
        int k1, long double &Rn_F1, long double &Rn_F2, long double &F1_F1, 
        long double &F2_F2, long double &F1_F2)
{
    Rn_F1 = 0;
    Rn_F2 = 0;
    F1_F1 = 0;
    F2_F2 = 0;
    F1_F2 = 0;

    int c = global_carg[i1][j1][k1];
    for(int oc=0; oc<c; ++oc)
    {
        indexes &a = arg[i1][j1][k1][oc];
        int i = a.i;
        int j = a.j;
        int k = a.k;
        if(!G[i][j][k]) continue;
        long double F1=A1(u,Z,i,j,k)+A1(Z,u,i,j,k)+A2(Z,i,j,k);
        long double F2=A1(Z,Z,i,j,k);
        long double RRn=R[i][j][k];
        long double hh=Hx[i]*Hy[j]*Hz[k];

        Rn_F1+=RRn*F1*hh;
        Rn_F2+=RRn*F2*hh;
        F1_F1+=F1*F1*hh;
        F2_F2+=F2*F2*hh;
        F1_F2+=F1*F2*hh;
    }

}


long double calc_alpha(long double Rn_F1, long double Rn_F2, long double F1_F1,
        long double F2_F2, long double F1_F2)
{
    // в этом случае получается квадратное уравнение
    if(fabs(F2_F2)<1e-15)
    {
        if(fabsl(F1_F1)<1e-15)
            return 0;
        return -Rn_F1/F1_F1;
    }

    long double alpha=0;
    long double x[3];
    int t1=cubic(x,6*F1_F2/4.0/F2_F2,(2*F1_F1+4*Rn_F2)/4.0/F2_F2,2*Rn_F1/4.0/F2_F2);
    alpha=x[0];
    long double m=2*alpha*Rn_F1+(2*Rn_F2+F1_F1)*alpha*alpha+2*alpha*alpha*alpha*F1_F2+alpha*alpha*alpha*alpha*F2_F2;
    for(int i=1; i<t1; ++i)
    {
        long double t;
        t=x[i];
        long double m1=2*t*Rn_F1+2*t*t*Rn_F2+t*t*F1_F1+2*t*t*t*F1_F2+t*t*t*t*F2_F2;
        if(m>m1)
            alpha=t,m=m1;
    }
    return alpha;
}



int cubic(long double *x,long double a,long double b,long double c)
{
    long double q,r,r2,q3;
    q=(a*a-3*b)/9;
    r=(2*a*a*a-9*a*b+27*c)/54;
    r2=r*r;
    q3=q*q*q;
    if(r2<q3)
    {
        long double t=acos(r/sqrt(q3))/3.0;
        x[0]=-2*sqrt(q)*cos(t)-a/3.0;
        x[1]=-2*sqrt(q)*cos(t+M_2PI/3)-a/3.0;
        x[2]=-2*sqrt(q)*cos(t-M_2PI/3)-a/3.0;
        return(3);
    }
    else
    {
        long double A=(r>0?-1:1)*powl(fabs(r)+sqrtl(r2-q3),1/3.0);
        long double B=(fabsl(A)>1e-15)?(q/A):0;
        x[0]=A+B-a/3.0;
        x[1]=-0.5*(A+B)-a/3.0;
        x[2]=sqrtl(3.0)*0.5*(A-B);
        if(fabsl(x[2])<1e-10) return(2);
        return(1);
    }
}

/////////////////////////////////////////////////////////////////////////////////
//Первый шаг
void tau_iter()
{
    long double Rn_F1=0,Rn_F2=0,F1_F1=0,F2_F2=0,F1_F2=0;
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                if(G[i][j][k])
                {
                    long double F1 = A1(U,R,i,j,k)+A1(R,U,i,j,k)+A2(R,i,j,k);
                    long double F2 = A1(R,R,i,j,k);
                    long double Rn = R[i][j][k];
                    long double hh = Hx[i]*Hy[j]*Hz[k];

                    Rn_F1+=Rn*F1*hh;
                    Rn_F2+=Rn*F2*hh;
                    F1_F1+=F1*F1*hh;
                    F2_F2+=F2*F2*hh;
                    F1_F2+=F1*F2*hh;
                }
            }
        }
    }
    long double tau=0;
    if(fabs(F2_F2)>1e-15)
    {
        long double x[3];
        int t1=cubic(x,6*F1_F2/4.0/F2_F2,(2*F1_F1+4*Rn_F2)/4.0/F2_F2,2*Rn_F1/4.0/F2_F2);
        tau=x[0];
        long double m=2*tau*Rn_F1+2*tau*tau*Rn_F2+tau*tau*F1_F1+2*tau*tau*tau*F1_F2+tau*tau*tau*tau*F2_F2;
        for(int i=1; i<t1; ++i)
        {
            long double t;
            t=x[i];
            long double m1=2*t*Rn_F1+2*t*t*Rn_F2+t*t*F1_F1+2*t*t*t*F1_F2+t*t*t*t*F2_F2;
            if(m>m1)
                tau=t,m=m1;
        }
    }
    // в этом случае получается квадратное уравнение
    // решать нужно по другому
    else
        tau=-Rn_F1/F1_F1;

    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                if(G[i][j][k])
                    U[i][j][k]+=tau*R[i][j][k];
    }
    // невязка
}

// Второй шаг схемы н/а
void alpha_iter()
{
    for(int i1=0; i1<Nx; ++i1)
    {
        for(int j1=0; j1<Ny; ++j1)
        {
            for(int k1=0; k1<Nz; ++k1)
            {
                int g=G[i1][j1][k1];
                if( (g) )
                {
                    Z[i1][j1][k1]=1.0;
                    long double Rn_F1=0,Rn_F2=0,F1_F1=0,F2_F2=0,F1_F2=0;
                    long double _Rn_F1=0,_Rn_F2=0,_F1_F1=0,_F2_F2=0,_F1_F2=0;
                    eval_scalars(U,R,i1,j1,k1,Rn_F1,Rn_F2,F1_F1,F2_F2,F1_F2);
                    long double alpha=calc_alpha(Rn_F1,Rn_F2,F1_F1,F2_F2,F1_F2);
                    U[i1][j1][k1]+=alpha;
                    Z[i1][j1][k1] = 0;
                    residual(i1,j1,k1);
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////
// ускорение
void speed_first()
{
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                U_2[i][j][k] = U_1[i][j][k];
                U_1[i][j][k] = U[i][j][k];
            }

        }
    }
    Rn_2=Rn_1;
    Rn_1=Rn;
}

void speed_work()
{
    long double Rn_F1=0,Rn_F2=0,F1_F1=0,F2_F2=0,F1_F2=0;
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                if(G[i][j][k])
                {

                    long double F1 = 2*A1(U_2,U_2,i,j,k) - A1(U,U_2,i,j,k) -
                        A1(U_2,U,i,j,k) + A2(U_2,i,j,k) - A2(U,i,j,k);

                    long double F2 = A1(U,U,i,j,k) - A1(U,U_2,i,j,k) -
                        A1(U_2,U,i,j,k) + A1(U_2,U_2,i,j,k);

                    long double Rn = A(U_2,U_2,i,j,k);
                    long double hh = Hx[i]*Hy[j]*Hz[k];

                    Rn_F1+=Rn*F1*hh;
                    Rn_F2+=Rn*F2*hh;
                    F1_F1+=F1*F1*hh;
                    F2_F2+=F2*F2*hh;
                    F1_F2+=F1*F2*hh;
                }
    }

    long double tau=0;
    if(fabs(F2_F2)>1e-15)
    {
        long double x[3];
        int t1=cubic(x,6*F1_F2/4.0/F2_F2,(2*F1_F1+4*Rn_F2)/4.0/F2_F2,2*Rn_F1/4.0/F2_F2);
        tau=x[0];
        long double m= 2*tau*Rn_F1+(2*Rn_F2+F1_F1)*tau*tau+2*tau*tau*tau*F1_F2+tau*tau*tau*tau*F2_F2;
        for(int i=1; i<t1; ++i)
        {
            long double t;
            t=x[i];
            long double m1=2*t*Rn_F1+2*t*t*Rn_F2+t*t*F1_F1+2*t*t*t*F1_F2+t*t*t*t*F2_F2;
            if(m>m1)
                tau=t,m=m1;
        }
    }
    // в этом случае получается квадратное уравнение
    // решать нужно по другому
    else if(F1_F2>1e-10)
    {
        long double a = 6.0*F1_F2;
        long double b = 2.0*(F1_F1 + 2.0*Rn_F2);
        long double c = 2.0*Rn_F1;

        long double D = b*b - 4*a*c;
        long double
            x1 = ( -b + sqrt(D) )/(2*a),
               x2 = ( -b - sqrt(D) )/(2*a);
        if( (a*x1*x1 + b*x1)<(a*x2*x2 + b*x2) )
            tau = x1;
        else
            tau = x2;
    }
    else
        tau=-Rn_F1/F1_F1;

    long double omega=tau;
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                if(G[i][j][k])
                {
                    long double tmp=U[i][j][k];    // сохранение для ускорения
                    U[i][j][k] = (1+omega)*U_2[i][j][k]-omega*U[i][j][k];
                }
    }
}

/*
 * Загрузка координат из файла
 */
void load_coords(const char *file_x_name,
        const char *file_y_name,
        const char *file_z_name)
{
    FILE *f = fopen(file_x_name,"r");

    if (!f)
    {
        printf("x coord file not found!\n");
        exit(1);
    }

    for(int i=0; i<Nx; ++i)
    {
        fscanf(f,"%LF\n", &Cx[i]);
    }

    fclose(f);

    f = fopen(file_y_name,"r");
    if (!f)
    {
        printf("y coord file not found!\n");
        exit(1);
    }

    for(int j=0; j<Ny; ++j)
    {
        fscanf(f,"%LF\n", &Cy[j]);
    }

    fclose(f);

    f = fopen(file_z_name,"r");
    if (!f)
    {
        printf("z coord file not found!\n");
        exit(1);
    }

    for(int k=0; k<dNz; ++k)
    {
        fscanf(f,"%LF\n", &Cz[k]);
    }

    fclose(f);

    printf("Coordinates has been loaded\n");
}

/*
 * Загрузка маски из файла
 */
void load_mask(const int type, const char *file_name)
{
    FILE *f = fopen(file_name,"r");
    if (!f)
    {
        printf("Mask file for type %d not found!\n", type);
        exit(1);
    }

    int level = 0;
    switch(type)
    {
        case VELOCITY_U:
            level = 0;
            break;
        case VELOCITY_V:
            level = 1;
            break;
        case VELOCITY_W:
            level = 2;
            break;
        case PRESSURE:
            level = 3;
            break;
    }

    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<dNz; ++k)
            {
                fscanf(f,"%d ", &G[i][j][k + level*dNz]);
            }
        }
    }
    fclose(f);
    printf("Mask for type %d has been loaded\n", type);
}


// Инициализирующие функции
void h_init()
{
    for(int i=1; i<Nx; ++i)
        Hx[i]=Cx[i]-Cx[i-1];
    Hx[0]=Hx[1];
    for(int j=1; j<Ny; ++j)
        Hy[j]=Cy[j]-Cy[j-1];
    Hy[0]=Hy[1];
    for(int k=1; k<dNz; ++k)
        Hz[k]=Cz[k]-Cz[k-1];
    Hz[0]=Hz[1];

    for(int i=0; i<dNz; ++i)
    {
        Hz[i+3*dNz] = Hz[i+2*dNz] = Hz[i+1*dNz] = Hz[i];
    }
}

void U_init()
{
    int vortex_inside_only = 1;
    // Первоначально все нули
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                U[i][j][k]=0;
    }

    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k< dNz; ++k)
            {
                long double p = p_left - (p_left-p_right)*i/(Nx-1);
                if(G[i][j][k + 3*dNz] == 1)
                {
                    U[i][j][k + 3*dNz] = p;
                }

                /* Pressure on in/out boundaries must be set and not calculated */
                if(G[i][j][k + 3*dNz] == 2)
                {
                    U[i][j][k + 3*dNz] = p;
                    G[i][j][k + 3*dNz] = 0;
                }

                if(G[i][j][k + 3*dNz] == 3)
                {
                    U[i][j][k + 3*dNz] = 0;
                    G[i][j][k + 3*dNz] = 0;
                }
            }
        }
    }
}

// Инициализация
void init()
{
    load_config();
    utils = new Utils(Nx, Ny, Nz);
    
    groupGenerator = new GroupsGenerator(Nx, Ny, Nz);
    groupGenerator->operator_nonlin = &A1;
    groupGenerator->operator_lin = &A2;

    U = utils->alloc_and_fill<long double>(Nx,Ny,Nz);
    R = utils->alloc_and_fill<long double>(Nx,Ny,Nz);
    G = utils->alloc_and_fill<int>(Nx,Ny,Nz);
    Z = utils->alloc_and_fill<long double>(Nx,Ny,Nz);

    U_1 = utils->alloc_and_fill<long double>(Nx,Ny,Nz);
    U_2 = utils->alloc_and_fill<long double>(Nx,Ny,Nz);

    Rn_1 = Rn_2 = 0;

    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                U_2[i][j][k] = U_1[i][j][k] = 0;
    }


    Hx = new long double [Nx];
    Hy = new long double [Ny];
    Hz = new long double [Nz];
    Cx = new long double [Nx];
    Cy = new long double [Ny];
    Cz = new long double [Nz];

    load_mask(VELOCITY_U, "u_area.mask");
    load_mask(VELOCITY_V, "v_area.mask");
    load_mask(VELOCITY_W, "w_area.mask");
    load_mask(PRESSURE, "pressure.mask");
    load_coords("area.x.coord","area.y.coord","area.z.coord");
    h_init();
    U_init();

    residual();
    Rn = norm(R);
    R0 = Rn;

    // вектор Z
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                Z[i][j][k]=0;
            }
        }
    }

    output = new Output(Nx, Ny, Nz, dNz, U, Cx, Cy, Cz, R, Hx, Hy, Hz, G);
    // группы
    if(generate_groups)
    {
        groupGenerator->init_gr();
        groupGenerator->set_gr();
        groupGenerator->print_gr();
    }
    else
    {
        groupGenerator->load_groups();
    }

    global_carg = groupGenerator->get_carg();
    arg = groupGenerator->get_arg();
}

// Основной цикл
void run()
{
    residual();
    iters = 0;
    time_t start_time = time(NULL);

    do
    {

#ifdef DEBUG
    output->print_vtk();
    printf("next step...");
    getchar();
#endif
         ++iters;
        speed_first();
        long double R1 = norm(R);

        tau_iter();
        residual();
        long double R2 = norm(R);

        alpha_iter();
        residual();
        long double R3 = norm(R);

        if( (iters>5)&&( (iters%100==0)||(iters%100==1)||(iters%100==2)||(iters%100==3)||(iters%100==4) ) )
        {
            speed_work();
            residual();
        }
        long double R4 = norm(R);

        if( (R4>R3) )
        {
            printf("Error on speed\n");
        }
        else if ( (R3>R2) )
        {
            printf("Error on alpha\n");
        }
        else if ( (R2>R1) )
        {
            printf("Error on tau\n");
        }

        Rn = norm(R);

        if(iters%100==0)
        {
            printf("%5d: %3.8LF %3.8LF\n", iters, Rn, Rn/R0);
        }
    }
    while (Rn/R0>eps);

    time_t end_time = time(NULL);
    printf("Time is %ld\n", end_time - start_time);
    output->print_vtk();
}

// Деструктор
void down()
{
    utils->del(U,Nx,Ny,Nz);
    utils->del(R,Nx,Ny,Nz);
    utils->del(G,Nx,Ny,Nz);
    utils->del(Z,Nx,Ny,Nz);

    delete [] Hx;
    delete [] Hy;
}

// Точка входа
int main()
{
    if(generate_groups)
    {
        return 0;
    }
    else
    {
        run();
        down();
    }
    return 0;
}

/*
 * Debugging functions
 */
void printm(long double ***m, int size_x, int size_y, int size_z)
{
    for (int i = 0; i < size_x; i++) {
        for (int j = 0; j < size_y; j++) {
            for (int k = 0; k < size_z; k++) {
                printf("%LF ", m[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

/*
 * Library initialization
 */
__attribute__((constructor)) void init_lib(void) {
    init();
}
