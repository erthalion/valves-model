#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
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

const int COORD_X = 0;
const int COORD_Y = 1;
const int COORD_Z = 2;

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

long double ***force_X, ***force_Y, ***force_Z;

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
long double force(int i, int j, int k);
void compute_boundary_forces(ImmersedBoundary *boundary);
void spread_force(ImmersedBoundary *boundary);
void interpolate(ImmersedBoundary *boundary);
void update_boundary_position(ImmersedBoundary *boundary);

void stop_handler();
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


#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

/*
 * @brief
 *      Return integer index of node by the coordinate.
 *      Assumed temporary that spatial steps is equals
 *      by Ox, Oy, Oz and grid is uniform.
 * @param coord
 *      Coordinate of node by Ox/Oy/Oz
 * @param type
 *      One of COORD_X, COORD_Y, COORD_Z
 * @return 
 *      Index of fluid node, that coordinate is lower than coord by selected axis
*/
int index(const long double coord, const int type)
{
    switch(type)
    {
        case COORD_X:
            return floor(coord / Hx[0]);
        case COORD_Y:
            return floor(coord / Hy[0]);
        case COORD_Z:
            return floor(coord / Hz[0]);
        default:
            throw "Incorrect type of axis";
    }
}

/*
 * @brief
 *      Return coordinate of node by the index.
 *      Assumed temporary that spatial steps is equals
 *      by Ox, Oy, Oz and grid is uniform.
 * @param index
 *      Index of node by Ox/Oy/Oz
 * @param type
 *      One of COORD_X, COORD_Y, COORD_Z
 * @return 
 *      Coordinate of fluid node
*/
long double coord(const int index, const int type)
{
    switch(type)
    {
        case COORD_X:
            return index * Hx[0];
        case COORD_Y:
            return index * Hy[0];
        case COORD_Z:
            return index * Hz[0];
        default:
            throw "Incorrect type of axis";
    }
}

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
// TODO: Add force terms to x, y, z equation accordingly
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
    long double force_term = 0;
    if(k==uk)
    {
        //printf("force_X[%d][%d][%d] = %LF\n", i, j, k, force_X[i][j][pk]);
        grad_p=(U[i][j][pk]-U[i-1][j][pk])/Hx[i];
        force_term = force_X[i][j][uk];
    }
    else if(k==vk)
    {
        //printf("force_Y[%d][%d][%d] = %LF\n", i, j, k, force_Y[i][j][pk]);
        grad_p=(U[i][j][pk]-U[i][j-1][pk])/Hy[j];
        force_term = force_Y[i][j][uk];
    }
    else if(k==wk)
    {
        //printf("force_Z[%d][%d][%d] = %LF\n", i, j, k, force_Z[i][j][pk]);
        grad_p=(U[i][j][pk]-U[i][j][pk-1])/Hz[k];
        force_term = force_Z[i][j][uk];
    }

    // TODO: Calculate force term like grad_p by k value
    return grad_p/rho - lap_u - force_term;
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
                //raise(SIGINT);
                long double p = 0;
                if (i == 0 || i == Nx-1)
                {
                    p = 0;
                }
                else
                {
                    p = p_left - (p_left-p_right)*(i-1)/(Nx-2);
                }

                if(G[i][j][k + 3*dNz] == 1 || G[i][j][k + 3*dNz] == 2 || G[i][j][k + 3*dNz] == 3)
                {
                    U[i][j][k + 3*dNz] = p;
                }

                /* Pressure on in/out boundaries must be set and not calculated */
                if(G[i][j][k + 3*dNz] == 2)
                {
                    U[i][j][k + 3*dNz] = p_left;
                    G[i][j][k + 3*dNz] = 0;
                }

                if(G[i][j][k + 3*dNz] == 3)
                {
                    U[i][j][k + 3*dNz] = 0;
                    G[i][j][k + 3*dNz] = 0;
                }

                //printf("%d ", G[i][j][k + 3*dNz]);
            }
            //printf("\n");
        }
        //printf("\n");
        //getchar();
    }
}

void check_boundary(const char* message, ImmersedBoundary *boundary)
{
    printf("%s\n", message);
    for(int n = 0; n < boundary->nodes_count; ++n)
    {
        printf("node %d: %lf, %lf\t %lf, %lf\t %lf, %lf\n", n,
                boundary->nodes[n].x,
                boundary->nodes[n].x_ref,
                boundary->nodes[n].y,
                boundary->nodes[n].y_ref,
                boundary->nodes[n].z,
                boundary->nodes[n].z_ref
              );
    }
    getchar();
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

    force_X = utils->alloc_and_fill<long double>(Nx, Ny, dNz);
    force_Y = utils->alloc_and_fill<long double>(Nx, Ny, dNz);
    force_Z = utils->alloc_and_fill<long double>(Nx, Ny, dNz);

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

    output = new Output(Nx, Ny, Nz, dNz, U, Cx, Cy, Cz, force_X, force_Y, force_Z, R, Hx, Hy, Hz, G);
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

void compute_fluid(int iteration)
{
    residual();
    iters = 0;
    do
    {
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

    //output->print_info(iters, R0, Rn);
    //output->print_vtk(iteration);
}

// Основной цикл
void run()
{
    ImmersedBoundary *boundary = new ImmersedBoundary();
    int iterations_count = 10;

    for (int i = 0; i < iterations_count; i++) {
        printf("Iteration %d\n", i);
        compute_boundary_forces(boundary);
        spread_force(boundary);
        compute_fluid(i);
        interpolate(boundary);
        update_boundary_position(boundary);
        output->print_vtk(i);
        output->print_boundary_vtk(i, boundary);
        //output->print_boundary(i, boundary, U, dNz);
    }

}

// Деструктор
void down()
{
    utils->del(U,Nx,Ny,Nz);
    utils->del(R,Nx,Ny,Nz);
    utils->del(G,Nx,Ny,Nz);
    utils->del(Z,Nx,Ny,Nz);

    utils->del(force_X, Nx, Ny, Nz);
    utils->del(force_Y, Nx, Ny, Nz);
    utils->del(force_Z, Nx, Ny, Nz);

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


void stop_handler(int s)
{
    printf("Caught signal %d\n",s);
    //output->print_info(iters, R0, Rn);
    output->print_vtk(99999999);
    down();
    exit(1); 
}


/*
 * Library initialization
 */
__attribute__((constructor)) void init_lib(void) {
    struct sigaction sigIntHandler;

    sigIntHandler.sa_handler = stop_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;

    // TODO: Need SIGTERM
    sigaction(SIGINT, &sigIntHandler, NULL);

    init();
}

void compute_boundary_forces(ImmersedBoundary *boundary)
{
    for(int n = 0; n < boundary->nodes_count; ++n) {
        boundary->nodes[n].x_force = 0;
        boundary->nodes[n].y_force = 0;
        boundary->nodes[n].z_force = 0;
    }

    const double area = 4 * M_PI * SQ(boundary->radius) / boundary->nodes_count;

    for(int n = 0; n < boundary->nodes_count; ++n) {
        boundary->nodes[n].x_force += -boundary->stiffness * (boundary->nodes[n].x - boundary->nodes[n].x_ref) * area;
        boundary->nodes[n].y_force += -boundary->stiffness * (boundary->nodes[n].y - boundary->nodes[n].y_ref) * area;
        boundary->nodes[n].z_force += -boundary->stiffness * (boundary->nodes[n].z - boundary->nodes[n].z_ref) * area;

    }

    return;
}


void spread_force(ImmersedBoundary *boundary)
{
    for(int i = 0; i < Nx; ++i) {
        for(int j = 1; j < Ny - 1; ++j) {
            for(int k = 1; k < dNz - 1; ++k) {
                force_X[i][j][k] = 0;
                force_Y[i][j][k] = 0;
                force_Z[i][j][k] = 0;
            }
        }
    }

    for(int n = 0; n < boundary->nodes_count; ++n) {
        int x_int = index(boundary->nodes[n].x, COORD_X);
        int y_int = index(boundary->nodes[n].y, COORD_Y);
        int z_int = index(boundary->nodes[n].z, COORD_Z);

        for(int i = x_int; i <= x_int + 1; ++i) {
            for(int j = y_int; j <= y_int + 1; ++j) {
                for(int k = z_int; k <= z_int + 1; ++k) {

                    const double dist_x = fabs(boundary->nodes[n].x - coord(i, COORD_X));
                    const double dist_y = fabs(boundary->nodes[n].y - coord(j, COORD_Y));
                    const double dist_z = fabs(boundary->nodes[n].z - coord(k, COORD_Z));

                    const double weight_x = 1 - dist_x;
                    const double weight_y = 1 - dist_y;
                    const double weight_z = 1 - dist_z;

                    force_X[i][j-1][k] += (boundary->nodes[n].x_force * weight_x * weight_y * weight_z)/boundary->nodes_count;
                    force_Y[i][j-1][k] += (boundary->nodes[n].y_force * weight_x * weight_y * weight_z)/boundary->nodes_count;
                    force_Z[i][j-1][k] += (boundary->nodes[n].z_force * weight_x * weight_y * weight_z)/boundary->nodes_count;
                }
            }
        }

    }

    return;
}


void interpolate(ImmersedBoundary *boundary)
{
    for(int n = 0; n < boundary->nodes_count; ++n)
    {
        boundary->nodes[n].x_vel = 0;
        boundary->nodes[n].y_vel = 0;
        boundary->nodes[n].z_vel = 0;

        int x_int = index(boundary->nodes[n].x, COORD_X);
        int y_int = index(boundary->nodes[n].y, COORD_Y);
        int z_int = index(boundary->nodes[n].z, COORD_Z);

        for(int i = x_int; i <= x_int + 1; ++i)
        {
            for(int j = y_int; j <= y_int + 1; ++j)
            {
                for(int k = z_int; k <= z_int + 1; ++k)
                {

                    const double dist_x = fabs(boundary->nodes[n].x - coord(i, COORD_X));
                    const double dist_y = fabs(boundary->nodes[n].y - coord(j, COORD_Y));
                    const double dist_z = fabs(boundary->nodes[n].z - coord(k, COORD_Z));

                    const double weight_x = 1 - abs(dist_x);
                    const double weight_y = 1 - abs(dist_y);
                    const double weight_z = 1 - abs(dist_z);

                    // interpolation from staggered grid before computation

                    // I don't know, why j-1 is needed..., may be it's related with fictive shapes?
                    long double velocity_U = U[i][j-1][k+VELOCITY_U*dNz];
                    long double velocity_V = U[i][j-1][k+VELOCITY_V*dNz];
                    long double velocity_W = U[i][j-1][k+VELOCITY_W*dNz];

                    // 1 is a density
                    //boundary->nodes[n].x_vel += (
                            //(
                             //velocity_U*SQ(Hx[i]) + 0.5 * force_X[i][j-1][k] / 1
                             //) * weight_x * weight_y * weight_z);

                    //boundary->nodes[n].y_vel += (
                            //(
                             //velocity_V*SQ(Hy[j]) + 0.5 * (force_Y[i][j-1][k]) / 1
                            //) * weight_x * weight_y * weight_z);
                    //boundary->nodes[n].z_vel += (
                            //(
                             //velocity_W*SQ(Hz[k]) + 0.5 * (force_Z[i][j-1][k]) / 1
                            //) * weight_x * weight_y * weight_z);

                    boundary->nodes[n].x_vel += (velocity_U*SQ(Hx[i]) * weight_x * weight_y * weight_z);
                    boundary->nodes[n].y_vel += (velocity_V*SQ(Hy[j]) * weight_x * weight_y * weight_z);
                    boundary->nodes[n].z_vel += (velocity_W*SQ(Hz[k]) * weight_x * weight_y * weight_z);
                }
            }
        }
    }

    return;
}


void update_boundary_position(ImmersedBoundary *boundary)
{
    for(int n = 0; n < boundary->nodes_count; ++n)
    {
        boundary->nodes[n].x += boundary->nodes[n].x_vel;
        boundary->nodes[n].y += boundary->nodes[n].y_vel;
        boundary->nodes[n].z += boundary->nodes[n].z_vel;
    }

    return;
}
