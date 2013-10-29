#include "output.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "utils.hpp"

using namespace std;

Output::Output(int Nx, int Ny, int Nz, int dNz
            , long double ***U
            , long double *Cx, long double *Cy, long double *Cz
            , long double ***R
            , long double *Hx, long double *Hy, long double *Hz
            , int ***G)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->dNz = dNz;

    this->Cx = Cx;
    this->Cy = Cy;
    this->Cz = Cz;

    this->Hx = Hx;
    this->Hy = Hy;
    this->Hz = Hz;

    this->U = U;
    this->R = R;
    this->G = G;
}


void Output::print_texplot_matrix(int iters)
{
    char name[40];
    sprintf(name, "Matrixdata%6d.dat", iters);
    FILE *f = fopen(name,"w");

    fprintf(f,"TITLE = Test\n");
    fprintf(f,"VARIABLES = X,Y,Z,U,V,W,P\n");
    fprintf(f,"ZONE T=Test,I=%d, J=%d, K=%d, F=POINT\n", Nx, Ny, dNz);
    for(int k=0; k<dNz; ++k)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%LF, %LF, %LF,", Cx[i], Cy[j], Cz[k]);

                fprintf(f,"%LF,", U[i][j][k]);

                fprintf(f,"%LF,", U[i][j][k+dNz]);

                fprintf(f,"%LF,", U[i][j][k+2*dNz]);

                fprintf(f,"%LF\n", U[i][j][k+3*dNz]);
            }
        }
    }

    fclose(f);

    strcpy(name,"");
    sprintf(name, "Res%6d.dat", iters);
    f = fopen(name,"w");

    fprintf(f,"TITLE = Test\n");
    fprintf(f,"VARIABLES = X,Y,Z,RU,RV,RW,P\n");
    fprintf(f,"ZONE T=Test,I=%d, J=%d, K=%d, F=POINT\n", Nx, Ny, dNz);
    for(int k=0; k<dNz; ++k)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%LF, %LF, %LF,", Cx[i], Cy[j], Cz[k]);

                fprintf(f,"%LF,", R[i][j][k]);

                fprintf(f,"%LF,", R[i][j][k+dNz]);

                fprintf(f,"%LF,", R[i][j][k+2*dNz]);

                fprintf(f,"%LF\n", R[i][j][k+3*dNz]);
            }
        }
    }

    fclose(f);

}

void Output::print_vtk()
{
    char output_path[] = "streamlines.vtk";
    //print_vtk_streamline_header(output_path, Nx, Ny-1, dNz-1);
    print_vtk_header(output_path, Nx, Ny-1, dNz-1);

    FILE *f = fopen(output_path,"a");

    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%LF %LF %LF\n", Cx[i], Cy[j], Cz[k]);
                //fprintf(f,"%LF\n",U[i][j][k+3*dNz]);
            }
        }
    }

    fclose(f);

    print_vtk_streamline_scalar_header(output_path, Nx, Ny-1, dNz-1);
    f = fopen(output_path,"a");

    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%LF\n", U[i][j][k+3*dNz]);
            }
        }
    }

    fclose(f);

    print_vtk_streamline_vector_header(output_path, Nx, Ny-1, dNz-1);

    f = fopen(output_path,"a");
    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                if(i==0)
                    fprintf(f,"%LF ",U[i+1][j][k]);
                else if(i==Nx-1)
                    fprintf(f,"%LF ",U[i][j][k]);
                else
                    fprintf(f,"%LF ", Hx[i]/(Hx[i+1]+Hx[i])*U[i+1][j][k]+Hx[i+1]/(Hx[i+1]+Hx[i])*U[i][j][k] );

                if(j==0)
                    fprintf(f,"%LF ",U[i][j][k+dNz]);
                else if(j==Ny-2)
                    fprintf(f,"%LF ",U[i][j+1][k+dNz]);
                else
                    fprintf(f,"%LF ", Hy[j]/(Hy[j]+Hy[j+1])*U[i][j+1][k+dNz]+
                            Hy[j+1]/(Hy[j]+Hy[j+1])*U[i][j][k+dNz]);

                if(k==0)
                    fprintf(f,"%LF ",U[i][j][k+2*dNz]);
                else if(k==dNz-2)
                    fprintf(f,"%LF ",U[i][j][k+1+2*dNz]);
                else
                    fprintf(f,"%LF ", Hz[k]/(Hz[k]+Hz[k+1])*U[i][j][k+1+2*dNz]+
                            Hz[k+1]/(Hz[k]+Hz[k+1])*U[i][j][k+2*dNz]);

                fprintf(f, "\n");
            }
        }
    }

    fclose(f);
}

void Output::print_vtk_streamline_vector_header(char *output_path, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"\nVECTORS uvw float\n");
    header.append(line);

    ofstream output_data;
    output_data.open(output_path,std::ios_base::app | std::ios_base::out);
    output_data << header;
    output_data.close();
}

void Output::print_vtk_streamline_scalar_header(char *output_path, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"\nPOINT_DATA %d\n", sizeX * sizeY * sizeZ);
    header.append(line);

    sprintf(line,"SCALARS p float\n");
    header.append(line);

    sprintf(line,"LOOKUP_TABLE default\n");
    header.append(line);

    ofstream output_data;
    output_data.open(output_path,std::ios_base::app | std::ios_base::out);
    output_data << header;
    output_data.close();
}

void Output::print_vtk_streamline_header(char *output_path, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"# vtk DataFile Version 1.0\n");
    header.append(line);

    sprintf(line,"Data file for valves model\n");
    header.append(line);

    sprintf(line,"ASCII\n");
    header.append(line);

    sprintf(line,"DATASET STRUCTURED_POINTS\n");
    header.append(line);

    sprintf(line,"DIMENSIONS %d %d %d\n",sizeX,sizeY,sizeZ);
    header.append(line);

    sprintf(line,"ORIGIN %f %f %f\n", 0.0, 0.0, 0.0);
    header.append(line);

    sprintf(line,"ASPECT_RATIO %f %f %f\n", 1.0, 1.0, 1.0);
    header.append(line);

    sprintf(line,"\nPOINT_DATA %d\n", sizeX*sizeY*sizeZ);
    header.append(line);

    sprintf(line,"SCALARS xyz double 3\n");
    header.append(line);

    sprintf(line,"LOOKUP_TABLE default\n");
    header.append(line);

    ofstream output_data;
    output_data.open(output_path);
    output_data << header;
    output_data.close();
}

void Output::print_texplot(int iters)
{
    char name[20];
    sprintf(name, "data%6d.dat", iters);
    FILE *f = fopen(name,"w");

    fprintf(f,"TITLE = Test\n");
    fprintf(f,"VARIABLES = X,Y,Z,U,V,W,P\n");
    fprintf(f,"ZONE T=Test,I=%d, J=%d, K=%d, F=POINT\n", Nx, Ny-1, dNz-1);
    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%LF, %LF, %LF,", Cx[i], Cy[j], Cz[k]);

                if(i==0)
                    fprintf(f,"%LF,",U[i+1][j][k]);
                else if(i==Nx-1)
                    fprintf(f,"%LF,",U[i][j][k]);
                else
                    fprintf(f,"%LF,", Hx[i]/(Hx[i+1]+Hx[i])*U[i+1][j][k]+Hx[i+1]/(Hx[i+1]+Hx[i])*U[i][j][k] );

                if(j==0)
                    fprintf(f,"%LF,",U[i][j][k+dNz]);
                else if(j==Ny-2)
                    fprintf(f,"%LF,",U[i][j+1][k+dNz]);
                else
                    fprintf(f,"%LF,", Hy[j]/(Hy[j]+Hy[j+1])*U[i][j+1][k+dNz]+
                            Hy[j+1]/(Hy[j]+Hy[j+1])*U[i][j][k+dNz]);

                if(k==0)
                    fprintf(f,"%LF,",U[i][j][k+2*dNz]);
                else if(k==dNz-2)
                    fprintf(f,"%LF,",U[i][j][k+1+2*dNz]);
                else
                    fprintf(f,"%LF,", Hz[k]/(Hz[k]+Hz[k+1])*U[i][j][k+1+2*dNz]+
                            Hz[k+1]/(Hz[k]+Hz[k+1])*U[i][j][k+2*dNz]);

                fprintf(f,"%LF\n",U[i][j][k+3*dNz]);
            }
        }
    }

    fclose(f);
}

/*
 * Визуализация расчетной области
 */
void Output::print_area()
{
    char output_path[] = "surface.vtk";
    print_vtk_header(output_path, Nx, Ny-1, dNz-1);

    FILE *f = fopen(output_path,"a");

    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%LF %LF %LF\n", Cx[i],Cy[j],Cz[k]);
            }
        }
    }

    fclose(f);

    print_vtk_data_header(output_path, Nx, Ny-1, dNz-1);

    f = fopen(output_path,"a");
    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%d\n", G[i][j][k]);
            }
        }
    }

    fclose(f);
}

/*
 * Визуализация расчетной области в STRUCTURED_POINTS
 */
void Output::print_area_points()
{
    char output_path[] = "surface.vtk";
    print_vtk_header_points(output_path, Nx, Ny-1, dNz-1);
    print_vtk_data_header(output_path, Nx, Ny-1, dNz-1);

    FILE *f = fopen(output_path,"a");
    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f,"%d\n", G[i][j][k]);
            }
        }
    }

    fclose(f);
}

/*
 * Запись в файл заголовка данных vtk
 */
void Output::print_vtk_data_header(char *output_path, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"POINT_DATA %d\n",sizeX*sizeY*sizeZ);
    header.append(line);

    sprintf(line,"SCALARS scalars int\n");
    header.append(line);

    sprintf(line,"LOOKUP_TABLE default\n");
    header.append(line);

    ofstream output_data;
    output_data.open(output_path,std::ios_base::app | std::ios_base::out);
    output_data << header;
    output_data.close();
}


/*
 * Запись в файл заголовка vtk
 */
void Output::print_vtk_header(char *output_path, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"# vtk DataFile Version 1.0\n");
    header.append(line);

    sprintf(line,"Data file for valves model\n");
    header.append(line);

    sprintf(line,"ASCII\n");
    header.append(line);

    sprintf(line,"DATASET STRUCTURED_GRID\n");
    header.append(line);

    sprintf(line,"DIMENSIONS %d %d %d\n",sizeX,sizeY,sizeZ);
    header.append(line);

    sprintf(line,"POINTS %d double\n",sizeX*sizeY*sizeZ);
    header.append(line);

    ofstream output_data;
    output_data.open(output_path);
    output_data << header;
    output_data.close();
}

/*
 * Запись в файл заголовка vtk STRUCTURED_POINTS
 */
void Output::print_vtk_header_points(char *output_path, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"# vtk DataFile Version 1.0\n");
    header.append(line);

    sprintf(line,"Data file for valves model\n");
    header.append(line);

    sprintf(line,"ASCII\n");
    header.append(line);

    sprintf(line,"DATASET STRUCTURED_POINTS\n");
    header.append(line);

    sprintf(line,"DIMENSIONS %d %d %d\n",sizeX,sizeY,sizeZ);
    header.append(line);

    sprintf(line,"ORIGIN %d %d %d\n",0,0,0);
    header.append(line);

    sprintf(line,"SPACING %d %d %d\n",1,1,1);
    header.append(line);

    ofstream output_data;
    output_data.open(output_path);
    output_data << header;
    output_data.close();
}

void Output::print_pressure()
{
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<dNz; ++k)
            {
                printf("%LF ", U[i][j][k+ 3*dNz]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
 }

void Output::print_info(int iters, long double R0, long double Rn)
{
    //print_texplot(iters);
    //print_texplot_matrix(iters);

    char fn [20];
    sprintf(fn, "info%6d.txt", iters);
    FILE *f = fopen(fn,"w");
    fprintf(f,"iters = %d\n\n", iters);
    fprintf(f,"r0 = %LF\n", R0);
    fprintf(f,"rn = %LF\n", Rn);
    fprintf(f,"rn/ro = %LF\n", Rn/R0);

    long double max_ri = fabs(R[0][0][0]);
    int im = 0, jm = 0, km = 0;
    for(int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
                if( max_ri < fabs(R[i][j][k]) )
                {
                    max_ri = fabs(R[i][j][k]);
                    im = i;
                    jm = j;
                    km = k;
                }
    fprintf(f,"max_ri = %LF at (%d,%d,%d)\n\n", max_ri, im, jm, km);

    long double s1=0, s2=0;

    /*for(int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j)
        {
            s1 += U[i][j][1]*Hx[i]*Hy[j];
            s2 += U[i][j][dNz-3]*Hx[i]*Hy[j];
        }*/


    /*for(int i=0; i<Nx; ++i)
        for(int k=0; k<dNz; ++k)
        {
            s1 += U[i][1][k]*Hx[i]*Hz[k];
            s2 += U[i][Ny-3][k]*Hx[i]*Hz[k];
        }*/

    for(int j=0; j<Ny; ++j)
        for(int k=0; k<dNz; ++k)
        {
            s1 += U[1][j][k]*Hy[j]*Hz[k];
            s2 += U[Nx-1][j][k]*Hy[j]*Hz[k];
        }

    /*for (int j = 0; j < Ny; j++)
    {
        for (int k = 0; k < dNz; k++)
        {
            printf("%LF ", U[1][j][k+3*dNz]);
        }
        printf("\n");
    }
    printf("\n");

    for (int j = 0; j < Ny; j++)
    {
        for (int k = 0; k < dNz; k++)
        {
            printf("%LF ", U[Nx-1][j][k+3*dNz]);
        }
        printf("\n");
    }
    printf("\n");
    getchar();*/

    /*for (int i = 0; i < Nx; i++)
    {
        for (int k = 0; k < dNz; k++)
        {
            printf("%LF ", U[i][1][k+3*dNz]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < Nx; i++)
    {
        for (int k = 0; k < dNz; k++)
        {
            printf("%LF ", U[i][Ny-3][k+3*dNz]);
        }
        printf("\n");
    }
    printf("\n");
    getchar();*/


    /*for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            printf("%LF ", U[i][j][1+3*dNz]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            printf("%LF ", U[i][j][dNz-3 + 3*dNz]);
        }
        printf("\n");
    }
    printf("\n");
    getchar();*/



    fprintf(f,"s1 = %0.18LF\n", s1);
    fprintf(f,"s2 = %0.18LF\n", s2);
    fprintf(f,"ds = %0.18LF\n", (long double)fabs(s1-s2));
    fclose(f);

}

void Output::print_boundary(int iter, ImmersedBoundary *boundary)
{
    char output_path[20];
    sprintf(output_path, "boundary%d.dat", iter);
    FILE *f = fopen(output_path,"a");

    for(int n = 0; n < boundary->nodes_count; ++n)
    {
        fprintf(f, "%1.8lf %1.8lf %1.8lf\n", boundary->nodes[n].x, boundary->nodes[n].y, boundary->nodes[n].z);
    }
    fclose(f);
}
