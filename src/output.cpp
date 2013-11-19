#include "output.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "utils.hpp"
#include <csignal>

using namespace std;

Output::Output(int Nx, int Ny, int Nz, int dNz
            , long double ***U
            , long double *Cx, long double *Cy, long double *Cz
            , long double ***force_X, long double ***force_Y, long double ***force_Z
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

    this->force_X = force_X;
    this->force_Y = force_Y;
    this->force_Z = force_Z;

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

void Output::print_vtk(int iteration)
{
    char output_path[20];
    sprintf(output_path, "streamlines%d.vtk", iteration);
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

    print_vtk_streamline_scalar_header(output_path, "p",  Nx, Ny-1, dNz-1);
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

    print_vtk_streamline_vector_header(output_path, "uvw", Nx, Ny-1, dNz-1);

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

    print_vtk_streamline_vector_header(output_path, "force", Nx, Ny-1, dNz-1);

    f = fopen(output_path,"a");
    for(int k=0; k<dNz-1; ++k)
    {
        for(int j=0; j<Ny-1; ++j)
        {
            for(int i=0; i<Nx; ++i)
            {
                fprintf(f, "%LF %LF %LF\n", force_X[i][j][k], force_Y[i][j][k], force_Z[i][j][k]);
            }
        }
    }

    fclose(f);

}

void Output::print_vtk_streamline_vector_header(char *output_path, const char *name, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"\nVECTORS %s float\n", name);
    header.append(line);

    ofstream output_data;
    output_data.open(output_path,std::ios_base::app | std::ios_base::out);
    output_data << header;
    output_data.close();
}

void Output::print_vtk_streamline_scalar_header(char *output_path, const char *name, int sizeX, int sizeY, int sizeZ)
{
    string header;
    char line [50];

    sprintf(line,"\nPOINT_DATA %d\n", sizeX * sizeY * sizeZ);
    header.append(line);

    sprintf(line,"SCALARS %s float\n", name);
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

// debug puposes
const int COORD_X = 0;
const int COORD_Y = 1;
const int COORD_Z = 2;

const int VELOCITY_U = 0;
const int VELOCITY_V = 1;
const int VELOCITY_W = 2;
const int PRESSURE = 3;

const long double hx = 0.2;
const long double hy = 0.1;
const long double hz = 0.1;

int debug_index(const long double coord, const int type)
{
    switch(type)
    {
        case COORD_X:
            return floor(coord / hx);
        case COORD_Y:
            return floor(coord / hy);
        case COORD_Z:
            return floor(coord / hz);
        default:
            throw "Incorrect type of axis";
    }
}

long double debug_coord(const int index, const int type)
{
    switch(type)
    {
        case COORD_X:
            return index * hx;
        case COORD_Y:
            return index * hy;
        case COORD_Z:
            return index * hz;
        default:
            throw "Incorrect type of axis";
    }
}

void Output::print_boundary(int iter, ImmersedBoundary *boundary, long double ***U, int dNz)
{
    char output_path[20];
    sprintf(output_path, "boundary%d.dat", iter);
    FILE *f = fopen(output_path,"a");

    // get middle circle of cylinder
    int circle = 2;
    for (int n = 0; n < 10; n++) {
        Node node = boundary->nodes[n + circle*10];

        fprintf(f, "%LF %LF %LF %LF %LF %LF", node.x, node.y, node.x_force, node.y_force, node.x_ref, node.y_ref);

        int x_int = debug_index(node.x, COORD_X);
        int y_int = debug_index(node.y, COORD_Y);
        int z_int = debug_index(node.z, COORD_Z);

        for(int i = x_int; i <= x_int + 1; ++i)
        {
            for(int j = y_int; j <= y_int + 1; ++j)
            {
                fprintf(f, " %LF %LF %LF %LF",
                        debug_coord(i, COORD_X),
                        debug_coord(j, COORD_Y),
                        U[i][j-1][z_int + VELOCITY_U*dNz],
                        U[i][j-1][z_int + VELOCITY_V*dNz]);
            }
        }
        fprintf(f, "\n");

    }
    fclose(f);


    sprintf(output_path, "velocity%d.dat", iter);
    f = fopen(output_path,"a");

    Node node = boundary->nodes[circle*10];
    int z_int = debug_index(node.z, COORD_Z);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny-1; j++) {
            int k = z_int;

            fprintf(f, "%LF %LF\n", U[i][j][k + VELOCITY_U*dNz], U[i][j][k + VELOCITY_V*dNz]);
        }
    }

    fclose(f);

    sprintf(output_path, "force%d.dat", iter);
    f = fopen(output_path,"a");

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny-1; j++) {
            int k = z_int;

            fprintf(f, "%LF %LF\n", force_X[i][j][k], force_Y[i][j][k]);
        }
    }


    fclose(f);
    sprintf(output_path, "boundary_force%d.dat", iter);
    f = fopen(output_path,"a");

    for (int n = 0; n < 10; n++) {
        Node node = boundary->nodes[n + circle*10];

        fprintf(f, "%LF %LF %LF %LF %LF %LF", node.x, node.y, node.x_force, node.y_force, node.x_ref, node.y_ref);

        int x_int = debug_index(node.x, COORD_X);
        int y_int = debug_index(node.y, COORD_Y);
        int z_int = debug_index(node.z, COORD_Z);

        for(int i = x_int; i <= x_int + 1; ++i)
        {
            for(int j = y_int; j <= y_int + 1; ++j)
            {
                fprintf(f, " %LF %LF %LF %LF",
                        debug_coord(i, COORD_X),
                        debug_coord(j, COORD_Y),
                        force_X[i][j-1][z_int],
                        force_Y[i][j-1][z_int]);
            }
        }
        fprintf(f, "\n");

    }
    fclose(f);

}

void Output::print_boundary_vtk(int iter, ImmersedBoundary *boundary)
{
    char output_path[20];
    sprintf(output_path, "boundary%d.vtk", iter);
    FILE *f = fopen(output_path,"a");

    print_vtk_unstructured_header(output_path, boundary->nodes_count);

    for(int n = 0; n < boundary->nodes_count; ++n)
    {
        fprintf(f, "%lf %lf %lf\n", boundary->nodes[n].x, boundary->nodes[n].y, boundary->nodes[n].z);
    }
    fclose(f);

    print_vtk_streamline_vector_header(output_path, "force", boundary->nodes_count, boundary->nodes_count, boundary->nodes_count);

    f = fopen(output_path,"a");
    for(int n = 0; n < boundary->nodes_count; ++n)
    {
        fprintf(f, "%lf %lf %lf\n", boundary->nodes[n].x_force, boundary->nodes[n].y_force, boundary->nodes[n].z_force);
    }
    fclose(f);
}
