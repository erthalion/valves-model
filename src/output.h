#ifndef OUTPUT_H
#define OUTPUT_H

#include "utils.hpp"
#include "ib.h"

class Output {
    int Nx, Ny, Nz, dNz;
    long double *Cx, *Cy, *Cz;
    long double ***force_X, ***force_Y, ***force_Z;
    long double *Hx, *Hy, *Hz;
    long double ***U;
    long double ***R;
    int ***G;

    public:
    Output();
    Output(int Nx, int Ny, int Nz, int dNz,
            long double ***U,
            long double *Cx, long double *Cy, long double *Cz,
            long double ***force_X, long double ***force_Y, long double ***force_Z,
            long double ***R,
            long double *Hx, long double *Hy, long double *Hz,
            int ***G);
    ~Output();

    void print_info(int iters, long double R0, long double Rn);
    void print_texplot(int iters);
    void print_texplot_matrix(int iters);
    void print_vtk_header(char *output_path, int sizeX, int sizeY, int sizeZ);
    void print_vtk_unstructured_header(char *output_path, int size);
    void print_vtk_data_header(char *output_path, int sizeX, int sizeY, int sizeZ);
    void print_area();
    void print_vtk_header_points(char *output_path, int sizeX, int sizeY, int sizeZ);
    void print_area_points();
    void print_vtk_streamline_header(char *output_path, int sizeX, int sizeY, int sizeZ);
    void print_vtk(int iteration);
    void print_vtk_streamline_vector_header(char *output_path, const char *name, int sizeX, int sizeY, int sizeZ);
    void print_vtk_streamline_scalar_header(char *output_path, const char *name,  int sizeX, int sizeY, int sizeZ);
    void print_pressure();
    void print_boundary(int iter, ImmersedBoundary *boundary, long double ***U, int dNz);
    void print_boundary_vtk(int iter, ImmersedBoundary *boundary);
    void dump_to_file(long double ***data, long data_size, const char *file_name);
    long double*** load_dump(long data_size, const char *file_name);
};

#endif
