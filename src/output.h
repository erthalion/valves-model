#ifndef OUTPUT_H
#define OUTPUT_H

#include "utils.hpp"
#include "ib.h"

#include <blitz/array.h>
using namespace blitz;

class Output {
    int Nx, Ny, Nz, dNz;
    Array<long double, 1> Cx, Cy, Cz;
    Array<long double, 3> force_X, force_Y, force_Z;
    Array<long double, 1> Hx, Hy, Hz;
    Array<long double, 3> U;
    Array<long double, 3> R;
    Array<int, 3> G;

    public:
    Output();
    Output(int Nx, int Ny, int Nz, int dNz,
            Array<long double, 3> U,
            Array<long double, 1> Cx, Array<long double, 1> Cy, Array<long double, 1> Cz,
            Array<long double, 3> force_X, Array<long double, 3> force_Y, Array<long double, 3> force_Z,
            Array<long double, 3> R,
            Array<long double, 1> Hx, Array<long double, 1> Hy, Array<long double, 1> Hz,
            Array<int, 3> G);
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
    void print_boundary(int iter, ImmersedBoundary *boundary, Array<long double, 3> U, int dNz);
    void print_boundary_vtk(int iter, ImmersedBoundary *boundary);
    void dump_to_file(Array<long double, 3> data, long data_size, const char *file_name);
    Array<long double, 3> load_dump(long data_size, const char *file_name);
};

#endif
