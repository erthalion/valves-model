#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "lib/inih/INIReader.h"
#define M_2PI 2*M_PI

using std::string;

/*
 * Kinematical viscosity
 */
const long double nu = 1e-1;

/*
 * Density
 */
const long double rho = 0.1;

const long double eps = 0.01;

/*
 * Boundary conditions for pressure
 */
const long double p_left = 0.01;
const long double p_right = 0;

/*
 * Grid dimensions
 */
int Nx, Ny, Nz, dNz;

/*
 * Boundary class name
 */
string BoundaryClass;

/*
 * Generate groups
 */
bool generate_groups;

void load_config()
{
    INIReader config("area.config");

    Nx = config.GetInteger("Area", "Nx", 10);
    Ny = config.GetInteger("Area", "Ny", 10);
    dNz = config.GetInteger("Area", "dNz", 10);
    Nz = config.GetInteger("Area", "Nz", 10);

    BoundaryClass = config.Get("Boundary", "Type", "RectangleBoundary");

    generate_groups = config.GetBoolean("Main", "GenerateGroups", false);
}

#endif
