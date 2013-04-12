#include "lib/inih/INIReader.h"
#define M_2PI 2*M_PI

/*
 * Variable type for U array
 */
const int VELOCITY_U = 0;
const int VELOCITY_V = 1;
const int VELOCITY_W = 2;
const int PRESSURE = 3;

/*
 * Kinematical viscosity
 */
const long double nu = 1e-2;

/*
 * Density
 */
const long double rho = 1;

const long double eps = 0.001;

/*
 * Boundary conditions for pressure
 */
const long double p_left = 1;
const long double p_right = 0;

/*
 * Grid dimensions
 */
int Nx, Ny, Nz, dNz;

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

    generate_groups = config.GetBoolean("Main", "GenerateGroups", false);
    printf("Generate groups is %s\n", (generate_groups) ? "true" : "false");

    printf("config has been loaded\n");
}
