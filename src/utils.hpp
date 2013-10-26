#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cmath>

// Типы
struct indexes
{
    int i,j,k;
};
typedef indexes matrix_ind[25];

class Node
{
    public:
        double x,y,z;
        double x_ref, y_ref, z_ref;
        double x_force, y_force, z_force;
        double x_vel, y_vel, z_vel;

        Node()
        {
            x = y = z = 0;
            x_ref = y_ref = z_ref = 0;
            x_force = y_force = z_force = 0;
            x_vel = y_vel = z_vel = 0;
        }
};

class ImmersedBoundary
{
    public:
        int nodes_count;
        long double radius;
        long double stiffness;
        Node *nodes;

        ImmersedBoundary()
        {
            this->nodes_count = 36;
            this->radius = 0.3;
            this->stiffness = 0.1;
            this->nodes = new Node[this->nodes_count];
        
            for(int n = 0; n < this->nodes_count; ++n) {
                this->nodes[n].x = 0.5 + this->radius * sin(2. * M_PI * (double) n / this->nodes_count);
                this->nodes[n].x_ref = 0.5 + this->radius * sin(2. * M_PI * (double) n / this->nodes_count);
                this->nodes[n].x_vel = 0;
                this->nodes[n].x_force = 0;
        
                this->nodes[n].y = 0.5 + this->radius * cos(2. * M_PI * (double) n / this->nodes_count);
                this->nodes[n].y_ref = 0.5 + this->radius * cos(2. * M_PI * (double) n / this->nodes_count);
                this->nodes[n].y_vel = 0;
                this->nodes[n].y_force = 0;
        
                this->nodes[n].z = 0.5 + this->radius * cos(2. * M_PI * (double) n / this->nodes_count);
                this->nodes[n].z_ref = 0.5 + this->radius * cos(2. * M_PI * (double) n / this->nodes_count);
                this->nodes[n].z_vel = 0;
                this->nodes[n].z_force = 0;
            }

        }
};

class Utils
{
    int Nx, Ny, Nz;

public:
	Utils()
    {
        this->Nx = 10;
        this->Ny = 10;
        this->Nz = 10;
    }

	Utils(int Nx, int Ny, int Nz)
    {
        this->Nx = Nx;
        this->Ny = Ny;
        this->Nz = Nz;
    }

    ~Utils()
    {
    }

    template<typename T>
    T*** alloc(const int x, const int y, const int z)
    {
        T*** a = new T **[x];
        for(int i=0; i<x; ++i)
        {
            a[i] = new T *[y];
            for(int j=0; j<y; ++j)
            {
                a[i][j] = new T [z];
            }
        }

        return a;
    }

    template<typename T>
    T*** alloc_and_fill(const int x, const int y, const int z)
    {
        T*** a = new T **[x];
        for(int i=0; i<x; ++i)
        {
            a[i] = new T *[y];
            for(int j=0; j<y; ++j)
            {
                a[i][j] = new T [z];
                for(int k=0; k<z; ++k)
                {
                    a[i][j][k] = 0;
                }
            }
        }

        return a;
    }

    template<typename T>
    void del(T ***&a, int x, int y, int z)
    {
        for(int i=0; i<x; ++i)
        {
            for(int j=0; j<y; ++j)
                delete [] a[i][j];
            delete [] a[i];
        }
        delete [] a;
    }
};


#endif
