#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <csignal>
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
        long double x,y,z;
        long double x_ref, y_ref, z_ref;
        long double x_force, y_force, z_force;
        long double x_vel, y_vel, z_vel;

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
            this->nodes_count = 50;
            this->radius = 0.3;
            this->stiffness = 30000;
            this->nodes = new Node[this->nodes_count];
        
            // In 5 half circles by 10 nodes
            for(int n = 0; n < 5; ++n) {
                for (int i = 0; i < 10; i++) {
                    long double x = 1.0 + this->radius * sin(-2. * M_PI * (double) i / 19);
                    this->nodes[i+n*10].x = x;
                    this->nodes[i+n*10].x_ref = x;
                    this->nodes[i+n*10].x_vel = 0;
                    this->nodes[i+n*10].x_force = 0;

                    long double y = 0.5 + this->radius * cos(-2. * M_PI * (double) i / 19);
                    this->nodes[i+n*10].y = y;
                    this->nodes[i+n*10].y_ref = y;
                    this->nodes[i+n*10].y_vel = 0;
                    this->nodes[i+n*10].y_force = 0;

                    long double z = 0.3 + 0.4 * n / 5 ;
                    this->nodes[i+n*10].z = z;
                    this->nodes[i+n*10].z_ref = z;
                    this->nodes[i+n*10].z_vel = 0;
                    this->nodes[i+n*10].z_force = 0;
                }
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
