#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <csignal>
#include <cmath>
#include <stdexcept>
#include <stdio.h>
#include <vector>

#include "ib.h"

using namespace std;

#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define CB(x) ((x) * (x) * (x))

struct indexes
{
    int i,j,k;
};
typedef indexes matrix_ind[25];

/*
 * Proxy must be used http://stackoverflow.com/a/2216055/1423473
 */
template <class T>
class matrix3D { 
    std::vector<T> data;
    int x_vector;
    int y_vector;
public:
    T &operator()(int x, int y, int z) {
        return this->data[x*x_vector*y_vector + y*y_vector + x];
    }

    matrix3D(int x, int y, int z) : data(x*y*z), x_vector(x), y_vector(y) {}
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

    /*
     * For contiguously memory allocation
     * T* buff = new T []x * y * z];
     * a[i][j] = (T*)(buff + i*y*z + j*z);
     */
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

    /*
     * For contiguously memory allocation
     * T* buff = new T []x * y * z];
     * a[i][j] = (T*)(buff + i*y*z + j*z);
     */
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

    ImmersedBoundary* get_immersed_boundary(string name)
    {
        if(name == "RectangleBoundary")
        {
            return new RectangleBoundary();
        }
        else if(name == "ValvesBoundary")
        {
            return new ValvesBoundary();
        }
        else if(name == "SphereBoundary")
        {
            return new SphereBoundary();
        }
        else if(name == "CylinderBoundary")
        {
            return new CylinderBoundary();
        }
        else
        {
            throw runtime_error("Incorrect boundary type");
        }
    }
};


#endif
