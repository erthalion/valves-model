#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <csignal>
#include <cmath>
#include <stdexcept>
#include <stdio.h>

using namespace std;

#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define CB(x) ((x) * (x) * (x))

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
        long double height;
        long double stiffness;
        Node *nodes;

        ImmersedBoundary() {};
        virtual long double get_area() = 0;
};


class RectangleBoundary: public ImmersedBoundary
{
    public:
        RectangleBoundary()
        {
            this->nodes_count = 36;
            this->stiffness = 900;
            this->nodes = new Node[this->nodes_count];

            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    this->nodes[j + i*6].x = 1;
                    this->nodes[j + i*6].x_ref = 1;
                    this->nodes[j + i*6].x_vel = 0;
                    this->nodes[j + i*6].x_force = 0;

                    this->nodes[j + i*6].y = 0.2 + i/5.0 * 0.4;
                    this->nodes[j + i*6].y_ref = 0.2 + i/5.0 * 0.4;
                    this->nodes[j + i*6].y_vel = 0;
                    this->nodes[j + i*6].y_force = 0;

                    this->nodes[j + i*6].z = 0.2 + j/5.0 * 0.4;
                    this->nodes[j + i*6].z_ref = 0.2 + j/5.0 * 0.4;
                    this->nodes[j + i*6].z_vel = 0;
                    this->nodes[j + i*6].z_force = 0;
                }
            }
        }

        long double get_area()
        {
            return 0.4/5.0 * 0.4/5.0;
        }
};


/*
 * Use 3d figure instead 2d surface as immersed boundary
 */

class CylinderBoundary: public ImmersedBoundary
{
    public:
        long double x_center;
        long double y_center;
        long double z_start;

        CylinderBoundary()
        {
            this->nodes_count = 220;
            this->stiffness = 2800;
            this->nodes = new Node[this->nodes_count];
            this->radius = 0.2;
            this->height = 0.4;
            this->x_center = 1.0;
            this->y_center = 0.4;
            this->z_start = 0.2;

            for (int i = 0; i < 10; i++) {
                for (int j = 0; j < 20; j++) {
                    long double x = this->x_center + this->radius * sin(double(j)/19.0 * 2 * M_PI);
                    long double x_init = this->x_center + this->radius * 0.9 * sin(double(j)/19.0 * 2 * M_PI);
                    this->nodes[j+i*20].x = x;
                    this->nodes[j+i*20].x_ref = x;
                    this->nodes[j+i*20].x_vel = 0;
                    this->nodes[j+i*20].x_force = 0;

                    long double y = this->y_center + this->radius * cos(double(j)/19.0 * 2 * M_PI);
                    long double y_init = this->y_center + this->radius * 0.9 * cos(double(j)/19.0 * 2 * M_PI);
                    this->nodes[j+i*20].y = y;
                    this->nodes[j+i*20].y_ref = y;
                    this->nodes[j+i*20].y_vel = 0;
                    this->nodes[j+i*20].y_force = 0;

                    long double z = this->z_start + this->height * double(i) / 9.0 ;
                    this->nodes[j+i*20].z = z;
                    this->nodes[j+i*20].z_ref = z;
                    this->nodes[j+i*20].z_vel = 0;
                    this->nodes[j+i*20].z_force = 0;
                }
            }

            // front/end
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 10; j++) {
                    long double x = this->x_center + this->radius/2.0 * sin(double(j)/8.0 * 2 * M_PI);
                    this->nodes[200+j+i*10].x = x;
                    this->nodes[200+j+i*10].x_ref = x;
                    this->nodes[200+j+i*10].x_vel = 0;
                    this->nodes[200+j+i*10].x_force = 0;

                    long double y = this->y_center + this->radius/2.0 * cos(double(j)/8.0 * 2 * M_PI);
                    this->nodes[200+j+i*10].y = y;
                    this->nodes[200+j+i*10].y_ref = y;
                    this->nodes[200+j+i*10].y_vel = 0;
                    this->nodes[200+j+i*10].y_force = 0;

                    long double z = this->z_start + this->height * double(i) / 1.0 ;
                    this->nodes[200+j+i*10].z = z;
                    this->nodes[200+j+i*10].z_ref = z;
                    this->nodes[200+j+i*10].z_vel = 0;
                    this->nodes[200+j+i*10].z_force = 0;
                }
            }


        }

        long double get_area()
        {
            return 2 * M_PI * this->radius * this->height / this->nodes_count;
        }
};


class SphereBoundary: public ImmersedBoundary
{
    public:
        long double x_center;
        long double y_center;
        long double z_center;

        SphereBoundary()
        {
            this->nodes_count = 100;
            this->stiffness = 1300;
            this->nodes = new Node[this->nodes_count];
            this->radius = 0.2;
            this->x_center = 1.0;
            this->y_center = 0.4;
            this->z_center = 0.4;

            for (int theta = 0; theta < 10; theta++) {
                for (int phi = 0; phi < 10; phi++) {
                    long double x = this->x_center + this->radius * cos(double(theta)/9.0 * 2 * M_PI) * sin(double(phi)/9.0 * M_PI);
                    this->nodes[phi + theta*10].x = x;
                    this->nodes[phi + theta*10].x_ref = x;
                    this->nodes[phi + theta*10].x_vel = 0;
                    this->nodes[phi + theta*10].x_force = 0;

                    long double y = this->y_center + this->radius * sin(double(theta)/9.0 * 2 * M_PI) * sin(double(phi)/9.0 * M_PI);
                    this->nodes[phi + theta*10].y = y;
                    this->nodes[phi + theta*10].y_ref = y;
                    this->nodes[phi + theta*10].y_vel = 0;
                    this->nodes[phi + theta*10].y_force = 0;

                    long double z = this->z_center + this->radius * cos(double(phi)/9.0 * M_PI);
                    this->nodes[phi + theta*10].z = z;
                    this->nodes[phi + theta*10].z_ref = z;
                    this->nodes[phi + theta*10].z_vel = 0;
                    this->nodes[phi + theta*10].z_force = 0;
                }
            }
        }

        long double get_area()
        {
            return 4 * M_PI * this->radius * this->radius/ this->nodes_count;
        }
};


class ValvesBoundary: public ImmersedBoundary
{
    public:
        ValvesBoundary()
        {
            this->nodes_count = 200;
            this->radius = 0.2;
            this->height = 0.8;
            this->stiffness = 400000;
            this->nodes = new Node[this->nodes_count];

            //In 5 half circles by 10 nodes
            for(int n = 0; n < 10; ++n) {
                for (int i = 0; i < 10; i++) {
                    long double x = 0.8 + this->radius * sin(-M_PI * 0.5 * (double)(9-i)/9);
                    this->nodes[i+n*20].x = x;
                    this->nodes[i+n*20].x_ref = x;
                    this->nodes[i+n*20].x_vel = 0;
                    this->nodes[i+n*20].x_force = 0;

                    long double y = 0.8 - this->radius * cos(-M_PI * 0.5 * (double)(9-i)/9);
                    this->nodes[i+n*20].y = y;
                    this->nodes[i+n*20].y_ref = y;
                    this->nodes[i+n*20].y_vel = 0;
                    this->nodes[i+n*20].y_force = 0;

                    long double z = 0 + this->height * n / 9 ;
                    this->nodes[i+n*20].z = z;
                    this->nodes[i+n*20].z_ref = z;
                    this->nodes[i+n*20].z_vel = 0;
                    this->nodes[i+n*20].z_force = 0;
                }

                for (int i = 0; i < 10; i++) {
                    long double x = 0.8 + this->radius * sin(-M_PI * 0.5 * (double)(9-i)/9);
                    this->nodes[i+n*20+10].x = x;
                    this->nodes[i+n*20+10].x_ref = x;
                    this->nodes[i+n*20+10].x_vel = 0;
                    this->nodes[i+n*20+10].x_force = 0;

                    long double y = 0 + this->radius * cos(-M_PI * 0.5 * (double)(9-i)/9);
                    this->nodes[i+n*20+10].y = y;
                    this->nodes[i+n*20+10].y_ref = y;
                    this->nodes[i+n*20+10].y_vel = 0;
                    this->nodes[i+n*20+10].y_force = 0;

                    long double z = 0 + this->height * n / 9 ;
                    this->nodes[i+n*20+10].z = z;
                    this->nodes[i+n*20+10].z_ref = z;
                    this->nodes[i+n*20+10].z_vel = 0;
                    this->nodes[i+n*20+10].z_force = 0;
                }
            }

        }

        long double get_area()
        {
            return 2 * M_PI * this->radius * this->height / this->nodes_count;
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
