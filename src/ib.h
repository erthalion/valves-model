#ifndef IB_H
#define IB_H

#include <cmath>

class Node
{
    public:
        long double x,y,z;
        long double x_ref, y_ref, z_ref;
        long double x_force, y_force, z_force;
        long double x_vel, y_vel, z_vel;

        Node();
};

class ImmersedBoundary
{
    public:
        int nodes_count;
        long double radius;
        long double height;
        long double stiffness;
        Node *nodes;

        ImmersedBoundary();
        virtual long double get_area() = 0;
};


class RectangleBoundary: public ImmersedBoundary
{
    public:
        RectangleBoundary();
        long double get_area();
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

        CylinderBoundary();
        long double get_area();
};


class SphereBoundary: public ImmersedBoundary
{
    public:
        long double x_center;
        long double y_center;
        long double z_center;

        SphereBoundary();
        long double get_area();
};


class ValvesBoundary: public ImmersedBoundary
{
    public:
        ValvesBoundary();
        long double get_area();
};

#endif
