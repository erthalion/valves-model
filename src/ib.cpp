#include "ib.h"

Node::Node()
{
    x = y = z = 0;
    x_ref = y_ref = z_ref = 0;
    x_force = y_force = z_force = 0;
    x_vel = y_vel = z_vel = 0;
}

ImmersedBoundary::ImmersedBoundary() {}

RectangleBoundary::RectangleBoundary()
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

long double RectangleBoundary::get_area()
{
    return 0.4/5.0 * 0.4/5.0;
}

CylinderBoundary::CylinderBoundary()
{
    this->nodes_count = 600;
    this->stiffness = 2800;
    this->nodes = new Node[this->nodes_count];
    this->radius = 0.3;
    this->height = 1.9;
    this->y_center = 0.4;
    this->z_center = 0.4;

    for (int i = 0; i < 30; i++) {
        for (int j = 0; j < 20; j++) {
            long double z = this->z_center + this->radius * sin(double(j)/19.0 * 2 * M_PI);
            this->nodes[j+i*20].z = z;
            this->nodes[j+i*20].z_ref = z;
            this->nodes[j+i*20].z_vel = 0;
            this->nodes[j+i*20].z_force = 0;

            long double y = this->y_center + this->radius * cos(double(j)/19.0 * 2 * M_PI);
            this->nodes[j+i*20].y = y;
            this->nodes[j+i*20].y_ref = y;
            this->nodes[j+i*20].y_vel = 0;
            this->nodes[j+i*20].y_force = 0;

            long double x = this->height * double(i) / 29.0 ;
            this->nodes[j+i*20].x = x;
            this->nodes[j+i*20].x_ref = x;
            this->nodes[j+i*20].x_vel = 0;
            this->nodes[j+i*20].x_force = 0;
        }
    }
}

long double CylinderBoundary::get_area()
{
    return 2 * M_PI * this->radius * this->height / this->nodes_count;
}

SphereBoundary::SphereBoundary()
{
    this->nodes_count = 100;
    this->stiffness = 1500;
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

long double SphereBoundary::get_area()
{
    return 4 * M_PI * this->radius * this->radius/ this->nodes_count;
}

ValvesBoundary::ValvesBoundary()
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

long double ValvesBoundary::get_area()
{
    return 2 * M_PI * this->radius * this->height / this->nodes_count;
}
