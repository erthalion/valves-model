#ifndef GROUPS_H
#define GROUPS_H

#include "utils.hpp"
#include <stdlib.h>
#include <stdio.h>

#include <blitz/array.h>

using namespace blitz;

class GroupsGenerator
{
   Array<int, 3> groups; 
   int num;
   Utils *utils;
   Array<matrix_ind, 4> args;
   Array<matrix_ind, 3> func;
   Array<int, 3> group_carg;
   int Nx, Ny, Nz;


   public:
       long double (*operator_nonlin)(Array<long double, 3> U1, Array<long double, 3> U2, int i, int j, int k);
       long double (*operator_lin)(Array<long double, 3> U, int i, int j, int k);

       GroupsGenerator(int Nx, int Ny, int Nz);
       ~GroupsGenerator();

       void set_gr();
       void print_gr();
       double random_gr();
       void add_gr(Array<matrix_ind, 4> &m, int &n, int i, int j, int k, int l, int s, int t);
       void init_gr();
       void load_groups();

       Array<int, 3> get_carg();// {return this->group_carg;}
       Array<matrix_ind, 4> get_args();// {return this->args;}
};

#endif
