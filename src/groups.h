#ifndef GROUPS_H
#define GROUPS_H

#include "utils.hpp"
#include <stdlib.h>
#include <stdio.h>

class GroupsGenerator
{
   int ***groups; 
   int num;
   Utils *utils;
   matrix_ind ***arg;
   matrix_ind ***func;
   int ***carg;
   int Nx, Ny, Nz;


   public:
       long double (*operator_nonlin)(long double ***U1, long double ***U2, int i, int j, int k);
       long double (*operator_lin)(long double ***U, int i, int j, int k);

       GroupsGenerator(int Nx, int Ny, int Nz);
       ~GroupsGenerator();

       void set_gr();
       void print_gr();
       double random_gr();
       void add_gr(matrix_ind ***&m, int &n, int i, int j, int k, int l, int s, int t);
       void init_gr();
       void load_groups();

       int*** get_carg() {return carg;}
       matrix_ind*** get_arg() {return arg;}
};

#endif
