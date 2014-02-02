#include "groups.h"
#include <math.h>
#include <csignal>
#include <omp.h>


GroupsGenerator::GroupsGenerator(int Nx, int Ny, int Nz)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;

    Utils *utils = new Utils(Nx, Ny, Nz);
    this->groups.resize(Nx,Ny,Nz);
    this->group_carg.resize(Nx,Ny,Nz);
    this->args.resize(Nx,Ny,Nz, 25);
    this->func.resize(Nx,Ny,Nz);
}

Array<int, 3> GroupsGenerator::get_carg()
{
    return this->group_carg.copy();
}

Array<matrix_ind, 4> GroupsGenerator::get_args()
{
    return this->args.copy();
}

GroupsGenerator::~GroupsGenerator()
{
    //utils->del(groups,Nx,Ny,Nz);
    //utils->del(arg,Nx,Ny,Nz);
    //utils->del(func,Nx,Ny,Nz);
    //utils->del(group_carg,Nx,Ny,Nz);
}

double GroupsGenerator::random_gr()
{
    int full = Nx*Ny*Nz;
    int number = rand();
    double z = (double)((unsigned int) ((double)number /
                (double) RAND_MAX * (double)full) + 1)/((double)full);
    return z;
}

void GroupsGenerator::add_gr(Array<matrix_ind, 4> &m, int &n, int i, int j, int k, int l, int s, int t)
{
    m(i, j, k, n).i = l;
    m(i, j, k, n).j = s;
    m(i, j, k, n).k = t;
    ++n;
}

void GroupsGenerator::init_gr()
{
    printf("init groups\n");
    num = -1;
    srand( (unsigned)time( NULL ) );

    Array<long double, 3> y(Nx, Ny, Nz);

    for(int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
            {
                group_carg(i, j, k) = 0;
                args(i, j, k, 0).i = 0;
                args(i, j, k, 0).j = 0;
                args(i, j, k, 0).k = 0;
                func(i, j, k, 0).i = 0;
                func(i, j, k, 0).j = 0;
                func(i, j, k, 0).k = 0;
                groups(i, j, k) = -1;
            };

    // basic
    int repeat = 5;
    for(int il=0; il<repeat; ++il)
    {
        printf("%4d\n",il);
        for(int a=0; a<Nx; ++a)
            for(int b=0; b<Ny; ++b)
                for(int c=0; c<Nz; ++c)
                    y(a, b, c) = random_gr();

        #pragma omp parallel
        {

        Array<long double, 3> e(Nx, Ny, Nz);

        #pragma omp for
        for(int i=0; i<Nx; ++i)
        {
            for(int j=0; j<Ny; ++j)
            {
                for(int k=0; k<Nz; ++k)
                {
                    e(i, j, k) = 1;
                    for(int l=0; l<Nx; ++l)
                    {
                        for(int s=0; s<Ny; ++s)
                        {
                            for(int t=0; t<Nz; ++t)
                            {
                                if ((fabs( (*operator_nonlin)(e,y,l,s,t)) > 1e-20)||
                                        (fabs( (*operator_nonlin)(y,e,l,s,t)) > 1e-20)||
                                        (fabs( (*operator_lin)(e,l,s,t)) > 1e-20))
                                {
                                    bool d = true;
                                    for(int c=0; (c<group_carg(i, j, k))&&d; ++c)
                                        if((args(i, j, k, c).i==l)&&
                                                (args(i, j, k, c).j==s)&&
                                                (args(i, j, k, c).k==t)) d=false;
                                    if(d)
                                    {
                                        add_gr(args,group_carg(i, j, k),i,j,k,l,s,t);
                                    }
                                }
                            }
                        }
                    }
                    e(i, j, k) = 0;
                }
            }
        }

        //utils->del(e,Nx,Ny,Nz);
        }
    }

    //utils->del(y,Nx,Ny,Nz);
    printf("groups has been initialized\n");
}

void GroupsGenerator::load_groups()
{
    printf("loading groups\n");
    FILE *f = fopen("oper_arg.txt","r");
    
    if (!f)
    {
        printf("oper_arg.txt not found!\n");
        exit(1);
    }

    for(int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
            {
                int ii,jj,kk;
                fscanf(f, "%3d %3d %3d", &ii, &jj, &kk);
                fscanf(f, "%3d\n", &group_carg(ii, jj, kk));
                for(int cc = 0; cc < group_carg(ii, jj, kk); ++cc)
                    fscanf(f, "%4d %4d %4d\n", &args(ii, jj, kk, cc).i, &args(ii, jj, kk, cc).j, &args(ii, jj, kk, cc).k);
            };
    fclose(f);

    FILE *fgr = fopen("gr.txt","r");
    if (!fgr)
    {
        printf("gr.txt not found!\n");
        exit(1);
    }

    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                fscanf(fgr,"%3d ", &groups(i, j, k));
            }
        }
        fprintf(fgr,"\n");
    }
    fclose(fgr);
    printf("groups have been loaded\n");
}
    //y = utils->alloc_and_fill<long double>(Nx,Ny,Nz);
    //st_id = None

void GroupsGenerator::print_gr()
{
    FILE *f = fopen("oper_arg.txt","w");
    for(int i=0; i<Nx; ++i)
        for(int j=0; j<Ny; ++j)
            for(int k=0; k<Nz; ++k)
            {
                // current indexes (i,j,k) and count of chanched indexes(group_carg)
                fprintf(f, "%3d %3d %3d %3d\n", i, j, k, group_carg(i, j, k));
                for(int c = 0; c<group_carg(i, j, k); ++c)
                    fprintf(f, "%4d %4d %4d\n", args(i, j, k, c).i, args(i, j, k, c).j, args(i, j, k, c).k);
            };
    fclose(f);

    FILE *fgr = fopen("gr.txt","w");
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                fprintf(fgr,"%3d ",groups(i, j, k));
            }
        }
        fprintf(fgr,"\n");
    }
    fclose(fgr);
}

void GroupsGenerator::set_gr()
{
    //int ***field;
    //field = utils->alloc_and_fill<int>(Nx,Ny,Nz);
    Array<int, 3> field(Nx, Ny, Nz);
    bool cells; // наличие компоненты, которой не назначена группа
    do
    {
        field = 0;
        //for(int i=0; i<Nx; ++i)
        //{
            //for(int j=0; j<Ny; ++j)
            //{
                //for(int k=0; k<Nz; ++k)
                //{
                    //field(i, j, k) = 0;
                //}
            //}
        //}

        cells = false;
        for(int i=0; i<Nx; ++i)
        {
            for(int j=0; j<Ny; ++j)
            {
                for(int k=0; k<Nz; ++k)
                {
                    if(groups(i, j, k)<0)
                    {
                        if(!cells)
                        {
                            cells = true;
                            ++num;
                        }
                        int count = group_carg(i, j, k);
                        Array<indexes, 1> t = args(i, j, k, Range::all());
                        bool push = true;
                        for(int c=0; (c<count)&&push; ++c)
                        {
                            if( field(t(c).i, t(c).j, t(c).k) ) push = false;
                        }
                        if(push)
                        {
                            groups(i, j, k) = num;
                            for(int c=0; c<count; ++c)
                                field(t(c).i, t(c).j, t(c).k) = 1;
                        }
                    }
                }
            }
        }
    }
    while(cells);
    ++num;
    //utils->del(field,Nx,Ny,Nz);
}
