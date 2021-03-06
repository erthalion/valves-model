# Hydrodynamic modeling valves

Automatically exported from code.google.com/p/valves-model . 
Initial (and mostly for testing purposes) version of stationary proflem of fluid flow in large blood vessels.
Significant part of this code is belongs to Geidarov N. A.

## Project structure:
* build
    this directory contains binary files and all that is need to run program (mask files, config files, group files etc.)

* CMakeFiles
    cmake directory for their files

* cmake_install.cmake, CMakeCache.txt, Makefile
    this files is needed for project build

* lib
    this directory contains all requirements libraries for project

* src
    code directory

* CMakeList.txt
    configuration file for cmake


## Source structure:

* area_builder
    this directory contains all python scripts, that is need to create mask files for calculation
    they have use numpy(http://www.numpy.org/) for creation of the arrays
    build_area.py is main script
    other subdirectories (cube, empty) contains script with same names for different area forms
    pressure_mask.py is used for creation pressure mask
    u_mask.py is used for creation velocity mask (x component of velocity vector)
    v_mask.py is used for creation velocity mask (y component of velocity vector)
    w_mask.py is used for creation velocity mask (z component of velocity vector)

* db
    this directory contains python script for upload calculation result into cloud mongodb storage
    at this moment there is mongolab & mongohq

* lib
    this directory contains source of libraries, which is used in project
    they used in source form for cross-platform
    at this moment there is inih library for parsing ini config files (http://code.google.com/p/inih/)
    structure of this library was slightly modified and it used as shared library

* visualize
    this directory contains all python scripts for visualization of the calculation results with mayavi2 (http://code.enthought.com/projects/mayavi/)
    mayavi is used as wrapper on VTK

* main.cpp
    main source file, which contains numerical algoryphm on c++

## Usefull commands:
```bash
Update on cluster with rsync (nusc it's short name from ssh config)
rsync -av valves-model/ dadolgov@nusc:~/valves-model --exclude .git --exclude results --exclude groups
```

## Requirements:
* inih http://code.google.com/p/inih/ New BSD License
* numpy http://www.numpy.org/ BSD License
* mayavi2 http://code.enthought.com/projects/mayavi/ BSD License
* VTK 5 http://www.vtk.org/ BSD License
* Cmake http://www.cmake.org/ New BSD License
