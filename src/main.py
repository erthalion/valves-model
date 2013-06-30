#!/usr/bin/python
from numpy.ctypeslib import load_library
from numpyctypes import c_ndarray

valves_lib = load_library(
    'libValveModel.Debug.so', '.')       # '.' is the directory of the C++ lib

def main():
    valves_lib.main()


if __name__ == '__main__':
    main()
