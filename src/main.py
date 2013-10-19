#!/usr/bin/python
from numpy.ctypeslib import load_library
from numpyctypes import c_ndarray


def main():
    with open('.build-info', 'r') as build_info:
        lib_name = build_info.read().strip()
        valves_lib = load_library(lib_name, '.')
        valves_lib.main()


if __name__ == '__main__':
    main()
