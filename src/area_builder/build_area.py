#!/usr/bin/python
# -*- coding: utf8 -*-

import coord
import mask
import empty_pressure_mask
import empty_u_mask
import empty_v_mask
import empty_w_mask

if __name__ == '__main__':
    coord.build_area()
    mask.build_area()
    empty_pressure_mask.build_area()
    empty_u_mask.build_area()
    empty_v_mask.build_area()
    empty_w_mask.build_area()
