#!/usr/bin/python
# -*- coding: utf8 -*-

import coord
import mask
from empty import empty_pressure_mask, empty_u_mask, empty_v_mask, empty_w_mask

if __name__ == '__main__':
    coord.build_area()
    mask.build_area()
    empty_pressure_mask.build_area()
    empty_u_mask.build_area()
    empty_v_mask.build_area()
    empty_w_mask.build_area()
