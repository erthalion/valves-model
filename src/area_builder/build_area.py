#!/usr/bin/python
# -*- coding: utf8 -*-

import mask
import coord
from empty import pressure_mask, u_mask, v_mask, w_mask

if __name__ == '__main__':
    coord.build_area()
    mask.build_area()
    pressure_mask.build_area()
    u_mask.build_area()
    v_mask.build_area()
    w_mask.build_area()
