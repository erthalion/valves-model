import numpy as np

def take_one_layer(array):
    rows = []
    for row in array[0,::]:
        row = np.logical_and(
                np.roll(row, 1),
                np.roll(row, -1)
                )
        rows.append(row)

    array = np.array(rows)
