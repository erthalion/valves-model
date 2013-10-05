import numpy as np

def take_off_first_last(array):
    rows = []
    for row in array:
        if np.sum(row) != 0:
            row[row == True] = False
            break
        rows.append(row)

    array = np.array(rows)

def take_one_layer_vert(array):
    rows = []
    for row in array[:0:]:
        row = np.logical_and(
                np.roll(row, 1, axis=2),
                np.roll(row, -1, axis=2)
                )
        rows.append(row)

    array = np.array(rows)
