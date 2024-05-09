'''
    A determinant of a VHM-7 matrix, that is, a matrix of form like this:

    0 0 0 0 0 #
    0 0 0 0 # #
    0 0 0 # # 0
    0 0 # # 0 0
    0 # # 0 0 0
    # # # # # #
    
    where # can be any sensible elements (s. t. the determinant can be found), and the others are strictly zero.

'''

import sys
import numpy as np
import random


n = int(sys.argv[-1])
matrix = np.zeros((n, n), dtype='int64')


for i in range(n):
    matrix[n - 1][i] = random.randint(1, 10)

    if i != n - 1:
        matrix[i][n - i - 1] = random.randint(1, 10)
        matrix[i + 1][n - i - 1] = random.randint(1, 10)


def det_vhm7_structure_preserving(matrix):
    # Calculate cumulative prods of a_i
    prods_a = [matrix[n - i - 1][i] for i in range(1, n)]
    for i in range(n-3, -1, -1):
        prods_a[i] *= prods_a[i + 1]
    
    prod_b = 1
    s = 0
    for i in range(n-1):
        s += ((-1) ** i) * matrix[n - 1][i] * prods_a[i] * prod_b
        prod_b *= matrix[n - i - 2][i]
    
    return ((-1) ** (n // 2)) * (s + ((-1) ** (n - 1)) * matrix[n - 1][n - 1] * prod_b)


def det_vhm7_structure_modifying(matrix):
    # Calculate cumulative prods of a_i
    for i in range(n-2, 0, -1):
        matrix[n - i - 1][i] *= matrix[n - i - 2][i + 1]

    prod_b = 1
    s = 0
    for i in range(n-1):
        s += ((-1) ** i) * matrix[n - 1][i] * matrix[n - i - 2][i + 1] * prod_b
        prod_b *= matrix[n - i - 2][i]
    
    return ((-1) ** (n // 2)) * (s + ((-1) ** (n - 1)) * matrix[n - 1][n - 1] * prod_b)


det_numpy = np.linalg.det(matrix)

det_sp = det_vhm7_structure_preserving(matrix)
det_sm = det_vhm7_structure_modifying(matrix)

print(f"VHM-7 regular det: {det_numpy}")
print(f"VHM-7 structure-preserving det: {det_sp}")
print(f"VHM-7 structure-modifying det : {det_sm}")
print(f"{'Equal' if det_sp == det_sm else 'not equal !!!!!!!!!!!!!!!!!!'}")