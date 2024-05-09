'''
    A determinant of an HM-7 matrix, that is, a matrix of form like this:

    # # # # # #
    # # 0 0 0 0
    0 # # 0 0 0
    0 0 # # 0 0
    0 0 0 # # 0
    0 0 0 0 # #
    
    where # can be any sensible elements (s. t. the determinant can be found), and the others are strictly zero.

'''

import sys
import numpy as np
import random


n = int(sys.argv[-1])
matrix = np.zeros((n, n), dtype='int64')


for i in range(n):
    matrix[0][i] = random.randint(1, 10)

    if i != 0:
        matrix[i][i] = random.randint(1, 10)
    
    if i != n - 1:
        matrix[i + 1][i] = random.randint(1, 10)


def det_hm7_structure_preserving(matrix):
    # Calculate cumulative prods of a_i
    prods_a = [matrix[i][i] for i in range(1, n)]
    for i in range(n-3, -1, -1):
        prods_a[i] *= prods_a[i + 1]
    
    prod_b = 1
    s = 0
    for i in range(n-1):
        s += ((-1) ** i) * matrix[0][i] * prods_a[i] * prod_b
        prod_b *= matrix[i + 1][i]
    
    return s + ((-1) ** (n - 1)) * matrix[0][n - 1] * prod_b


def det_hm7_structure_modifying(matrix):
    # Calculate cumulative prods of a_i
    for i in range(n-2, 0, -1):
        matrix[i][i] *= matrix[i + 1][i + 1]

    prod_b = 1
    s = 0
    for i in range(n-1):
        s += ((-1) ** i) * matrix[0][i] * matrix[i + 1][i + 1] * prod_b
        prod_b *= matrix[i + 1][i]
    
    return s + ((-1) ** (n - 1)) * matrix[0][n - 1] * prod_b


det_numpy = np.linalg.det(matrix)

det_sp = det_hm7_structure_preserving(matrix)
det_sm = det_hm7_structure_modifying(matrix)



print(f"HM-7 regular det: {det_numpy}")
print(f"HM-7 structure-preserving det: {det_sp}")
print(f"HM-7 structure-modifying det : {det_sm}")
print(f"{'Equal' if det_sp == det_sm else 'not equal !!!!!!!!!!!!!!!!!!'}")