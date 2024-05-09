'''
    A characteristic polynomial of an HM-7 matrix, that is, a matrix of form like this:

    # # # # # #
    # # 0 0 0 0
    0 # # 0 0 0
    0 0 # # 0 0
    0 0 0 # # 0
    0 0 0 0 # #
    
    where # can be any sensible elements (s. t. the characteristic polynomial can be found), and the others are strictly zero.

'''

import sys
import random
from sympy import Matrix, Symbol

if __name__ == '__main__':
    N = int(sys.argv[-1]) if __name__ == '__main__' else 0


    matrix = Matrix([[0 for _ in range(N)] for _ in range(N)])
    for i in range(N):
        matrix[0, i] = random.random() * 10

        if i != 0:
            matrix[i, i] = random.random() * 10
        
        if i != N - 1:
            matrix[i + 1, i] = random.random() * 10


def charpoly_hm7_structure_preserving(matrix):
    n = matrix.rows
    Lambda = Symbol('L')

    if n == 1:
        return Lambda - matrix[0, 0]
    
    # Calculate cumulative prods of a_i
    prods_a = [Lambda - matrix[i, i] for i in range(1, n)]
    for i in range(n-3, -1, -1):
        prods_a[i] *= prods_a[i + 1]
    
    matrix[0, 0] = Lambda - matrix[0, 0]
    prod_b = -matrix[1, 0]
    s = matrix[0, 0] * prods_a[0]

    for i in range(1, n-1):
        s += ((-1) ** i) * (-matrix[0, i]) * prods_a[i] * prod_b
        prod_b *= -matrix[i + 1, i]

    matrix[0, 0] = Lambda - matrix[0, 0]

    return s + ((-1) ** (n - 1)) * (-matrix[0, n-1]) * prod_b


def charpoly_hm7_structure_modifying(matrix):
    n = matrix.rows
    Lambda = Symbol('L')

    if n == 1:
        return Lambda - matrix[0, 0]

    matrix[0, 0] = Lambda - matrix[0, 0]
    matrix[n - 1, n - 1] = Lambda - matrix[n - 1, n - 1]
    
    # Calculate cumulative prods of a_i
    for i in range(n-2, 0, -1):
        matrix[i, i] = (Lambda - matrix[i, i]) * matrix[i + 1, i + 1]

    prod_b = -matrix[1, 0]
    s = matrix[0, 0] * matrix[1, 1]
    for i in range(1, n-1):
        s += ((-1) ** i) * (-matrix[0, i]) * matrix[i + 1, i + 1] * prod_b
        prod_b *= -matrix[i + 1, i]

    return s + ((-1) ** (n - 1)) * (-matrix[0, n - 1]) * prod_b


if __name__ == '__main__':
    charpoly_sympy = matrix.charpoly('L')
    charpoly_sp = charpoly_hm7_structure_preserving(matrix)
    charpoly_sm = charpoly_hm7_structure_modifying(matrix)


    print(f"HM-7 regular characteristic polynomial             : {charpoly_sympy.as_expr().expand()}")
    print(f"HM-7 structure-preserving characteristic polynomial:     {charpoly_sp.expand()}")
    print(f"HM-7 structure-modifying characteristic polynomial :     {charpoly_sm.expand()}")