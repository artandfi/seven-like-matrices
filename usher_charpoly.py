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

from sympy import Matrix, Symbol, Eq, solve
from hm7_charpoly import charpoly_hm7_structure_preserving, charpoly_hm7_structure_modifying


L = Symbol('L')
n = 6

matrix = Matrix([
    [0.72, 0,    0,    3.6 * (L - 1), 5.1 * (L - 1), 7.5 * L],
    [0.28, 0.69, 0,    0,             0,             0      ],
    [0,    0.31, 0.75, 0,             0,             0      ],
    [0,    0,    0.25, 0.77,          0,             0      ],
    [0,    0,    0,    0.23,          0.63,          0      ],
    [0,    0,    0,    0,             0.37,          0      ],
])


charpoly_sympy = matrix.charpoly(L)
charpoly_sympy = charpoly_sympy.subs("_L", L)
char_eq = Eq(charpoly_sympy, 0)
eigenvalues_sympy = solve(char_eq, L)

print(eigenvalues_sympy)


charpoly_sp = charpoly_hm7_structure_preserving(matrix)
charpoly_sm = charpoly_hm7_structure_modifying(matrix)


print(f"HM-7 regular characteristic polynomial: {charpoly_sympy.as_expr().expand()}")
print(f"HM-7 structure-preserving characteristic polynomial: {charpoly_sp.expand()}")
print(f"HM-7 structure-modifying characteristic polynomial: {charpoly_sm.expand()}")

char_eq_sympy = Eq(charpoly_sympy, 0)
eigenvalues_sympy = solve(char_eq_sympy, L)

optimal_lambda_sympy = max(e for e in eigenvalues_sympy if e.is_real and e > 1)
print(f"Optimal eigenvalue: {optimal_lambda_sympy}")

char_eq_sp = Eq(charpoly_sp, 0)
eigenvalues_sp = solve(char_eq_sp, L)

optimal_lambda_sp = max(e for e in eigenvalues_sp if e.is_real and e > 1)
print(f"Optimal eigenvalue: {optimal_lambda_sp}")

char_eq_sm = Eq(charpoly_sm, 0)
eigenvalues_sm = solve(char_eq_sm, L)

optimal_lambda_sm = max(e for e in eigenvalues_sm if e.is_real and e > 1)
print(f"Optimal eigenvalue: {optimal_lambda_sm}")