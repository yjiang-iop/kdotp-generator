# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
# Adapted from the kdotp-symmetry package by: Dominik Gresch <greschd@gmx.ch>  © 2017-2018, ETH Zurich, Institut für Theoretische Physik
"""
Utilities for handling algebraic expressions, such as turning them to vector or matrix form.
"""

import random
import operator
from functools import reduce
from itertools import combinations_with_replacement, product
import sympy as sp


def expr_to_vector(
    expr,
    basis,
    K_VEC = sp.symbols('kx, ky, kz'),
    *,
    random_fct=lambda: random.randint(-100, 100),
    **kwargs  # pylint: disable=unused-argument
):
    """
    Converts an algebraic (sympy) expression into vector form.

    :param expr: Algebraic expression   # R*f(k) 
    :type expr: sympy.Expr            

    :param basis: Basis of the vector space, w.r.t. which the vector will be expressed. # f(k)
    :type basis: list[sympy.Expr]

    :param random_fct: Function creating random numbers on which the expression will be evaluated.
    # YJ: change random number range from -100:100 to -1000:1000
    """
    dim = len(basis)
    # create random values for the coordinates and evaluate
    # both the basis functions and the expression to generate
    # the linear equation to be solved
    A = []
    b = []  # pylint: disable=invalid-name
    for _ in range(2 * dim):
        if sp.Matrix(A).rank() >= len(basis):
            break
        vals = [(k, random_fct()) for k in K_VEC]
        A.append([bs.subs(vals) for bs in basis])
        b.append(expr.subs(vals))
    else:
        # ============= added by YJ: generate again the random number ================ #
        A = []
        b = []
        for _ in range(2 * dim):
            if sp.Matrix(A).rank() >= len(basis):
                break
            vals = [(k, random_fct()) for k in K_VEC]
            A.append([bs.subs(vals) for bs in basis])
            b.append(expr.subs(vals))
        else:
        # =========================================================================== #
        # this could happen if the random_fct is bad, or the 'basis' is not linearly independent
            raise ValueError(
                'Could not find a sufficient number of linearly independent vectors'
            )
    res = sp.linsolve((sp.Matrix(A), sp.Matrix(b)))
    if len(res) != 1:
        raise ValueError(
            'Invalid result {res} when trying to match expression {expr} to basis {basis}.'
            .format(res=res, expr=expr, basis=basis)
        )
    vec = next(iter(res))
   #vec = tuple(v.nsimplify() for v in vec)
    vec = tuple(0 if abs(sp.nsimplify(v))<1e-7 else sp.nsimplify(v) for v in vec) # modified by YJ
    # check consistency
    if not expr.equals(sum(v * b for v, b in zip(vec, basis))):
        raise ValueError(
            "Vector {vec} in basis {basis} does not match expression {expr}".
            format(vec=vec, basis=basis, expr=expr)
        )
    return vec


def monomial_basis(degree_list, kp_variable = 'k'):
    """
    Returns the product basis of (kx, ky, kz), with monomials of the given degrees.

    :param degrees: Degree of the monomials. # Modified by YJ: degree_list=[d1,d2..], with each di denotes the degree of one type of kp_variable 
    :type degrees: list

    :param kp_variable: The variable of the hamiltonian, can be anyone of 'k', 'E', 'B', 'e', 'k E', 'k B', 'E B', 'k E B'
    :type kp_variable: :py:class:str

    Example:
        >>> import kdotp_symmetry as kp
        >>> kp.monomial_basis([1,2], kp_variable='k E')
        [Ex*kx**2, Ex*kx*ky, Ex*kx*kz, Ex*ky**2, Ex*ky*kz, Ex*kz**2, Ey*kx**2, Ey*kx*ky, Ey*kx*kz, \
            Ey*ky**2, Ey*ky*kz, Ey*kz**2, Ez*kx**2, Ez*kx*ky, Ez*kx*kz, Ez*ky**2, Ez*ky*kz, Ez*kz**2]

    """
    assert all(deg >= 0 for deg in degree_list), ('Degrees must be non-negative integers', degree_list)
    assert len(degree_list) == len(kp_variable.split()), ('Input degree_list and kp_variable do not have the same length!',degree_list, kp_variable)
    
    num_input = len(degree_list)
    kp_variable = kp_variable.split()
    Base_vec = ''
    for t in kp_variable:
        if t == 'k':
            Base_vec += 'kx ky kz '
        elif t == 'E':
            Base_vec += 'Ex Ey Ez '
        elif t == 'B':
            Base_vec += 'Bx By Bz '
        elif t == 'e':
            Base_vec += 'ex ey ez '
    Base_vec = sp.symbols(Base_vec)

    # for strain tensors, the degree should time 2
    degree = degree_list.copy() 
    for i in range(num_input):
        if Base_vec[3*i] == sp.symbols('ex'):
            degree[i] *= 2
    
    basis = []
    if num_input == 1:  
        monomial_tuples = combinations_with_replacement(Base_vec, degree[0])
        #e.g, deg=2: [(kx, kx), (kx, ky), (kx, kz), (ky, ky), (ky, kz), (kz, kz)]
        basis.extend(reduce(operator.mul, m, sp.Integer(1)) for m in monomial_tuples)
    else:
        basis = []
        for i, deg in enumerate(degree):
            variable = Base_vec[3*i: 3*(i+1)]
            monomial_tuples = combinations_with_replacement(variable, deg)
            basis.append([reduce(operator.mul, m, sp.Integer(1)) for m in monomial_tuples])
        basis = [ reduce(operator.mul, m, sp.Integer(1)) for m in product(*basis) ]
    return basis


def matrix_to_expr_operator(matrix_form, repr_has_cc=False, K_VEC=sp.symbols('kx, ky, kz')):
    """Returns a function that operates on expression, corresponding to the given ``matrix_form`` 
    which operates on a vector in real space. ``repr_has_cc`` determines whether the symmetry contains time reversal."""
    # k-form and r-form of the matrix are related by A -> A^-1^T (For Orthogonal matrix A, A^-1^T=A)
    # => matrix^T gives g^-1 in k-space coordinates (-1 in g^-1 gives A.T)
    # Change sign if the representation has complex conjugation

    k_matrix_form = sp.Matrix(matrix_form).T       # g^-1 acts on k. same for E,B,epsilon.
   #k_matrix_form = sp.Matrix(matrix_form).T.inv() # if g acts on k. Bilbao convertion. 

    def transform_matrix(matrix, K_VEC):
        if K_VEC == sp.symbols('kx, ky, kz'): # P: k->-k, T: k->-k
            if repr_has_cc:
                matrix *= -1
        elif K_VEC == sp.symbols('Bx, By, Bz'): # P: B->B, T: B->-B
            matrix *= sp.det(k_matrix_form)
            if repr_has_cc:
                matrix *= -1
        elif K_VEC == sp.symbols('ex, ey, ez'): # P: e->e, T: e->e
            matrix *= sp.det(k_matrix_form)
        elif K_VEC == sp.symbols('Ex, Ey, Ez'): # P: E->-E, T: E->E
             pass
        else:
            raise ValueError('Wrong K_VEC! K_VEC='+str(K_VEC))
        return matrix

    if len(K_VEC) == 3:
        k_matrix_form = transform_matrix(k_matrix_form, K_VEC)
    elif len(K_VEC) == 6: # in this case, k_mat is the direct sum of two mat
        mat1 = transform_matrix(k_matrix_form, K_VEC[0:3])
        mat2 = transform_matrix(k_matrix_form, K_VEC[3:6])
        k_matrix_form = sp.zeros(6)
        k_matrix_form[0:3,0:3] = mat1
        k_matrix_form[3:6,3:6] = mat2
    elif len(K_VEC) == 9: # in this case, k_mat is the direct sum of three mat
        mat1 = transform_matrix(k_matrix_form, K_VEC[0:3])
        mat2 = transform_matrix(k_matrix_form, K_VEC[3:6])
        mat3 = transform_matrix(k_matrix_form, K_VEC[6:9])
        k_matrix_form = sp.zeros(9)
        k_matrix_form[0:3,0:3] = mat1
        k_matrix_form[3:6,3:6] = mat2
        k_matrix_form[6:9,6:9] = mat3
    else:
        raise ValueError('Wrong K_VEC! K_VEC='+str(K_VEC))

    substitution = list(zip(K_VEC, k_matrix_form @ sp.Matrix(K_VEC))) # k_mat @ k

    def expr_operator(expr):
        return expr.subs(substitution, simultaneous=True)

    return expr_operator
