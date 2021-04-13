# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
# Adapted from the kdotp-symmetry package by: Dominik Gresch <greschd@gmx.ch>  © 2017-2018, ETH Zurich, Institut für Theoretische Physik
"""
Defines functions to construct the basis of the symmetry-constrained Hamiltonian.
"""

import sympy as sp
from sympy.physics.quantum import TensorProduct
import numpy as np
from functools import reduce
import scipy.linalg as la

from ._expr_utils import monomial_basis, expr_to_vector, matrix_to_expr_operator
from ._repr_utils import hermitian_to_vector, hermitian_basis, repr_to_matrix_operator, check_orthogonal, frobenius_product, solve_linear_system_numpy
from ._repr_utils import hermitian_pauli_basis, hermitian_pauli_basis_symbols 
from ._linalg import intersection_basis, nullspace_blocked
from ._to_matrix import to_matrix
from ._logging_setup import LOGGER
from ._decompose_kp import decompose_kp


def symmetric_hamiltonian(
    symmetry_operations,   
    kp_variable = 'k',
    order = [0],
    repr_basis = 'pauli',
    msg_num = None,
    kvec = None,
):
    r"""
    Calculates the basis of the symmetric Hamiltonian for a given set of symmetry operations.

    :param symmetry_operations: The symmetry operations that the Hamiltonian should respect.
    :type symmetry_operations: :py:class: `dict` with keys 'rotation_matrix', 'repr_matrix', 'repr_has_cc'. 

    :param kp_variable: The variable of the hamiltonian, can be anyone of 'k', 'E', 'B', 'e', 'k E', 'k B', 'E B', 'k E B'
    :type kp_variable: :py:class:str

    :param order: The list of orders of the monomials. Each number in the list specifies the order of a variable.
    :type order: :py:class:`list` of :py:class:`int`

    :param repr_basis: The basis for the hermitian matrices, with the same size as the representations. 
                       By default, the :py:func:`.hermitian_pauli_basis` of the appropriate size is used.
    :type repr_basis: :py:class:`list` of :py:mod:`sympy` matrices
    
    :param msg_num & kvec: two string used to denote the magnetic space group and little group k, 
                             used to locate linear representations in order to decompose kp hamiltonian.
    :type msg_num & kvec: :py:class:str

    :returns: Basis for the symmetric Hamiltonian, as a :py:class:`list` of :py:mod:`sympy` matrix expressions.
    # Modified by YJ: if msg_num and kvec is specified, also return lists of decomposed repr and expr basis, otherwise return empty lists.
    """
    # for sympy or numpy matrices
    try:
        repr_matrix_size = symmetry_operations[0]['repr_matrix'].shape[0]
    # for plain lists -- this doesn't work for sympy matrices because
    # their 'len' is the total number of elements
    except AttributeError:
        repr_matrix_size = len(symmetry_operations[0]['repr_matrix'])

    repr_basis_type = 'pauli' if repr_basis == 'pauli' else None
    if repr_basis == 'auto':
        repr_basis = hermitian_basis(repr_matrix_size)
    elif repr_basis == 'pauli':
        repr_basis = hermitian_pauli_basis(repr_matrix_size)
        repr_basis_symbols = hermitian_pauli_basis_symbols(repr_matrix_size)

    if repr_basis not in ['auto', 'pauli']:
        check_orthogonal(repr_basis)
    
    Base_vec = ''
    for t in kp_variable.split():
        if t == 'k':
            Base_vec += 'kx ky kz '
        elif t == 'E':
            Base_vec += 'Ex Ey Ez '
        elif t == 'B':
            Base_vec += 'Bx By Bz '
        elif t == 'e':
            Base_vec += 'ex ey ez '
    Base_vec = sp.symbols(Base_vec)
    expr_basis = monomial_basis(order, kp_variable) 

    expr_dim = len(expr_basis)
    repr_dim = len(repr_basis)
    repr_basis_norm_squares = [frobenius_product(b, b) for b in repr_basis]
    full_dim = expr_dim * repr_dim
    full_basis = [
        sp.Matrix(x) for x in np.outer(expr_basis, repr_basis).
        reshape(full_dim, repr_matrix_size, repr_matrix_size).tolist()
    ]

    invariant_bases = []
    expr_mat_collection = []
    repr_mat_collection = []
    for isym_op, sym_op in enumerate(symmetry_operations):
        LOGGER.info('Calculating matrix form of expression.')
        expr_mat = to_matrix(
            operator=matrix_to_expr_operator(
                sym_op['rotation_matrix'], repr_has_cc = sym_op['repr_has_cc'],
                K_VEC = Base_vec
            ),
            basis=expr_basis,
            to_vector_fct=expr_to_vector,
            K_VEC = Base_vec
        )
        expr_mat_collection.append(expr_mat)

        LOGGER.info('Calculating matrix form of representation.')
        repr_mat = to_matrix(
            operator=repr_to_matrix_operator(
                sym_op['repr_matrix'], complex_conjugate = sym_op['repr_has_cc']
            ),
            basis=repr_basis,
            to_vector_fct=hermitian_to_vector,
            to_vector_kwargs=dict(basis_norm_squares=repr_basis_norm_squares)
        )
        repr_mat_collection.append(repr_mat)

        # outer product
        LOGGER.info('Calculating outer product.')
        full_mat = TensorProduct(expr_mat, repr_mat)

        # get Eig(F \ocross G, 1) basis
        mat = full_mat - sp.eye(full_dim)
        LOGGER.info('Calculating nullspace.')
        nullspace_basis = nullspace_blocked(mat, simplify=sp.nsimplify)
        # Modified by YJ: reshape here is necessary. The original np.array(nullspace_basis).tolist() will run into bugs for python>3.8
        curr_basis = [ bs.reshape(1, expr_dim*repr_dim) for bs in nullspace_basis ]

        if len(curr_basis) != _numeric_nullspace_dim(mat):
            raise ValueError(
                'Analytic and numeric dimensions of the nullspace of the matrix {mat} do not match'
                .format(mat=mat)
            )
        invariant_bases.append(curr_basis)

    LOGGER.info('Calculating basis intersection.')
    basis_vectors = intersection_basis(*invariant_bases)

    # ===== Added by YJ: decompose the kp model into symmetric basis  ===== #
    decomposed_repr_vec, decomposed_repr_mat, decomposed_expr, ir_str_list = [], [], [], []
    for basis_vector in basis_vectors:
        tmp_repr_vec, tmp_repr_mat, tmp_expr, linear_ir_str = decompose_kp(basis_vector, repr_basis, expr_basis, symmetry_operations, Base_vec, msg_num, kvec)
        decomposed_repr_vec.append(tmp_repr_vec)
        decomposed_repr_mat.append(tmp_repr_mat)
        decomposed_expr.append(tmp_expr)
        ir_str_list.append(linear_ir_str)
    
    LOGGER.info('Expanding basis vectors.')
    basis_vectors_expanded, decomposed_repr_symbols = [], []
    for full_vec, repr_vec in zip(basis_vectors, decomposed_repr_vec):
        basis_vectors_expanded.append( sum((v * b for v, b in zip(full_vec, full_basis)), sp.zeros(repr_matrix_size)) )
        decomposed_repr_symbols.append([ reduce(lambda x, y : x+' + '+y, [str(sp.nsimplify(v)) + '* ' + b if v != 1 else b\
                                        for v, b in zip(tmp, repr_basis_symbols) if v != 0]) for tmp in repr_vec ]) \
                                        if repr_basis_type == 'pauli' else [None] * len(repr_vec) 

    _print_result(basis_vectors_expanded, basis_vectors, decomposed_expr, decomposed_repr_mat, decomposed_repr_symbols, ir_str_list)
    return basis_vectors_expanded, decomposed_expr, decomposed_repr_mat


def _numeric_nullspace_dim(mat):
    """Numerically computes the nullspace dimension of a matrix."""
    mat_numeric = np.array(mat.evalf().tolist(), dtype=complex)
    eigenvals = la.eigvals(mat_numeric)
    return np.sum(np.isclose(eigenvals, np.zeros_like(eigenvals)))


def _print_result(kpmodels, basis_vecs, expr_basis_vecs, repr_basis_mats, repr_basis_symbols, ir_str_list):
    """ Print the result of kp models and decompoed basis"""
    if len(kpmodels) == 0:
        print('No symmetry-allowed kp models.')
    else:
        print('Number of independent kp models:', len(kpmodels))
        for ith, kp, base_vec, rep, exp, rep_sym, ir in zip(range(len(kpmodels)), kpmodels, basis_vecs, repr_basis_mats, expr_basis_vecs, repr_basis_symbols, ir_str_list):
            print('-----------------------------------------------------')
            print('%d-th kp model:'%(ith+1))
            print(kp)
            print('Basis vector:', base_vec)
            if exp == None:
                print('Fail to decompose kp.')
            else:
                if ir:
                    print('\nDecomposed basis using linear IR:', ir)
                else:
                    print('\nDecomposed basis (not symmetric):')
                print('Coefficient basis:')
                for ie in exp:
                    print(ie)
                print('\nMatrix basis:')
                for isym, ib in zip(rep_sym, rep):
                    print('Symbol:',isym, '  Expression:', ib, '\n')
        print('-----------------------------------------------------')



