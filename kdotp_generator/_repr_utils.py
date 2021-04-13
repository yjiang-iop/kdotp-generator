# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
# Adapted from the kdotp-symmetry package by: Dominik Gresch <greschd@gmx.ch>  © 2017-2018, ETH Zurich, Institut für Theoretische Physik
"""
Utilities for handling symmetry representations, such as converting them to matrix form.
"""

import numpy as np
import sympy as sp
from sympy.physics.quantum import TensorProduct
from sympy.core.numbers import I
from sympy import sqrt


def frobenius_product(A, B):
    r"""
    Returns the Frobenius scalar product <A, B> = Tr(A^\dagger B) for two matrices.
    """
    return (sp.simplify(A.H) @ sp.simplify(B)).trace().nsimplify()


def hermitian_basis(dim):
    """
    Returns a basis of the hermitian matrices of size ``dim`` that is orthogonal w.r.t. the Frobenius scalar product.

    :param dim: size of the matrices
    :type dim:  int

    Example:

        >>> import kdotp_symmetry as kp
        >>> kp.hermitian_basis(2)
        [Matrix([
        [1, 0],
        [0, 0]]), Matrix([
        [0, 0],
        [0, 1]]), Matrix([
        [0, 1],
        [1, 0]]), Matrix([
        [0, -I],
        [I,  0]])]
    """
    basis = []
    # diagonal entries
    for i in range(dim):
        basis.append(sp.SparseMatrix(dim, dim, {(i, i): 1}))

    # off-diagonal entries
    for i in range(dim):
        for j in range(i + 1, dim):
            # real
            basis.append(sp.SparseMatrix(dim, dim, {(i, j): 1, (j, i): 1}))

            # imag
            basis.append(
                sp.SparseMatrix(dim, dim, {(i, j): -sp.I,
                                           (j, i): sp.I})
            )

    assert len(basis) == dim**2
    return basis

def hermitian_pauli_basis(dim):
    '''
    Returns a basis of the hermitian pauli matrices (or Gellmann matrics) of size ``dim`` that is orthogonal w.r.t. the Frobenius scalar product.

    :param dim: size of the matrices
    :type dim:  int

    '''
    s0 = sp.eye(2)
    s1 = sp.Matrix([[0, 1], [1, 0]])
    s2 = sp.Matrix([[0, -I], [I, 0]])
    s3 = sp.Matrix([[1, 0], [0, -1]])
    pauli_basis = [s0, s1, s2, s3]

    g0 = sp.eye(3)
    g1 = sp.SparseMatrix(3, 3, {(0, 1): 1, (1, 0): 1})
    g2 = sp.SparseMatrix(3, 3, {(0, 1): -I, (1, 0): I})
    g3 = sp.diag(1,-1,0)
    g4 = sp.SparseMatrix(3, 3, {(0, 2): 1, (2, 0): 1})
    g5 = sp.SparseMatrix(3, 3, {(0, 2): -I, (2, 0): I})
    g6 = sp.SparseMatrix(3, 3, {(1, 2): 1, (2, 1): 1})
    g7 = sp.SparseMatrix(3, 3, {(1, 2): -I, (2, 1): I})
    g8 = sp.diag(1,1,-2)/sp.sqrt(3)
    gellmann_basis = [g0, g1, g2, g3, g4, g5, g6, g7, g8]

    if dim == 1:
        return [sp.eye(1)]
    elif dim == 2:
        return pauli_basis
    elif dim == 4:
        return [TensorProduct(si, sj) for si in pauli_basis for sj in pauli_basis]
    elif dim == 8:
        return [TensorProduct(TensorProduct(si, sj), sk) for si in pauli_basis for sj in pauli_basis for sk in pauli_basis]
    elif dim == 16:
        return [TensorProduct(TensorProduct(TensorProduct(si, sj), sk), sl) for si in pauli_basis for sj in pauli_basis \
                    for sk in pauli_basis for sl in pauli_basis]

    elif dim == 3: # Gellmann matrix
        return gellmann_basis
    elif dim == 6:
        return [TensorProduct(si, sj) for si in pauli_basis for sj in gellmann_basis]
    elif dim == 9:
        return [TensorProduct(si, sj) for si in gellmann_basis for sj in gellmann_basis]
    elif dim == 12:
        return [TensorProduct(TensorProduct(si, sj), sk) for si in pauli_basis for sj in pauli_basis for sk in gellmann_basis]
    else:   
       #raise ValueError('input dim not supported!')
        return hermitian_basis(dim)

def hermitian_pauli_basis_symbols(dim):
    '''
    Returns a list of the symbols of hermitian pauli matrices (or Gellmann matrics) of size ``dim``, e.g 's0 * s1 '

    :param dim: size of the matrices
    :type dim:  int

    '''
    pauli_basis = [str(i) for i in range(4)]

    gellmann_basis = [str(i) for i in range(9)]

    if dim == 1:
        return ['Identity']
    elif dim == 2:
        return [ 's_'+ s for s in pauli_basis]
    elif dim == 4:
        return [ 'Gamma_'+si+sj for si in pauli_basis for sj in pauli_basis]
    elif dim == 8:
        return [ 'Gamma_'+si+sj+sk for si in pauli_basis for sj in pauli_basis for sk in pauli_basis]
    elif dim == 16:
        return [ 'Gamma_'+si+sj+sk+sl for si in pauli_basis for sj in pauli_basis \
                    for sk in pauli_basis for sl in pauli_basis]
    elif dim == 3: # Gellmann matrix
        return [ 'G_'+g for g in gellmann_basis]
    elif dim == 6:
        return [ 'G_'+si+sj for si in pauli_basis for sj in gellmann_basis]
    elif dim == 9:
        return [ 'G_'+si+sj for si in gellmann_basis for sj in gellmann_basis]
    elif dim == 12:
        return [ 'G_'+si+sj+sk for si in pauli_basis for sj in pauli_basis for sk in gellmann_basis]
    else:   
        return [ 'H_'+str(i) for i in range(dim*dim) ]


def check_orthogonal(basis):
    """Check orthogonality for a given ``basis``."""
    for i, b_i in enumerate(basis):
        for offset, b_j in enumerate(basis[i:]):
            j = i + offset
            frob_product = frobenius_product(b_i, b_j)
            if i == j:
                pass # Modified by YJ: norm need not be 1, only orthogonal is required
               #if frob_product != 1:
               #    raise ValueError(
               #        'Basis element {} has norm {}, not one.'.format(
               #            i, frob_product
               #        )
               #    )
            else:
                if frob_product != 0:
                    return False
                   #raise ValueError(
                   #    'Basis elements {}, {} are not orthogonal.'.format(
                   #        i, j
                   #    )
                   #)
    return True


def hermitian_to_vector(mat, basis, basis_norm_squares=None):
    """
    Returns a the vector representing the ``mat`` w.r.t. the given *orthogonal* ``basis``.
    """
    vec = tuple(
        frobenius_product(mat, b) / norm_sq for b, norm_sq in zip(
            basis, basis_norm_squares
            or [frobenius_product(b, b) for b in basis]
        )
    )
    vec = tuple(v.nsimplify() for v in vec)
    # check consistency
    if not mat.equals(
        sum((v * b for v, b in zip(vec, basis)), sp.zeros(*mat.shape))
    ):
        raise ValueError(
            'Vector {vec} in basis {basis} does not match mat {mat}'.
            format(vec=vec, basis=basis, mat=mat)
        )
    return vec


def repr_to_matrix_operator(matrix_representation, complex_conjugate=False):
    """
    Converts a symmetry representation into the corresponding matrix operator.

    :param matrix_representation: Real-space matrix form of the symmetry representation.
    :type matrix_representation: sympy.Matrix

    :param complex_conjugate: Specifies whether the representation contains complex conjugation.
    :type complex_conjugate: bool
    """
    matrix_representation = sp.Matrix(matrix_representation)

    def operator(matrix):
        if complex_conjugate:
            return sp.simplify(matrix_representation @ matrix.conjugate() @ matrix_representation.H)
        return sp.simplify(matrix_representation @ matrix @ matrix_representation.H) # simplify() add by YJ

    return operator


def solve_linear_system_numpy(B, basis):
    """
    Added by YJ: solve a linear system using numpy lstsq
    :param: B = [n1, n2, ...], basis = [A1, A2, ...], solve a linear system A1*x1 + A2*x2 + ... = B

    Method: transform basis into A, and B into b, s.t. A*x = b, then use numpy to solve.
    """
    assert basis[0].shape == B.shape and basis[0].shape[0] == basis[0].shape[1]
    nA, dim = len(basis), basis[0].shape[0]
    A, b = np.zeros((dim**2, nA), dtype=complex), np.zeros((dim**2, 1), dtype=complex)
    for ia, Ai in enumerate(basis):
        A[:,ia] = np.array(Ai).reshape(dim**2).astype(complex)
    b[:,0] = np.array(B).reshape(dim**2).astype(complex)
    sol, res = np.linalg.lstsq(A, b, rcond=None)[0:2]
    assert np.linalg.norm(res) < 1e-4, (res, sol)
    return tuple(sol.reshape(nA))

