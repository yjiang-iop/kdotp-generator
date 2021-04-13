# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
"""
Defines functions to decompose the kp model into symmetric basis according to MSG linear coirreps
"""

import numpy as np
import sympy as sp
from sympy.core.numbers import I
from functools import reduce
from scipy.linalg import null_space
import os

from ._expr_utils import expr_to_vector, matrix_to_expr_operator
from ._repr_utils import hermitian_to_vector, hermitian_basis, repr_to_matrix_operator, check_orthogonal, frobenius_product, solve_linear_system_numpy
from ._to_matrix import to_matrix


def float_to_sympy(a):
    # convert a float number to a sympy number
    def is_int(num):
        # return True if the input number is close enough to an integer
        return True if abs(num - round(num)) < 1e-8 else False

    if abs(a) < 1e-8:
        return 0
    sign = 1 if abs(a) == a else -1
    a = abs(a)
    sqrt_dict = {1:1, 2:1.4142135624, 3:1.7320508075, 4:2, 5:2.23606797749, 6:2.449489742783178}
    for num, sqrt_num in sqrt_dict.items():
        if is_int(sqrt_num / a):
            time = int(round(sqrt_num / a))
            return sp.sqrt(num) / time * sign
        elif is_int(a / sqrt_num):
            time = int(round(a / sqrt_num))
            return sp.sqrt(num) * time * sign

    if abs(a - 0.6830127019) < 1e-8:
        return (sp.sqrt(3) + 1)/4 * sign
    elif abs(a - 0.1830127019) < 1e-8:
        return (sp.sqrt(3) - 1)/4 * sign
    else:
        raise ValueError('Fail to identify float,a='+str(a))

def numpy_to_sympy(A):
    # convert a numpy square matrix to sympy matrix
    shape = np.shape(A)
    return np.array([[ float_to_sympy(np.real(A[i,j])) + I * float_to_sympy(np.imag(A[i,j])) for j in range(shape[1]) ] for i in range(shape[0]) ])

def round_mat(A):
    # round the small decimal part of matrix A
    def round_num(a):
        # round the small decimal part of float a 
        real, imag = np.real(a), np.imag(a)
        real = np.round(real) if abs(real - np.round(real)) < 1e-8 else real
        imag = np.round(imag) if abs(imag - np.round(imag)) < 1e-8 else imag
        return real + imag * 1j

    shape = np.shape(A)
    return np.array([[ round_num(A[i,j]) for j in range(shape[1]) ] for i in range(shape[0]) ])

def phase_correction(C):
    # try to remove the phase ambiguity of C and turn C into an integer matrix if possible
    for i in range(np.shape(C)[0]):
        for j in range(np.shape(C)[1]):
            if abs(C[i,j]) < 1e-8: #c=0
                C[i,j] = 0
                continue
            if abs(abs(C[i,j]) - 1) < 1e-8: #c=exp(i*theta)
                C = C/C[i,j]
    return round_mat(C)
    

def find_similarity_transformation(A,B):
    # find a unitary matrix C, s.t A = CBC^-1 ==> AC=CB
    # Method: tranform into a linear eq of C, and then solve C
    """
    Example:
    >>> A = [np.diag((1,-1)), np.array([[0,1],[1,0]])]
    >>> B = [np.diag((1,-1)),np.array([[0,1j],[-1j,0]])]
    >>> find_similarity_transformation(A,B)
    [[1.+0.j 0.+0.j],[0.+0.j 0.+1.j]]
    """
    def transform_AC_eq_CB(A,B):
        # transform AC-CB=0 into a linear eq of C, i.e, mat*C'=0, where C' is the flattened column vector of C.
        dim = np.shape(A[0])[0]
        num_mat = len(A)
        mat = np.zeros((num_mat * dim**2, dim**2), dtype=complex)
        for imat in range(num_mat):
            a, b = A[imat], B[imat]
            for i1 in range(dim):
                for i2 in range(dim):
                    for j1 in range(dim):
                        for j2 in range(dim):
                            if j2 == i2:
                                mat[imat*dim**2 + i1*dim+i2, j1*dim+j2] += a[i1, j1]
                            if j1 == i1:
                                mat[imat*dim**2 + i1*dim+i2, j1*dim+j2] -= b[j2,i2]
        return mat

    def transform_colC_to_matC(colC, dim):
        # transfrom N^2 dim column vec C to N*N unitary matrix C
        # shape(C) =(dim^2,m), where m is the number of indepent null vec
        assert np.shape(colC)[0] == dim**2 and np.shape(colC)[1] > 0, ('Wrong colC!',colC)
        trial_C_list = [ colC[:,ii].reshape((dim,dim)) for ii in range(np.shape(colC)[1])]
        trial_C_list.append(reduce(lambda x, y : x + y, trial_C_list))
        trial_C_list.append((colC[:,0] + colC[:,-1]).reshape((dim,dim)))

        for tmp in trial_C_list:
            CCdagger = np.dot(tmp, np.conjugate(tmp).T)
            if CCdagger[0,0].real != 0 and np.linalg.norm(CCdagger/CCdagger[0,0].real - np.eye(dim)) < 1e-6:
                return tmp / np.sqrt(CCdagger[0,0].real)
        else:
            raise ValueError('Fail to find unitary C! null vec='+str(colC))

    dim = np.shape(A[0])[0]
    trans_mat = transform_AC_eq_CB(A, B)
    tmpC = null_space(trans_mat)
    C = transform_colC_to_matC(tmpC, dim)
    C = phase_correction(phase_correction(C))
    assert sum([np.linalg.norm(np.dot(A[i],C)-np.dot(C,B[i]))>1e-6 for i in range(len(A))])<1e-6,'Wrong C! C='+str(C)
    return C

def rep_similarity_transformation(A, B, block_dim = None, nir_list = None):
    """
    for two equivalent representations A, B, find unitary C s.t C^-1AC = B
    If B is irreducible, block_dim=None, while if B is reducible, block_dim is a list to denote the dim of each block irrep.
    Algorithm:
    I. If B is irreducible: 
         M := sum(A_i \otimes B_i^*)/g, where g=len(B)
         Calculate the eigenvector v with eigenvalue 1 of M (assert there exists one such v, and only one)
         C = v.reshape(n,n)/sqrt(g), where (n,n) is the shape of B
    II. If B is reducible, each B_i=(B_i^a1, B_i^a2, ...) must be block-diagonal, with the dim of each block given in block_dim=[d1,d2..]
         M^a := sum(A_i \otimes B_i^a*)/g, calculate eigenvector v^a with eigenvalue 1 of M^a
         C^a = v.reshape(n,d^a), 
         C=(C^a1, C^a2, ...)

    Limitations: This algorithm takes only rep matrix of unitary symmetries.
    For anti-unitary symmetries, the transformation rule is C^-1 A C^*=B, which have not been considered.
    """
    def mat_eq(M1, M2):
        return True if np.max(np.abs(M1 - M2)) < 1e-8 else False

    def eigvec_of_eigval_1(M):
        eig_vals, eig_vecs = np.linalg.eig(M)
        locate_1 = [ np.abs(d - 1) < 1e-8 for d in eig_vals ]
        assert sum(locate_1) > 0, ('M does not have eigenvalue 1', locate_1,eig_vals)
        if sum(locate_1) == 1:
            v = eig_vecs[:, locate_1.index(1)]
            return 1, [v]
        else:
            num_1 = sum(locate_1)
            v_list = [ eig_vecs[:,icol] for icol, col in enumerate(locate_1) if col == True ]
            return num_1, v_list

    assert len(A) == len(B) and np.shape(A[0]) == np.shape(B[0])
    nop, dim = len(A), np.shape(A[0])[0]
    if not block_dim:
        M = sum([ np.kron(Ai, np.conjugate(Bi)) for Ai, Bi in zip(A, B) ]) / nop
        num_1, v = eigvec_of_eigval_1(M)
        assert num_1 == 1
        C = v[0].reshape((dim, dim)) * np.sqrt(dim)
    else:
        C = np.zeros((dim,dim),dtype=complex)
        for ith, IRdim in enumerate(block_dim):
            start = sum([ d * n for d, n in zip(block_dim[:ith], nir_list[:ith]) ])
            end = start + IRdim
            Ma = sum([ np.kron(Ai, np.conjugate(Bi[start:end, start:end])) for Ai, Bi in zip(A, B) ]) / nop
            num_1, va = eigvec_of_eigval_1(Ma)
            assert num_1 == nir_list[ith], (num_1, nir_list, ith)
            for nth, v in enumerate(va):
                Ca = v.reshape((dim, IRdim))
                C[:, start + nth * IRdim : end + nth * IRdim] = Ca * np.sqrt(IRdim)

    C = phase_correction(C)
    assert mat_eq(C @ np.conjugate(C).T, np.eye(dim)), ('C not unitary!', C)
    assert all([ mat_eq(Ai @ C, C @ Bi) for Ai, Bi in zip(A, B) ]), 'C fails to tranform A to B!'        
    return C  



def irrep_to_linear_irrep(irmat_gen, k_dict, symmetry_operations):
    # for an input irrep irmat, identify which linear irrep corresponds to it, using characters
    # Return: linear irrep, and the similarity transformation matrix
    def generate_all_irmat(irmat, k_dict, symmetry_operations):
        # For input irmat, generate the irmat for all unitary ops
        def mat_in_array(mat, array):
            # return True/False, and the index of mat in array if True
            tmp = [ np.max(np.abs(mat - ii)) < 1e-8 for ii in array]
            if sum(tmp) == 0:
                return False, None
            elif sum(tmp) == 1:
                return True, tmp.index(1)
            else:
                raise ValueError('mat appear in array more than once!', mat, array)
        
        rotC = k_dict['rotC'][0:len(k_dict['rotC'])]
        un_rot_set, un_rep_set, au_rot_set, au_rep_set = [], [], [], []
        for isym, sym_op in enumerate(symmetry_operations):
            if sym_op['repr_has_cc']:
                au_rot_set.append(np.array(sym_op['rotation_matrix']).astype(int))
                au_rep_set.append(irmat[isym])
            else:
                un_rot_set.append(np.array(sym_op['rotation_matrix']).astype(int))
                un_rep_set.append(irmat[isym])

        for iop1, op1 in enumerate(un_rot_set):
            for iop2, op2 in enumerate(un_rot_set):
                op3 = op1 @ op2
                if not mat_in_array(op3, un_rot_set)[0]:
                    un_rot_set.append(op3)
                    un_rep_set.append(un_rep_set[iop1] @ un_rep_set[iop2])

        # sort irmat_full according to rotC
        un_rep_sorted = []
        for irot, rot in enumerate(rotC):
            inarray, index = mat_in_array(rot, un_rot_set)
            assert inarray, ('Input rotations cannot generate the little group!', rotC, '\n', un_rot_set)
            un_rep_sorted.append(un_rep_set[index])
        
        return un_rep_sorted
                    
    
    def identify_linear_irrep(irmat, k_dict, toln=1e-8):
        def ortho_thm_factor(nir, ir, k_dict):
            # for MSGs with anti-unitary ops, need to divide torison number I
            if not k_dict['au_op']['exist']:
                return nir
            else:
                if ir['reality'] == -1:
                    nir /= 4
                elif ir['reality'] == 0:
                    nir /= 2
                return nir

        def mat_direct_sum(mats1, mats2):
            # for two input lists of mats, return their direct sum [[mat1, 0],[0,mat2]].
            if np.sum(np.abs(mats1[0])) == 0: # if mats1 are empty:
                return mats2
            dim1r, dim1c = np.shape(mats1[0])
            dim2r, dim2c = np.shape(mats2[0])
            nop = len(mats1)

            direct_sum = [] 
            for mat1, mat2 in zip(mats1, mats2):
                sum_mat = np.zeros((dim1r+dim2r, dim1c+dim2c), dtype=complex)
                sum_mat[0:dim1r, 0:dim1c] = mat1
                sum_mat[dim1r:,dim1c:] = mat2
                direct_sum.append(sum_mat)
            return direct_sum

        # for an input irrep irmat, identify which linear irrep corresponds to it, using characters
        nop = len(k_dict['rotC'])
        irmat_characters = [ np.trace(mat) for mat in irmat ]
        
        coirs = k_dict['coirreps']
        BR = [0] * len(coirs)
        for ir, irrep in enumerate(coirs):
            assert irrep['label'][-1] != 'd'
            nir = sum([ irmat_characters[iop].conjugate() * irrep['characters'][iop] for iop in range(nop) ]) / nop
            nir = ortho_thm_factor(nir, irrep, k_dict)
            assert abs(nir.imag) < toln and abs(nir-np.round(nir)) < toln, ('Non-integer nir!', nir, irrep['label'], irmat_characters, irrep['characters'])
            nir = int(round(nir.real))
            BR[ir] = nir

        # BR may contain more than one coirrep
        assert np.shape(irmat[0])[0] == sum([ coirs[ir]['dim'] * nir for ir, nir in enumerate(BR) ]), ('Wrong decomposed BR dim!', BR)
        if sum(BR) == 1:
            return BR, coirs[BR.index(1)]['matrices']
        else:
            direct_sum = [[]]*nop
            for ir, nir in enumerate(BR):
                if nir == 1:
                    direct_sum = mat_direct_sum(direct_sum, coirs[ir]['matrices'])
                elif nir > 1:
                    for i in range(nir):
                        direct_sum = mat_direct_sum(direct_sum, coirs[ir]['matrices'])
            return BR, direct_sum
         
    irmat_all = generate_all_irmat(irmat_gen, k_dict, symmetry_operations)
    linear_ir_vec, linear_irmat = identify_linear_irrep(irmat_all, k_dict)
    linear_ir_str = reduce(lambda x,y: x+y, [ str(nir) + ' ' + k_dict['coirreps'][ir]['label'] + '  ' for ir, nir in enumerate(linear_ir_vec) if nir > 0])
    try:
        similarity_mat = find_similarity_transformation(irmat_all, linear_irmat)
    except:
        linear_ir_dim = [ k_dict['coirreps'][ir]['dim'] for ir, nir in enumerate(linear_ir_vec) if nir > 0 ]
        nir_vec = [ nir for nir in linear_ir_vec if nir > 0 ]
        similarity_mat = rep_similarity_transformation(irmat_all, linear_irmat, linear_ir_dim, nir_vec)

    similarity_mat = numpy_to_sympy(similarity_mat)
    return linear_ir_str, similarity_mat


def transform_basis(kp_ir_mats, ir_repr_vec, ir_repr_mat, ir_expr, symmetry_operations, msg_num=None, kvec=None):
    # For input expr, repr basis, find correspond Bilbao linear coirrep and similarity matrix, and then transform them. 
    kp_ir_mats = [ round_mat(m) for m in kp_ir_mats ] 
    Locate = os.path.abspath(os.path.dirname(__file__))
    if msg_num != None and '.' not in str(msg_num):
        # if input msg_num is an sg num, find corresponding type2 msg_num (msg_num always has an '.' in it)
        sg_to_type2msg = np.load(Locate+'/MSG_linear_coir_data/sg_to_type2msg.npy', allow_pickle=True)[0]
        msg_num = sg_to_type2msg[str(msg_num)]

    msg_linear_coir_data = np.load(Locate+'/MSG_linear_coir_data/'+msg_num+'.npy', allow_pickle=True)
    linear_ir_data = [ kdict for kdict in msg_linear_coir_data if np.max(np.abs(kdict['kvecC'] - kvec)) < 1e-8 ][0]

    linear_ir_str, sim_mat = irrep_to_linear_irrep(kp_ir_mats, linear_ir_data, symmetry_operations)

    # D(g) = C * D'(g) * C^-1, g|psi> = |psi> D(g) ==> g|psi>C = |psi>C * D'(g)  ==> |psi> -> |psi>C
    repr_vec_transformed = [ reduce(lambda x,y: x+y, [ sp.Matrix(v) * b for v, b in zip(ir_repr_vec, col) ]) for col in sim_mat.T ]
    repr_mat_transformed = [ sum([ v * b for v, b in zip(ir_repr_mat, col)], sp.zeros(np.shape(ir_repr_mat[0])[0])) for col in sim_mat.T ]
    expr_basis_transformed = [ sum([ v * b for v, b in zip(ir_expr, col)]) for col in sim_mat.T ]

    return repr_vec_transformed, repr_mat_transformed, expr_basis_transformed, linear_ir_str



def decompose_kp(basis_vector, repr_basis, expr_basis, symmetry_operations, Base_vec, msg_num=None, kvec=None):
    """ decompose the kp hamiltonian into symmetric repr and expr basis, using linear coirreps of the little group. 
    
    :param basis_vector: The kp hamiltonian written in the full basis
    :type basis_vector: :py:class:`list` of :py:class:`int`
    
    :param repr_basis: The basis for the hermitian matrices, with the same size as the representations. 
                       By default, the :py:func:`.hermitian_pauli_basis` of the appropriate size is used.
    :type repr_basis: :py:class:`list` of :py:mod:`sympy` matrices

    :param expr_basis: The basis for the monomial functions that are considered.
    :type expr_basis: :py:class:`list` of :py:mod:`sympy` expressions

    :param symmetry_operations: The symmetry operations that the Hamiltonian should respect.
    :type symmetry_operations: :py:class: `list` of `symmetry_representation.SymmetryOperation`

    :param Base_vec: The variable of the hamiltonian
    :type kp_variable: :py:class: `list` of :py:mod:`sympy` expressions
    
    :param msg_num & kvec: two string used to denote the magnetic space group and little group k, 
                             used to locate linear representations in order to decompose kp hamiltonian.
    :type msg_num & kvec: :py:class:str

    :returns: decomposed repr and expr basis, as well as the linear coirrep label.
    """

    def vec_in_array(vec, array):
        # return whether vec is in array(can time a coefficient), and if yes, its location in array, and the coeff
        array = [ np.array(a.copy()).astype(complex) for a in array ]
        array2 = array.copy()
        array2.append(np.array(vec).astype(complex))
        if np.linalg.matrix_rank(array) == np.linalg.matrix_rank(array2) - 1:
            return False, None, None
        elif np.linalg.matrix_rank(array) == np.linalg.matrix_rank(array2):
            for irow, row in enumerate(array):
                tmp = [row, np.array(vec).astype(complex)]
                if np.linalg.matrix_rank(tmp) == 1:
                    for ii, iv in enumerate(vec):
                        if abs(iv) > 1e-9:
                            coeff = iv / row[ii]
                            coeff = int(round(coeff)) if abs(coeff - round(coeff)) < 1e-9 else coeff
                            return True, irow, coeff
            else:
               #raise Exception('Fail to find linearly-dependent vec!',vec, array)
                return False, None, None # vec is an linear combination of rows in array, need checking.
        else:
            raise Exception('Wrong array rank!', array)

    def _extract_basis_using_expr(basis_vector, repr_basis, expr_basis):
        # input basis_vectors is a list of vectors written in full basis, eg, [k1*X1, k1*X2, k1*X3, k2*X1, k2*X2, k2*X3, k3*X1, k3*X2, k3*X3]
        repr_dim, expr_dim = len(repr_basis), len(expr_basis)
        repr_vec, expr_vec = [], [] 
        for i1 in range(expr_dim): # i1'th k basis, eg, [k1*X1, k1*X2, k1*X3], extract [nX1, nX2, nX3]
            if sum(np.abs(basis_vector[i1*repr_dim : (i1+1)*repr_dim])) != 0:
                repr_ = basis_vector[i1*repr_dim : (i1+1)*repr_dim]
                in_array, index, coeff = vec_in_array(repr_, repr_vec)
                if not in_array:
                    repr_vec.append(repr_)
                    expr_ = [0] * expr_dim
                    expr_[i1] = 1 
                    expr_vec.append(expr_)
                else:
                    expr_vec[index][i1] += coeff

        repr_vec_expression = [ sp.Matrix(reduce(lambda x, y : x + y, [v * b for v, b in zip(tmp, repr_basis)])) for tmp in repr_vec ]
        expr_vec_expression = [ sp.simplify(sum([sp.nsimplify(v) * b for v, b in zip(tmp, expr_basis)])) for tmp in expr_vec ]
        return repr_vec, expr_vec, repr_vec_expression, expr_vec_expression
                        
    def _extract_basis_using_repr(basis_vector, repr_basis, expr_basis):
        def _extract_expr_vec(vec, n_repr, repr_dim, expr_dim):
            # extract the expr_ corresponding to nth repr from vec
            expr_ = []
            for jj in range(expr_dim):
                expr_.append(vec[jj * repr_dim + n_repr])
            return expr_

        repr_dim, expr_dim = len(repr_basis), len(expr_basis)
        repr_vec, expr_vec = [], [] # eg [ [[1,0,0],[0,1,0]], [...] ]
        for i in range(repr_dim): 
            # For i'th repr basis, if its corresponding expr_ not zero, gives one term in decomposition, ki*Xi (i.e, repr_, expr_)
            expr_ = _extract_expr_vec(basis_vector, i, repr_dim, expr_dim)
            if np.sum(np.abs(expr_)) > 1e-6:
                repr_ = np.zeros(repr_dim)
                repr_[i] = 1

                in_array, index, coeff = vec_in_array(expr_, expr_vec)
                if not in_array:
                    repr_vec.append(repr_)
                    expr_vec.append(expr_)
                else:
                    repr_vec[index][i] += coeff
        
        repr_vec_expression = [ sp.simplify(sp.Matrix(reduce(lambda x, y : x + y, [sp.nsimplify(v) * b for v, b in zip(tmp, repr_basis)]))) for tmp in repr_vec ]
        expr_vec_expression = [ sp.simplify(sum([sp.nsimplify(v) * b for v, b in zip(tmp, expr_basis)])) for tmp in expr_vec ]
        return repr_vec, expr_vec, repr_vec_expression, expr_vec_expression
                        
    def _extract_irrep(basis_expression, symmetry_operations, Base_vec, in_type='expr'):
        # eg, basis_expression = [ [1,0,0],[0,1,0] ], a set of invariant vectors
        # return the irrep mat of these vectors
        irrep = []
        if in_type == 'expr':
            for isym_op, sym_op in enumerate(symmetry_operations):
                expr_mat = to_matrix(
                    operator=matrix_to_expr_operator(
                        sym_op['rotation_matrix'], repr_has_cc = sym_op['repr_has_cc'],
                        K_VEC = Base_vec
                    ),
                    basis=basis_expression,
                    to_vector_fct=expr_to_vector,
                    K_VEC = Base_vec
                )
                irrep.append(expr_mat)
        else: 
            if check_orthogonal(basis_expression): # basis orthogonal, use frobenius_produce to solve
                repr_basis_norm_squares = [frobenius_product(b, b) for b in basis_expression]
                for isym_op, sym_op in enumerate(symmetry_operations):
                    repr_mat = to_matrix(
                        operator=repr_to_matrix_operator(
                            sym_op['repr_matrix'], complex_conjugate = sym_op['repr_has_cc']
                        ),
                        basis=basis_expression,
                        to_vector_fct=hermitian_to_vector,
                        to_vector_kwargs=dict(basis_norm_squares=repr_basis_norm_squares)
                    )
                    irrep.append(repr_mat)
            else: # repr basis not orthogonal, use numpy to directly solve
                for isym_op, sym_op in enumerate(symmetry_operations):
                    repr_mat = to_matrix(
                        operator=repr_to_matrix_operator(
                            sym_op['repr_matrix'], complex_conjugate = sym_op['repr_has_cc']
                        ),
                        basis=basis_expression,
                        to_vector_fct=solve_linear_system_numpy
                    )
                    irrep.append(repr_mat)
        return irrep

    def _check_irrep_expr(basis_expression, symmetry_operations, irrep, Base_vec):
        # check if R*f(k) == f(k) * D (=D.T * f), where R is the symmetry_operation, D is the irmat
        for irmat, sym_op in zip(irrep, symmetry_operations):
            operator = matrix_to_expr_operator(sym_op['rotation_matrix'], repr_has_cc = sym_op['repr_has_cc'], K_VEC=Base_vec)
            basis_expression = sp.Matrix(basis_expression)
            Rf = operator(basis_expression)
            Df = irmat.transpose() @ basis_expression
            assert Rf.simplify() == Df.simplify(), (Rf, Df, basis_expression, irmat, sym_op.rotation_matrix)
    
    
    repr_dim, expr_dim = len(repr_basis), len(expr_basis)
    # First decompose kp using expression basis
    repr_vec, expr_vec, repr_expression, expr_expression = _extract_basis_using_expr(basis_vector, repr_basis, expr_basis) 
    repr_irrep = _extract_irrep(repr_expression, symmetry_operations, Base_vec, in_type = 'repr')
    _check_irrep_expr(expr_expression, symmetry_operations, repr_irrep, Base_vec)
    repr_irrep = [ np.array(mat).astype(np.complex) for mat in repr_irrep ]

    linear_ir_str = None
    if msg_num:
        try:
            repr_vec, repr_expression, expr_expression, linear_ir_str = transform_basis(repr_irrep, repr_vec, repr_expression, expr_expression, \
                symmetry_operations, msg_num, kvec)
        except:
            # If fail to decompose using expression basis, redo the decomposion using representation basis
            repr_vec, expr_vec, repr_expression, expr_expression = _extract_basis_using_repr(basis_vector, repr_basis, expr_basis) 
            repr_irrep = _extract_irrep(repr_expression, symmetry_operations, Base_vec, in_type = 'repr')
            _check_irrep_expr(expr_expression, symmetry_operations, repr_irrep, Base_vec)
            repr_irrep = [ np.array(mat).astype(np.complex) for mat in repr_irrep ]
            try:
                repr_vec, repr_expression, expr_expression, linear_ir_str = transform_basis(repr_irrep, repr_vec, repr_expression, expr_expression,\
                    symmetry_operations, msg_num, kvec)
            except:
                print('Warning: fail to transform basis according to linear coirrep! The decomposed basis is not symmetric.')

    return repr_vec, repr_expression, expr_expression, linear_ir_str 


