#!/usr/bin/env python
# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences

import numpy as np
import sympy as sp
from sympy.core.numbers import I
import sympy.physics.matrices as sm
from sympy.physics.quantum import TensorProduct
import kdotp_generator as kp

# This example is of non-magnetic SG 123, GM, GM6d + GM9d 
# creating the symmetry operations
# The rotation_matrix and repr_matrix are recommended to use sp.Matrix(), 
# Avoid float numbers, use sp.Rational() and sp.sqrt() instead.
c4z = {
    'rotation_matrix': [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    'repr_matrix': TensorProduct(sp.diag(1,-1), sp.Matrix([[-sp.sqrt(2)+sp.sqrt(2)*I, 0], [0, -sp.sqrt(2)-sp.sqrt(2)*I]])/2),
    'repr_has_cc': False 
}

c2 = {
    'rotation_matrix': [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
    'repr_matrix': TensorProduct(sp.eye(2), sp.diag([[0, I], [I, 0]])),
    'repr_has_cc': False 
}

P = {
    'rotation_matrix': -sp.eye(3),
    'repr_matrix': sp.diag(1, 1, -1, -1),
    'repr_has_cc': False
}

time_reversal = {
    'rotation_matrix': sp.eye(3),
    'repr_matrix': TensorProduct(sp.eye(2), sp.Matrix([[0, -1], [1, 0]])),
    'repr_has_cc': True 
}

def print_result(sym_ops, kp_variable, order, msg_num, kvec):
    """prints the basis for a given order of k"""

    num_input = len(kp_variable.split())
    if num_input == 1:
        order_list = [ [i] for i in range(order+1) ]
    elif num_input == 2:
        order_list = [ [i,j] for i, j in zip(range(order+1), np.arange(order,-1,-1)) ]
    elif num_input == 3:
        order_list = [ [i,j,k] for i in range(order+1) for j in range(order+1) for k in range(order+1) if i+j+k==order ] 

    print('\nPrint kp results:')
    for order in order_list:
        print('\n==========  Result of msg %s  order = %s  =========='%(msg_num, str(order)))
        kpmodel, irrep_expr_basis, irrep_repr_basis = kp.symmetric_hamiltonian(
            sym_ops,
            kp_variable = kp_variable,
            order = order,
            repr_basis = 'pauli',
            msg_num = msg_num,
            kvec = kvec
        )


if __name__ == '__main__':
    kp_variable = 'k' # kp_variable can be either 'k', 'E', 'B', 'e', 'k E', 'k B', 'E B', 'k E B' (use space to seperate different variable)
    order = 2         # Generate all kp models with order <= 2
    msg_num = 123     # Optional parameter. Gives the msg number. For non-magnetic SG (i.e., type-2 msg), either the SG number or the msg number is applicable.
    kvec = [0,0,0]    # Optional parameter. Gives the kvec of the corresponding kpoint (convertional lattice).
    sym_ops = [c4z, c2, P, time_reversal]
    print_result(sym_ops, kp_variable, order, msg_num, kvec)

