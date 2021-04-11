#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
# Adapted from the kdotp-symmetry package by: Dominik Gresch <greschd@gmx.ch>  © 2017-2018, ETH Zurich, Institut für Theoretische Physik

import sympy as sp
from sympy.core.numbers import I
import sympy.physics.matrices as sm
from sympy.physics.quantum import TensorProduct
import symmetry_representation as sr
import kdotp_symmetry as kp

# This example is from TaAs2, SG 12, kpoint M
# creating the symmetry operations
c2y = sr.SymmetryOperation(
    rotation_matrix=[[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    repr_matrix=sp.diag(I, -I, I, -I),
    repr_has_cc=False,
)

parity = sr.SymmetryOperation(
    rotation_matrix=-sp.eye(3),
    repr_matrix=sp.diag(1, 1, -1, -1),
    repr_has_cc=False,
)

time_reversal = sr.SymmetryOperation(
    rotation_matrix=sp.eye(3),
    repr_matrix=TensorProduct(sp.eye(2), sp.Matrix([[0, -1], [1, 0]])),
    repr_has_cc=True,
)


def print_result(sym_ops, order, kp_variable, msg_num, kvec):
    """prints the basis for a given order of k"""

    num_input = len(kp_variable.split())
    if num_input == 1:
        order_list = [ [i] for i in range(order+1) ]
    elif num_input == 2:
        order_list = [ [i,j] for i, j in zip(range(order+1), np.arange(order,-1,-1)) ]
    elif num_input == 3:
        order_list = [ [i,j,k] for i in range(order+1) for j in range(order+1) for k in range(order+1) if i+j+k==order ] 

    print('Print kp results:')
    for order in order_list:
        print('\n\n==========  Result of msg %s  order = %s  =========='%(msg_num, str(order)))
        kpmodel, irrep_expr_basis, irrep_repr_basis = kp.symmetric_hamiltonian(
            sym_ops,
            expr_basis = kp.monomial_basis(order, kp_variable),
            repr_basis = 'pauli',
            kp_variable = kp_variable,
            msg_num = msg_num,
            kvec = kvec
        )


if __name__ == '__main__':
    kp_variable = 'k' # kp_variable can be either 'k', 'E', 'B', 'k E', 'k B', 'E B', 'e'
    order = 2  # Generate all kp models with order <= 2
    msg_num = '12.59'
    kvec = [0,1,0.5]
    sym_ops = [c2y, parity, time_reversal]
    print_result(sym_ops, order, kp_variable, msg_num, kvec)

