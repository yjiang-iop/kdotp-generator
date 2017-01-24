#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    24.01.2017 10:19:47 CET
# File:    linalg.py

from collections import namedtuple

import sympy as sp
import numpy as np

ZassenhausResult = namedtuple('ZassenhausResult', ['sum', 'intersection'])

def zassenhaus(basis_a, basis_b):
    r"""
    Given two bases ``basis_a`` and ``basis_b`` of vector spaces :math:`U, W \subseteq V`, computes bases of :math:`U + W` and :math:`U \cap W` using the Zassenhaus algorithm.
    """
    # handle the case where one of the bases is empty
    if len(basis_a) == 0 or len(basis_b) == 0:
        return ZassenhausResult(sum=basis_a + basis_b, intersection=[])
    
    # Set up the Zassenhaus matrix from the given bases.
    A = sp.Matrix(basis_a)
    B = sp.Matrix(basis_b)
    dim = A.shape[1]
    if B.shape[1] != dim:
        raise ValueError('Inconsistent dimensions of the two bases given.')
    zassenhaus_mat = A.row_join(A).col_join(B.row_join(sp.zeros(*B.shape)))
    mat, pivot = zassenhaus_mat.rref()
    
    # idx is the row index of the first row belonging to the intersection basis
    idx = np.searchsorted(pivot, dim)
    plus_basis = mat[:idx, :dim].tolist()
    # The length of the pivot table is used to get rid of all-zero rows
    int_basis = mat[idx:len(pivot), dim:].tolist()
    return ZassenhausResult(sum=plus_basis, intersection=int_basis)
    
def intersection_basis(*bases):
    bases_sorted = sorted(bases, key=lambda x: len(x))
    current_basis = bases_sorted.pop(0)
    for b in bases_sorted:
        if len(current_basis) == 0:
            break
        current_basis = zassenhaus(current_basis, b).intersection
    return current_basis
    
if __name__ == '__main__':
    print(intersection_basis([[1, 1, 0], [0, 1, 0]], [[0, 1, 1], [0, 0, 1]], [[0, 2, 1]], []))
    print(zassenhaus([], [[1, 2, 3]]))