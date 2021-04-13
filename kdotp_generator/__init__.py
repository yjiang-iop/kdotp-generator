# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
# Adapted from the kdotp-symmetry package by: Dominik Gresch <greschd@gmx.ch>  © 2017-2018, ETH Zurich, Institut für Theoretische Physik
"""A tool for computing k.p effective Hamiltonians with couplings to external fields including E, B, and epsilon, under given symmetry constraints."""

__version__ = '1.0'

from ._expr_utils import monomial_basis
from ._repr_utils import hermitian_basis, hermitian_pauli_basis, hermitian_pauli_basis_symbols
from ._symmetric_hamiltonian import symmetric_hamiltonian
from ._decompose_kp import decompose_kp

__all__ = ['symmetric_hamiltonian', 'monomial_basis', 'hermitian_basis',  
           'hermitian_pauli_basis', 'hermitian_pauli_basis_symbols', 'decompose_kp']
