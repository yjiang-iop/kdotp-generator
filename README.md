# *kdotp-generator* 
This package is used to automatically generate kdotp effective Hamiltonians, developed based on the [*kdotp-symmetry*](https://kdotp-symmetry.greschd.ch/en/latest/) package by D. Gresch. 

To use this package, users can either input their symmetry data and compute the effective Hamiltonians, or directly refer to the pre-computed results for magnetic space groups (MSGs) in the `MSG_kp_results` folder.

If you are using this package, please cite the following works:
- [*A kp effective Hamiltonian generator*](http://arxiv.org/abs/2104.08493), Yi Jiang, Zhong Fang, and Chen Fang.
- [*Identifying Topological Semimetals*](https://www.research-collection.ethz.ch/handle/20.500.11850/308602), D. Gresch, PhD Thesis.

## Installation
Pip is used to install the package:
> pip install kdotp-generator 

The pre-computed results for MSGs are not included in the pip installation. They can be downloaded from the `MSG_kp_results` folder in the [Github page of this package](https://github.com/yjiang-iop/kdotp-generator).

Recently, we also add a new folder `MSG_kp_results_pair`, which contains the kp results for any pair of irreducible representations.

## Example
Two examples are given in the *examples* folder, i.e., *input1.py* and *input2.py*, which calculate the effective Hamiltonians for k and (E, B), respectively, for non-magnetic space group 123 P4/mmm, representation $\Gamma_6, \Gamma_9$, up to the second order. 


## Inputs
Three parameters need to be specified to run the main function `symmetric_hamiltonian`:
1. `symmetry_operations`
2. `kp_variable`
3. `order`

#### symmetry_operations
The symmetry operations of the group generators. They are generated using python `dictionary`, with three keys:
- `rotation_matrix`: the O(3) (real-space) rotation matrix of the operation.
- `repr_matrix`: the representation matrix of the operation, either reducible or irreducible. 
- `repr_has_cc`: a boolean flag to determine whether the operation is unitary or anti-unitary. If `repr_has_cc=True`, the operation becomes g*T, with g being unitary and specified in `rotation_matrix`, and `repr_matrix` contains a complex conjugation.

Both `rotation_matrix` and `repr_matrix` are recommended to use `sympy.Matrix()`, and non-integer numbers should be typed by `sympy.Rational()` and `sympy.sqrt()`. If there are exponential numbers, converting them to square numbers is more efficient to compute, e.g., `sp.exp(I*Pi/4)=(sp.sqrt(2)+sp.sqrt(2)*I)/2`.

#### kp_variable
The variable of the effective Hamiltonians. When using the default value `'k'`, the kp Hamiltonian is calculated. Other variables include `'E'`, `'B'`, and `'e'` (the strain tensor epsilon), as well as combinations of them like `'k E'`, `'k B'`, and `'E B'` (use space to separate different variables). Three or more variables may also be used, but the computation could be very slow.

#### order
The order of the Hamiltonian to be computed. If `kp_variable` contains only one variable, `order=[n]`, while if `kp_variable` contains two variables, `order=[n1,n2]`. 

The monomial basis is automatically generated from `order` and `kp_variable`:
```python
>>>import kdotp_generator as kp
>>> kp.monomial_basis([1], kp_variable='k')
[kx, ky, kz]
>>> kp.monomial_basis([2], kp_variable='k')
[kx**2, kx*ky, kx*kz, ky**2, ky*kz, kz**2]
>>> kp.monomial_basis([1,1], kp_variable='E B')
[Bx*Ex, By*Ex, Bz*Ex, Bx*Ey, By*Ey, Bz*Ey, Bx*Ez, By*Ez, Bz*Ez]
```

#### Optional parameters
There are also three optional parameters `repr_basis`, `msg_num`, and `kvec`. The last two are used to specify the MSG number and the conventional coordinate of the little group. When specified, the Hamiltonian is decomposed into symmetric bases of monomial functions f(k) and hermitian matrices, while if not, the Hamiltonian is decomposed simply by monomial functions.

The `repr_basis` parameter specifies what hermitian matrix basis is adopted, which can be `'pauli'` (default value), `'auto'`, or other user-specified bases. The dimension of the hermitian matrices is determined by the representation matrices.

When `repr_basis='auto'`, the natural basis of n-dimensional hermitian matrices is used.
When `repr_basis='pauli'`:
- If dim=2, the 4 Pauli matrices are used, denoted by *s_i, i=0~3*.
- If dim=4, the 16 Gamma matrices are used, denoted by *Gamma_ij, i,j=0~3*.  
  The cases of dim=2^n are similar.
- If dim=3, the 9 Gellmann matrices are used, denoted by *G_i, i=0~8*.
- If dim=6, the 4*9 Pauli matrices direct product with Gellmann matrices are used, denoted by *G_ij, i=0-3, j=0-8*. The cases of dim=3 *2^n are similar.
- For other dimensions, the natural basis of hermitian matrices is used.

Users can use the `hermitian_pauli_basis` function check the matrix form of basis:
```python
>>>import kdotp_generator as kp
>>>kp.hermitian_pauli_basis(dim=2)
[Matrix([
[1, 0],
[0, 1]]), Matrix([
[0, 1],
[1, 0]]), Matrix([
[0, -I],
[I,  0]]), Matrix([
[1,  0],
[0, -1]])]
>>> kp.hermitian_pauli_basis_symbols(dim=2)
['s_0', 's_1', 's_2', 's_3']
```


## Outputs
An example output is listed below:
```python
# ===  Result of msg 10.44  GM3d_GM4d  order = [1]  ===
Matrix([[kz, 0], [0, kz]])
Basis vector: [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

Decomposed basis using linear IR: 1 GM1  
Coefficient basis:
kz

Matrix basis:
Symbol: s_0   Expression: Matrix([[1, 0], [0, 1]]) 
```
For each independent kp basis, there are 5 parts:
- The matrix form of the Hamiltonian.
- The basis vector, written in the direct product basis, i.e., monomial basis direct product with the hermitian matrix basis. For example, when `order=1` and `dim=2`, the basis is:
  `[kx*S_0, kx*S_1, kx*S_2, kx*S_3, ky*S_0, ky*S_1, ky*S_2, ky*S_3, kz*S_0, kz*S_1, kz*S_2, kz*S_3]`
- Linear IR: the linear irreducible corepresentation (coirrep) label of the decomposed basis. The definition of linear coirreps can be found in the `MSG_kp_results` folder.
- Coefficient basis and Matrix basis: the decomposed symmetric basis. 
When fails to decompose into symmetric basis, the output will be marked as `not symmetric`. This message only means the decomposition fails, and the kp results are still correct.


## Tips
There are a few tips on the usage of this package:
- For `order>3` and `dim>8`, the computation could be very slow (~ a few hours). 
- In the package, the strain tensor epsilon $e_{uv}$ is treated as two independent polar vectors $e_u e_v$.

## Modifications from the *kdotp-symmetry* package
The original *kdotp-symmetry* package is well-written and the main program structure is maintained in the *kdop-generator* package. We remark here the major modifications:
- Generalize the input so that it can compute the effective Hamiltonians of external fields and their couplings to k.
- Add a post-processing step to decompose the Hamiltonian into symmetrical monomial function and Hermitian matrices bases using linear coirreps.
- Pre-compute effective Hamiltonians in magnetic space groups up to the third order.
- Slightly improve the efficiency of the code.

## Bug fixing
- 2022.5: We update the kp results for some type-IV MSGs, which have problematic representation matrices for some antiunitary operators.
