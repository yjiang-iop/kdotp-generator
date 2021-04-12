# *kdotp-generator* 
This package is used to automatically generate $k\cdot p$ effective Hamilonians, developed based on the [*kdotp-symmetry*](https://kdotp-symmetry.greschd.ch/en/latest/) package by D. Gresch. 

To use this package, users can either input their symmetry data and compute the effective Hamiltonians, or directly refer to our pre-computed results for magnetic space groups (MSGs).

If you are using this package, please cite the following works:
- [*A $k\cdot p$ effective Hamiltonian generator*](), Yi Jiang, Chen Fang.
- [*Identifying Topological Semimetals*](https://www.research-collection.ethz.ch/handle/20.500.11850/308602), D. Gresch, PhD Thesis.

## Installation
Pip is used to install the package:
> pip install kdotp-generator 

The pre-computed results for MSGs are not included in the pip installation. They can be downloaded from the Github page of [the package](https://github.com/yjiang-iop/kdotp-generator).

## Example
Two examples are given in the *examples* folder, i.e., *input1.py* and *input2.py*, which calculate the effective Hamiltonians for $k$ and $k, E$, respectively, for non-magnetic space group 123 $P4/mmm$, representation $\bar{\Gamma}_6, \bar{\Gamma}_6$, up to the second order. 




## Tips


​
## Modifications from the *kdotp-symmetry* package
The original *kdotp-symmetry* package is well-written. We follow its main program structure and make the following modifications:
- Generalize the input so that it can compute the effective Hamiltonians of external fields and their couplings to $k$.
- Add a post-processing step to decompose the effective Hamiltonian symmetrically using linear coirreps.
- Slightly improve the efficiency of the code.
- Pre-compute effective Hamiltonians in magnetic space groups up to the third order.
