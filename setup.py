# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
# Adapted from the kdotp-symmetry package by: Dominik Gresch <greschd@gmx.ch>  © 2017-2018, ETH Zurich, Institut für Theoretische Physik
"""
Usage: pip install -e .[dev]
"""

import re
import sys
if sys.version_info < (3, 6):
    raise 'must use Python version 3.6 or higher'

from setuptools import setup

README = """A tool for computing k.p effective Hamiltonians with couplings to external fields including E, B, and epsilon, under given symmetry constraints."""

with open('./kdotp_generator/__init__.py', 'r') as f:
    MATCH_EXPR = "__version__[^'\"]+(['\"])([^'\"]+)"
    VERSION = re.search(MATCH_EXPR, f.read()).group(2)

setup(
    name='kdotp-generator',
    version=VERSION,
    url='https://github.com/yjiang-iop/kdotp-generator',
    author='Yi Jiang, Dominik Gresch',
    author_email='jiangyi14@mails.ucas.ac.cn, greschd@gmx.ch',
    description=
    'A tool for computing k.p effective Hamiltonians.',
    install_requires=[
        'sympy', 'numpy', 'scipy', 'networkx>=2'
    ],
    python_requires=">=3.6",
    extras_require={
        'dev': [
            'pytest', 'yapf==0.29', 'pre-commit', 'prospector==1.2.0',
            'pylint==2.4.4', 'sphinx', 'sphinx_rtd_theme'
        ]
    },
    long_description=README,
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English', 'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Development Status :: 4 - Beta'
    ],
    license='Apache 2.0',
    packages=['kdotp_generator', 'kdotp_generator.MSG_linear_coir_data'],
    include_package_data=True
)
