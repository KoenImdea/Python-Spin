"""
Based on the code by M. Ternes
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


Defenition of the Hamiltonians goes here
"""
import numpy as np

"""
////////////////////////////////////////////////////////////////////////////////
// Basic spin operators
////////////////////////////////////////////////////////////////////////////////
"""

def s_plus(Spin):
    """
    prepare the S+ matrix
    S.S must be the spin of the system
    """
    a = np.flip(np.arange(1,((2*Spin.S)+1)))
    e = np.sqrt(np.outer(np.arange(1,((2*Spin.S)+1)), a))
    return np.diagonal(e,-1)
