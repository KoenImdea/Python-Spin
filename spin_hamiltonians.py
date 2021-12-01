"""
Based on the code by M. Ternes
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


Defenition of the Hamiltonians goes here
"""
import numpy as np
j = complex(0,1)
from scipy.linalg import expm

class Spin():
    def __init__(self):
        pass


"""
////////////////////////////////////////////////////////////////////////////////
// Basic spin operators
////////////////////////////////////////////////////////////////////////////////
"""

def s_plus(Spin):
    """
    prepare the S+ matrix
    Spin.S must be the spin of the system
    """
    a = np.flip(np.arange(1,((2*Spin.S)+1)))
    e = np.sqrt((np.arange(1,((2*Spin.S)+1))* a))
    return np.diag(e, +1)

def s_minus(Spin):
    """
    prepare the S- matrix
    Spin.S must be the spin of the system
    """
    a = np.flip(np.arange(1,((2*Spin.S)+1)))
    e = np.sqrt((np.arange(1,((2*Spin.S)+1))* a))
    return np.diag(e, -1)

def s_x(Spin):
    """
    prepare the Sx matrix
    Spin.S must be the spin of the system
    """
    return 0.5*(s_plus(Spin)+s_minus(Spin))

def s_y(Spin):
    """
    prepare the Sy matrix
    Spin.S must be the spin of the system
    """
    return 0.5*j*(s_plus(Spin)-s_minus(Spin))

def s_z(Spin):
    """
    prepare the Sz matrix with diagonal elements -S, -S+1,...,+S
    Spin.S must be the spin of the system
    """
    e = np.arange(-Spin.S, Spin.S+1)
    return np.diag(np.flip(e))

def s_arb(Spin, varargin):
    """
    prepare a S matrix with arbitary directions
    S.S must be the spin of the system
    varargin must contain either the rotation angles theta and phi
    or the direction (x,y,z)
    """
    if len(varangin) > 2 then:
        x = varangin[0]
        y = varangin[0]
        z = varangin[0]
        theta = np.arccos(z/np.sqrt(x**2+y**2+z**2))
        if np.isnan(theta):
            theta = 0
        phi = np.arctan2(y, x)
    else:
        theta = varangin[0]
        phi = varangin[1]
    e=expm(-j*s_y(Spin)*theta)*s_z(Spin)*expm(i*s_y(Spin)*theta)
    return expm(-j*s_z(Spin)*phi)*e*expm(j*s_z(Spin)*phi)

def s_1(Spin):
    """
    prepare the unity matrix for system S
    S.S must be the spin of the system
    """
    return np.diag(np.ones((Spin.S*2)+1))

def s_2(Spin):
    """
    prepare the S^2 matrix
    S.S must be the spin of the system
    """
    return s_x(Spin)*s_x(Spin)+s_y(Spin)*s_y(Spin)+s_z(Spin)*s_z(Spin)

def s_rot(Spin, phi, n):
    """
    prepare the rotational matrix
    of the spin system S
    rotation angle phi
    rotational axis n (must have length 1)
    """
    e = expm(-j*phi*(s_x(Spin)*n[0]+s_y(Spin)*n[1]+s_z(Spin)*n[2]))
    e[np.abs(e) < 1E-10] = 0
    return e

def dmatrix(n):
    """
    create the Density Matrix for spin 1/2 electrons which have a polarization
    in the direction n and a polarization strength -1<|n|<1.
    """
    S = Spin()
    S.S=1/2
    return (1/2*s_1(S)+n[0]*s_x(S)+n[1]*s_y(S)+n[2]*s_z(S))
