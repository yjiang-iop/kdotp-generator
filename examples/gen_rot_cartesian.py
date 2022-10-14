import sympy as sp
from sympy import pi, sin, cos, sqrt

def get_rotation(angle, axis, det):
    '''
    This function is used to generate O(3) matrix in Cartesian coordinates.
    Inputs: rotation angle (-2*pi to 2*pi), rotation axis (length-3 array), and determinate (+1 or -1)
    
    Formula:
    R[i,j] = cos(w) delta_ij + (1 - cos(w)) * ni * nj - sin(w) sum(epsilon_ijk nk)
    '''
    assert -2 * pi <= angle <= 2 * pi, ('angle should be in radian', angle)
    assert det == 1 or det == -1, det
    n = sp.Matrix(axis) / sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)

    R = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            if i == j:
                R[i, j] += cos(angle)
            R[i, j] += (1 - cos(angle)) * n[i] * n[j]
    R += -sin(angle) * sp.Matrix([[0, n[2], -n[1]], [-n[2], 0, n[0]], [n[1], -n[0], 0]])
    R *= det
    return R
            

if __name__ == '__main__':
    # two examples:

    # C6z
    angle = pi/3
    axis = [0, 0, 1]
    det = 1
    R = get_rotation(angle, axis, det)
    print(R)
    
    # M_{1,-1,0}
    angle = pi
    axis = [1, -1, 0]
    det = -1
    R = get_rotation(angle, axis, det)
    print(R)
