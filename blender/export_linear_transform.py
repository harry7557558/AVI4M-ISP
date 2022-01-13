import bpy
import numpy as np
import scipy.optimize
from math import sin, cos


EPSILON = 1e-6
PI = 3.1415926535897932384626


def clean_round_arr(arr, th=0.):
    """Set the element of the array to @th if the difference is less than EPSILON """
    for i in range(len(arr)):
        if abs(arr[i]-th) < EPSILON:
            arr[i] = th
    return arr


def string_join_array(arr):
    """array to comma-separated string"""
    return ','.join(["{:.2f}".format(t) for t in arr])


def get_euler_angles(rotmat, rotation_mode: str):
    """Get the Euler angles of a 3x3 rotation matrix
       I'm too lazy to write a function for each rotation mode so I just call optimization
    """
    assert rotmat.shape == (3, 3)
    assert rotation_mode.upper() in ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX']

    def loss_function(rxyz):
        rx, ry, rz = rxyz
        matx = [[1., 0., 0.],
                [0., cos(rx), -sin(rx)],
                [0., sin(rx), cos(rx)]]
        maty = [[cos(ry), 0., sin(ry)],
                [0., 1., 0.],
                [-sin(ry), 0., cos(ry)]]
        matz = [[cos(rz), -sin(rz), 0.],
                [sin(rz), cos(rz), 0.],
                [0., 0., 1.]]
        mat = np.identity(3)
        for axis in rotation_mode:
            if axis == 'X':
                mat = np.matmul(matx, mat)
            if axis == 'Y':
                mat = np.matmul(maty, mat)
            if axis == 'Z':
                mat = np.matmul(matz, mat)
        return np.sum((mat-rotmat)**2.0)

    a0 = [0., 0., 0.]
    a = scipy.optimize.minimize(loss_function, a0, method='Nelder-Mead',
                                options={'xatol': 1e-2*EPSILON, 'fatol': (1e-2*EPSILON)**2})
    a = np.average(a['final_simplex'][0], (0))
    loss = loss_function(a)
    if not loss < EPSILON:
        print(rotmat)
        print(a, loss)
        assert False

    a = (a+PI) % (2.0*PI)-PI
    return clean_round_arr(a)


def get_transform_code(obj):
    """Return the transformation information of an object as code string"""

    # get name
    name = obj.name
    if '.' in name:
        name = name[:name.find('.')]

    # get matrix
    matrix = np.array(obj.matrix_world)
    matrix /= matrix[3][3]   # normalize scaling
    assert not np.any(matrix[3][0:3])  # no perspective
    assert matrix[3][3] == 1.0

    # get translation components
    translate = matrix.T[3][0:3]

    # get scaling components
    matrix3 = matrix[0:3, 0:3]
    scale = np.zeros((4))  # |sx|, |sy|, |sz|, reflection
    scale[0] = np.linalg.norm(matrix3.T[0])
    scale[1] = np.linalg.norm(matrix3.T[1])
    scale[2] = np.linalg.norm(matrix3.T[2])
    determinant = np.linalg.det(matrix3)
    assert abs(determinant) > EPSILON
    scale[3] = determinant / np.prod(scale[0:3])
    assert 1.0-EPSILON < abs(scale[3]) < 1.0+EPSILON
    clean_round_arr(scale, 1.)
    clean_round_arr(scale, -1.)

    # get rotation components
    rotmat = np.matmul(matrix3, np.diag(scale[3]/scale[0:3]))
    assert np.linalg.norm(np.matmul(rotmat, rotmat.T) -
                          np.identity(3)) < EPSILON  # orthogonal matrix
    rotate_mode = 'XYZ'
    angles = get_euler_angles(rotmat, rotate_mode)

    # export code
    code_t = "(p-vec3(" + string_join_array(translate) + "))"
    code_r = ""
    for i in range(3):
        if abs(angles[i]) == 0.0:
            continue
        axis = rotate_mode[i]
        code_r += "rot"+axis.lower() + \
            "(" + string_join_array([-angles[i]/PI]) + "*PI)*"
    code_s = "/vec3(" + string_join_array(scale[3]*scale[0:3]) + ")"
    code_s1 = "*vec4(1,1,1,"+string_join_array([min(scale[0:3])])+")"
    code = name+"(("+code_r+code_t+")"+code_s+", col_required)"+code_s1
    code_full = "d = smin(d, " + code + ", 0.01);"
    print(code_full)


if __name__ == "__main__":

    print("======== Get Transform ========")

    for obj in bpy.data.objects:
        if obj.name in ["Camera", "Light"]:
            continue
        if not obj.visible_get():
            continue
        get_transform_code(obj)

    print(end='\n')
