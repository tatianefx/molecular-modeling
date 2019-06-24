#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from math import sqrt, acos, pi
from scipy.linalg import expm
import numpy as np
import math


def calculates_unit_norm(v):
    x = v[0]**2
    y = v[1]**2
    z = v[2]**2

    return sqrt(x + y + z)


def normalize_vector(v):
    c = calculates_unit_norm(v)

    x = v[0]/c
    y = v[1]/c
    z = v[2]/c

    return [x, y, z]


def calculates_distance_a_b(a, b):
    x = (a[0] - b[0])**2
    y = (a[1] - b[1])**2
    z = (a[2] - b[2])**2

    return sqrt(x + y + z)


def calculates_angle_a_b_c(a, b, c):
    ab = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]

    ab_vec = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2])
    bc_vec = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2])

    ab_norm = [ab[0] / ab_vec, ab[1] / ab_vec, ab[2] / ab_vec]
    bc_norm = [bc[0] / bc_vec, bc[1] / bc_vec, bc[2] / bc_vec]

    res = ab_norm[0] * bc_norm[0] + ab_norm[1] * bc_norm[1] + ab_norm[2] * bc_norm[2]

    return acos(res) * 180.0 / pi


# p1 = [-1.4105, 1.1507, 0.1821]
# p2 = [-0.7085, -0.1136, 0.3937]
# p3 = [0.7470, 0.0903, 0.0308]
# p4 = [1.2492, 1.1165, -0.4047]
# torsion_angle = calculates_torsion_angle(p1, p2, p3, p4)
# print(torsion_angle)
def calculates_torsion_angle(p1, p2, p3, p4):
    """ Function to calculate the torsion angle """

    # Call calc_q_vectors(p1,p2,p3,p4) function
    q1, q2, q3 = calc_q_vectors(p1, p2, p3, p4)
    # Call calc_cross_vectors(q1,q2,q3) function
    q1_x_q2, q2_x_q3 = calc_cross_vectors(q1, q2, q3)
    # Call calc_nornals(q1_x_q2,q2_x_q3) function
    n1, n2 = calc_nornals(q1_x_q2, q2_x_q3)
    # Call calc_orthogonal_unit_vectors(n2,q2) function
    u1, u2, u3 = calc_orthogonal_unit_vectors(n2, q2)
    # Call calc_dihedral_angle(u1,u2,u3) function
    return calc_dihedral_angle(n1, u1, u2, u3)


def calc_q_vectors(p1, p2, p3, p4):
    """Function to calculate q vectors"""

    # Calculate coordinates for vectors q1, q2 and q3
    q1 = np.subtract(p2, p1) # b - a
    q2 = np.subtract(p3, p2) # c - b
    q3 = np.subtract(p4, p3) # d - c
    return q1, q2, q3


def calc_cross_vectors(q1, q2, q3):
    """Function to calculate cross vectors"""

    # Calculate cross vectors
    q1_x_q2 = np.cross(q1, q2)
    q2_x_q3 = np.cross(q2, q3)
    return q1_x_q2, q2_x_q3


def calc_nornals(q1_x_q2, q2_x_q3):
    """Function to calculate normal vectors to planes"""

    # Calculate normal vectors
    n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2, q1_x_q2))
    n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3, q2_x_q3))
    return n1, n2


def calc_orthogonal_unit_vectors(n2, q2):
    """Function to calculate orthogonal unit vectors"""

    # Calculate unit vectors
    u1 = n2
    u3 = q2/(np.sqrt(np.dot(q2, q2)))
    u2 = np.cross(u3, u1)
    return u1, u2, u3


def calc_dihedral_angle(n1, u1, u2, u3):
    """Function to calculate dihedral angle"""

    # Calculate cosine and sine
    cos_theta = np.dot(n1,u1)
    sin_theta = np.dot(n1,u2)
    # Calculate theta
    theta = -math.atan2(sin_theta,cos_theta) # it is different from Fortran math.atan2(y,x)
    theta_deg = np.degrees(theta)
    # Show results
    print("theta (rad) = %8.3f"%theta)
    print("theta (deg) = %8.3f"%theta_deg)

    return theta_deg


def calculates_xyz_to_rotate(a, b):
    ab = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

    ab_vec = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2])

    ab_norm = [ab[0] / ab_vec, ab[1] / ab_vec, ab[2] / ab_vec]

    return ab_norm


def rotation_euler(v, xyz):
    """
    # https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    Rotate vector v (or array of vectors) by the euler angles xyz

    :param v: vector or array of vectors
    :param xyz: vector euler angles xyz
    :return: vector or array of vectors
    """

    for theta, axis in zip(xyz, np.eye(3)):
        v = np.dot(np.array(v), expm(np.cross(np.eye(3), axis*-theta)))
    return v.tolist()


def translation(m, xyz):
    for v in m:
        v[0] = v[0] + xyz[0]
        v[1] = v[1] + xyz[1]
        v[2] = v[2] + xyz[2]

    return m


def resulting_vector(a, b):
    ab = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    return ab


def resulting_vector_with_distance(a, b, distance):
    # m1 = [
    #     [-1, -1, 0],
    #     [0, 0, 0],
    #     [1, 1, 0]
    # ]
    #
    # m2 = [
    #     [0, -2, 0],
    #     [1, -1, 0],
    #     [2, 0, 0]
    # ]
    #
    # v = vector_result_with_distance(m1[2], m2[0], 1)
    # print(v)
    #
    # m = translation(m2, v)
    # print(m)
    #
    # >> [2, 4, 0]
    # >> [[2, 2, 0], [3, 3, 0], [4, 4, 0]]

    v = resulting_vector(a, b)

    v[0] = (a[0] - v[0]) + distance
    v[1] = (a[1] - v[1]) + distance
    v[2] = (a[2] - v[2]) + distance

    return v





