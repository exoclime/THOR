import math
import numpy as np


def HSV_to_RGB(H, S, V):
    """
    Computes an rgb value with components between [0,1] from
    H [0:360], S [0:1], V [0:1] values
    """
    C = V*S
    H_p = H/60.0
    X = C*(1.0-abs((H_p % 2) - 1.0))

    R1G1B1 = [0.0, 0.0, 0.0]
    if 0.0 <= H_p and H_p <= 1.0:
        R1G1B1 = [C, X, 0.0]
    elif 1.0 <= H_p and H_p <= 2.0:
        R1G1B1 = (X, C, 0.0)
    elif 2.0 <= H_p and H_p <= 3.0:
        R1G1B1 = (0.0, C, X)
    elif 3.0 <= H_p and H_p <= 4.0:
        R1G1B1 = (0.0, X, C)
    elif 4.0 <= H_p and H_p <= 5.0:
        R1G1B1 = (X, 0.0, C)
    elif 5.0 <= H_p and H_p < 6.0:
        R1G1B1 = (C, 0.0, X)

    m = V-C

    RGB = (R1G1B1[0] + m,
           R1G1B1[1] + m,
           R1G1B1[2] + m)

    return RGB


def spherical(r, theta, phi):
    return np.array((
        r*math.cos(theta)*math.cos(phi),
        r*math.sin(theta)*math.cos(phi),
        r*math.sin(phi)),
        dtype=np.float32)
