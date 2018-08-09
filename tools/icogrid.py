
import numpy as np
from utilities import spherical


class ico:
    def pt_in_quad(quad, p):
        # computes if a point is in the quadrilateral

        return False

    def __init__(self, g, lonlat):
        """computes icosahedron subdivided by g levels"""
        self.g = g
        self.n_s = pow(2, g)
        self.n_r = self.n_s*self.n_s
        self.lonlat = lonlat

        print("Subdivision level:", self.g)
        print("number of points on side of rhomboid, n_s:", self.n_s)
        print("number of points in rhomboid, n_r:", self.n_r, )
        print("number of points in grid, 10*n_r+2:", 10*self.n_r+2)
        print("halos: ", 10*4*(self.n_s+2))

        print("maps: ", 10*pow(self.n_s+2, 2))

        # indexing function through rhombis
        # level

        g = int(pow((num_points - 2)/10, 1/4)) - 2
        num_rhombi = 10

        num_points_side_region = int(pow(2, 4))
        nl_reg = num_points_side_region
        nl2 = int(pow(num_points_side_region, 2))
        # kxl = int(sqrt(num_subrhombi))
        kxl = int(pow((num_points - 2)/10, 1/2))//num_points_side_region

        # nfaces
        num_subrhombi = kxl*kxl
        nfaces = num_subrhombi

        def idx(fc, kx, ky, i, j):
            return nl2*(fc*nfaces + ky*kxl + kx) + j*nl_reg + i
