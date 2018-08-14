import math
import numpy as np
from utilities import spherical


class barycentric_coordinates:
    def prepare_barycentric_coordinates(self, c, vertices):
        """Return the barycentric coordinates of a ray from the origin to triangles
        defined by three points.

        arrays index
        c = [plane, point]
        vertex = [vertex, coordinate]
        plane: plane index
        point: point index, 3 points per plane, CCW
        coord: coord per point, x,y,z
        """
        # c : [triangle idx, vertex idx]
        # vertices: [ vertex idx, coord idx ]
        num_planes = c.shape[0]

        # ordered triangle corners
        # dim: [vertex_idx, coordinate] list of vectors
        self.c0 = vertices[c[:, 0]]
        self.c1 = vertices[c[:, 1]]
        self.c2 = vertices[c[:, 2]]
        # plane equations
        # compute vectors in plane
        # dim: [vector idx, coordinate] list of vectors
        self.v0 = self.c1 - self.c0
        self.v1 = self.c2 - self.c0
        # compute normals
        print(self.v0.shape, self.v1.shape)
        # dim: [vector idx, coordinate] list of vectors
        self.n = np.cross(self.v0, self.v1)
        print(self.n.shape, self.c0.shape)
        # compute distance to origin / not used now
        # self.D = -np.tensordot(self.n, self.c0, axes=(1, 1))

    def get_barycentric_coordinates(self, r):
        """Get barycentric coordinates of a vector compared to all triangles 
        in the list of triangles. 

        Returns
         dot: dot == 0: no intersection, r is parallel to plane
         dot > 0:  plane normal and ray are in same direction, we are on
                   the correct side of plane
         dot < 0:  plane normal and ray are opposed direction, ray goes away
                   from triangle

         t: r multiplier to find P
         u, v: barycentric coordinates

        t, u, v need to be normalised by 1/dot. Not done to avoid NaNs
        """
        # r: [radial vector idx, coordinates] list of vectors
        # n: [normal vector idx, coordinates] li8st of vectors
        # should return a matrix of r's compared to all n's
        # return of dim [normal, radial]

        dot = np.tensordot(self.n, r, axes=(-1, -1))
        print("dot:", dot.shape)
        #    it = dot != 0
        # intersection_t = np.zeros(dot.shape)
        # intersection_t[it] = D[it]/dot[it]

        # distance to origin: [vector idx] list of distances of one plane to its origin
        d = np.cross(self.c0, self.v0)
        print("d", d.shape)
        t_i = -np.sum(d*self.v1, axis=1)
        print("t_i:", t_i.shape)
        # u and v barycentric coordinates of each r to each triangles: [vector idx, radial]
        u_i = np.tensordot(np.cross(self.c0, self.v1), r, axes=(-1, -1))
        print("u_i:", u_i.shape)

        v_i = -np.tensordot(d, r, axes=(-1, -1))
        print("v_i:", v_i.shape)

        return dot, t_i, u_i, v_i


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

        num_points = lonlat.shape[0]

        # indexing function through rhombis
        # level
        vertices = np.zeros((num_points, 3),
                            dtype=np.float32)
        r = 1.0
        for i in range(vertices.shape[0]):
            # spherical -> theta vertical mvt, declination-> latitude
            # -> phi -> horizontal mvt, azimuth -> longitude
            vertices[i, :] = spherical(r,
                                       lonlat[i, 1],
                                       lonlat[i, 0])
        self.num_rhombi = 10

        num_points_side_region = int(pow(2, 4))
        # number of points in one subrhombi side
        self.nl_reg = num_points_side_region
        # number of points in one subrhombi region
        self.nl2 = int(pow(num_points_side_region, 2))
        # kxl = int(sqrt(num_subrhombi))
        # number of subrhombi regions on one rhombi side
        self.kxl = int(pow((num_points - 2)/10, 1/2))//num_points_side_region
        self.num_points = num_points
        print("Subdivision level:", self.g)
        print("number of points on side of rhomboid:", self.nl_reg)
        print("number of points in rhomboid:", self.nl2)
        print("number of points in grid", num_points)

        # nfaces
        # number of subrhombi in one rhombi
        self.num_subrhombi_per_face = self.kxl*self.kxl

        self.nfaces = self.num_subrhombi_per_face

        # list of triangles
        # CCW points
        self.num_triangles = self.num_rhombi * \
            self.num_subrhombi_per_face*pow(self.nl_reg, 2)*2

        self.triangles = np.ones((self.num_triangles, 3), dtype=np.int32)*-1
        self.triangle_build_idx = 0
        # triangles neighbours
        # CCW neighbours matching points above
        # edge from 0 -> 1, neighbour 0
        # edge from 1 -> 2, neighbour 1
        # edge from 2 -> 0, neighbour 2
        self.neighbours = np.ones((self.num_triangles, 3), dtype=np.int32)*-1

        self.rhombi_neighbours = [
            [4, 1, 5, 9],  # idx: 0
            [0, 2, 6, 5],  # idx: 1
            [1, 3, 7, 6],  # idx: 2
            [2, 4, 8, 7],  # idx: 3
            [3, 0, 9, 8],  # idx: 4

            [0, 1, 6, 9],  # idx: 5
            [1, 2, 7, 5],  # idx: 6
            [2, 3, 8, 6],  # idx: 7
            [3, 4, 9, 7],  # idx: 8
            [4, 0, 5, 8],  # idx: 9
        ]

        print("preparing triangle mesh")
        self.build_triangle_mesh()
        print("triangle mesh done", self.triangle_build_idx)

        print("prepare coordinates seracher")
        self.barycentric_coordinates = barycentric_coordinates()
        self.barycentric_coordinates.prepare_barycentric_coordinates(
            self.triangles, vertices)
        l = 100
        rs = np.zeros((l, 3))

        for i in range(l):
            rs[i, :] = spherical(1.0, 45.0, ((i/(l-1))*2.0 - 1)*math.pi)

        dot, t_i, u_i, v_i = self.barycentric_coordinates.get_barycentric_coordinates(
            rs)
        print("done")

    def idx_p(self, fc, kx, ky, i, j):
        """point index in subrhombus, in face fc, for subrhombus kx, ky,
        for point coordinates i,j. This does not take into account halos and poles"""
        return self.nl2*(fc*self.nfaces + ky*self.kxl + kx) + j*self.nl_reg + i

    def idx_t(self, f, kx, ky, i_, j_, n):
        """triangle points indices for face f, subrhombus kx, ky, square coord i, j
        triangle 0 or 1 (two triangles per square coord)
        for 16x16x2 triangles per subrhombus)"""
        # offset indices for triangles overlaping neighbours
        i = i_ - 1
        j = j_ - 1

        fc_t_l = self.rhombi_neighbours[f][0]
        fc_b_l = self.rhombi_neighbours[f][3]
        # index for border point
        cn = self.nl_reg - 1
        # TODO check orders of edges
        if kx == 0 and ky == 0 and i < 0 and j < 0:
            if f < 5:
                if n == 0:
                    i_c = self.idx_p(f, 0, 0, 0, 0)
                    i_c_t = self.idx_p(fc_b_l, 0, self.kxl-1, 0, cn)
                    i_c_b = self.idx_p(fc_t_l,
                                       self.kxl - 1, self.kxl - 1, cn, cn)
                    return np.array((i_c, i_c_t, i_c_b), dtype=np.int32)
                else:

                    i1 = self.idx_p(f, 0, self.kxl - 1, 0, cn)

                    i2 = self.num_points - 2

                    i3 = self.idx_p(fc_t_l,
                                    0,
                                    self.kxl - 1,
                                    0,
                                    cn)

                    return np.array((i1, i2, i3), dtype=np.int32)

            else:
                if n == 0:
                    i_c = self.idx_p(f, 0, 0, 0, 0)
                    i_c_t = self.idx_p(fc_t_l, self.kxl-1, 0, cn, 0)
                    i_c_b = self.idx_p(fc_b_l,
                                       self.kxl - 1, self.kxl - 1, cn, cn)
                    return np.array((i_c, i_c_t, i_c_b), dtype=np.int32)
                else:
                    i1 = self.idx_p(f, self.kxl - 1, 0, cn, 0)

                    i3 = self.num_points - 1

                    i2 = self.idx_p(fc_b_l,
                                    self.kxl - 1,
                                    0,
                                    cn,
                                    0)

                    return np.array((i1, i2, i3), dtype=np.int32)

        elif kx == 0 and i < 0 and j < 0:
            if f < 5:
                i1 = self.idx_p(f, kx, ky, 0, 0)
                i4 = self.idx_p(f, kx, ky-1, 0, cn)
                i3 = self.idx_p(fc_t_l, self.kxl-1 - ky+1, self.kxl-1,
                                0, cn)
                i2 = self.idx_p(fc_t_l, self.kxl-1 - ky, self.kxl-1,
                                cn, cn)

                if n == 0:
                    return np.array((i1, i2, i3), dtype=np.int32)
                else:
                    return np.array((i1, i3, i4), dtype=np.int32)
            else:
                i1 = self.idx_p(f, kx, ky, 0, 0)
                i3 = self.idx_p(fc_t_l, self.kxl-1, ky-1, cn, cn)
                i4 = self.idx_p(f, kx, ky-1, 0, cn)
                i2 = self.idx_p(fc_t_l, self.kxl-1, ky, cn, 0)

                if n == 0:
                    return np.array((i1, i2, i3), dtype=np.int32)
                else:
                    return np.array((i1, i3, i4), dtype=np.int32)
        elif kx == 0 and i < 0:

            if f < 5:
                i1 = self.idx_p(f, kx, ky, 0, j+1)
                i4 = self.idx_p(f, kx, ky, 0, j)

                i3 = self.idx_p(fc_t_l,
                                self.kxl - 1 - ky,
                                self.kxl - 1,
                                self.nl_reg - 1 - (j),
                                cn)

                i2 = self.idx_p(fc_t_l,
                                self.kxl - 1 - ky,
                                self.kxl - 1,
                                self.nl_reg - 1 - (j+1),
                                cn)

                if n == 0:
                    return np.array((i1, i2, i3), dtype=np.int32)
                else:
                    return np.array((i3, i4, i1), dtype=np.int32)
            else:
                i1 = self.idx_p(f, kx, ky, 0, j+1)
                i4 = self.idx_p(f, kx, ky, 0, j)

                i3 = self.idx_p(fc_t_l,
                                self.kxl - 1,
                                ky,
                                cn,
                                j)
                i2 = self.idx_p(fc_t_l,
                                self.kxl - 1,
                                ky,
                                cn,
                                j+1)

                if n == 0:
                    return np.array((i1, i2, i3), dtype=np.int32)
                else:
                    return np.array((i3, i4, i1), dtype=np.int32)

        elif ky == 0 and j < 0 and i < 0:
            # done up to here CCW check
            if f < 5:
                i1 = self.idx_p(f, kx, ky, 0, 0)
                i2 = self.idx_p(f, kx-1, ky, cn, 0)
                i3 = self.idx_p(fc_b_l, kx-1, self.kxl-1, cn, cn)
                i4 = self.idx_p(fc_b_l, kx, self.kxl-1, 0, cn)

                if n == 0:
                    return np.array((i1, i2, i3), dtype=np.int32)
                else:
                    return np.array((i3, i4, i1), dtype=np.int32)
            else:
                i1 = self.idx_p(f, kx, ky, 0, 0)
                i2 = self.idx_p(f, kx-1, ky, cn, 0)
                i3 = self.idx_p(fc_b_l, self.kxl-1, self.kxl-1 - kx + 1,
                                cn, 0)
                i4 = self.idx_p(fc_b_l, self.kxl-1, self.kxl-1 - kx,
                                cn, cn)

                if n == 0:
                    return np.array((i3, i1, i2), dtype=np.int32)
                else:
                    return np.array((i1, i3, i4), dtype=np.int32)

        elif ky == 0 and j < 0:
            if f < 5:
                i1 = self.idx_p(f, kx, ky, i+1, 0)

                i3 = self.idx_p(fc_b_l,
                                kx,
                                self.kxl - 1,
                                i,
                                cn)

                i4 = self.idx_p(fc_b_l,
                                kx,
                                self.kxl - 1,
                                i+1,
                                cn)
                i2 = self.idx_p(f, kx, ky, i, 0)
                if n == 0:
                    return np.array((i1, i2, i3), dtype=np.int32)
                else:
                    return np.array((i3, i4, i1), dtype=np.int32)
            else:
                i1 = self.idx_p(f, kx, ky, i+1, 0)
                i2 = self.idx_p(f, kx, ky, i, 0)

                i3 = self.idx_p(fc_b_l,
                                self.kxl - 1,
                                self.kxl - 1 - kx,
                                cn,
                                cn - i)
                i4 = self.idx_p(fc_b_l,
                                self.kxl - 1,
                                self.kxl - 1 - kx,
                                cn, cn - (i+1))

                if n == 0:
                    return np.array((i3, i1, i2), dtype=np.int32)
                else:
                    return np.array((i1, i3, i4), dtype=np.int32)
        else:
            # non border rhombus
            i1 = 0
            i2 = 0
            i3 = 0
            i4 = 0
            if i < 0 and j < 0:
                # corner
                i1 = self.idx_p(f, kx-1, ky-1, cn, cn)
                i2 = self.idx_p(f, kx, ky-1, 0, cn)
                i3 = self.idx_p(f, kx, ky, 0, 0)
                i4 = self.idx_p(f, kx-1, ky, cn, 0)
            elif i < 0:
                # x axis
                i1 = self.idx_p(f, kx-1, ky, cn, j)
                i2 = self.idx_p(f, kx, ky, 0, j)
                i3 = self.idx_p(f, kx, ky, 0, j+1)
                i4 = self.idx_p(f, kx-1, ky, cn, j+1)
                #i1, i2, i3, i4 = (0, 0, 0, 0)
            elif j < 0:
                # y axis
                i1 = self.idx_p(f, kx, ky-1, i, cn)
                i2 = self.idx_p(f, kx, ky-1, i+1, cn)
                i3 = self.idx_p(f, kx, ky, i+1, 0)
                i4 = self.idx_p(f, kx, ky, i, 0)
            else:
                # inside rhombus
                i1 = self.idx_p(f, kx, ky, i, j)
                i2 = self.idx_p(f, kx, ky, i+1, j)
                i3 = self.idx_p(f, kx, ky, i+1, j+1)
                i4 = self.idx_p(f, kx, ky, i, j+1)
                # if f == 9:
                #     i1, i2, i3, i4 = (0, 0, 0, 0)
            # TODO: check that they are CCW
            if n == 0:
                return np.array((i1, i2, i3), dtype=np.int32)
            else:
                return np.array((i3, i4, i1), dtype=np.int32)

    def build_triangle_mesh_for_subrhombus(self, face, kx, ky):
        for i in range(self.nl_reg):
            for j in range(self.nl_reg):
                # make two triangles
                idx = self.triangle_build_idx

                self.triangles[idx+0, :] = self.idx_t(face, kx, ky, i, j, 0)
                self.triangles[idx+1, :] = self.idx_t(face, kx, ky, i, j, 1)
                # TODO: make neighbours list
                self.triangle_build_idx += 2

    def build_triangle_mesh_for_rhombus(self, face):
        for kx in range(self.kxl):
            for ky in range(self.kxl):
                self.build_triangle_mesh_for_subrhombus(face, kx, ky)

    def build_triangle_mesh(self):
        self.triangle_build_idx = 0
        """build a list of triangles, with links to corners and neighbours"""
        for face in range(10):
            self.build_triangle_mesh_for_rhombus(face)
