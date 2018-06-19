
from PyQt5.QtGui import (QOpenGLShader,
                         QOpenGLShaderProgram,
                         QMatrix4x4,
                         QVector4D,
                         QPainterPath,
                         QFont,
                         QSurfaceFormat)

import OpenGL.GL as gl
import numpy as np

import math
from math import cos, sin, pi, pow, sqrt

import matplotlib.pyplot as plt


class ico:
    def __init__(self, g):
        """computes icosahedron subdivided by g levels"""
        self.g = g
        self.n_s = pow(2, g)
        self.n_r = self.n_s*self.n_s

        print("Subdivision level:", self.g)
        print("number of points on side of rhomboid, n_s:", self.n_s)
        print("number of points in rhomboid, n_r:", self.n_r, )
        print("number of points in grid, 10*n_r+2:", 10*self.n_r+2)
        print("halos: ", 10*4*(self.n_s+2))

        print("maps: ", 10*pow(self.n_s+2, 2))


def spherical(r, theta, phi):
    return np.array((
        r*math.cos(theta)*math.cos(phi),
        r*math.sin(theta)*math.cos(phi),
        r*math.sin(phi)),
        dtype=np.float32)


class IcoGridPainter:
    def __init__(self):
        self.level = 0

    def set_grid_data(self, lonlat, colors):
        self.lonlat = lonlat
        self.colors = colors

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def initializeGL(self):

        # create

        # num_vertices = 12

        # vertices = np.zeros((num_vertices, 3),
        #                     dtype=np.float32)
        num_points = int(len(self.lonlat)/2)
        vertices = np.zeros((num_points, 3),
                            dtype=np.float32)

        # r = 1.2
        # w = 2.0*math.acos(1.0/(2.0*sin(pi/5.0)))

        # angles = [
        #     [0.0, pi/2.0],
        #     [0.0, -pi/2.0],
        #     [-pi/5.0, pi/2.0-w],
        #     [pi/5.0, pi/2.0-w],
        #     [3.0*pi/5.0, pi/2.0-w],
        #     [pi, pi/2.0-w],
        #     [-3.0*pi/5.0, pi/2.0-w],
        #     [0.0, -(pi/2.0-w)],
        #     [2.0*pi/5.0, -(pi/2.0-w)],
        #     [4.0*pi/5.0, -(pi/2.0-w)],
        #     [-4.0*pi/5.0, -(pi/2.0-w)],
        #     [-2.0*pi/5.0, -(pi/2.0-w)]
        # ]

        # for i, a in enumerate(angles):
        #     vertices[i, :] = spherical(r, a[0], a[1])

        vertices_colors = np.zeros((vertices.shape[0], 3), dtype=np.float32)
        elements = np.zeros(vertices.shape[0], dtype=np.uint32)

        # indexing function through rhombis
        # level
        g = int(pow((num_points - 2)/10, 1/4))
        num_rhombi = 10
        # nfaces
        num_subrhombi = int(pow(4.0, g - 4))
        nfaces = num_subrhombi
        num_points_side_region = int(pow(2, 4))
        nl_reg = num_points_side_region
        nl2 = int(pow(num_points_side_region, 2))
        kxl = int(sqrt(num_subrhombi))

        def idx(fc, kx, ky, i, j):
            return nl2*(fc*nfaces + ky*kxl + kx) + j*nl_reg + i

        # num segments in rhombis
        num_segments = int(2*(pow(2, g)*(pow(2, g)+1)))*num_rhombi
        # num segments in halo

        print("number of segments:", num_segments)

        # build line array for rhombis
        lines = np.zeros((num_segments, 2), dtype=np.uint32)

        lines_idx = 0
        for fc in range(num_rhombi):
            # sub rhombis
            for kx in range(kxl):
                for ky in range(kxl):
                    # inside one rombi
                    # horizontal
                    for i in range(nl_reg):
                        for j in range(nl_reg-1):
                            i1 = idx(fc, kx, ky, i, j)
                            i2 = idx(fc, kx, ky, i, j+1)

                            lines[lines_idx][0] = i1
                            lines[lines_idx][1] = i2

                            lines_idx += 1
                    # vertical
                    for i in range(nl_reg-1):
                        for j in range(nl_reg):
                            i1 = idx(fc, kx, ky, i, j)
                            i2 = idx(fc, kx, ky, i+1, j)

                            lines[lines_idx][0] = i1
                            lines[lines_idx][1] = i2

                            lines_idx += 1
        print("number of lines", lines_idx)

        # build triangles
        # triangles in rhombi and subrhombi
        num_triangles = int(pow(nl_reg - 1, 2))*2*num_subrhombi*num_rhombi
        # halos
        # each subrhombus takes care of two of its borders,
        # nl_reg-1 squares, * 2 for triangles, * 2 for two borders

        # neighbours of major rhombi face indexes
        # use clockwhise addressing, top left, top right, bottom right, bottom left
        rhombi_neighbours = [
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
        num_triangles += 2*(nl_reg-1)*2*num_subrhombi*num_rhombi
        # TODO: corners
        print("number of triangles", num_triangles)

        triangles = np.zeros((num_triangles, 3), dtype=np.uint32)

        draw_rhomb = [0, 1, 2, 3, 4, 5, 6, 7, 8,  9]
        triangle_idx = 0
        for fc in range(num_rhombi):
            if fc not in draw_rhomb:
                continue
            # sub rhombis
            for kx in range(kxl):
                for ky in range(kxl):
                    # inside one sub-rombi
                    for i in range(nl_reg-1):
                        for j in range(nl_reg-1):
                            i1 = idx(fc, kx, ky, i, j)
                            i2 = idx(fc, kx, ky, i, j+1)
                            i3 = idx(fc, kx, ky, i+1, j)
                            i4 = idx(fc, kx, ky, i+1, j+1)
                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i2
                            triangles[triangle_idx][2] = i4

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i4
                            triangles[triangle_idx][2] = i3
                            triangle_idx += 1
                    # sub rhombi halo
                    # top left halo
                    fc_top = fc
                    kx_top = kx - 1
                    if kx == 0 and fc < 5:
                        fc_top = rhombi_neighbours[fc][0]
                        kx_top = kxl - 1
                        print(kx, ky, kx_top)

                        for i in range(nl_reg-1):
                            i1 = idx(fc, kx, ky, 0, i)
                            i2 = idx(fc, kx, ky, 0, i+1)
                            i3 = idx(fc_top, kxl-1-ky, kxl-1,
                                     nl_reg - 1 - (i), nl_reg-1)
                            i4 = idx(fc_top, kxl-1-ky, kxl-1,
                                     nl_reg - 1 - (i+1), nl_reg-1)
                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i2

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i2
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i4
                            triangle_idx += 1
                    elif kx == 0 and fc > 4:
                        fc_top = rhombi_neighbours[fc][0]
                        for i in range(nl_reg-1):
                            i1 = idx(fc, kx, ky, 0, i)
                            i2 = idx(fc, kx, ky, 0, i+1)
                            i3 = idx(fc_top, kxl-1, ky, nl_reg - 1, i)

                            i4 = idx(fc_top, kxl-1, ky, nl_reg - 1, i+1)

                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i2

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i2
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i4
                            triangle_idx += 1
                    else:
                        for i in range(nl_reg-1):
                            i1 = idx(fc, kx, ky, 0, i)
                            i2 = idx(fc, kx, ky, 0, i+1)
                            i3 = idx(fc, kx-1, ky, nl_reg-1, i)
                            i4 = idx(fc, kx-1, ky, nl_reg-1, i+1)

                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i2

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i2
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i4
                            triangle_idx += 1
                    if ky == 0 and fc > 4:
                        fc_bot = rhombi_neighbours[fc][3]
                        for i in range(nl_reg-1):
                            i1 = idx(fc, kx, ky,  i, 0)
                            i2 = idx(fc, kx, ky,  i+1, 0)
                            i3 = idx(fc_bot, kxl-1, kxl-1-kx,
                                     nl_reg-1, nl_reg-1 - i)
                            i4 = idx(fc_bot, kxl-1, kxl-1-kx,
                                     nl_reg-1, nl_reg-1-(i+1))
                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i2

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i2
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i4
                            triangle_idx += 1
                    elif ky == 0 and fc < 5:
                        fc_bot = rhombi_neighbours[fc][3]
                        ky_bot = kxl - 1
                        print(kx, ky, ky_bot)
                        for i in range(nl_reg-1):
                            i1 = idx(fc, kx, ky,  i, 0)
                            i2 = idx(fc, kx, ky,  i+1, 0)
                            i3 = idx(fc_bot, kx, ky_bot, i, nl_reg-1)
                            i4 = idx(fc_bot, kx, ky_bot, i+1, nl_reg-1)
                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i2

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i2
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i4
                            triangle_idx += 1
                    else:
                        for i in range(nl_reg-1):
                            i1 = idx(fc, kx, ky,  i, 0)
                            i2 = idx(fc, kx, ky,  i+1, 0)
                            i3 = idx(fc, kx, ky-1, i, nl_reg-1)
                            i4 = idx(fc, kx, ky-1, i+1, nl_reg-1)
                            triangles[triangle_idx][0] = i1
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i2

                            triangle_idx += 1
                            triangles[triangle_idx][0] = i2
                            triangles[triangle_idx][1] = i3
                            triangles[triangle_idx][2] = i4
                            triangle_idx += 1

        print("number of created triangles:", triangle_idx)

        # R = 6371000.0
        # for i in range(vertices.shape[0]):
        #     vertices[i][0] = self.xyz[i*3 + 0]/R
        #     vertices[i][1] = self.xyz[i*3 + 1]/R
        #     vertices[i][2] = self.xyz[i*3 + 2]/R
        #     elements[i] = i

        #     vertices_colors[i][0] = self.colors[i][0]
        #     vertices_colors[i][1] = self.colors[i][1]
        #     vertices_colors[i][2] = self.colors[i][2]

        r = 1.0
        for i in range(vertices.shape[0]):
            vertices[i, :] = spherical(r,
                                       self.lonlat[2*i+0],
                                       self.lonlat[2*i+1])

            vertices_colors[i][0] = self.colors[i][0]
            vertices_colors[i][1] = self.colors[i][1]
            vertices_colors[i][2] = self.colors[i][2]

        # self.grid_elements_count = elements.size
        # self.grid_elements_count = lines.size
        self.grid_elements_count = triangles.size

        self.vao_grid = gl.glGenVertexArrays(1)
        gl.glBindVertexArray(self.vao_grid)
        gl.glEnableVertexAttribArray(0)

        self.grid_vertex_count = vertices.size
        self.grid_vbo = self.create_vbo(vertices)

        # self.grid_normals_vbo = self.create_normals_vbo(vertices)
        self.grid_colors_vbo = self.create_colors_vbo(vertices_colors)
        # elements = np.zeros((len(self.neighbours), 2), dtype=np.uint32)

        #        fig = plt.figure()
        #        ax = fig.add_subplot(1, 1, 1)

        # self.grid_elements_vbo = self.create_elements_vbo(lines)
        self.grid_elements_vbo = self.create_elements_vbo(triangles)

    def paint_grid(self):
        # display grid
        #gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
        print("paint grid", self.vao_grid, self.grid_vertex_count)
        gl.glBindVertexArray(self.vao_grid)

        c = QVector4D(1.0, 1.0, 1.0, 1.0)

        self.shader_manager.set_colour(c)
        # gl.glDrawArrays(gl.GL_POINTS, 0,
        #                self.grid_vertex_count)

        gl.glDrawElements(gl.GL_TRIANGLES,
                          self.grid_elements_count,
                          gl.GL_UNSIGNED_INT,
                          None)
        #gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)

    def create_elements_vbo(self, elements):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, vbo)

        gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER,
                        elements.nbytes,
                        elements,
                        gl.GL_STATIC_DRAW)

        # gl.glVertexAttribPointer(0, 1,
        #                         gl.GL_INT, False,
        #                         0, None)

        # gl.glEnableVertexAttribArray(0)

        return vbo

    def create_vbo(self, vertices):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        gl.glBufferData(gl.GL_ARRAY_BUFFER,
                        vertices.nbytes,
                        vertices,
                        gl.GL_STATIC_DRAW)

        gl.glVertexAttribPointer(0, 3,
                                 gl.GL_FLOAT, False,
                                 0, None)

        gl.glEnableVertexAttribArray(0)

        return vbo

    def create_normals_vbo(self, vertices):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        gl.glBufferData(gl.GL_ARRAY_BUFFER,
                        vertices.nbytes,
                        vertices,
                        gl.GL_STATIC_DRAW)

        gl.glVertexAttribPointer(2, 3,
                                 gl.GL_FLOAT, False,
                                 0, None)

        gl.glEnableVertexAttribArray(2)

        return vbo

    def create_colors_vbo(self, colors):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        gl.glBufferData(gl.GL_ARRAY_BUFFER,
                        colors.nbytes,
                        colors,
                        gl.GL_STATIC_DRAW)

        gl.glVertexAttribPointer(1, 3,
                                 gl.GL_FLOAT, False,
                                 0, None)
        gl.glEnableVertexAttribArray(1)

        return vbo


if __name__ == '__main__':
    ico = ico(4)
