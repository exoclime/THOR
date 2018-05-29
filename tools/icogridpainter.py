
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
from math import cos, sin, pi

import matplotlib.pyplot as plt


def spherical(r, theta, phi):
    return np.array((
        r*math.cos(theta)*math.cos(phi),
        r*math.sin(theta)*math.cos(phi),
        r*math.sin(phi)),
        dtype=np.float32)


class IcoGridPainter:
    def __init__(self):
        self.level = 0

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def initializeGL(self):

        # create

        num_vertices = 12

        vertices = np.zeros((num_vertices, 3),
                            dtype=np.float32)

        r = 1.2
        w = 2.0*math.acos(1.0/(2.0*sin(pi/5.0)))

        angles = [
            [0.0, pi/2.0],
            [0.0, -pi/2.0],
            [-pi/5.0, pi/2.0-w],
            [pi/5.0, pi/2.0-w],
            [3.0*pi/5.0, pi/2.0-w],
            [pi, pi/2.0-w],
            [-3.0*pi/5.0, pi/2.0-w],
            [0.0, -(pi/2.0-w)],
            [2.0*pi/5.0, -(pi/2.0-w)],
            [4.0*pi/5.0, -(pi/2.0-w)],
            [-4.0*pi/5.0, -(pi/2.0-w)],
            [-2.0*pi/5.0, -(pi/2.0-w)]
        ]

        for i, a in enumerate(angles):
            vertices[i, :] = spherical(r, a[0], a[1])

        elements = np.zeros(vertices.shape[0], dtype=np.uint32)

        for i in range(vertices.shape[0]):
            elements[i] = i

        self.grid_elements_count = elements.size

        self.vao_grid = []

        gl.glEnableVertexAttribArray(0)
        self.vao_grid = gl.glGenVertexArrays(1)
        gl.glBindVertexArray(self.vao_grid)

        self.grid_vertex_count = vertices.size
        self.grid_vbo = self.create_vbo(vertices)

        #self.grid_normals_vbo = self.create_normals_vbo(vertices)
        # self.grid_colors_vbo = self.create_colors_vbo(
        #    self.vertices_colors_arrays[c])
        # elements = np.zeros((len(self.neighbours), 2), dtype=np.uint32)

        #        fig = plt.figure()
        #        ax = fig.add_subplot(1, 1, 1)

        self.grid_elements_vbo = self.create_elements_vbo(elements)

    def paint_grid(self):
        # display grid
        print("paint grid", self.vao_grid, self.grid_vertex_count)
        gl.glBindVertexArray(self.vao_grid)

        c = QVector4D(1.0, 1.0, 1.0, 1.0)

        # self.shader_manager.set_colour(c)
        # gl.glDrawArrays(gl.GL_POINTS, 0,
        #                self.grid_vertex_count)

        gl.glDrawElements(gl.GL_POINTS,
                          self.grid_elements_count,
                          gl.GL_UNSIGNED_INT,
                          None)

    def create_elements_vbo(self, elements):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, vbo)

        gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER,
                        elements.nbytes,
                        elements,
                        gl.GL_STATIC_DRAW)

        # gl.glVertexAttribPointer(0, 3,
        #                          gl.GL_FLOAT, False,
        #                          0, None)

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
