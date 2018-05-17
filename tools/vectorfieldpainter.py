
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

import matplotlib.pyplot as plt


class FieldPainter:
    def __init__(self):
        pass

    def set_field(self, grid, moments):
        self.grid = grid
        self.moments = moments
        # (image, vertex, color) 3 dim

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def initializeGL(self):

        # create
        num_images = self.moments.shape[0]
        print("num fields", num_images)
        vertices = np.zeros((num_images, len(self.grid)//2, 2, 3),
                            dtype=np.float32)
        print("vertices files", vertices.shape)
        r = 1.2
        #        for i in range(vertices.shape[0]):
        #            vertices[i, 0] = i/vertices.shape[0]
        #            vertices[i, 1] = i/vertices.shape[0]
        #            vertices[i, 2] = 0.0

        # k = 1.0/(np.max(self.moments)*1.0)
        k = 0.0
        for m in range(num_images):
            for j in range(len(self.grid)//2):
                theta = self.grid[2*j] + math.pi/2.0     # -pi:pi -> 0:2*pi
                # -pi/2: pi/2 -> 0:pi
                phi = self.grid[2*j+1] + math.pi
                r_0 = r*math.cos(phi)*math.cos(theta)
                r_1 = r*math.cos(phi)*math.sin(theta)
                r_2 = r*math.sin(phi)

                vertices[m, j, 0, 0] = 0
                vertices[m, j, 0, 1] = 0
                vertices[m, j, 0, 2] = 0

                vertices[m, j, 1, 0] = r_0 + k*r_0
                vertices[m, j, 1, 1] = r_1 + k*r_1
                vertices[m, j, 1, 2] = r_2 + k*r_2

                # vertices[m, j, 1, 0] = r_0 + k*self.moments[m, j, 0]
                # vertices[m, j, 1, 1] = r_1 + k*self.moments[m, j, 1]
                # vertices[m, j, 1, 2] = r_2 + k*self.moments[m, j, 2]

        self.vao_vector_field = []
        self.field_vbos = []
        for c in range(num_images):
            gl.glEnableVertexAttribArray(0)
            vao_field = gl.glGenVertexArrays(1)
            gl.glBindVertexArray(vao_field)

            self.field_vertex_count = vertices[m, :, :, :].size
            vbo = self.create_vbo(vertices[m, :, :, :])
            self.field_vbos.append(vbo)
            print("create field", vao_field)
            self.vao_vector_field.append(vao_field)

    def paint_field(self, idx):
        # display grid
        print("paint field", idx, self.field_vertex_count,
              self.vao_vector_field[idx])

        gl.glBindVertexArray(self.vao_vector_field[idx])

        c = QVector4D(1.0, 1.0, 1.0, 1.0)

        self.shader_manager.set_colour_field(c)
        gl.glDrawArrays(gl.GL_LINES, 0,
                        self.field_vertex_count)

        # gl.glDrawElements(gl.GL_TRIANGLES,
        #                  self.grid_elements_count,
        #                  gl.GL_UNSIGNED_INT,
        #                  None)

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
