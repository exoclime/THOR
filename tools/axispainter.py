
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


class AxisPainter:
    def __init__(self):
        pass

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def initializeGL(self):

        # create
        gl.glEnableVertexAttribArray(0)
        self.vao_axis = gl.glGenVertexArrays(1)
        gl.glBindVertexArray(self.vao_axis)

        vertices = np.zeros((6, 3),
                            dtype=np.float32)

        colors = np.zeros((6, 3),
                          dtype=np.float32)

        vertices[0, :] = np.array((0.0, 0.0, 0.0), dtype=np.float32)
        vertices[1, :] = np.array((1.0, 0.0, 0.0), dtype=np.float32)
        vertices[2, :] = np.array((0.0, 0.0, 0.0), dtype=np.float32)
        vertices[3, :] = np.array((0.0, 1.0, 0.0), dtype=np.float32)
        vertices[4, :] = np.array((0.0, 0.0, 0.0), dtype=np.float32)
        vertices[5, :] = np.array((0.0, 0.0, 1.0), dtype=np.float32)

        colors[0, :] = np.array((1.0, 0.0, 0.0), dtype=np.float32)
        colors[1, :] = np.array((1.0, 0.0, 0.0), dtype=np.float32)
        colors[2, :] = np.array((0.0, 1.0, 0.0), dtype=np.float32)
        colors[3, :] = np.array((0.0, 1.0, 0.0), dtype=np.float32)
        colors[4, :] = np.array((0.0, 0.0, 1.0), dtype=np.float32)
        colors[5, :] = np.array((0.0, 0.0, 1.0), dtype=np.float32)

        self.grid_vertex_count = vertices.size
        self.grid_vbo = self.create_vbo(vertices)

        self.grid_colors_count = colors.size

        self.grid_colors_vbo = self.create_colors_vbo(colors)

    def paint_axis(self):
        # display grid
        print("paint axis", self.vao_axis, self.grid_vertex_count)
        gl.glBindVertexArray(self.vao_axis)

        #c = QVector4D(1.0, 1.0, 1.0, 1.0)

        # self.shader_manager.set_colour(c)
        gl.glDrawArrays(gl.GL_LINES, 0,
                        self.grid_vertex_count)

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
