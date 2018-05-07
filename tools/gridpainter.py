
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


class GridPainter:
    def __init__(self):
        pass

    def set_grid(self, grid, neighbours):
        self.grid = grid
        self.neighbours = neighbours

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def initializeGL(self):

        # create
        gl.glEnableVertexAttribArray(0)
        self.vao_grid = gl.glGenVertexArrays(1)
        gl.glBindVertexArray(self.vao_grid)

        vertices = np.zeros((len(self.grid)//2, 3),
                            dtype=np.float32)

        r = 1.0
#        for i in range(vertices.shape[0]):
#            vertices[i, 0] = i/vertices.shape[0]
#            vertices[i, 1] = i/vertices.shape[0]
#            vertices[i, 2] = 0.0
        for i in range(len(self.grid)//2):
            lat = self.grid[2*i]
            lon = self.grid[2*i+1]
            vertices[i, 0] = r*math.sin(lat)*math.cos(lon)
            vertices[i, 1] = r*math.sin(lat)*math.sin(lon)
            vertices[i, 2] = r*math.cos(lat)

        print(vertices)
        self.grid_vertex_count = vertices.size
        self.grid_vbo = self.create_vbo(vertices)

        elements = np.zeros((len(self.neighbours), 2), dtype=np.uint32)

        for i in range(vertices.shape[0]):
            for j in range(6):
                elements[i*6 + j, 0] = i
                elements[i*6 + j, 1] = self.neighbours[i*6 + j]

        self.grid_elements_count = elements.size

        self.grid_elements_vbo = self.create_elements_vbo(elements)

    def paint_grid(self):
        # display grid
        print("paint grid", self.vao_grid, self.grid_vertex_count)
        gl.glBindVertexArray(self.vao_grid)

        c = QVector4D(1.0, 1.0, 1.0, 1.0)

        self.shader_manager.set_colour(c)
        # gl.glDrawArrays(gl.GL_POINTS, 0,
        #                self.grid_vertex_count)

        gl.glDrawElements(gl.GL_LINES,
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
