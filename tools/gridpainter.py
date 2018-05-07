
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


class GridPainter:
    def __init__(self):
        pass

    def set_grid(self, grid, neighbours, vertices_colors_arrays):
        self.grid = grid
        self.neighbours = neighbours
        self.vertices_colors_arrays = vertices_colors_arrays
        # (image, vertex, color) 3 dim

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def initializeGL(self):

        # create

        vertices = np.zeros((len(self.grid)//2, 3),
                            dtype=np.float32)

        r = 1.0
        #        for i in range(vertices.shape[0]):
        #            vertices[i, 0] = i/vertices.shape[0]
        #            vertices[i, 1] = i/vertices.shape[0]
        #            vertices[i, 2] = 0.0

        for j in range(len(self.grid)//2):
            theta = self.grid[2*j] + math.pi/2.0     # -pi:pi -> 0:2*pi
            # -pi/2: pi/2 -> 0:pi
            phi = self.grid[2*j+1] + math.pi
            vertices[j, 0] = r*math.cos(phi)*math.cos(theta)
            vertices[j, 1] = r*math.cos(phi)*math.sin(theta)
            vertices[j, 2] = r*math.sin(phi)

        print(vertices)

        elements = np.zeros((6*vertices.shape[0], 3), dtype=np.uint32)

        for i in range(vertices.shape[0]):
            elements[6*i+0, :] = np.array((i,
                                           self.neighbours[i*6 + 0],
                                           self.neighbours[i*6 + 1]), dtype=np.uint32)

            elements[6*i+1, :] = np.array((i,
                                           self.neighbours[i*6 + 1],
                                           self.neighbours[i*6 + 2]), dtype=np.uint32)
            elements[6*i+2, :] = np.array((i,
                                           self.neighbours[i*6 + 2],
                                           self.neighbours[i*6 + 3]), dtype=np.uint32)
            elements[6*i+3, :] = np.array((i,
                                           self.neighbours[i*6 + 3],
                                           self.neighbours[i*6 + 4]), dtype=np.uint32)
            elements[6*i+4, :] = np.array((i,
                                           self.neighbours[i*6 + 4],
                                           self.neighbours[i*6 + 5]), dtype=np.uint32)
            elements[6*i+5, :] = np.array((i,
                                           self.neighbours[i*6 + 5],
                                           self.neighbours[i*6 + 0]), dtype=np.uint32)

        print("number of triangles", elements.shape[0])
        sorted_elements = np.unique(np.sort(elements, axis=1), axis=0)
        elements = sorted_elements
        print("number of triangles after sort", elements.shape[0])
        self.grid_elements_count = elements.size

        self.vao_grid = []

        for c in range(self.vertices_colors_arrays.shape[0]):
            gl.glEnableVertexAttribArray(0)
            vao_grid = gl.glGenVertexArrays(1)
            gl.glBindVertexArray(vao_grid)

            self.grid_vertex_count = vertices.size
            self.grid_vbo = self.create_vbo(vertices)

            self.grid_normals_vbo = self.create_normals_vbo(vertices)
            self.grid_colors_vbo = self.create_colors_vbo(
                self.vertices_colors_arrays[c])
            # elements = np.zeros((len(self.neighbours), 2), dtype=np.uint32)

            #        fig = plt.figure()
            #        ax = fig.add_subplot(1, 1, 1)

            self.grid_elements_vbo = self.create_elements_vbo(elements)

            self.vao_grid.append(vao_grid)

    def paint_grid(self, idx):
        # display grid
        print("paint grid", self.vao_grid[idx], self.grid_vertex_count)
        gl.glBindVertexArray(self.vao_grid[idx])

        c = QVector4D(1.0, 1.0, 1.0, 1.0)

        # self.shader_manager.set_colour(c)
        # gl.glDrawArrays(gl.GL_POINTS, 0,
        #                self.grid_vertex_count)

        gl.glDrawElements(gl.GL_TRIANGLES,
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
