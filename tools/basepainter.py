
import OpenGL.GL as gl
import numpy as np
import ctypes
import math


class BasePainter:
    def __init__(self):
        self.shader_manager = None

    def set_shader_manager(self, shader_manager):
        self.shader_manager = shader_manager

    def create_elements_vbo(self, elements, dynamic=False):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, vbo)
        if dynamic:
            gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER,
                            elements.nbytes,
                            elements,
                            gl.GL_DYNAMIC_DRAW)
        else:
            gl.glBufferData(gl.GL_ELEMENT_ARRAY_BUFFER,
                            elements.nbytes,
                            elements,
                            gl.GL_STATIC_DRAW)
        return vbo

    def create_vbo(self, vertices, dynamic=False):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        if dynamic:
            gl.glBufferData(gl.GL_ARRAY_BUFFER,
                            vertices.nbytes,
                            vertices,
                            gl.GL_DYNAMIC_DRAW)
        else:
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

    def create_colors_vbo(self, colors, dynamic=False, num_colors=3):
        vbo = gl.glGenBuffers(1)

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        if dynamic:
            gl.glBufferData(gl.GL_ARRAY_BUFFER,
                            colors.nbytes,
                            colors,
                            gl.GL_DYNAMIC_DRAW)
        else:
            gl.glBufferData(gl.GL_ARRAY_BUFFER,
                            colors.nbytes,
                            colors,
                            gl.GL_STATIC_DRAW)

        gl.glVertexAttribPointer(1, num_colors,
                                 gl.GL_FLOAT, False,
                                 0, None)
        gl.glEnableVertexAttribArray(1)

        return vbo

    def update_vbo(self, vao, vbo, data, offset=0):
        gl.glBindVertexArray(vao)
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, vbo)

        gl.glBufferSubData(gl.GL_ARRAY_BUFFER,
                           offset,
                           data.nbytes,
                           data)
