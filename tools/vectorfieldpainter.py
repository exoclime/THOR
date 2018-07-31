
from PyQt5.QtGui import (QOpenGLShader,
                         QOpenGLShaderProgram,
                         QMatrix4x4,
                         QVector4D,
                         QPainterPath,
                         QFont,
                         QSurfaceFormat)

import OpenGL.GL as gl
import numpy as np
import ctypes
import math
from math import cos, sin, pi, pow, sqrt

import matplotlib.pyplot as plt

from utilities import HSV_to_RGB, spherical

from basepainter import BasePainter


vertex_shader_field = """
#version 400\n
layout(location = 0) in vec3 vp;

uniform highp mat4 view;
uniform highp mat4 projection;
uniform highp mat4 model;
void main() {
   gl_Position = projection * view * model *  vec4(vp, 1.0);
}"""

fragment_shader_field = """
#version 400\n

out vec4 colour_out;
uniform highp vec4 colour;
void main() {
   colour_out = colour;
}"""


class VectorFieldPainter(BasePainter):
    def __init__(self):
        BasePainter.__init__(self)
        self.level = 0
        self.draw_idx = 0
        self.altitude = 0
        self.display_data = False

    def set_grid_data(self, dataset):
        self.dataset = dataset

    def set_shader_manager(self, shader_manager):
        BasePainter.set_shader_manager(self, shader_manager)

        self.shader_manager.add_shaders("icos_field",
                                        vertex=vertex_shader_field,
                                        fragment=fragment_shader_field,
                                        uniforms=["model", "view", "projection", "colour"])

    def initializeGL(self):

        # create

        # num_vertices = 12

        # vertices = np.zeros((num_vertices, 3),
        #                     dtype=np.float32)
        lonlat = self.dataset.get_grid()
        num_points, num_levels = self.dataset.get_dim()

        vertices = np.zeros((num_points, 3),
                            dtype=np.float32)
        print("number of points: ", num_points)

        r = 1.0
        for i in range(vertices.shape[0]):
            vertices[i, :] = spherical(r,
                                       lonlat[i, 0],
                                       lonlat[i, 1])

        # build VAOs for all input data
        self.vao_list = []

        self.num_levels = num_levels
        self.num_points = num_points
        self.last_loaded = -1

        self.display_vector_field = True
        vector_data = np.zeros(
            (self.num_levels, self.num_points, 2, 3), dtype=np.float32)
        vector_elements = np.zeros((num_levels, self.num_points, 2),
                                   dtype=np.uint32)

        for l in range(num_levels):
            for i in range(self.num_points):
                vector_elements[l, i, 0] = 2*l*self.num_points + 2*i + 0
                vector_elements[l, i, 1] = 2*l*self.num_points + 2*i + 1

        self.vao_vector_field = self.create_vector_field_vao(
            vector_data, vector_elements)
        self.field_element_count = num_points*2
        self.last_field_loaded = -1

    def update_field(self):
        if self.draw_idx == self.last_field_loaded:
            return

        data = np.array(self.dataset.get_field_data(self.draw_idx), copy=True)

        self.update_field_vbo(data)
        self.last_field_loaded = self.draw_idx

    def create_vector_field_vao(self, vertices, elements):
        vao = gl.glGenVertexArrays(1)

        gl.glBindVertexArray(vao)
        gl.glEnableVertexAttribArray(0)

        self.vector_field_vbo = self.create_vbo(vertices, dynamic=True)

        self.create_elements_vbo(elements)

        return vao

    def paint_vector_field(self):
        if not self.display_vector_field and not self.display_data:
            return

        self.update_field()
        gl.glBindVertexArray(self.vao_vector_field)
        idx = 4*self.field_element_count * self.altitude
        gl.glDrawElements(gl.GL_LINES,
                          int(self.field_element_count),
                          gl.GL_UNSIGNED_INT,
                          ctypes.c_void_p(idx))

    def update_field_vbo(self, field):
        print("update field")
        self.update_vbo(self.vao_vector_field,
                        self.vector_field_vbo,
                        field)
