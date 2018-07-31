
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
from basepainter import BasePainter


vertex_shader_axis = """
#version 400\n
layout(location = 0) in vec3 vp;
layout(location = 1) in vec3 vertex_colour;
out vec3 colour_;

uniform highp mat4 view;
uniform highp mat4 projection;
uniform highp mat4 model;
void main() {
   gl_Position = projection * view * model *  vec4(vp, 1.0);
   colour_ = vertex_colour;
}"""

fragment_shader_axis = """
#version 400\n
in vec3 colour_;
out vec4 colour_out;

void main() {
   colour_out = vec4(colour_, 1.0);
}"""


class AxisPainter(BasePainter):
    def __init__(self):
        BasePainter.__init__(self)

    def set_shader_manager(self, shader_manager):
        BasePainter.set_shader_manager(self, shader_manager)
        self.shader_manager.add_shaders("axis",
                                        vertex=vertex_shader_axis,
                                        fragment=fragment_shader_axis,
                                        uniforms=["model", "view", "projection", "colour"])

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

        gl.glDrawArrays(gl.GL_LINES, 0,
                        self.grid_vertex_count)
