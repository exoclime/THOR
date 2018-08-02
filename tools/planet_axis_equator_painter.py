
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

from utilities import HSV_to_RGB, spherical, cylindrical

from basepainter import BasePainter

vertex_shader = """#version 400\n
                       layout(location = 0) in vec3 vp;
                       layout(location = 1) in vec3 vertex_colour;


                       out VS_OUT {
                           out vec3 colour;
                       } vs_out;


                       uniform highp mat4 view;
                       uniform highp mat4 projection;
                       uniform highp mat4 model;
                       void main() {
                         vs_out.colour = vertex_colour;
                          //colour = vec3(1.0,1.0,0.0);
                          gl_Position = projection * \
                              view * model *  vec4(vp, 1.0);
                       }"""

geometry_shader = """
#version 400\n


layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in VS_OUT {
   vec3 colour;
} gs_in[];

out GS_OUT {
   vec3 fcolor;
   vec3 dist;
} gs_out;

void main(void)
{
   float WIN_SCALE = 800.0;
   // taken from 'Single-Pass Wireframe Rendering'
   vec2 p0 = WIN_SCALE * gl_in[0].gl_Position.xy/gl_in[0].gl_Position.w;
   vec2 p1 = WIN_SCALE * gl_in[1].gl_Position.xy/gl_in[1].gl_Position.w;
   vec2 p2 = WIN_SCALE * gl_in[2].gl_Position.xy/gl_in[2].gl_Position.w;
                           vec2 v0 = p2-p1;
   vec2 v1 = p2-p0;
   vec2 v2 = p1-p0;
   float area = abs(v1.x*v2.y - v1.y * v2.x);

   gs_out.fcolor = gs_in[0].colour;

   gs_out.dist = vec3(area/length(v0),0,0);
   gl_Position = 1.5*gl_in[0].gl_Position;
   EmitVertex();
   gs_out.fcolor = gs_in[1].colour;
   gs_out.dist = vec3(0,area/length(v1),0);
   gl_Position = gl_in[1].gl_Position;
   EmitVertex();
   gs_out.fcolor = gs_in[2].colour;
   gs_out.dist = vec3(0,0,area/length(v2));
   gl_Position = gl_in[2].gl_Position;
   EmitVertex();
   EndPrimitive();
}"""

fragment_shader = """
#version 400\n
in GS_OUT {
   vec3 fcolor;
   vec3 dist;
} fs_in;

out vec4 frag_colour;
uniform float wire_limit;
void main() {
   float nearD = min(min(fs_in.dist[0],fs_in.dist[1]),fs_in.dist[2]);
   float edgeIntensity = exp2(-1.0*nearD*nearD);
   if (nearD < wire_limit)
      frag_colour = vec4(1.0, 0.1, 0.1, 0.3 );
   else
      frag_colour = vec4(fs_in.fcolor, 0.3);
      //frag_colour = vec4(fs_in.fcolor, 1.0);
      //frag_colour = vec4(1.0,0.0,0.0, 1.0);
      //                                   frag_colour = (edgeIntensity * vec4( 0.1, 0.1, 0.1, 1.0 )) +
      //                                                 (1.0-edgeIntensity)*vec4(colour, 1.0);
}"""


class PlanetAxisEquatorPainter(BasePainter):
    def __init__(self):
        BasePainter.__init__(self)
        self.equator_radius_inner = 0.9
        self.equator_radius_outer = 1.2
        self.pole_length = 1.3
        self.pole_radius = 0.02
        self.wireframe = False

    def set_shader_manager(self, shader_manager):
        BasePainter.set_shader_manager(self, shader_manager)

        self.shader_manager.add_shaders("axis_equator",
                                        vertex=vertex_shader,
                                        geometry=geometry_shader,
                                        fragment=fragment_shader,
                                        uniforms=["model", "view", "projection", "wire_limit"])

    def initializeGL(self):

        pole_num_sides = 60

        pole_vertices = np.zeros((pole_num_sides, 4,  3),
                                 dtype=np.float32)
        pole_elements = np.zeros((pole_num_sides, 2, 3), dtype=np.uint32)

        pole_colors = np.ones((pole_num_sides, 4,  3),
                              dtype=np.float32)

        for i in range(pole_num_sides):
            pole_vertices[i, 0, :] = cylindrical(self.pole_radius,
                                                 self.pole_length,
                                                 i*2.0*math.pi/pole_num_sides)
            pole_vertices[i, 1, :] = cylindrical(self.pole_radius,
                                                 self.pole_length,
                                                 (i+1)*2.0*math.pi/pole_num_sides)

            pole_vertices[i, 2, :] = cylindrical(self.pole_radius,
                                                 -self.pole_length,
                                                 (i+1)*2.0*math.pi/pole_num_sides)

            pole_vertices[i, 3, :] = cylindrical(self.pole_radius,
                                                 -self.pole_length,
                                                 i*2.0*math.pi/pole_num_sides)

            pole_elements[i, 0, 0] = 4*i + 0
            pole_elements[i, 0, 1] = 4*i + 1
            pole_elements[i, 0, 2] = 4*i + 2
            pole_elements[i, 1, 0] = 4*i + 0
            pole_elements[i, 1, 1] = 4*i + 2
            pole_elements[i, 1, 2] = 4*i + 3

        self.pole_elements_count = pole_elements.size

        self.vao_pole = self.create_vao(pole_vertices,
                                        pole_elements,
                                        pole_colors)

        equator_num_sides = 60

        equator_vertices = np.zeros((equator_num_sides, 4,  3),
                                    dtype=np.float32)
        equator_elements = np.zeros((equator_num_sides, 2, 3), dtype=np.uint32)

        equator_colors = np.ones((equator_num_sides, 4,  3),
                                 dtype=np.float32)

        for i in range(equator_num_sides):
            equator_vertices[i, 0, :] = cylindrical(self.equator_radius_outer,
                                                    0.0,
                                                    i*2.0*math.pi/equator_num_sides)
            equator_vertices[i, 1, :] = cylindrical(self.equator_radius_outer,
                                                    0.0,
                                                    (i+1)*2.0*math.pi/equator_num_sides)

            equator_vertices[i, 2, :] = cylindrical(self.equator_radius_inner,
                                                    0.0,
                                                    (i+1)*2.0*math.pi/equator_num_sides)

            equator_vertices[i, 3, :] = cylindrical(self.equator_radius_inner,
                                                    0.0,
                                                    i*2.0*math.pi/equator_num_sides)

            equator_elements[i, 0, 0] = 4*i + 0
            equator_elements[i, 0, 1] = 4*i + 1
            equator_elements[i, 0, 2] = 4*i + 2
            equator_elements[i, 1, 0] = 4*i + 0
            equator_elements[i, 1, 1] = 4*i + 2
            equator_elements[i, 1, 2] = 4*i + 3

        self.equator_elements_count = equator_elements.size

        self.vao_equator = self.create_vao(equator_vertices,
                                           equator_elements,
                                           equator_colors)

    def create_vao(self, vertices, triangles, colors):
        vao = gl.glGenVertexArrays(1)

        gl.glBindVertexArray(vao)
        gl.glEnableVertexAttribArray(0)

        self.create_vbo(vertices)

        self.create_colors_vbo(colors)

        self.create_elements_vbo(triangles)

        return vao

    def paint(self):
        # display grid
        if self.wireframe:
            self.shader_manager.set_uniform(
                "axis_equator", "wire_limit", 0.00001)
        else:
            self.shader_manager.set_uniform(
                "axis_equator", "wire_limit", -10.0)

        gl.glBindVertexArray(self.vao_pole)

        gl.glDrawElements(gl.GL_TRIANGLES,
                          int(self.pole_elements_count),
                          gl.GL_UNSIGNED_INT,
                          None)

        gl.glBindVertexArray(self.vao_equator)

        gl.glDrawElements(gl.GL_TRIANGLES,
                          int(self.equator_elements_count),
                          gl.GL_UNSIGNED_INT,
                          None)
