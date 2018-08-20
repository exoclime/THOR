
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
from icogrid import ico
import time

vertex_shader_marker = """
# version 400\n
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

fragment_shader_marker = """
# version 400\n
in vec3 colour_;
out vec4 colour_out;

void main() {
   colour_out = vec4(colour_, 1.0);
}"""


vertex_shader_point = """
# version 400\n
layout(location = 0) in vec3 vp;
//layout(location = 1) in vec3 vertex_colour;


out GS_OUT {
   vec3 fcolor;
} gs_out;

uniform highp mat4 view;
uniform highp mat4 projection;
uniform highp mat4 model;

void main() {
   gl_Position = projection * view * model *  vec4(vp, 1.0);
   //gl_Position = vec4(vp, 1.0);
   gs_out.fcolor = vec3(0.0,1.0,0.0);
}"""


fragment_shader_point = """
# version 400\n

in GS_OUT {
   vec3 fcolor;
} fs_in;


out vec4 colour_out;

void main() {
   colour_out = vec4(fs_in.fcolor, 1.0);
}"""


geometry_shader_point = """
#version 400\n


layout (points) in;
layout (triangle_strip, max_vertices = 3) out;
//layout (points, max_vertices = 3) out;

in VS_OUT {
   vec3 colour;
} gs_in[];

out GS_OUT {
   vec3 fcolor;
} gs_out;


uniform highp mat4 view;
uniform highp mat4 projection;
uniform highp mat4 model;


void main(void)
{
float a = 0.05;
float w = gl_in[0].gl_Position.w;
vec4 p = gl_in[0].gl_Position;

gl_Position = projection * view * model * (p + vec4(a, 0.0,a,0.0));
gs_out.fcolor = gs_in[0].colour;

EmitVertex();
gl_Position = projection * view * model * (p + vec4( -a,-a,0.0,0.0));
gs_out.fcolor = gs_in[0].colour;

EmitVertex();

gl_Position = projection * view * model * (p + vec4( -a,a,0.0,0.0));
gs_out.fcolor = gs_in[0].colour;

EmitVertex();
EndPrimitive();
}
"""


class MarkerGridPainter(BasePainter):
    def __init__(self):
        BasePainter.__init__(self)
        self.level = 0
        self.wireframe = False
        self.draw_idx = 0
        self.altitude = 0
        self.show_data = True

    def set_grid_data(self, dataset):
        self.dataset = dataset

    def set_shader_manager(self, shader_manager):
        BasePainter.set_shader_manager(self, shader_manager)

        self.shader_manager.add_shaders("marker_grid",
                                        vertex=vertex_shader_marker,
                                        fragment=fragment_shader_marker,
                                        uniforms=["model", "view", "projection"])

        self.shader_manager.add_shaders("uv_point_shader",
                                        vertex=vertex_shader_point,
                                        # geometry=geometry_shader_point,
                                        fragment=fragment_shader_point,
                                        uniforms=["model", "view", "projection"])

    def initializeGL(self):

        # create

        # num_vertices = 12

        # vertices = np.zeros((num_vertices, 3),
        #                     dtype=np.float32)
        lonlat = self.dataset.get_grid()
        num_points, num_levels = self.dataset.get_dim()

        g = int(pow((num_points - 2)/10, 1/4)) - 2

        print("Start computation for barycentric coordinates")
        start = time.time()
        self.ico = ico(g, lonlat)

        # get some triangles
        l = 4
        rs = np.zeros((l, 3))

        theta = 25.0*math.pi/(2.0*90.0)
        phi = 10.0*math.pi/(2.0*90.0)
        r = spherical(1.0, theta, phi)
        p_z = 0.0
        p_y = math.sqrt(1.0/(1+math.pow(r[1]/r[0], 2.0)))
        p_x = -r[1]/r[0]*p_y
        p1 = np.array((p_x, p_y, p_z))
        p2 = np.cross(p1, r)

        for i in range(l):
            alpha = (i/(l))*math.pi*2.0
            rs[i, :] = cos(alpha)*p1+sin(alpha)*p2

#            rs[i, :] = spherical(1.0, 25.0*math.pi/(2.0*90.0),
#                                 ((i/(l))*2.0 - 1)*math.pi)

        dot, t_i, u_i, v_i = self.ico.barycentric_coordinates.get_barycentric_coordinates(
            rs)

        print("dot", dot.shape)
        print("t_i", t_i.shape)
        print("u_i", u_i.shape)
        print("v_i", v_i.shape)
        # for i in range(num_points):
        #    print(dot[i], t_i[i]/dot[i], u_i[i]/dot[i], v_i[i]/dot[i])

        # get indexes that actually match:

        # t = t_i/dot
        u = u_i/dot
        v = v_i/dot
        print(u.shape, v.shape, u_i.shape, v_i.shape, dot.shape)
        condition = np.logical_and.reduce((dot > 0.0,
                                           u > 0.0,
                                           u < 1.0,
                                           v > 0.0,
                                           v < 1.0,
                                           u+v < 1.0))
        w = np.where(condition)
        stop = time.time()
        print("Barycentric coordinates done: {} s".format(stop-start))

        a = 0.01
        rr = 1.002
        # build UV coordinates

        print("w0", w[0].shape)
        uv_vertices = np.zeros((3*w[0].shape[0], 3), dtype=np.float32)
        print("uv", u.shape, v.shape)
        for i in range(w[0].shape[0]):
            triangle_idx = w[0][i]
            triangle = self.ico.triangles[triangle_idx]
            print("tri", triangle)
            u_p = u[w[0][i], w[1][i]]
            v_p = v[w[0][i], w[1][i]]
            print("u:\t{}\tv:{}".format(u_p, v_p))

            i1 = triangle[0]
            i2 = triangle[1]
            i3 = triangle[2]
            p1 = spherical(rr, lonlat[i1, 1], lonlat[i1, 0])
            print("p1", p1)
            p2 = spherical(rr, lonlat[i2, 1], lonlat[i2, 0])
            p3 = spherical(rr, lonlat[i3, 1], lonlat[i3, 0])
            v0 = p2 - p1
            v1 = p3 - p1

            pk = p1 + u_p*v1 + v_p * v0
            #uv_vertices[i] = p1 + u_p*v0 + v_p * v1
            uv_vertices[3*i+0] = pk + \
                np.array((a, 0.0, 0.0), dtype=np.float32)
            uv_vertices[3*i+1] = pk + np.array((-a, -a, 0.0), dtype=np.float32)
            uv_vertices[3*i+2] = pk + np.array((-a,  a, 0.0), dtype=np.float32)
        self.uv_count = uv_vertices.size

        self.vao_uv = self.create_uv_vao(uv_vertices)

        # for r in range(l):
        #     print(rs[r, :])

        # for i in range(w[0].shape[0]):
        #     idx = w[0][i], w[1][i]
        #     t = t_i[idx[0]]/dot[idx]
        #     print(idx, dot[idx], t, u[idx], v[idx])

        self.triangles = False
        r = 1.002
        if self.triangles:
            triangle_mesh = self.ico.triangles[w[0]]
            elements = np.zeros((triangle_mesh.shape[0], 3), np.int32)
            vertices = np.zeros((triangle_mesh.shape[0] * 3, 3),
                                dtype=np.float32)
            colors = np.zeros((triangle_mesh.shape[0] * 3, 3),
                              dtype=np.float32)
            cnt = 0

            c1 = np.array((1.0, 0.0, 0.0), dtype=np.float32)
            c2 = np.array((0.0, 1.0, 0.0), dtype=np.float32)
            c3 = np.array((0.0, 0.0, 1.0), dtype=np.float32)

            for t in range(triangle_mesh.shape[0]):
                triangle = triangle_mesh[t]
                i1 = triangle[0]
                i2 = triangle[1]
                i3 = triangle[2]

                p1 = spherical(r, lonlat[i1, 1], lonlat[i1, 0])
                p2 = spherical(r, lonlat[i2, 1], lonlat[i2, 0])
                p3 = spherical(r, lonlat[i3, 1], lonlat[i3, 0])

                vertices[cnt, :] = p1
                colors[cnt, :] = c1
                idx1 = cnt
                cnt += 1
                vertices[cnt, :] = p2
                colors[cnt, :] = c2
                idx2 = cnt
                cnt += 1

                vertices[cnt, :] = p3
                colors[cnt, :] = c3
                idx3 = cnt
                cnt += 1

                elements[t, 0] = idx1
                elements[t, 1] = idx2
                elements[t, 2] = idx3

            self.elements_count = elements.size
            self.vao_list = self.create_marker_vao(vertices,
                                                   elements,
                                                   colors)
        else:
            vertices = np.zeros((num_points, 3),
                                dtype=np.float32)
            colors = np.zeros((num_points, 3),
                              dtype=np.float32)

            triangle_mesh = self.ico.triangles[w[0]]

            elements = np.zeros((triangle_mesh.shape[0], 3, 2), np.int32)
            print(triangle_mesh.shape, elements.shape)
            cnt = 0
            for t in range(triangle_mesh.shape[0]):
                triangle = triangle_mesh[t]

                elements[t, 0, 0] = triangle[0]
                elements[t, 0, 1] = triangle[1]
                elements[t, 1, 0] = triangle[1]
                elements[t, 1, 1] = triangle[2]
                elements[t, 2, 0] = triangle[2]
                elements[t, 2, 1] = triangle[0]
            self.elements_count = elements.size
            for i in range(vertices.shape[0]):
                # spherical -> theta vertical mvt, declination-> latitude
                # -> phi -> horizontal mvt, azimuth -> longitude
                vertices[i, :] = spherical(r,
                                           lonlat[i, 1],
                                           lonlat[i, 0])
                colors[i, :] = np.array((1.0, 1.0, 1.0), dtype=np.float32)

            self.vao_list = self.create_marker_vao(vertices,
                                                   elements,
                                                   colors)

            num_radii = rs.shape[0]
            radii = np.zeros((rs.shape[0], 2, 3), dtype=np.float32)
            radii_elements = np.zeros((rs.shape[0], 2), dtype=np.int32)
            radii_colors = np.ones(radii.shape, dtype=np.float32)
            r2 = 1.1

            for i in range(num_radii):
                radii[i, 0, :] = np.zeros(3, dtype=np.float32)
                radii[i, 1, :] = r2*rs[i]
                radii_elements[i, 0] = 2*i
                radii_elements[i, 0] = 2*i+1
            self.radii_elements_count = radii_elements.size
            self.vao_radii = self.create_marker_vao(
                radii, radii_elements, radii_colors)

    # def update_colors(self):
    #     if not self.show_data:

    #         self.update_colors_vbo(self.grid_color_data)
    #         self.update_values_vbo(self.grid_scalar_data)
    #     else:
    #         if self.draw_idx == self.last_loaded:
    #             return
    #         data, data_scalar = self.dataset.get_color_data(self.draw_idx)

    #         data = np.array(data, copy=True)
    #         data_scalar = np.array(data_scalar, copy=True)
    #         self.update_colors_vbo(data)
    #         self.update_values_vbo(data_scalar)
    #         self.last_loaded = self.draw_idx

    def create_uv_vao(self, vertices):
        vao = gl.glGenVertexArrays(1)

        gl.glBindVertexArray(vao)
        gl.glEnableVertexAttribArray(0)

        self.create_vbo(vertices)

        return vao

    def create_marker_vao(self, vertices, elements, colors):
        vao = gl.glGenVertexArrays(1)

        gl.glBindVertexArray(vao)
        gl.glEnableVertexAttribArray(0)

        self.create_vbo(vertices)

        self.color_vbo = self.create_colors_vbo(colors)
        self.create_elements_vbo(elements)

        return vao

    def paint_uv(self):
        gl.glBindVertexArray(self.vao_uv)
        gl.glDrawArrays(gl.GL_TRIANGLES,
                        0,
                        int(self.uv_count))

    def paint_grid(self):
        # display grid
        # self.update_colors()

        gl.glBindVertexArray(self.vao_list)
        if self.triangles:
            gl.glDrawElements(gl.GL_TRIANGLES,
                              int(self.elements_count),
                              gl.GL_UNSIGNED_INT,
                              None)

        else:
            gl.glDrawElements(gl.GL_LINES,
                              int(self.elements_count),
                              gl.GL_UNSIGNED_INT,
                              None)

        gl.glBindVertexArray(self.vao_radii)
        gl.glDrawElements(gl.GL_LINES,
                          int(self.radii_elements_count),
                          gl.GL_UNSIGNED_INT,
                          None)
        # ctypes.c_void_p(idx))

    def update_colors_vbo(self, colors):
        # print("update colors")
        self.update_vbo(self.vao_list,
                        self.color_vbo,
                        colors)
