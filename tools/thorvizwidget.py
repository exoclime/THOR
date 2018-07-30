"""
Visualisation widget for THOR datasets
"""

from PyQt5.QtWidgets import (QApplication,
                             QGridLayout,
                             QLabel,
                             QOpenGLWidget,
                             QWidget)

from PyQt5.QtGui import (QOpenGLShader,
                         QOpenGLShaderProgram,
                         QMatrix4x4,
                         QVector4D,
                         QPainterPath,
                         QFont,
                         QSurfaceFormat)

from PyQt5.QtCore import (QPointF)

from PyQt5 import QtCore
from PyQt5 import QtGui
import OpenGL.GL as gl
import numpy as np
import math
import time

from shadermanager import ShaderManager
from gridpainter import GridPainter
from axispainter import AxisPainter
from vectorfieldpainter import FieldPainter
from icogridpainter import IcoGridPainter


class ThorVizWidget(QOpenGLWidget):

    def __init__(self, parent=None):
        super(ThorVizWidget, self).__init__(parent)
        print("initialise ThorVizWidget")

        fmt = QSurfaceFormat()
        fmt.setVersion(4, 6)
        fmt.setSamples(4)
        fmt.setProfile(QSurfaceFormat.NoProfile)
        fmt.setSwapBehavior(QSurfaceFormat.TripleBuffer)
        self.setFormat(fmt)

        #self.grid_painter = GridPainter()
        self.ico_grid_painter = IcoGridPainter()
        self.axis_painter = AxisPainter()
        #self.field_painter = FieldPainter()

        self.m_matrixUniform = 0

        self.cam_theta = 0.0
        self.cam_phi = 0.0

        self.cam_r = -5.0
        self.cam_orientation = QMatrix4x4()

        self.right_button_pressed = False
        self.left_button_pressed = False
        self.image = 0

    def set_grid_data(self, *args):
        self.ico_grid_painter.set_grid_data(*args)

    def set_grid(self, grid, neighbours, vertices_color_arrays, moments):
        self.grid_painter.set_grid(
            grid, neighbours, vertices_color_arrays)
        self.field_painter.set_field(grid, moments)

    def set_level(self, level):
        self.ico_grid_painter.altitude = level
        self.update()

    def set_image(self, img):
        self.image = img
        self.ico_grid_painter.draw_idx = self.image
        self.update()

    def show_data(self, b):
        if b:
            self.ico_grid_painter.draw_idx = self.image + 1

            self.update()

    def show_grid(self, b):
        if b:
            self.ico_grid_painter.draw_idx = 0
            self.update()

    def show_wireframe(self, b):
        self.ico_grid_painter.wireframe = b
        self.update()

    def initialize(self):
        self.shader_manager = ShaderManager(self)
        # self.grid_painter.set_shader_manager(self.shader_manager)
        self.ico_grid_painter.set_shader_manager(self.shader_manager)
        self.axis_painter.set_shader_manager(self.shader_manager)
        # self.field_painter.set_shader_manager(self.shader_manager)
        self.projection = QMatrix4x4()
        self.projection.perspective(35.0, 1.0, 0.01, 15.0)

    def initializeGL(self):
        print("OpenGL: " + str(gl.glGetString(gl.GL_VERSION)))
        print("GLSL: " + str(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)))

        version = QtGui.QOpenGLVersionProfile()
        version.setVersion(4, 5)
 #       QtGui.QSurfaceFormat.DoubleBuffer
        version.setProfile(QtGui.QSurfaceFormat.CoreProfile)

        fmt = self.format()
        print("surface format: ", self.format())
        print("major {}, minor {}".format(fmt.majorVersion(),
                                          fmt.minorVersion()))
        print("OpenGL: " + str(gl.glGetString(gl.GL_VERSION)))
        print("GLSL: " + str(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)))

        self.initialize()
        # self.shader_manager.start_paint()
        # self.grid_painter.initializeGL()

        # self.shader_manager.end_paint()
        # self.shader_manager.start_paint_field()
        # self.field_painter.initializeGL()
        self.ico_grid_painter.initializeGL()
        # self.shader_manager.end_paint_field()
        # self.shader_manager.start_paint_pc()
        self.axis_painter.initializeGL()

        # self.shader_manager.end_paint_pc()
        gl.glEnable(gl.GL_DEPTH_TEST)

        gl.glEnable(gl.GL_BLEND)

        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        # gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
        self.idx = 0

    def paintGL(self):
        # print(time.time())
        gl.glClear(gl.GL_COLOR_BUFFER_BIT
                   | gl.GL_DEPTH_BUFFER_BIT)
 #        gl.glLoadIdentity()

        self.resize()
        ##############################################
        self.shader_manager.start_paint_pc()

        colour = QVector4D(1.0, 1.0, 1.0, 1.0)
        # self.shader_manager.set_colour_pc(colour)
        # view transformation
        view = QMatrix4x4()
        view.translate(0.0,  0.0, self.cam_r)
        view *= self.cam_orientation
        print(self.cam_theta, self.cam_phi, self.cam_r)
        #view.rotate(-self.cam_theta, 0.0, 1.0, 0.0)
        #view.rotate(-self.cam_phi, 1.0, 0.0, 0.0)
        self.shader_manager.set_view_pc(view)

        # projection transformation
        self.shader_manager.set_projection_pc(self.projection)

        model = QMatrix4x4()
        #        model.translate(-0.4, -0.8, 0.0)
        #        model.scale(0.05, 0.05, 1.0)
        self.shader_manager.set_model_pc(model)

        # self.ico_grid_painter.paint_grid()
        #         self.grid_painter.paint_grid(self.image)
        #         self.idx = (self.idx+1) % 3
        self.shader_manager.end_paint_pc()

        #################################
        self.shader_manager.start_paint_field()

        colour = QVector4D(1.0, 1.0, 1.0, 1.0)
        self.shader_manager.set_colour_field(colour)
        # view transformation
        self.shader_manager.set_view_field(view)

        # projection transformation
        self.shader_manager.set_projection_field(self.projection)

        model = QMatrix4x4()
        self.shader_manager.set_model_field(model)

        # self.field_painter.paint_field(self.image)
        self.ico_grid_painter.paint_vector_field()
        self.shader_manager.end_paint_field()
        ############################################
        self.shader_manager.start_paint_pc()
        model = QMatrix4x4()
        self.shader_manager.set_model(model)
        gl.glViewport(0, 0, self.width//10, self.height//10)
        #view = QMatrix4x4()
        #view.translate(0.0,  0.0, -4.0)
        #view.rotate(self.cam_theta, 1.0, 0.0, 0.0)
        #view.rotate(self.cam_phi, 0.0, 1.0, 0.0)

        self.shader_manager.set_view_pc(view)
        # projection transformation

        self.shader_manager.set_projection_pc(self.projection)

        model = QMatrix4x4()
#        model.translate(-0.4, -0.8, 0.0)
#        model.scale(0.05, 0.05, 1.0)
        self.shader_manager.set_model_pc(model)

        # self.axis_painter.paint_axis()

        self.shader_manager.end_paint_pc()

    def resizeGL(self, width, height):
        self.width = width
        self.height = height

    def resize(self):
        gl.glViewport(0, 0, self.width, self.height)

        aspectratio = self.width/self.height

        self.projection.setToIdentity()
        scale = 1.0
        self.projection.perspective(35.0, aspectratio, 0.01, 15.0)

    def mousePressEvent(self, e):

        self.x = e.x()
        self.y = e.y()
        if e.button() == QtCore.Qt.LeftButton:
            self.left_button_pressed = True
        elif e.button() == QtCore.Qt.RightButton:
            self.right_button_pressed = True

    def mouseReleaseEvent(self, e):
        self.x = e.x()
        self.y = e.y()

        if e.button() == QtCore.Qt.LeftButton:
            self.left_button_pressed = False
        elif e.button() == QtCore.Qt.RightButton:
            self.right_button_pressed = False

    def mouseMoveEvent(self, e):
        dx = e.x() - self.x
        dy = e.y() - self.y
        self.x = e.x()
        self.y = e.y()

        k = 50
        if self.left_button_pressed:

            drx = k*dx/self.width
            dry = k*dy/self.height
            c = self.cam_orientation.inverted()[0]

            v_x = c.column(0)
            v_y = c.column(1)
            v = (drx*v_y + dry*v_x).normalized()
            alpha = math.sqrt(math.pow(drx, 2.0)+math.pow(dry, 2.0))
            self.cam_orientation.rotate(alpha, v.x(), v.y(), v.z())

        if self.right_button_pressed:
            dr = dy/100.0
            self.cam_r += dr
            if self.cam_r > -1.0:
                self.cam_r = 0.0
            if self.cam_r < -20.0:
                self.cam_r = -20.0
        self.update()
