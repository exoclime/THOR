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
from axispainter import AxisPainter
from vectorfieldpainter import VectorFieldPainter
from icogridpainter import IcoGridPainter
from planet_axis_equator_painter import PlanetAxisEquatorPainter
from markergrid import MarkerGridPainter


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

        self.ico_grid_painter = IcoGridPainter()
        self.axis_painter = AxisPainter()
        self.field_painter = VectorFieldPainter()
        self.axis_equator_painter = PlanetAxisEquatorPainter()
        #self.marker_grid_painter = MarkerGridPainter()

        self.cam_theta = 0.0
        self.cam_phi = 0.0

        self.cam_r = -5.0
        self.cam_orientation = QMatrix4x4()
        self.cam_orientation.rotate(-90, 1.0, 0.0, 0.0)

        self.right_button_pressed = False
        self.left_button_pressed = False
        self.mid_button_pressed = False

        self.image = 0

    def set_grid_data(self, *args):
        self.ico_grid_painter.set_grid_data(*args)
        self.field_painter.set_grid_data(*args)
        # self.marker_grid_painter.set_grid_data(*args)

    def set_level(self, level):
        self.ico_grid_painter.altitude = level
        self.field_painter.altitude = level
        self.update()

    def set_image(self, img):
        self.image = img
        self.ico_grid_painter.draw_idx = self.image
        self.field_painter.draw_idx = self.image
        self.update()

    def show_data(self, b):
        if b:
            self.ico_grid_painter.show_data = True
            self.ico_grid_painter.draw_idx = self.image
            self.field_painter.show_data = True
            self.update()

    def show_grid(self, b):
        if b:
            self.ico_grid_painter.show_data = False
            self.ico_grid_painter.draw_idx = 0
            self.field_painter.show_data = False
            self.field_painter.draw_idx = 0
            self.update()

    def show_wireframe(self, b):
        self.ico_grid_painter.wireframe = b
        self.update()

    def initialize(self):
        self.shader_manager = ShaderManager(self)
        self.ico_grid_painter.set_shader_manager(self.shader_manager)
        self.axis_painter.set_shader_manager(self.shader_manager)
        self.field_painter.set_shader_manager(self.shader_manager)
        self.axis_equator_painter.set_shader_manager(self.shader_manager)
        # self.marker_grid_painter.set_shader_manager(self.shader_manager)

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
        self.field_painter.initializeGL()
        self.ico_grid_painter.initializeGL()
        self.axis_painter.initializeGL()
        self.axis_equator_painter.initializeGL()
        # self.marker_grid_painter.initializeGL()
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
        self.shader_manager.load("icos_grid")

        colour = QVector4D(1.0, 1.0, 1.0, 1.0)
        # view transformation
        view = QMatrix4x4()
        view.translate(0.0,  0.0, self.cam_r)
        view *= self.cam_orientation
        self.shader_manager.set_uniform("icos_grid", "view", view)

        # projection transformation
        self.shader_manager.set_uniform(
            "icos_grid", "projection", self.projection)

        model = QMatrix4x4()
        self.shader_manager.set_uniform("icos_grid", "model", model)

        self.ico_grid_painter.paint_grid()
        self.shader_manager.release("icos_grid")
        # ##############################################
        # self.shader_manager.load("marker_grid")

        # colour = QVector4D(1.0, 1.0, 1.0, 1.0)
        # # view transformation
        # view = QMatrix4x4()
        # view.translate(0.0,  0.0, self.cam_r)
        # view *= self.cam_orientation
        # self.shader_manager.set_uniform("marker_grid", "view", view)

        # # projection transformation
        # self.shader_manager.set_uniform(
        #     "marker_grid", "projection", self.projection)

        # model = QMatrix4x4()
        # self.shader_manager.set_uniform("marker_grid", "model", model)

        # self.marker_grid_painter.paint_grid()
        # self.shader_manager.release("marker_grid")

        #################################
        self.shader_manager.load("icos_field")
        view = QMatrix4x4()
        view.translate(0.0,  0.0, self.cam_r)
        view *= self.cam_orientation
        colour = QVector4D(1.0, 1.0, 1.0, 1.0)
        self.shader_manager.set_uniform("icos_field", "colour", colour)
        # view transformation
        self.shader_manager.set_uniform("icos_field", "view", view)

        # projection transformation
        self.shader_manager.set_uniform(
            "icos_field", "projection", self.projection)

        model = QMatrix4x4()
        self.shader_manager.set_uniform("icos_field", "model", model)
        self.field_painter.paint_vector_field()
        self.shader_manager.release("icos_field")
        ############################################
        self.shader_manager.load("axis_equator")
        # view transformation
        view = QMatrix4x4()
        view.translate(0.0,  0.0, self.cam_r)
        view *= self.cam_orientation
        self.shader_manager.set_uniform("axis_equator", "view", view)

        # projection transformation
        self.shader_manager.set_uniform(
            "axis_equator", "projection", self.projection)

        model = QMatrix4x4()
        self.shader_manager.set_uniform("axis_equator", "model", model)
        self.axis_equator_painter.paint()

        self.shader_manager.release("axis_equator")

        ############################################
        self.shader_manager.load("axis")

        gl.glViewport(0, 0, self.width//8, self.height//8)
        #         #view.translate(0.0,  0.0, -4.0)
        #         #view.rotate(self.cam_theta, 1.0, 0.0, 0.0)
        #         #view.rotate(self.cam_phi, 0.0, 1.0, 0.0)
        view = QMatrix4x4()
        view.translate(0.0,  0.0, self.cam_r)
        view *= self.cam_orientation
        self.shader_manager.set_uniform("axis", "view", view)
        # projection transformation

        self.shader_manager.set_uniform("axis", "projection", self.projection)

        model = QMatrix4x4()
        # #        model.translate(-0.4, -0.8, 0.0)
        # #        model.scale(0.05, 0.05, 1.0)

        self.shader_manager.set_uniform("axis", "model", model)

        self.axis_painter.paint_axis()

        self.shader_manager.release("axis")

    def resizeGL(self, width, height):
        self.width = width
        self.height = height

    def resize(self):
        gl.glViewport(0, 0, self.width, self.height)

        aspectratio = self.width/self.height

        self.projection.setToIdentity()
        scale = 1.0
        self.projection.perspective(35.0, aspectratio, 0.01, 15.0)

    def reset_position(self):
        self.cam_orientation = QMatrix4x4()
        self.cam_orientation.rotate(90, 1.0, 0.0, 0.0)
        self.cam_r = -5.0
        self.update()

    def mousePressEvent(self, e):

        self.x = e.x()
        self.y = e.y()
        if e.button() == QtCore.Qt.LeftButton:
            self.left_button_pressed = True
        elif e.button() == QtCore.Qt.RightButton:
            self.right_button_pressed = True
        elif e.button() == QtCore.Qt.MidButton:
            self.mid_button_pressed = True

    def mouseReleaseEvent(self, e):
        self.x = e.x()
        self.y = e.y()

        if e.button() == QtCore.Qt.LeftButton:
            self.left_button_pressed = False
        elif e.button() == QtCore.Qt.RightButton:
            self.right_button_pressed = False
        elif e.button() == QtCore.Qt.MidButton:
            self.mid_button_pressed = False

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

        if self.mid_button_pressed:
            drx = k*dx/self.width
            c = self.cam_orientation.inverted()[0]

            v = c.column(2)

            self.cam_orientation.rotate(drx, v.x(), v.y(), v.z())

        if self.right_button_pressed:
            dr = dy/100.0
            self.cam_r += dr
            if self.cam_r > -1.0:
                self.cam_r = 0.0
            if self.cam_r < -20.0:
                self.cam_r = -20.0
        self.update()
