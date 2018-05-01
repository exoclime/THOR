"""
Manage Shader programs and associated variables

TODO: might not be the best way to manage multiple different objects.
      Need to learn more about how to structure rendering and what is optimal
"""

from PyQt5.QtGui import (QOpenGLShader,
                         QOpenGLShaderProgram,
                         QMatrix4x4,
                         QVector4D)


class ShaderManager:

    vertex_shader = """#version 400\n
                       in vec3 vp;
                       uniform highp mat4 view;
                       uniform highp mat4 projection;
                       uniform highp mat4 model;
                       void main() {
                          gl_Position = projection * \
                              view * model *  vec4(vp, 1.0);
                       }"""

    fragment_shader = """#version 400\n
                        out vec4 frag_colour;
                        uniform highp vec4 colour;
                        void main() {

                           frag_colour = colour;
                        }"""

    def __init__(self, parent):
        self.m_program = QOpenGLShaderProgram(parent)

        if self.m_program.addShaderFromSourceCode(QOpenGLShader.Vertex,
                                                  self.vertex_shader):
            print("initialised vertex shader")
        else:
            print("vertex shader failed")

        self.m_program.addShaderFromSourceCode(QOpenGLShader.Fragment,
                                               self.fragment_shader)

        self.m_program.link()

        self.m_viewUniform = self.m_program.uniformLocation('view')
        self.m_projectionUniform = self.m_program.uniformLocation('projection')
        self.m_modelUniform = self.m_program.uniformLocation('model')

        self.m_colourUniform = self.m_program.uniformLocation('colour')

    def start_paint(self):
        self.m_program.bind()

    def end_paint(self):
        self.m_program.release()

    def set_colour(self, colour):
        self.m_program.setUniformValue(self.m_colourUniform, colour)

    def set_view(self, view):
        self.m_program.setUniformValue(
            self.m_viewUniform, view)

    def set_projection(self, projection):
        self.m_program.setUniformValue(
            self.m_projectionUniform, projection)

    def set_model(self, model):
        self.m_program.setUniformValue(
            self.m_modelUniform, model)
