#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Simulation einer 2D Normalverteilung. Gezeigt werden neben zuf.
Realisierungen auch die Höhenlinie der Dichte für K=9.
"""

from __future__ import unicode_literals, division, print_function
from functools import partial

import numpy as np
import matplotlib as mp
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from PyQt4 import QtGui, QtCore


class Canvas(FigureCanvasQTAgg):
    """Ultimately, this is a QWidget"""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        """Creates a figure and axes and draws periodically on it."""
        # Create a figure and axes
        fig = mp.figure.Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_axes((0.1, 0.1, 0.85, 0.85))
        self.axes.set_aspect('equal')

        # Initialize widget and update timer
        super(Canvas, self).__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)
        self.updateGeometry()
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.plot_gaussian_samples)

        # plot parameter defaults
        self.sigma1 = self.sigma2 = 1.0
        self.rho = 0
        self.points_per_update = 100

        # plot data
        self.draw_contour_line()
        self.plot_gaussian_samples()
        timer.start(100)

    def set_param(self, name, value):
        """Update a plot parameter, clear axes"""
        setattr(self, name, value)
        self.axes.clear()
        self.draw_contour_line()

    def draw_contour_line(self):
        """Draw the theoretical contour line"""
        self.axes.set_xlim(-6, 6)
        self.axes.set_ylim(-5, 5)
        # get parameters and math functions
        o1, o2, r = self.sigma1, self.sigma2, self.rho
        sqrt, sin, cos, pi = np.sqrt, np.sin, np.cos, np.pi
        # calculate ellipse parameters
        K = 3 ** 2
        g = 0.5 * np.arctan(2 * r * o1 * o2 / (o2**2 - o1**2)) \
            if o1 != o2 else pi/4
        a = o1 * o2 * sqrt(K * (1 - r**2) / (
            (o1 * sin(g)) ** 2 +
            (o2 * cos(g)) ** 2 +
            2 * r * o1 * o2 * sin(g) * cos(g)
        ))
        b = o1 * o2 * sqrt(K * (1 - r**2) / (
            (o1 * cos(g)) ** 2 +
            (o2 * sin(g)) ** 2 -
            2 * r * o1 * o2 * sin(g) * cos(g)
        ))
        # add contour line (ellipse)
        self.axes.add_artist(mp.patches.Ellipse(
            xy=(0, 0), width=2 * a, height=2 * b, angle=180 / pi * g,
            facecolor='none', edgecolor='r', zorder=2, linewidth=2
        ))

    def plot_gaussian_samples(self):
        """Put some samples of the current distribution on the axes"""
        o1, o2, r = self.sigma1, self.sigma2, self.rho
        # get two std norm distributed vectors
        x, y = np.random.normal(0, 1, (2, self.points_per_update))
        # scaling parameters
        r1, r2 = np.sqrt((1 - r) / 2), np.sqrt((1 + r) / 2)
        # mix the random vectors to get desired correlation
        x, y = o1 * (x * r1 + y * r2), o2 * (x * r1 - y * r2)
        # plot and draw
        self.axes.plot(x, y, 'ko', zorder=1, alpha=0.5, ms=2)
        self.draw()


class FigureCanvasWithControls(QtGui.QWidget):

    def __init__(self):
        super(FigureCanvasWithControls, self).__init__()
        layout = QtGui.QVBoxLayout(self)

        canvas = Canvas(self, width=5, height=4, dpi=100)
        params = (('sigma1', 0.1, 2.0, 1.0),
                  ('sigma2', 0.1, 2.0, 1.0),
                  ('rho', -0.99, 0.99, 0.0))

        # create a control for each figure parameter
        for name, lo, hi, default in params:
            row = QtGui.QHBoxLayout()
            layout.addLayout(row)

            # label
            label = QtGui.QLabel(name)
            label.setMinimumWidth(50)
            label.setAlignment(QtCore.Qt.AlignRight)
            row.addWidget(label, 0)

            # value slider
            slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
            slider.setRange(0, 200)
            slider.setSingleStep(2)
            slider.setPageStep(10)
            row.addWidget(slider, 1)

            # value display
            text = QtGui.QLineEdit()
            text.setReadOnly(True)
            text.setMaximumWidth(50)
            text.setFocusPolicy(QtCore.Qt.NoFocus)
            text.setAlignment(QtCore.Qt.AlignRight)
            row.addWidget(text, 0)

            def update(name, lo, hi, text, value):
                """Convert int slider value to target range"""
                value = value / slider.maximum() * (hi - lo) + lo
                canvas.set_param(name, value)
                text.setText("{:.2f}".format(value))

            # update figure canvas on value change
            slider.valueChanged.connect(partial(update, name, lo, hi, text))
            slider.setValue(round((default - lo) / (hi - lo) * slider.maximum()))

        layout.addWidget(canvas)


if __name__ == '__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    win = FigureCanvasWithControls()
    win.show()
    sys.exit(app.exec_())

