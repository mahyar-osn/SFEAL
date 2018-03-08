#!/usr/bin/env python -W ignore::DeprecationWarning
import warnings
warnings.filterwarnings("ignore")
import os
os.environ['ETS_TOOLKIT'] = 'wx'
import matplotlib
matplotlib.use('wx')
from mayavi import mlab
import numpy


class FIGURE:

    def __init__(self, figure='SFEAL', bgcolor=(1., 1., 1.), res=20):
        self.figure = mlab.figure(figure, bgcolor=bgcolor)
        self.plots = {}
        self.pmesh = None
        self.mode = None
        self.sigma_1 = None
        self.sigma_2 = None
        self.interval = None
        self.resolution = res

    def clear(self, label=None):
        if label == None:
            labels = self.plots.keys()
        else:
            labels = [label]

        mlab.figure(self.figure.name)

        for label in labels:
            mlab_obj = self.plots.get(label)
            if mlab_obj != None:
                if mlab_obj.name == 'Surface':
                    mlab_obj.parent.parent.parent.remove()
                else:
                    mlab_obj.parent.parent.remove()
                self.plots.pop(label)

    def plot_surfaces(self, label, X, T, scalars=None, color=None, rep='surface', opacity=1.0):

        mlab.figure(self.figure.name)

        if color == None:
            color = (1, 0, 0)

        mlab_obj = self.plots.get(label)
        if mlab_obj == None:
            if scalars == None:
                self.plots[label] = mlab.triangular_mesh(X[:, 0], X[:, 1], X[:, 2], T, color=color, opacity=opacity,
                                                         representation=rep)
            else:
                self.plots[label] = mlab.triangular_mesh(X[:, 0], X[:, 1], X[:, 2], T, scalars=scalars, opacity=opacity)

        else:
            self.figure.scene.disable_render = True
            view = mlab.view()
            roll = mlab.roll()

            if X.shape[0] == mlab_obj.mlab_source.x.shape[0]:
                if scalars == None:
                    mlab_obj.mlab_source.set(x=X[:, 0], y=X[:, 1], z=X[:, 2])
                    mlab_obj.actor.property.color = color
                    mlab_obj.actor.property.opacity = opacity
                else:
                    mlab_obj.mlab_source.set(x=X[:, 0], y=X[:, 1], z=X[:, 2], scalars=scalars, opacity=opacity)

            else:
                self.clear(label)
                if scalars == None:
                    self.plots[label] = mlab.triangular_mesh(X[:, 0], X[:, 1], X[:, 2], T, color=color, opacity=opacity,
                                                             representation=rep)
                else:
                    self.plots[label] = mlab.triangular_mesh(X[:, 0], X[:, 1], X[:, 2], T, scalars=scalars,
                                                             opacity=opacity)
            mlab.view(*view)
            mlab.roll(roll)
            self.figure.scene.disable_render = False

    def plot_lines(self, label, X, color=None, size=0):

        nPoints = 0
        for x in X:
            nPoints += x.shape[0]

        Xl = numpy.zeros((nPoints, 3))
        connections = []

        ind = 0
        for x in X:
            Xl[ind:ind + x.shape[0], :] = x
            for l in range(x.shape[0] - 1):
                connections.append([ind + l, ind + l + 1])
            ind += x.shape[0]
        connections = numpy.array(connections)

        mlab.figure(self.figure.name)

        if color == None:
            color = (1, 0, 0)
        if size == None:
            size = 1

        mlab_obj = self.plots.get(label)
        if mlab_obj == None:
            self.plots[label] = mlab.points3d(Xl[:, 0], Xl[:, 1], Xl[:, 2], color=color, scale_factor=0)
            self.plots[label].mlab_source.dataset.lines = connections
            mlab.pipeline.surface(self.plots[label], color=(1, 1, 1),
                                  representation='wireframe',
                                  line_width=size,
                                  name='Connections')
        else:
            self.figure.scene.disable_render = True
            self.clear(label)
            self.plots[label] = mlab.points3d(Xl[:, 0], Xl[:, 1], Xl[:, 2], color=color, scale_factor=0)
            self.plots[label].mlab_source.dataset.lines = connections
            # ~ self.plots[label].mlab_source.update()
            mlab.pipeline.surface(self.plots[label], color=color,
                                  representation='wireframe',
                                  line_width=size,
                                  name='Connections')
            self.figure.scene.disable_render = False

    def spectrum(self, pmesh, mode, s1, s2, fissure=False):
        self.mode = mode
        self.sigma_1 = s1
        self.sigma_2 = s2

        print "\n\t=========================================\n"
        print "\t   SURFACE VARIATION ALONG MODE %d AXIS" % mode
        print "\n\t=========================================\n"

        pmesh.nodes['weights'].values[1:] = 0  # Need to reset weights to zero
        pmesh.update_pca_nodes()

        pmesh.nodes['weights'].values[mode] = self.sigma_1
        pmesh.update_pca_nodes()

        X, T = pmesh.get_surfaces(res=self.resolution)

        pmesh.nodes['weights'].values[1:] = 0  # Reset weights to zero again
        pmesh.update_pca_nodes()

        pmesh.nodes['weights'].values[mode] = self.sigma_2
        pmesh.update_pca_nodes()
        X2, T = pmesh.get_surfaces(res=self.resolution)

        dx = X2 - X
        dx = numpy.sqrt(numpy.sum(dx * dx, 1))

        self.plot_surfaces('SPECTRUM', X2, T, scalars=dx, opacity=1.0)

        if fissure:
            lines = []

            ### Left lung fissure
            line_index = [3]
            lines = pmesh.append_lines(lines, [2, 7, 12, 16, 31, 34, 35], line_index)
            line_index = [4]
            lines2 = pmesh.append_lines(lines, [21, 24], line_index)

            ## Righ lung fissure
            line_index = [1]
            lines3 = pmesh.append_lines(lines, [57, 58, 59, 60, 61, 62, 63, 64, 65, 66], line_index)
            line_index = [2]
            lines4 = pmesh.append_lines(lines, [47, 59, 64, 85, 91, 97], line_index)

            self.plot_lines('lines', lines, size=3, color=(0 / 255.0, 0 / 255.0, 0 / 255.0))

    def animation(self, pmesh, mode, s1, s2, t=20, fissure=False):
        import scipy

        self.mode = mode
        self.sigma_1 = s1
        self.sigma_2 = s2
        self.interval = t

        pmesh.nodes['weights'].values[1:] = 0  # Need to reset weights to zero
        pmesh.update_pca_nodes()

        i = 0
        print "\n\t=========================================\n"
        print "\t   CHANGING WEIGHTS ALONG MODE %d AXIS\n" % mode

        for w in scipy.linspace(self.sigma_1, self.sigma_2, self.interval):
            pmesh.nodes['weights'].values[mode] = w
            print "\t   sigma = %.2f" %w
            pmesh.update_pca_nodes()

            if fissure:
                lines = []
                ### Left lung fissure
                line_index = [3]
                lines = pmesh.append_lines(lines, [2, 7, 12, 16, 31, 34, 35], line_index)
                line_index = [4]
                lines2 = pmesh.append_lines(lines, [21, 24], line_index)

                ### Righ lung fissure
                line_index = [1]
                lines3 = pmesh.append_lines(lines, [57, 58, 59, 60, 61, 62, 63, 64, 65, 66], line_index)
                line_index = [2]
                lines4 = pmesh.append_lines(lines, [47, 59, 64, 85, 91, 97], line_index)

                self.plot_lines('lines', lines, size=3, color=(0 / 255.0, 0 / 255.0, 0 / 255.0), )

                Xz, Tz = pmesh.get_surfaces(res=self.resolution)
                self.plot_surfaces('ANIMATION', Xz, Tz, color=(244 / 255.0, 164 / 255.0, 96 / 255.0), opacity=1.0)
            else:

                Xz, Tz = pmesh.get_surfaces(res=self.resolution)
                self.plot_surfaces('ANIMATION', Xz, Tz, color=(244 / 255.0, 164 / 255.0, 96 / 255.0), opacity=1.0)
        print "\n\t=========================================\n"

    def show_lung(self, mesh, fissure=False):
        import morphic
        path_to_mesh = mesh
        mesh = morphic.Mesh(path_to_mesh)
        X, T = mesh.get_surfaces(res=self.resolution)

        print "\n\t=========================================\n"
        print "\t   LUNG SURFACE FROM SUBJECT %s" % path_to_mesh
        print "\n\t=========================================\n"

        self.plot_surfaces('A SUBJECTS MESH', X, T,
                            color=(238 / 255.0, 213 / 255.0, 183 / 255.0))

        if fissure:
            lines = []
            ### Left lung fissure
            line_index = [3]
            lines = mesh.append_lines(lines, [2, 7, 12, 16, 31, 34, 35], line_index)
            line_index = [4]
            lines2 = mesh.append_lines(lines, [21, 24], line_index)

            ### Righ lung fissure
            line_index = [1]
            lines3 = mesh.append_lines(lines, [57, 58, 59, 60, 61, 62, 63, 64, 65, 66], line_index)
            line_index = [2]
            lines4 = mesh.append_lines(lines, [47, 59, 64, 85, 91, 97], line_index)

            self.plot_lines('lines', lines, size=3, color=(0 / 255.0, 0 / 255.0, 0 / 255.0), )