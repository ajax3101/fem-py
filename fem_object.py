#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                       Описание объекта расчета
###################################################################

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from math import fabs, floor
from fem_mesh import TMesh
from fem_fem import TFEM
from fem_params import TFEMParams
from fem_static import TFEMStatic
from fem_dynamic import TFEMDynamic
from fem_defs import eps
from fem_error import TFEMException


# Вывод сообщения об ошибке
def error(err_msg):
    print('\033[1;31m%s\033[1;m' % err_msg)


class TObject:
    def __init__(self):
        self.__params__ = TFEMParams()  # Параметры расчета
        self.__mesh__ = TMesh()         # КЭ-модель
        self.__results__ = []           # Список результатов расчета для перемещений, деформаций, ...

    def set_mesh(self, name):
        try:
            self.__mesh__.load(name)
            print('Object: %s' % self.object_name())
            print('Points: %d' % len(self.__mesh__.x))
            print('FE: %d - %s' % (len(self.__mesh__.fe), self.__mesh__.fe_name()))
        except TFEMException as err:
            err.print_error()
            return False
        return True

    # Название объекта
    def object_name(self):
        return os.path.splitext(os.path.basename(self.__mesh__.mesh_file))[0]

    def set_problem_type(self, problem_type):
        self.__params__.problem_type = problem_type

    def set_solve_method(self, solve_method):
        self.__params__.solve_method = solve_method

    def set_eps(self, e):
        self.__params__.eps = e

    def set_width(self, width):
        self.__params__.width = width

    def set_precision(self, precision):
        self.__params__.precision = precision

    def set_elasticity(self, e, m):
        self.__params__.e = e
        self.__params__.m = m

    def set_density(self, density):
        self.__params__.density = density

    def set_time(self, t0, t1, th):
        self.__params__.t0 = t0
        self.__params__.t1 = t1
        self.__params__.th = th

    def set_damping(self, damping):
        self.__params__.damping = damping

    def set_names(self, names):
        self.__params__.names = names

    def add_boundary_condition(self, e, p, d):
        self.__params__.add_boundary_condition(e, p, d)

    def add_initial_condition(self, e, d):
        self.__params__.add_initial_condition(e, '', d)

    def add_volume_load(self, e, p, d):
        self.__params__.add_volume_load(e, p, d)

    def add_surface_load(self, e, p, d):
        self.__params__.add_surface_load(e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.__params__.add_concentrated_load(e, p, d)

    def add_variable(self, var, val):
        self.__params__.add_variable(var, val)

    def calc(self):
        fem = TFEM()
        if self.__params__.problem_type == 'static':
            fem = TFEMStatic()
        elif self.__params__.problem_type == 'dynamic':
            fem = TFEMDynamic()
        fem.set_mesh(self.__mesh__)
        fem.set_params(self.__params__)
        ret = fem.calc()
        if ret:
            self.__results__ = fem.get_result()
        return ret

    # Вывод результатов расчета
    def print_result(self, *argv):
        file = sys.stdout
        try:
            if len(argv) == 1:
                file = open(argv[0], 'w')
        except IOError:
            error('Error: unable to open file %s' % argv[0])
            return
        if self.__params__.problem_type == 'static':
            self.__print__(file)
        else:
            t = self.__params__.t0
            while t <= self.__params__.t1:
                file.write('t = %5.2f\n' % t)
                self.__print__(file, t)
                t += self.__params__.th
                if fabs(t - self.__params__.t1) < self.__params__.eps:
                    t = self.__params__.t1
        file.close()

    # Вывод результатов расчета для одного момента времени
    def __print__(self, file, t=0):
        # Определение ширины позиции
        len1 = len('%+*.*E' % (self.__params__.width, self.__params__.precision, 3.14159))
        len2 = len('%d' % len(self.__mesh__.x))
        # Вывод заголовка
        file.write('| %*s  (' % (len2, 'N'))
        for i in range(0, self.__mesh__.freedom):
            file.write(' %*s' % (len1, self.__params__.names[i]))
            if i < self.__mesh__.freedom - 1:
                file.write(',')
        file.write(') |')
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %*s |' % (len1, self.__results__[i].name))
        file.write('\n')
        for i in range(0, len(self.__mesh__.x)):
            file.write('| %*d  (' % (len2, i + 1))
            file.write(' %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.x[i]))
            if len(self.__mesh__.y):
                file.write(', %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.y[i]))
            if len(self.__mesh__.z):
                file.write(', %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.z[i]))
            file.write(') | ')
            for k in range(0, len(self.__results__)):
                if self.__results__[k].t == t:
                    file.write('%+*.*E' %
                               (self.__params__.width, self.__params__.precision, self.__results__[k].results[i]))
                    file.write(' | ')
            file.write('\n')
        file.write('\n')
        # Печать итогов
        file.write('|  %*s  ' % (len2, ' '))
        for i in range(0, self.__mesh__.freedom):
            file.write(' %*s' % (len1, ' '))
            if i < self.__mesh__.freedom - 1:
                file.write(' ')
        file.write('  |')
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %*s |' % (len1, self.__results__[i].name))
        file.write('\n')
        file.write('|   %*s  |' % (self.__mesh__.freedom*(len1 + 1) + self.__mesh__.freedom + len2, 'min:'))
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %+*.*E ' % (self.__params__.width, self.__params__.precision, self.__results__[i].min()))
                file.write('|')
        file.write('\n')
        file.write('|   %*s  |' % (self.__mesh__.freedom*(len1 + 1) + self.__mesh__.freedom + len2, 'max:'))
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %+*.*E ' % (self.__params__.width, self.__params__.precision, self.__results__[i].max()))
                file.write('|')
        file.write('\n\n\n')

    # Визуализация заданной функции
    def plot(self, fun_name, t=0):
        # Проверка корректности задания времени
        if self.__params__.problem_type == 'dynamic' and \
                ((t < self.__params__.t0 or t > self.__params__.t1) or t % self.__params__.th > eps):
            error('Error: incorrectly specified the time: %5.2f' % t)
            return
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.__results__)):
            if self.__results__[i].name == fun_name and self.__results__[i].t == t:
                index = i
                break
        if index == -1:
            error('Error: \'%s\' is not a recognized function name' % fun_name)
            return
        # Визуализация результата
        if self.__mesh__.fe_type == 'fe_1d_2':
            self.__plot_1d_linear__(index)
        elif self.__mesh__.fe_type == 'fe_2d_3':
            self.__plot_2d_tri__(index)
        elif self.__mesh__.fe_type == 'fe_2d_4':
            self.__plot_2d_quad__(index)
        elif self.__mesh__.fe_type == 'fe_3d_4':
            self.__plot_3d_tet__(index)
        elif self.__mesh__.fe_type == 'fe_3d_8':
            self.__plot_3d_hex__(index)
        # Задание заголовка
        if self.__params__.problem_type == 'dynamic':
            fun_name += ' (t = %5.2f)' % t

        plt.gcf().canvas.set_window_title('Result image')
        plt.title(fun_name)
        plt.show()

    # Визуализация заданной функции в случае одномерного линейного КЭ
    def __plot_1d_linear__(self, index):
        plt.plot(self.__mesh__.x, self.__results__[index].results, '-', linewidth=2)

    # Визуализация заданной функции в случае плоской треугольной сетки
    def __plot_2d_tri__(self, index):
        plt.figure()
        plt.gca().set_aspect('equal')
#        c_map = cm.get_cmap(name='terrain', lut=None)
        c_map = cm.get_cmap(name='spectral', lut=None)
        plt.triplot(self.__mesh__.x, self.__mesh__.y, self.__mesh__.fe, lw=0.5, color='white')
        plt.tricontourf(self.__mesh__.x, self.__mesh__.y, self.__mesh__.fe, self.__results__[index].results, cmap=c_map)
        plt.colorbar()

    # Визуализация заданной функции в случае плоской четырехугольной сетки
    def __plot_2d_quad__(self, index):
        plt.figure()
        plt.gca().set_aspect('equal')
        c_map = cm.get_cmap(name='spectral', lut=None)
        tri = np.array([np.array([T[0], T[1], T[2]]) for T in self.__mesh__.fe])
        plt.tricontourf(self.__mesh__.x, self.__mesh__.y, tri, self.__results__[index].results, cmap=c_map)
        tri = np.array([np.array([T[0], T[2], T[3]]) for T in self.__mesh__.fe])
        plt.tricontourf(self.__mesh__.x, self.__mesh__.y, tri, self.__results__[index].results, cmap=c_map)
        plt.colorbar()

    # Визуализация заданной функции в случае КЭ в форме тетраэдра
    def __plot_3d_tet__(self, index):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        triangle_vertices = np.array([np.array([[self.__mesh__.x[T[0]], self.__mesh__.y[T[0]], self.__mesh__.z[T[0]]],
                                                [self.__mesh__.x[T[1]], self.__mesh__.y[T[1]], self.__mesh__.z[T[1]]],
                                                [self.__mesh__.x[T[2]], self.__mesh__.y[T[2]], self.__mesh__.z[T[2]]]])
                                      for T in self.__mesh__.surface])

        c_map = cm.ScalarMappable()
        c_map.set_array([self.__results__[index].min(), self.__results__[index].max()])
        face_colors = self.get_surface_color(self.__results__[index].results)
        coll = Poly3DCollection(triangle_vertices, facecolors=face_colors, edgecolors='none')
        ax.add_collection(coll)
        ax.set_xlim(min(self.__mesh__.x), max(self.__mesh__.x))
        ax.set_ylim(min(self.__mesh__.y), max(self.__mesh__.y))
        ax.set_zlim(min(self.__mesh__.z), max(self.__mesh__.z))
        plt.colorbar(c_map)

    # Визуализация заданной функции в случае кубического КЭ
    def __plot_3d_hex__(self, index):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect("auto")
        ax.set_autoscale_on(True)

        # for i in range(0, len(self.__mesh__.surface)):
        #     for j in range(0, len(self.__mesh__.surface[i])):
        #         ind1 = j
        #         ind2 = j + 1 if j < len(self.__mesh__.surface[i]) - 1 else 0
        #         x = [self.__mesh__.x[self.__mesh__.surface[i][ind1]], self.__mesh__.x[self.__mesh__.surface[i][ind2]]]
        #         y = [self.__mesh__.y[self.__mesh__.surface[i][ind1]], self.__mesh__.y[self.__mesh__.surface[i][ind2]]]
        #         z = [self.__mesh__.z[self.__mesh__.surface[i][ind1]], self.__mesh__.z[self.__mesh__.surface[i][ind2]]]
        #         ax.plot3D(x, y, z, color="w")

        triangle_vertices1 = np.array([np.array([[self.__mesh__.x[T[0]], self.__mesh__.y[T[0]], self.__mesh__.z[T[0]]],
                                                 [self.__mesh__.x[T[1]], self.__mesh__.y[T[1]], self.__mesh__.z[T[1]]],
                                                 [self.__mesh__.x[T[2]], self.__mesh__.y[T[2]], self.__mesh__.z[T[2]]],
                                                 [self.__mesh__.x[T[0]], self.__mesh__.y[T[0]], self.__mesh__.z[T[0]]],
                                                 [self.__mesh__.x[T[2]], self.__mesh__.y[T[2]], self.__mesh__.z[T[2]]],
                                                 [self.__mesh__.x[T[3]], self.__mesh__.y[T[3]], self.__mesh__.z[T[3]]]])
                                       for T in self.__mesh__.surface])

        c_map = cm.ScalarMappable()
        c_map.set_array([self.__results__[index].min(), self.__results__[index].max()])
        face_colors = self.get_surface_color(self.__results__[index].results)
        coll = Poly3DCollection(triangle_vertices1, facecolors=face_colors, edgecolors='none')
        ax.add_collection(coll)

        ax.set_xlim(min(self.__mesh__.x), max(self.__mesh__.x))
        ax.set_ylim(min(self.__mesh__.y), max(self.__mesh__.y))
        ax.set_zlim(min(self.__mesh__.z), max(self.__mesh__.z))
        plt.colorbar(c_map)

    # Определене цвета поверхностной грани
    def get_surface_color(self, res):
        # Градации цвета (спектр)
        colors = [  # красный - желтый
                    [1.00, 0.00, 0.00], [1.00, 0.25, 0.00], [1.00, 0.50, 0.00], [1.00, 0.75, 0.00],
                    # желтый - зеленый
                    [1.00, 1.00, 0.00], [0.75, 1.00, 0.00], [0.50, 1.00, 0.00], [0.25, 1.00, 0.00],
                    # зеленый - фиолетовый
                    [0.00, 1.00, 0.00], [0.00, 1.00, 0.25], [0.00, 1.00, 0.50], [0.00, 1.00, 0.75],
                    # фиолетовый - синий
                    [0.00, 1.00, 1.00], [0.00, 0.75, 1.00], [0.00, 0.50, 1.00], [0.00, 0.00, 1.00]]
        u_min = min(res)
        u_max = max(res)
        face_colors = []
        for i in range(0, len(self.__mesh__.surface)):
            u = 0
            for j in range(0, len(self.__mesh__.surface[i])):
                u += res[self.__mesh__.surface[i][j]]
            u /= len(self.__mesh__.surface[i])
            index = floor((u - u_min)/((u_max - u_min)/16.0))
            if index > 15:
                index = 15
            face_colors.append(colors[index])
        return face_colors
