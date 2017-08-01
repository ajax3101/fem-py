#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Конечно-элементная модель объекта расчета
###################################################################

import math
from fem_error import TFEMException

# Типы конечных элементов
FEType = [
    'fe_1d_2',
    'fe_2d_3',
    'fe_2d_4',
    'fe_3d_4',
    'fe_3d_8'
]


class TMesh:
    def __init__(self):
        self.mesh_file = ''     # Имя файла с данными о геометрической модели
        self.fe_type = ''       # Тип КЭ
        self.x = []             # Координаты узлов
        self.y = []
        self.z = []
        self.surface = []       # Связи граничных элементов
        self.fe = []            # Связи в КЭ
        self.freedom = 0        # Кол-во степеней свободы

    @staticmethod
    def get_fe_data(t):
        if t == 3:
            return 'fe_2d_3', 2, 3, 2
        elif t == 4:
            return 'fe_3d_4', 3, 4, 3
        elif t == 8:
            return 'fe_3d_8', 4, 8, 3
        elif t == 24:
            return 'fe_2d_4', 2, 4, 2
        elif t == 34:
            return 'fe_1d_2', 0, 2, 1
        else:
            raise TFEMException('unknown_fe_err')

    def load(self, name):
        try:
            self.mesh_file = name
            file = open(self.mesh_file)
            lines = file.readlines()
            file.close()
        except IOError:
            raise TFEMException('read_file_err')
        self.fe_type, size_surface, size_fe, self.freedom = self.get_fe_data(int(lines[0]))
        # Кол-во узлов
        n = int(lines[1])
        # Считываем узлы
        index = 2
        for i in range(0, n):
            self.x.append(float(lines[2 + i].split()[0]))
            if self.freedom > 1:
                self.y.append(float(lines[2 + i].split()[1]))
            if self.freedom > 2:
                self.z.append(float(lines[2 + i].split()[2]))
            index += 1
        # Кол-во КЭ
        n = int(lines[index])
        index += 1
        # Считываем КЭ
        for i in range(0, n):
            row = []
            for j in range(0, size_fe):
                row.append(int(lines[index].split()[j]))
            self.fe.append(row)
            index += 1
        # Кол-во ГЭ
        n = int(lines[index])
        index += 1
        # Считываем ГЭ
        for i in range(0, n):
            row = []
            for j in range(0, size_surface):
                row.append(int(lines[index].split()[j]))
            self.surface.append(row)
            index += 1

    def fe_name(self):
        if self.fe_type == 'fe_1d_2':
            return 'one-dimensional linear element (2 nodes)'
        elif self.fe_type == 'fe_2d_3':
            return 'linear triangular element (3 nodes)'
        elif self.fe_type == 'fe_2d_4':
            return 'quadrilateral element (4 nodes)'
        elif self.fe_type == 'fe_3d_4':
            return 'linear tetrahedron (4 nodes)'
        elif self.fe_type == 'fe_3d_8':
            return 'cube element (8 nodes)'

    def get_coord(self, i):
        return self.x[i], self.y[i] if (len(self.y)) else 0, self.z[i] if (len(self.z)) else 0

    # Вычисление длины (площади) заданного граничного элемента
    def square(self, index):
        x = [0]*len(self.surface[index])
        y = [0]*len(self.surface[index])
        z = [0]*len(self.surface[index])
        for i in range(0, len(self.surface[index])):
            x[i], y[i], z[i] = self.get_coord(self.surface[index][i])
        s = 0
        if len(x) == 2:     # Граничный элемент - отрезок
            s = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2)
        elif len(x) == 3:   # Граничный элемент - треугольник
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif len(x) == 4:   # Граничный элемент - четырехугольник
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))

            a = math.sqrt((x[0] - x[3])**2 + (y[0] - y[3])**2 + (z[0] - z[3])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[3])**2 + (y[2] - y[3])**2 + (z[2] - z[3])**2)
            p = 0.5*(a + b + c)
            s += math.sqrt(p*(p - a)*(p - b)*(p - c))
        return s

    # Вычисление объема (длины, площади) заданного конечного элемента
    def volume(self, index):
        x = [0]*len(self.fe[index])
        y = [0]*len(self.fe[index])
        z = [0]*len(self.fe[index])
        for i in range(0, len(self.fe[index])):
            x[i], y[i], z[i] = self.get_coord(self.fe[index][i])
        v = 0
        if self.fe_type == 'fe_1d_2':
            v = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2)
        elif self.fe_type == 'fe_2d_3':
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            v = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.fe_type == 'fe_2d_4':
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            v = math.sqrt(p*(p - a)*(p - b)*(p - c))
            a = math.sqrt((x[0] - x[3])**2 + (y[0] - y[3])**2 + (z[0] - z[3])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[3])**2 + (y[2] - y[3])**2 + (z[2] - z[3])**2)
            p = 0.5*(a + b + c)
            v += math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.fe_type == 'fe_3d_4':
            a = (x[1] - x[0])*(y[2] - y[0])*(z[3] - z[0]) + (x[3] - x[0])*(y[1] - y[0])*(z[2] - z[0]) + \
                (x[2] - x[0])*(y[3] - y[0])*(z[1] - z[0])
            b = (x[3] - x[0])*(y[2] - y[0])*(z[1] - z[0]) + (x[2] - x[0])*(y[1] - y[0])*(z[3] - z[0]) + \
                (x[1] - x[0])*(y[3] - y[0])*(z[2] - z[0])
            v = math.fabs(a - b)/6.0
        elif self.fe_type == 'fe_3d_8':
            ref = [[0, 1, 4, 7], [4, 1, 5, 7], [1, 2, 6, 7], [1, 5, 6, 7], [1, 2, 3, 7], [0, 3, 1, 7]]
            for i in range(0, 6):
                a = (x[ref[i][1]] - x[ref[i][0]])*(y[ref[i][2]] - y[ref[i][0]])*(z[ref[i][3]] - z[ref[i][0]]) + \
                    (x[ref[i][3]] - x[ref[i][0]])*(y[ref[i][1]] - y[ref[i][0]])*(z[ref[i][2]] - z[ref[i][0]]) + \
                    (x[ref[i][2]] - x[ref[i][0]])*(y[ref[i][3]] - y[ref[i][0]])*(z[ref[i][1]] - z[ref[i][0]])
                b = (x[ref[i][3]] - x[ref[i][0]])*(y[ref[i][2]] - y[ref[i][0]])*(z[ref[i][1]] - z[ref[i][0]]) + \
                    (x[ref[i][2]] - x[ref[i][0]])*(y[ref[i][1]] - y[ref[i][0]])*(z[ref[i][3]] - z[ref[i][0]]) + \
                    (x[ref[i][1]] - x[ref[i][0]])*(y[ref[i][3]] - y[ref[i][0]])*(z[ref[i][2]] - z[ref[i][0]])
                v += math.fabs(a - b)/6.0
        return v
