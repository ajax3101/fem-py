#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Реализация функций форм конечных элементов
###################################################################

import math
from abc import abstractmethod
from numpy.linalg import solve, LinAlgError
from numpy import array
from numpy import zeros
from numpy.linalg import det
from numpy.linalg import inv
from fem_error import TFEMException
from fem_defs import eps


# Абстрактный базовый класс, описывающий конечный элемент (КЭ)
class TFE:
    def __init__(self):
        self.size = 0
        self.density = self.damping = 0.0
        self.e = []  # Модуль Юнга и коэффициент Пуассона
        self.m = []
        self.x = []  # Координаты вершин КЭ
        self.y = []
        self.z = []
        self.vx = []  # Компоненты вектора объемной нагрузки
        self.vy = []
        self.vz = []
        self.K = [[]]  # Матрица жесткости
        self.M = [[]]  # ... масс
        self.D = [[]]  # ... демпфирования
        self.c = [[]]  # Коэффициенты функций форм

    # Задание координат
    def set_coord(self, *args):
        if len(args) == 1:
            self.x = args[0]
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
        elif len(args) == 3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        self.__create__()

    # Задание объемной нагрузки
    def set_volume_load(self, *args):
        if len(args) == 1:
            self.vx = args[0]
        elif len(args) == 2:
            self.vx = args[0]
            self.vy = args[1]
        elif len(args) == 3:
            self.vx = args[0]
            self.vy = args[1]
            self.vz = args[2]

    # Задание параметров упругости
    def set_elasticity(self, p1, p2):
        self.e = p1
        self.m = p2

    # Задание плотности
    def set_density(self, d):
        self.density = d

    # Задание параметра демпфирования
    def set_damping(self, d):
        self.damping = d

    @abstractmethod
    # Вычисление функций форм КЭ
    def __create__(self):
        raise NotImplementedError('Method TFE.__create__ is pure virtual')

    # Формирование матриц жесткости, масс и демпфирования
    @abstractmethod
    def generate(self, is_static=True):
        raise NotImplementedError('Method TFE.generate is pure virtual')

    # Вычисления стандартных результатов КЭ
    @abstractmethod
    def calc(self, u):
        return [[]]


# Линейный (двухузловой) одномерный КЭ
class TFE1D2(TFE):
    def __init__(self):
        super().__init__()
        self.size = 2
        self.K = [
            [0, 0, 0],
            [0, 0, 0]
        ]
        self.M = [
            [0, 0, 0],
            [0, 0, 0]
        ]
        self.D = [
            [0, 0, 0],
            [0, 0, 0]
        ]
        self.c = [
            [0, 0],
            [0, 0]
        ]

    def __length__(self):
        return math.fabs(self.x[1] - self.x[0])

    def __create__(self):
        if self.__length__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        self.c[0][0] = self.x[1]/(self.x[1] - self.x[0])
        self.c[0][1] = -1.0/(self.x[1] - self.x[0])
        self.c[1][0] = self.x[0]/(self.x[0] - self.x[1])
        self.c[1][1] = -1.0/(self.x[0] - self.x[1])

    def calc(self, u):
        res = zeros((2, 2))
        res[0][0] = res[0][1] = u[0]*self.c[0][1] + u[1]*self.c[1][1]
        res[1][0] = res[1][1] = self.e[0]*(u[0]*self.c[0][1] + u[1]*self.c[1][1])
        return res

    def generate(self, is_static=True):
        self.K[0][0] = 2.0*self.__length__()*self.e[0]*self.c[0][1]**2
        self.K[0][1] = 2.0*self.__length__()*self.e[0]*self.c[0][1]*self.c[1][1]
        self.K[1][1] = 2.0*self.__length__()*self.e[0]*self.c[1][1]**2

        # Вычисление интеграла для объемных сил
        if len(self.vx):
            self.K[0][2] += self.vx[0]*(self.c[0][1]*(self.x[1]**2 - self.x[0]**2)*0.5 +
                                        self.c[0][0]*(self.x[1] - self.x[0]))
            self.K[1][2] += self.vx[1]*(self.c[1][1]*(self.x[1]**2 - self.x[0]**2)*0.5 +
                                        self.c[1][0]*(self.x[1] - self.x[0]))
        if not is_static:
            # Формирование матрицы массы
            k00 = (self.c[0][0]**2*(self.x[1] - self.x[0])) + \
                  (self.c[0][0]*self.c[0][1]*(self.x[1]**2 - self.x[0]**2)) + \
                  (1.0/3.0*self.c[0][1]**2*(self.x[1]**3 - self.x[0]**3))
            k01 = (self.c[0][0]*self.c[1][0]*(self.x[1] - self.x[0])) + \
                  (0.5*(self.c[0][0]*self.c[1][1] + self.c[0][1]*self.c[1][0])*(self.x[1]**2 - self.x[0]**2)) + \
                  (1.0/3.0*self.c[0][1]*self.c[1][1]*(self.x[1]**3 - self.x[0]**3))
            k11 = (self.c[1][0]**2*(self.x[1] - self.x[0])) + \
                  (self.c[1][0]*self.c[1][1]*(self.x[1]**2 - self.x[0]**2)) + \
                  (1.0/3.0*self.c[1][1]**2*(self.x[1]**3 - self.x[0]**3))

            self.M = [[0, 0], [0, 0]]
            self.M[0][0] = self.density*k00
            self.M[0][1] = self.density*k01
            self.M[1][0] = self.density*k01
            self.M[1][1] = self.density*k11

            # Формирование матрицы демпфирования
            self.D = [[0, 0], [0, 0]]
            self.D[0][0] = self.damping*k00
            self.D[0][1] = self.damping*k01
            self.D[1][0] = self.damping*k01
            self.D[1][1] = self.damping*k11


# Линейный (трехузловой) треугольный КЭ
class TFE2D3(TFE):
    def __init__(self):
        super().__init__()
        self.size = 3
        self.K = zeros((6, 7))
        self.M = zeros((6, 6))
        self.D = zeros((6, 6))
        self.c = zeros((3, 3))

    def __square__(self):
        a = math.sqrt((self.x[0] - self.x[1])*(self.x[0] - self.x[1]) + (self.y[0] - self.y[1])*(self.y[0] - self.y[1]))
        b = math.sqrt((self.x[0] - self.x[2])*(self.x[0] - self.x[2]) + (self.y[0] - self.y[2])*(self.y[0] - self.y[2]))
        c = math.sqrt((self.x[2] - self.x[1])*(self.x[2] - self.x[1]) + (self.y[2] - self.y[1])*(self.y[2] - self.y[1]))
        p = 0.5*(a + b + c)
        return math.sqrt(p*(p - a)*(p - b)*(p - c))

    def __create__(self):
        det0 = self.y[2]*self.x[1] - self.y[2]*self.x[0] - self.y[0]*self.x[1] - self.y[1]*self.x[2] + \
               self.y[1]*self.x[0] + self.y[0]*self.x[2]
        if math.fabs(det0) < eps:
            raise TFEMException('incorrect_fe_err')

        index = [[2, 1], [0, 2], [1, 0]]

        for i in range(0, 3):
            det1 = self.y[index[i][0]]*self.x[index[i][1]] - self.y[index[i][1]]*self.x[index[i][0]]
            det2 = self.y[index[i][1]] - self.y[index[i][0]]
            det3 = self.x[index[i][0]] - self.x[index[i][1]]
            self.c[i][0] = det1/det0
            self.c[i][1] = det2/det0
            self.c[i][2] = det3/det0

    def calc(self, u):
        m = self.m[0]
        g = self.e[0]/(2.0 + 2.0*m)
        k = self.e[0]/(1.0 - m**2)
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.c[j][1]
                dy[i][j] = self.c[j][2]
                res[0][i] += u[2*j]*dx[i][j]
                res[1][i] += u[2*j + 1]*dy[i][j]
                res[2][i] += u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j]
                res[3][i] += k*(u[2*j]*dx[i][j] + m*u[2*j + 1]*dy[i][j])
                res[4][i] += k*(u[2*j + 1]*dy[i][j] + m*u[2*j]*dx[i][j])
                res[5][i] += g*(u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j])
        return res

    def generate(self, is_static=True):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*self.e[0]/(1.0 - self.m[0]**2)
        # Производные функций формы
        shape_dx = array([self.c[0][1], self.c[1][1], self.c[2][1]])
        shape_dy = array([self.c[0][2], self.c[1][2], self.c[2][2]])
        # Матрица градиентов
        b = array([
            [shape_dx[0], 0.0, shape_dx[1], 0.0, shape_dx[2], 0.0],
            [0.0, shape_dy[0], 0.0, shape_dy[1], 0.0, shape_dy[2]],
            [shape_dy[0], shape_dx[0], shape_dy[1], shape_dx[1], shape_dy[2], shape_dx[2]]
        ])
        # Формирование локальных матриц жесткости, масс и демпфирования
        local_k = b.conj().transpose().dot(d).dot(b)*self.__square__()
        for i in range(0, 6):
            for j in range(i, 6):
                self.K[i][j] = local_k[i][j]
                if len(self.vx) or len(self.vy):
                    self.K[i][6] = self.vy[0]*self.__square__()/3.0 if j % 2 == 0 else self.vx[0]*self.__square__()/3.0
                if not is_static:
                    self.M[i][j] = self.density*self.__square__()/12.0 if i == j else self.density*self.__square__()/6.0
                    self.D[i][j] = self.damping*self.__square__()/12.0 if i == j else self.density*self.__square__()/6.0


# Линейный (четырехузловой) тетраэдральный КЭ
class TFE3D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        self.K = zeros((12, 13))
        self.c = zeros((4, 4))

    def __volume__(self):
        a = (self.x[1] - self.x[0])*(self.y[2] - self.y[0])*(self.z[3] - self.z[0]) + \
            (self.x[3] - self.x[0])*(self.y[1] - self.y[0])*(self.z[2] - self.z[0]) + \
            (self.x[2] - self.x[0])*(self.y[3] - self.y[0])*(self.z[1] - self.z[0])
        b = (self.x[3] - self.x[0])*(self.y[2] - self.y[0])*(self.z[1] - self.z[0]) + \
            (self.x[2] - self.x[0])*(self.y[1] - self.y[0])*(self.z[3] - self.z[0]) + \
            (self.x[1] - self.x[0])*(self.y[3] - self.y[0])*(self.z[2] - self.z[0])
        return math.fabs(a - b)/6.0

    def __create__(self):
        if self.__volume__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a = array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.z[i]
            b[j] = 1.0
            x = solve(a, b)
            self.c[j] = list(x)

    def calc(self, u):
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        l = 2.0*self.m[0]*g/(1.0 - 2.0*self.m[0])
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        dz = zeros((self.size, self.size))
        res = zeros((12, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.c[j][1]
                dy[i][j] = self.c[j][2]
                dz[i][j] = self.c[j][3]
                res[0][i] += u[3*j]*dx[i][j]
                res[1][i] += u[3*j + 1]*dy[i][j]
                res[2][i] += u[3*j + 2]*dz[i][j]
                res[3][i] += u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j]
                res[4][i] += u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j]
                res[5][i] += u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j]
                res[6][i] += 2.0*g*u[3*j]*dx[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[7][i] += 2.0*g*u[3*j + 1]*dy[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[8][i] += 2.0*g*u[3*j + 2]*dz[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[9][i] += g*(u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j])
                res[10][i] += g*(u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j])
                res[11][i] += g*(u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j])
        return res

    def generate(self, is_static=True):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
            ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])
        # Производные функций формы
        shape_dx = array([self.c[0][1], self.c[1][1], self.c[2][1], self.c[3][1]])
        shape_dy = array([self.c[0][2], self.c[1][2], self.c[2][2], self.c[3][2]])
        shape_dz = array([self.c[0][3], self.c[1][3], self.c[2][3], self.c[3][3]])
        # Матрица градиентов
        b = array([
            [shape_dx[0], 0.0, 0.0, shape_dx[1], 0.0, 0.0, shape_dx[2], 0.0, 0.0, shape_dx[3], 0.0, 0.0],
            [0.0, shape_dy[0], 0.0, 0.0, shape_dy[1], 0.0, 0.0, shape_dy[2], 0.0, 0.0, shape_dy[3], 0.0],
            [0.0, 0.0, shape_dz[0], 0.0, 0.0, shape_dz[1], 0.0, 0.0, shape_dz[2], 0.0, 0.0, shape_dz[3]],
            [shape_dy[0], shape_dx[0], 0.0, shape_dy[1], shape_dx[1], 0.0, shape_dy[2], shape_dx[2], 0.0, shape_dy[3],
             shape_dx[3], 0.0],
            [0.0, shape_dz[0], shape_dy[0], 0.0, shape_dz[1], shape_dy[1], 0.0, shape_dz[2], shape_dy[2], 0.0,
             shape_dz[3], shape_dy[3]],
            [shape_dz[0], 0.0, shape_dx[0], shape_dz[1], 0.0, shape_dx[1], shape_dz[2], 0.0, shape_dx[2], shape_dz[3],
             0.0, shape_dx[3]]
        ])
        # Формирование локальных матриц жесткости, масс и демпфирования
        local_k = b.conj().transpose().dot(d).dot(b)*self.__volume__()
        for i in range(0, 12):
            if len(self.vx) or len(self.vy) or len(self.vz):
                if i % 3 == 0:
                    self.K[i][12] = 0.25*self.vx[0]*self.__volume__()
                elif i % 3 == 1:
                    self.K[i][12] = 0.25*self.vy[0]*self.__volume__()
                elif i % 3 == 2:
                    self.K[i][12] = 0.25*self.vz[0]*self.__volume__()
            for j in range(i, 12):
                self.K[i][j] = local_k[i][j]
                if not is_static:
                    self.M[i][j] = 0.1*self.density*self.__volume__() if i == j else 0.05*self.density*self.__volume__()
                    self.D[i][j] = 0.1*self.damping*self.__volume__() if i == j else 0.05*self.density*self.__volume__()


# Билинейный четырехузловой двумерный КЭ
class TFE2D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        self.K = zeros((8, 9))
        self.c = zeros((4, 4))

    def __square__(self):
        return math.sqrt((self.x[0] - self.x[1])**2 + (self.y[0] - self.y[1])**2)

    def __create__(self):
        if self.__square__() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a = array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.x[i]*self.y[i]
            b[j] = 1.0
            x = solve(a, b)
            self.c[j] = list(x)

    def generate(self, is_static=True):
        # Параметры квадратур Гаусса
        xi = [-0.57735027, -0.57735027, 0.57735027, 0.57735027]
        eta = [-0.57735027, 0.57735027, -0.57735027, 0.57735027]
        w = [1.0, 1.0, 1.0, 1.0]
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
            ])*self.e[0]/(1.0 - self.m[0]**2)
        # Формирование локальных матриц жесткости, масс и демпфирования
        local_k = zeros((8, 8))
        local_m = zeros((8, 8))
        volume_load = zeros(8)
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(w)):
            # Изопараметрические функции формы и их производные
            shape = array([
                0.25*(1.0 - xi[i])*(1.0 - eta[i]),
                0.25*(1.0 + xi[i])*(1.0 - eta[i]),
                0.25*(1.0 + xi[i])*(1.0 + eta[i]),
                0.25*(1.0 - xi[i])*(1.0 + eta[i])
                ])
            shape_dxi = array([
                -0.25*(1.0 - eta[i]),
                0.25*(1.0 - eta[i]),
                0.25*(1.0 + eta[i]),
                -0.25*(1.0 + eta[i])
                ])
            shape_deta = array([
                -0.25*(1.0 - xi[i]),
                -0.25*(1.0 + xi[i]),
                0.25*(1.0 + xi[i]),
                0.25*(1.0 - xi[i])
                ])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x), sum(shape_dxi*self.y)],
                [sum(shape_deta*self.x), sum(shape_deta*self.y)]
                ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Матрица градиентов
            b = array([
                [shape_dx[0], 0.0, shape_dx[1], 0.0, shape_dx[2], 0.0, shape_dx[3], 0.0],
                [0.0, shape_dy[0], 0.0, shape_dy[1], 0.0, shape_dy[2], 0.0, shape_dy[3]],
                [shape_dy[0], shape_dx[0], shape_dy[1], shape_dx[1], shape_dy[2], shape_dx[2], shape_dy[3], shape_dx[3]]
                ])
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                [shape[0], 0.0, shape[1], 0.0, shape[2], 0.0, shape[3], 0.0],
                [0.0, shape[0], 0.0, shape[1], 0.0, shape[2], 0.0, shape[3]]
                ])
            local_k += b.conj().transpose().dot(d).dot(b)*jacobian*w[i]
            if not is_static:
                local_m += c.conj().transpose().dot(c)*jacobian*w[i]
            # Учет объемной нагрузки
            if len(self.vx) or len(self.vy):
                for j in range(0, 4):
                    volume_load[2*j] += self.vx[j]*shape[j]*jacobian*w[i]
                    volume_load[2*j + 1] += self.vy[j]*shape[j]*jacobian*w[i]

        for i in range(0, 8):
            for j in range(i, 8):
                self.K[i][j] = local_k[i][j]
                if not is_static:
                    self.M[i][j] = self.density*local_m[i][j]
                    self.D[i][j] = self.damping*local_m[i][j]
            self.K[i][8] = volume_load[i]
        import sys
#        print('*******************************')
#        for i in range(0, len(self.K)):
#            for j in range(0, len(self.K[0])):
#                sys.stdout.write('%f' % self.K[i][j])
#                sys.stdout.write(' ')
#            sys.stdout.write('\n')
#        print('*******************************')

    def calc(self, u):
        m = self.m[0]
        g = self.e[0]/(2.0 + 2.0*m)
        k = self.e[0]/(1.0 - m**2)
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        res = zeros((6, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = self.c[j][1] + self.c[j][3]*self.y[i]
                dy[i][j] = self.c[j][2] + self.c[j][3]*self.x[i]
                res[0][i] += u[2*j]*dx[i][j]
                res[1][i] += u[2*j + 1]*dy[i][j]
                res[2][i] += u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j]
                res[3][i] += k*(u[2*j]*dx[i][j] + m*u[2*j + 1]*dy[i][j])
                res[4][i] += k*(u[2*j + 1]*dy[i][j] + m*u[2*j]*dx[i][j])
                res[5][i] += g*(u[2*j]*dy[i][j] + u[2*j + 1]*dx[i][j])
        return res


# Восьмиузловой призматический КЭ
class TFE3D8(TFE):
    def __init__(self):
        super().__init__()
        self.size = 8
        self.K = zeros((24, 25))
        self.c = zeros((8, 8))

    def __create__(self):
        a = zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0]*self.size)
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i]
                a[i][2] = self.y[i]
                a[i][3] = self.z[i]
                a[i][4] = self.x[i]*self.y[i]
                a[i][5] = self.x[i]*self.z[i]
                a[i][6] = self.y[i]*self.z[i]
                a[i][7] = self.x[i]*self.y[i]*self.z[i]
            b[j] = 1.0
            try:
                x = solve(a, b)
            except LinAlgError:
                raise TFEMException('incorrect_fe_err')
            self.c[j] = list(x)

    def generate(self, is_static=True):
        # Параметры квадратур Гаусса
        xi = [-0.57735027, -0.57735027, -0.57735027, -0.57735027,  0.57735027, 0.57735027, 0.57735027, 0.57735027]
        eta = [-0.57735027, -0.57735027, 0.57735027, 0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027]
        psi = [-0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027, 0.57735027]
        w = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
            ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])
        # Формирование локальных матриц жесткости, масс и демпфирования
        local_k = zeros((24, 24))
        local_m = zeros((24, 24))
        volume_load = zeros(24)
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(w)):
            # Изопараметрические функции формы и их производные
            shape = array([
                0.125*(1.0 - xi[i])*(1.0 - eta[i])*(1.0 - psi[i]),
                0.125*(1.0 + xi[i])*(1.0 - eta[i])*(1.0 - psi[i]),
                0.125*(1.0 + xi[i])*(1.0 + eta[i])*(1.0 - psi[i]),
                0.125*(1.0 - xi[i])*(1.0 + eta[i])*(1.0 - psi[i]),
                0.125*(1.0 - xi[i])*(1.0 - eta[i])*(1.0 + psi[i]),
                0.125*(1.0 + xi[i])*(1.0 - eta[i])*(1.0 + psi[i]),
                0.125*(1.0 + xi[i])*(1.0 + eta[i])*(1.0 + psi[i]),
                0.125*(1.0 - xi[i])*(1.0 + eta[i])*(1.0 + psi[i])
                ])
            shape_dxi = array([
                -0.125*(1.0 - eta[i])*(1.0 - psi[i]),
                0.125*(1.0 - eta[i])*(1.0 - psi[i]),
                0.125*(1.0 + eta[i])*(1.0 - psi[i]),
                -0.125*(1.0 + eta[i])*(1.0 - psi[i]),
                -0.125*(1.0 - eta[i])*(1.0 + psi[i]),
                0.125*(1.0 - eta[i])*(1.0 + psi[i]),
                0.125*(1.0 + eta[i])*(1.0 + psi[i]),
                -0.125*(1.0 + eta[i])*(1.0 + psi[i])
                ])
            shape_deta = array([
                -0.125*(1.0 - xi[i])*(1.0 - psi[i]),
                -0.125*(1.0 + xi[i])*(1.0 - psi[i]),
                0.125*(1.0 + xi[i])*(1.0 - psi[i]),
                0.125*(1.0 - xi[i])*(1.0 - psi[i]),
                -0.125*(1.0 - xi[i])*(1.0 + psi[i]),
                -0.125*(1.0 + xi[i])*(1.0 + psi[i]),
                0.125*(1.0 + xi[i])*(1.0 + psi[i]),
                0.125*(1.0 - xi[i])*(1.0 + psi[i])
                ])
            shape_dpsi = array([
                -0.125*(1.0 - xi[i])*(1.0 - eta[i]),
                -0.125*(1.0 + xi[i])*(1.0 - eta[i]),
                -0.125*(1.0 + xi[i])*(1.0 + eta[i]),
                -0.125*(1.0 - xi[i])*(1.0 + eta[i]),
                0.125*(1.0 - xi[i])*(1.0 - eta[i]),
                0.125*(1.0 + xi[i])*(1.0 - eta[i]),
                0.125*(1.0 + xi[i])*(1.0 + eta[i]),
                0.125*(1.0 - xi[i])*(1.0 + eta[i])
                ])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x), sum(shape_dxi*self.y), sum(shape_dxi*self.z)],
                [sum(shape_deta*self.x), sum(shape_deta*self.y), sum(shape_deta*self.z)],
                [sum(shape_dpsi*self.x), sum(shape_dpsi*self.y), sum(shape_dpsi*self.z)]
                ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = (inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta) + \
                       (inverted_jacobi[0, 2]*shape_dpsi)
            shape_dy = (inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta) + \
                       (inverted_jacobi[1, 2]*shape_dpsi)
            shape_dz = (inverted_jacobi[2, 0]*shape_dxi + inverted_jacobi[2, 1]*shape_deta) + \
                       (inverted_jacobi[2, 2]*shape_dpsi)
            # Матрица градиентов
            b = array([
                [shape_dx[0], 0.0, 0.0, shape_dx[1], 0.0, 0.0, shape_dx[2], 0.0, 0.0, shape_dx[3], 0.0, 0.0,
                 shape_dx[4], 0.0, 0.0, shape_dx[5], 0.0, 0.0, shape_dx[6], 0.0, 0.0, shape_dx[7], 0.0, 0.0],
                [0.0, shape_dy[0], 0.0, 0.0, shape_dy[1], 0.0, 0.0, shape_dy[2], 0.0, 0.0, shape_dy[3], 0.0, 0.0,
                 shape_dy[4], 0.0, 0.0, shape_dy[5], 0.0, 0.0, shape_dy[6], 0.0, 0.0, shape_dy[7], 0.0],
                [0.0, 0.0, shape_dz[0], 0.0, 0.0, shape_dz[1], 0.0, 0.0, shape_dz[2], 0.0, 0.0, shape_dz[3], 0.0, 0.0,
                 shape_dz[4], 0.0, 0.0, shape_dz[5], 0.0, 0.0, shape_dz[6], 0.0, 0.0, shape_dz[7]],
                [shape_dy[0], shape_dx[0], 0.0, shape_dy[1], shape_dx[1], 0.0, shape_dy[2], shape_dx[2], 0.0,
                 shape_dy[3], shape_dx[3], 0.0, shape_dy[4], shape_dx[4], 0.0, shape_dy[5], shape_dx[5], 0.0,
                 shape_dy[6], shape_dx[6], 0.0, shape_dy[7], shape_dx[7], 0.0],
                [0.0, shape_dz[0], shape_dy[0], 0.0, shape_dz[1], shape_dy[1], 0.0, shape_dz[2], shape_dy[2], 0.0,
                 shape_dz[3], shape_dy[3], 0.0, shape_dz[4], shape_dy[4], 0.0, shape_dz[5], shape_dy[5], 0.0,
                 shape_dz[6], shape_dy[6], 0.0, shape_dz[7], shape_dy[7]],
                [shape_dz[0], 0.0, shape_dx[0], shape_dz[1], 0.0, shape_dx[1], shape_dz[2], 0.0, shape_dx[2],
                 shape_dz[3], 0.0, shape_dx[3], shape_dz[4], 0.0, shape_dx[4], shape_dz[5], 0.0, shape_dx[5],
                 shape_dz[6], 0.0, shape_dx[6], shape_dz[7], 0.0, shape_dx[7]],
                ])
            local_k += b.conj().transpose().dot(d).dot(b)*jacobian*w[i]
            if not is_static:
                # Вспомогательная матрица для построения матриц масс и демпфирования
                c = array([
                    [shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3], 0.0, 0.0,
                     shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7], 0.0, 0.0],
                    [0.0, shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3], 0.0,
                     0.0, shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7], 0.0],
                    [0.0, 0.0, shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3],
                     0.0, 0.0, shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7]]
                    ])
                local_m += c.conj().transpose().dot(c)*jacobian*w[i]
            # Учет объемной нагрузки
            if len(self.vx) or len(self.vy) or len(self.vz):
                for j in range(0, self.size):
                    volume_load[3*j + 0] += self.vx[j]*shape[j]*jacobian*w[j]
                    volume_load[3*j + 1] += self.vy[j]*shape[j]*jacobian*w[j]
                    volume_load[3*j + 2] += self.vz[j]*shape[j]*jacobian*w[j]
        for i in range(0, 24):
            for j in range(i, 24):
                self.K[i][j] = local_k[i][j]
                if not is_static:
                    self.M[i][j] = self.density*local_m[i][j]
                    self.D[i][j] = self.damping*local_m[i][j]
            self.K[i][24] = volume_load[i]


#        print('******************************************')
#        import sys
#        sz = 24
#        for i in range(0, sz):
#            for j in range(0, sz):
#                sys.stdout.write('%+E\t' % local_k[i][j])
#            sys.stdout.write('\n')
#        print('******************************************')

    def calc(self, u):
        g = self.e[0]/(2.0 + 2.0*self.m[0])
        l = 2.0*self.m[0]*g/(1.0 - 2.0*self.m[0])
        dx = zeros((self.size, self.size))
        dy = zeros((self.size, self.size))
        dz = zeros((self.size, self.size))
        res = zeros((12, self.size))
        for i in range(0, self.size):
            for j in range(0, self.size):
                dx[i][j] = (self.c[j][1] + self.c[j][4]*self.y[i]) + \
                           (self.c[j][5]*self.z[i] + self.c[j][7]*self.y[i]*self.z[i])
                dy[i][j] = (self.c[j][2] + self.c[j][4]*self.x[i]) + \
                           (self.c[j][6]*self.z[i] + self.c[j][7]*self.x[i]*self.z[i])
                dz[i][j] = (self.c[j][3] + self.c[j][5]*self.x[i]) + \
                           (self.c[j][6]*self.y[i] + self.c[j][7]*self.x[i]*self.y[i])
                res[0][i] += u[3*j]*dx[i][j]
                res[1][i] += u[3*j + 1]*dy[i][j]
                res[2][i] += u[3*j + 2]*dz[i][j]
                res[3][i] += u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j]
                res[4][i] += u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j]
                res[5][i] += u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j]
                res[6][i] += 2.0*g*u[3*j]*dx[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[7][i] += 2.0*g*u[3*j + 1]*dy[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[8][i] += 2.0*g*u[3*j + 2]*dz[i][j] + l*(u[3*j]*dx[i][j] + u[3*j + 1]*dy[i][j] + u[3*j + 2]*dz[i][j])
                res[9][i] += g*(u[3*j]*dy[i][j] + u[3*j + 1]*dx[i][j])
                res[10][i] += g*(u[3*j]*dz[i][j] + u[3*j + 2]*dx[i][j])
                res[11][i] += g*(u[3*j + 1]*dz[i][j] + u[3*j + 2]*dy[i][j])
        return res
