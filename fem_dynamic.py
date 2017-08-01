#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

import math
from scipy.sparse import lil_matrix, coo_matrix
from numpy import zeros, savez, load
from fem_defs import INIT_U, INIT_V, INIT_W, INIT_U_T, INIT_V_T, INIT_W_T, INIT_U_T_T, INIT_V_T_T, INIT_W_T_T
from fem_static import TFEMStatic


# Сохранение разреженной матрицы в файл
def save_matrix(file_name, matrix):
    matrix_coo = matrix.tocoo()
    row = matrix_coo.row
    col = matrix_coo.col
    data = matrix_coo.data
    shape = matrix_coo.shape
    savez(file_name, row=row, col=col, data=data, shape=shape)


# Загрузка разреженной матрицы из файла
def load_matrix(file_name):
    matrix = load(file_name)
    return coo_matrix((matrix['data'], (matrix['row'], matrix['col'])), shape=matrix['shape']).tolil()


class TFEMDynamic(TFEMStatic):
    def __init__(self):
        super().__init__()
        self.__global_matrix_mass__ = lil_matrix((0, 0))        # Глобальная матрица масс (ГММ)
        self.__global_matrix_damping__ = lil_matrix((0, 0))     # Глобальная матрица демпфирования (ГМД)

    # Расчет динамической задачи методом конечных элементов
    def __calc_problem__(self):
        size = len(self.__mesh__.x)*self.__mesh__.freedom
        self.__global_matrix_stiffness__ = lil_matrix((size, size))
        self.__global_matrix_mass__ = lil_matrix((size, size))
        self.__global_matrix_damping__ = lil_matrix((size, size))
        self.__global_load__ = [0]*size
        fe = self.__create_fe__()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        fe.set_damping(self.__params__.damping)
        fe.set_density(self.__params__.density)

        # Создание глобальных матриц жесткости, масс и демпфирования
        self.__progress__.set_process('Assembling global stiffness, mass and damping matrix...', 1,
                                      len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = [0]*len(self.__mesh__.fe[i])
            y = [0]*len(self.__mesh__.fe[i])
            z = [0]*len(self.__mesh__.fe[i])
            # Настройка КЭ
            for j in range(len(self.__mesh__.fe[i])):
                x[j], y[j], z[j] = self.__mesh__.get_coord(self.__mesh__.fe[i][j])
            fe.set_coord(x, y, z)
            fe.generate(False)
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly__(fe, i)
        # Формирование левой части СЛАУ
        self.__create_dynamic_matrix__()
        # Учет начальных условий
        u0, ut0, utt0 = self.__prepare_initial_condition__()
        # Итерационный процесс по времени
        t = self.__params__.t0
        while t <= self.__params__.t1:
            print('t = %5.2f' % t)
            # Формирование правой части СЛАУ
            self.__create_dynamic_vector__(u0, ut0, utt0, t)
            # Учет краевых условий
            self.__use_boundary_condition__()
            # Решение СЛАУ
            if not self.__solve__():
                print('The system of equations is not solved!')
                return False
            u0, ut0, utt0 = self.__calc_dynamic_results__(u0, ut0, utt0, t)
            t += self.__params__.th
            if math.fabs(t - self.__params__.t1) < self.__params__.eps:
                t = self.__params__.t1
        print('**************** Success! ****************')
        return True

    # Извлечение начальных условий
    def __prepare_initial_condition__(self):
        u0 = zeros(len(self.__mesh__.x)*self.__mesh__.freedom)
        ut0 = zeros(len(self.__mesh__.x)*self.__mesh__.freedom)
        utt0 = zeros(len(self.__mesh__.x)*self.__mesh__.freedom)
        parser = self.__create_parser__()
        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'initial':
                counter += 1
        self.__progress__.set_process('Using initial conditions...', 1, counter*len(self.__mesh__.x))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'initial':
                parser.set_code(self.__params__.bc_list[i].expression)
                value = parser.run()
                direct = self.__params__.bc_list[i].direct
                for j in range(0, len(self.__mesh__.x)):
                    self.__progress__.set_progress(counter)
                    counter += 1
                    if direct & INIT_U:
                        u0[j*self.__mesh__.freedom + 0] = value
                    if direct & INIT_V:
                        u0[j*self.__mesh__.freedom + 1] = value
                    if direct & INIT_W:
                        u0[j*self.__mesh__.freedom + 2] = value
                    if direct & INIT_U_T:
                        ut0[j*self.__mesh__.freedom + 0] = value
                    if direct & INIT_V_T:
                        ut0[j*self.__mesh__.freedom + 1] = value
                    if direct & INIT_W_T:
                        ut0[j*self.__mesh__.freedom + 2] = value
                    if direct & INIT_U_T_T:
                        utt0[j*self.__mesh__.freedom + 0] = value
                    if direct & INIT_V_T_T:
                        utt0[j*self.__mesh__.freedom + 1] = value
                    if direct & INIT_W_T_T:
                        utt0[j*self.__mesh__.freedom + 2] = value
        return u0, ut0, utt0

    # Вычисление напряжений, деформаций, скоростей и ускорений
    def __calc_dynamic_results__(self, u0, ut0, utt0, t):
        # Вычисление деформаций и напряжений
        super().__calc_results__(t)
        # Вычисление скоростей и ускорений (конечными разностями)
        th = self.__params__.th
        for i in range(0, len(self.__mesh__.x)):
            for j in range(0, self.__mesh__.freedom):
                u_0 = u0[i*self.__mesh__.freedom + j]
                u_1 = self.__global_load__[i*self.__mesh__.freedom + j]
                u_t_0 = ut0[i*self.__mesh__.freedom + j]
                u_t_1 = (u_1 - u_0)/th
                u_t_t_1 = (u_t_1 - u_t_0)/th
                u0[i*self.__mesh__.freedom + j] = u_1
                ut0[i*self.__mesh__.freedom + j] = u_t_1
                utt0[i*self.__mesh__.freedom + j] = u_t_t_1
                self.__result__[len(self.__result__) - 2*self.__mesh__.freedom + j].results[i] = u_t_1
                self.__result__[len(self.__result__) - self.__mesh__.freedom + j].results[i] = u_t_t_1
        return u0, ut0, utt0

    # Добавление ЛМЖ, ЛММ и ЛМД к ГМЖ
    def __assembly__(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self.__mesh__.fe[index][i//self.__mesh__.freedom]*self.__mesh__.freedom + i % self.__mesh__.freedom
            for j in range(i, len(fe.K)):
                l = self.__mesh__.fe[index][j//self.__mesh__.freedom]*self.__mesh__.freedom + j % self.__mesh__.freedom
                self.__global_matrix_stiffness__[k, l] += fe.K[i][j]
                self.__global_matrix_mass__[k, l] += fe.M[i][j]
                self.__global_matrix_damping__[k, l] += fe.D[i][j]
                if k != l:
                    self.__global_matrix_stiffness__[l, k] += fe.K[i][j]
                    self.__global_matrix_mass__[l, k] += fe.M[i][j]
                    self.__global_matrix_damping__[l, k] += fe.D[i][j]
            self.__global_load__[k] += fe.K[i][len(fe.K)]

    # Формирование левой части (матрицы) уравнения квазистатического равновесия
    def __create_dynamic_matrix__(self):
        theta = 1.37
        self.__progress__.set_process('Creating static part of global matrix...', 1, 2)
        self.__global_matrix_stiffness__ += \
            self.__global_matrix_mass__.dot(6.0/(theta**2**self.__params__.th**2)) + \
            self.__global_matrix_damping__.dot(6.0/(3.0/(theta*self.__params__.th)))
        self.__progress__.set_progress(2)

    # Формирование правой части (вектора) уравнения квазистатического равновесия
    def __create_dynamic_vector__(self, u0, ut0, utt0, t):
        theta = 1.37
        k1 = 6.0/(theta**2**self.__params__.th**2)
        k2 = 3.0/(theta*self.__params__.th)
        k3 = 0.5*(theta*self.__params__.th)
        self.__global_load__ = zeros(len(self.__mesh__.x)*self.__mesh__.freedom)
        # Вычисление компонент нагрузки для текущего момента времени
        self.__prepare_concentrated_load__(t)
        self.__prepare_surface_load__(t)
        self.__prepare_volume_load__(t)
        self.__global_load__ += \
            self.__global_matrix_mass__.dot(u0.dot(k1) + ut0.dot(2.0*k2) + utt0.dot(2.0)) + \
            self.__global_matrix_damping__.dot(u0.dot(k2) + ut0.dot(2.0) + utt0.dot(k3))

    # Определение кол-ва результатов в зависимости от размерности задачи
    def __num_result__(self):
        res = 0
        if self.__mesh__.freedom == 1:
            # u, Exx, Sxx, ut, utt
            res = 5
        elif self.__mesh__.freedom == 2:
            # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
            res = 12
        elif self.__mesh__.freedom == 3:
            # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
            res = 21
        return res

    # Индекс функции в зависимости от размерности задачи
    def __index_result__(self, i):
        ret = 0
        # u, Exx, Sxx, ut, utt
        index4 = [4, 7, 13, 19, 22]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
        index5 = [4, 5, 7, 8, 10, 13, 14, 16, 19, 20, 22, 23]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
        index6 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        if self.__mesh__.freedom == 1:
            ret = index4[i]
        elif self.__mesh__.freedom == 2:
            ret = index5[i]
        elif self.__mesh__.freedom == 3:
            ret = index6[i]
        return ret
