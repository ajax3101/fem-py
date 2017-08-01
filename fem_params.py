#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#              Параметры решения задачи с помощью МКЭ
###################################################################

from fem_error import TFEMException


# Тип задачи: статика или динамика
ProblemType = [
    'static',
    'dynamic'
]

# Метод решения СЛАУ (точный, приближенный
SolveMethod = [
    'direct',
    'iterative'
]

# Стандартные имена функций (перемещения, деформации и напряжения) и их агрументов
StdName = [
    'x',    # 0  - идентификатор первого аргумента (x)
    'y',    # 1  - ... второго аргумента
    'z',    # 2  - ... третьего аргумента
    't',    # 3  - ... времени
    'U',    # 4  - компонента вектора перемещений по первому направлению (x)
    'V',    # 5  - ... по y
    'W',    # 6  - ... по z
    'Exx',  # 7  - компонента тензора нормальных деформаций по xx
    'Eyy',  # 8  - ... по yy
    'Ezz',  # 9  - ... по zz
    'Exy',  # 10 - компонента тензора тангенциальных деформаций по xу
    'Exz',  # 11 - ... по xz
    'Eyz',  # 12 - ... по yz
    'Sxx',  # 13 - компонента тензора нормальных напряжений по xx
    'Syy',  # 14 - ... по yy
    'Szz',  # 15 - ... по zz
    'Sxy',  # 16 - компонента тензора тангенциальных напряжений по xу
    'Sxz',  # 17 - ... по xz
    'Syz',  # 18 - ... по yz
    'Ut',   # 19 - скорость по x
    'Vt',   # 20 - ... по y
    'Wt',   # 21 - ... по z
    'Utt',  # 22 - ускорение по x
    'Vtt',  # 23 - ... по y
    'Wtt'   # 24 - ... по z
]

# Типы краевых условий
BoundaryCondition = [
    'initial',
    'boundary',
    'volume',
    'surface',
    'concentrated'
]


# Описание краевого условия
class TBoundaryCondition:
    def __init__(self):
        self.direct = 0         # Направление (номер функции, для которой задается условие): 0 - по x и т.д.
        self.expression = ''    # Функциональное выражение, определяющее значение условия (например, 10^5)
        self.predicate = ''     # Предикат отбора узлов
        self.type = ''          # Тип краевого условия


# Базовые параметры расчета задачи теории упругости с помощью МКЭ
class TFEMParams:
    def __init__(self):
        self.problem_type = ''  # Тип задачи
        self.solve_method = ''  # Метод решения СЛАУ
        self.width = 12         # Формат вывода результатов
        self.precision = 5
        self.eps = 1.0E-6       # Точность вычислений
        self.density = 0        # Плотность материала
        self.damping = 0        # Параметр демпфирования
        self.t0 = 0             # Начальный момент времени расчета
        self.t1 = 0             # Конечный момент времени расчета
        self.th = 0             # Шаг по времени
        self.e = []             # Коэффициент упругости (модуль Юнга)
        self.m = []             # Коэффициент Пуассона
        self.names = StdName    # Список имен функций и их аргументов
        self.bc_list = []       # Список краевых условий
        self.var_list = {}      # Список вспомогательных переменных и их значений

    def __add_condition__(self, t, e, p, d):
        c = TBoundaryCondition()
        c.type = t
        c.direct = d
        c.expression = e
        c.predicate = p
        self.bc_list.append(c)

    def add_boundary_condition(self, e, p, d):
        self.__add_condition__('boundary', e, p, d)

    def add_initial_condition(self, e, p, d):
        self.__add_condition__('initial', e, p, d)

    def add_volume_load(self, e, p, d):
        self.__add_condition__('volume', e, p, d)

    def add_surface_load(self, e, p, d):
        self.__add_condition__('surface', e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.__add_condition__('concentrated', e, p, d)

    def add_variable(self, var, val):
        self.var_list.setdefault(var, val)

    def check_params(self):
        if self.solve_method == '':
            raise TFEMException('solve_method_err')
        if self.problem_type == '':
            raise TFEMException('problem_type_err')
        if not len(self.e) or self.e[0] == 0:
            raise TFEMException('elasticity_err')
        if self.problem_type == 'dynamic' and (self.t0 == self.t1 or self.th <= 0):
            raise TFEMException('time_err')
