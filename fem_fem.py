#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
# Реализация вычислений задач заданного типа методо конечных элементов
#######################################################################

from abc import abstractmethod
from fem_mesh import TMesh
from fem_params import TFEMParams
from fem_progress import TProgress
from fem_fe import TFE, TFE1D2, TFE2D3, TFE2D4, TFE3D4, TFE3D8
from fem_parser import TParser
from fem_error import TFEMException


# Абстрактный базовый класс, реализующий МКЭ
class TFEM:
    def __init__(self):
        self.__mesh__ = TMesh()                                 # Дискретная модель объекта
        self.__params__ = TFEMParams()                          # Параметры расчета
        self.__progress__ = TProgress()                         # Индикатор прогресса расчета
        self.__result__ = []                                    # Список результатов расчета

    @abstractmethod
    def __calc_problem__(self):
        raise NotImplementedError('Method TFEM.__calc_problem__ is pure virtual')

    # Добавление локальной матрицы жесткости (масс, демпфирования) к глобальной
    @abstractmethod
    def __assembly__(self, fe, index):
        raise NotImplementedError('Method TFEM.__assembly__ is pure virtual')

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...) 
    @abstractmethod
    def __calc_results__(self):
        raise NotImplementedError('Method TFEM.calc_results is pure virtual')

    # Прямое решение СЛАУ
    @abstractmethod
    def __solve_direct__(self):
        raise NotImplementedError('Method TFEM.__solve_direct__ is pure virtual')

    # Приближенное решение СЛАУ
    @abstractmethod
    def __solve_iterative__(self):
        raise NotImplementedError('Method TFEM.__solve_iterative__ is pure virtual')

    # Решение СЛАУ
    def __solve__(self):
        ret = False
        if self.__params__.solve_method == 'direct':
            ret = self.__solve_direct__()
        elif self.__params__.solve_method == 'iterative':
            ret = self.__solve_iterative__()
        return ret

    # Создание нужного типа КЭ
    def __create_fe__(self):
        fe = TFE()
        if self.__mesh__.fe_type == 'fe_1d_2':
            fe = TFE1D2()
        elif self.__mesh__.fe_type == 'fe_2d_3':
            fe = TFE2D3()
        elif self.__mesh__.fe_type == 'fe_2d_4':
            fe = TFE2D4()
        elif self.__mesh__.fe_type == 'fe_3d_4':
            fe = TFE3D4()
        elif self.__mesh__.fe_type == 'fe_3d_8':
            fe = TFE3D8()
        return fe

    # Настройка парсера
    def __create_parser__(self):
        parser = TParser()
        for i in range(0, len(self.__params__.names)):
            parser.add_variable(self.__params__.names[i])
        for key, value in self.__params__.var_list.items():
            parser.add_variable(key, value)
        return parser

    # Запуск процедуры расчета
    def calc(self):
        try:
            # Проверка наличия и соответствия необходимых параметров расчета
            self.__params__.check_params()
            ret = self.__calc_problem__()
        except TFEMException as err:
            ret = False
            err.print_error()
        return ret

    # Задание сетки
    def set_mesh(self, mesh):
        self.__mesh__ = mesh

    # Задание параметров расчета
    def set_params(self, params):
        self.__params__ = params

    # Возврат результатов расчета
    def get_result(self):
        return self.__result__
