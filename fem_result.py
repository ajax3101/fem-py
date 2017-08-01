#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#      Реализация контейнера для хранения результатов расчета
###################################################################


class TResult:
    def __init__(self):
        self.name = ''      # Имя функции
        self.results = []   # Узловые значения
        self.t = 0          # Значение времени, для которого выполнен расчет

    def min(self):
        return min(self.results)

    def max(self):
        return max(self.results)

