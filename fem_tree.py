#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
# Реализация дерева разбора арифметических и логических выражений
###################################################################

import math
from abc import abstractmethod


# Абстрактный базовый класс значения выражения
class TNode:
    @abstractmethod
    def value(self):
        raise NotImplementedError('Method TNode.value is pure virtual')


# Класс, реализующий дерево разбора арифметических выражений
class TTree:
    def __init__(self, *args):
        if len(args) == 0:
            self.node = TRealNode(0)
        elif len(args) == 1:
            self.node = TRealNode(args[0])
        elif len(args) == 2:
            self.node = TUnaryNode(args[0], args[1])
        elif len(args) == 3:
            self.node = TBinaryNode(args[0], args[1], args[2])

    def value(self):
        return self.node.value()


# Вещественная переменная
class TRealNode(TNode):
    def __init__(self, val):
        self.__val__ = val

    def value(self):
        return self.__val__


# Унарная операция
class TUnaryNode(TNode):
    def __init__(self, op, val):
        self.__op__ = op
        self.__val__ = val

    def value(self):
        if self.__op__ == '-':
            return -self.__val__.value()
        elif self.__op__ == '+':
            return +self.__val__.value()
        elif self.__op__ == 'abs':
            return math.fabs(self.__val__.value())
        elif self.__op__ == 'sin':
            return math.sin(self.__val__.value())
        elif self.__op__ == 'cos':
            return math.cos(self.__val__.value())
        elif self.__op__ == 'tan':
            return math.tan(self.__val__.value())
        elif self.__op__ == 'exp':
            return math.exp(self.__val__.value())
        elif self.__op__ == 'asin':
            return math.asin(self.__val__.value())
        elif self.__op__ == 'acos':
            return math.acos(self.__val__.value())
        elif self.__op__ == 'atan':
            return math.atan(self.__val__.value())
        elif self.__op__ == 'sinh':
            return math.sinh(self.__val__.value())
        elif self.__op__ == 'cosh':
            return math.cosh(self.__val__.value())
        elif self.__op__ == 'not':
            return 1 if self.__val__.value() == 0 else 0


# Бинарная операция
class TBinaryNode(TNode):
    def __init__(self, left, op, right):
        self.left = left
        self.__op__ = op
        self.right = right

    def value(self):
        if self.__op__ == '+':
            return self.left.value() + self.right.value()
        elif self.__op__ == '-':
            return self.left.value() - self.right.value()
        elif self.__op__ == '*':
            return self.left.value()*self.right.value()
        elif self.__op__ == '/':
            return self.left.value()/self.right.value()
        elif self.__op__ == '^':
            return math.pow(self.left.value(), self.right.value())
        elif self.__op__ == '=':
            return 1 if self.left.value() == self.right.value() else 0
        elif self.__op__ == '<>':
            return 0 if self.left.value() == self.right.value() else 1
        elif self.__op__ == '<':
            return 1 if self.left.value() < self.right.value() else 0
        elif self.__op__ == '<=':
            return 1 if self.left.value() <= self.right.value() else 0
        elif self.__op__ == '>':
            return 1 if self.left.value() > self.right.value() else 0
        elif self.__op__ == '>=':
            return 1 if self.left.value() >= self.right.value() else 0
        elif self.__op__ == 'or':
            return 0 if (self.left.value() or self.right.value()) == 0 else 1
        elif self.__op__ == 'and':
            return 0 if (self.left.value() and self.right.value()) == 0 else 1
        elif self.__op__ == 'atan2':
            return math.atan2(self.left.value(), self.right.value())
