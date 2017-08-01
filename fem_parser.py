#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
# Реализация интерпретатора арифметических и логических выражений
###################################################################

from fem_error import TFEMException
from fem_tree import TTree

# Типы лексем
token = [
    'delimiter',
    'digit',
    'function',
    'variable',
    'end'
]

# Функции
function = [
    'sin',
    'cos',
    'tan',
    'exp',
    'asin',
    'acos',
    'atan',
    'atan2',
    'sinh',
    'cosh',
    'tanh',
    'abs'
]

# Логические операции
boolean = [
    'not',
    'and',
    'or'
]

# Арифметические операции
operation = [
    '+',
    '-',
    '*',
    '/',
    '^',
    '>',
    '>=',
    '=',
    '<>',
    '<',
    '<='
]


# Класс, реализующий разбор и выполнение арифметических и логических выражений
class TParser:
    def __init__(self):
        self.error = self.code = self.token = self.token_type = ''
        self.result = TTree()
        self.variables = {}

    def add_variable(self, var, val=0.0):
        if var in self.variables:
            self.say_error('redefinition_err')
        self.variables.setdefault(var, val)

    def set_variable(self, var, val):
        if var not in self.variables:
            self.say_error('undef_err')
        self.variables[var] = val

    def set_code(self, c):
        # Удаление пробелов ...
        # self.code = re.sub(r'\s', '', c)
        self.code = c
        self.compile()
#        try:
#            self.compile()
#        except TFEMException as e:
#            e.print_error()

    def get_exp(self, result):
        self.get_token()
        if len(self.token) == 0:
            self.say_error('syntax_err')
        result = self.token_or(result)
#        result = self.token_add(result)
        return result

    def token_or(self, result):
        result = self.token_and(result)
        hold = TTree()
        while not self.token == 'end' and self.token == 'or':
            self.get_token()
            hold = self.token_and(hold)
            result = TTree(result, 'or', hold)
        return result

    def token_and(self, result):
        result = self.token_not(result)
        hold = TTree()
        while not self.token == 'end' and self.token == 'and':
            self.get_token()
            hold = self.token_not(hold)
            result = TTree(result, 'and', hold)
        return result

    def token_not(self, result):
        sign = self.token
        if (self.token_type == 'delimiter') and sign == 'not':
            self.get_token()
        result = self.token_add(result)
        if sign == 'not':
            result = TTree(sign, result)
        return result

    def token_add(self, result):
        result = self.token_mul(result)
        hold = TTree()
        while not self.token == 'end' and (self.token == '+' or self.token == '-' or self.token == '>' or
                                           self.token == '<' or self.token == '>=' or self.token == '<=' or
                                           self.token == '<>' or self.token == '='):
            sign = self.token
            self.get_token()
            hold = self.token_mul(hold)
            result = TTree(result, sign, hold)
        return result

    def token_mul(self, result):
        result = self.token_pow(result)
        hold = TTree()
        while self.token != 'end' and (self.token == '*' or self.token == '/'):
            sign = self.token
            self.get_token()
            hold = self.token_pow(hold)
            result = TTree(result, sign, hold)
        return result

    def token_pow(self, result):
        result = self.token_un(result)
        hold = TTree()
        while self.token != 'end' and (self.token == '^'):
            self.get_token()
            hold = self.token_brackets(hold)
            result = TTree(result, '^', hold)
        return result

    def token_un(self, result):
        sign = self.token
        if (self.token_type == 'delimiter') and (sign == '+' or sign == '-'):
            self.get_token()
        result = self.token_brackets(result)
        if sign == '+' or sign == '-':
            result = TTree(sign, result)
        return result

    def token_brackets(self, result):
        if self.token != 'end' and self.token == '(' and self.token_type == 'delimiter':
            self.get_token()
            result = self.token_or(result)
            if self.token != ')':
                self.say_error('brackets_err')
            self.get_token()
        else:
            result = self.token_prim(result)
        return result

    def token_prim(self, result):
        if self.token_type == 'digit':
            val = float(self.token)
            result = TTree(val)
        elif self.token_type == 'function':
            result = self.token_func(result)
        elif self.token_type == 'variable':
            result = TTree(self.variables[self.token])
        else:
            self.say_error('syntax_err')
        self.get_token()
        return result

    def token_func(self, result):
        fun_token = self.token
        hold = TTree()
        self.get_token()
        if self.token == 'end' or self.token != '(':
            self.say_error('syntax_err')
        self.get_token()
        result = self.token_add(result)
        if fun_token == 'atan2':
            if self.token != ',':
                self.say_error('syntax_err')
            self.get_token()
            hold = self.token_add(hold)
            result = TTree(result, 'atan2', hold)
        else:
            result = TTree(fun_token, result)
        if self.token != ')':
            self.say_error('syntax_err')
        return result

    def get_token(self):
        self.token_type = self.token = ''
        # Обработка пустой строки
        if len(self.code) == 0:
            self.token = 'end'
            self.token_type = 'delimiter'
            return
        # Удаление ведущих пробелов ...
        i = 0
        while i < len(self.code) and (self.code[i] == ' ' or self.code[i] == '\t'):
            i += 1
        self.code = self.code[i:]
        if len(self.code) == 0:
            self.token = 'end'
            self.token_type = 'delimiter'
            return
        # Обработка разделителей
        if '+-*/()^=><,'.find(self.code[0]) != -1:
            self.token = self.code[0]
            self.code = self.code[1:len(self.code)]
            # Проверка на наличие двойного разделителя
            if len(self.code) and '+-*/()^=><,'.find(self.code[0]) != -1:
                if self.token + self.code[0] in operation:
                    self.token += self.code[0]
                    self.code = self.code[1:len(self.code)]
            self.token_type = 'delimiter'
            return
        # Обработка чисел
        if self.code[0].isdigit():
            i = 0
            while i < len(self.code) and self.code[i].isdigit():
                i += 1
            self.token = self.code[0:i] if (i < len(self.code)) else self.code[0:]
            self.code = self.code[i:]
            if len(self.code) > 0 and self.code[0] == '.':
                self.token += '.'
                self.code = self.code[1:]
                i = 0
                while i < len(self.code) and self.code[i].isdigit():
                    i += 1
                self.token += self.code[0:i] if (i < len(self.code)) else self.code[0:]
                self.code = self.code[i:]
            if len(self.code) > 0 and (self.code[0] == 'E' or self.code[0] == 'e'):
                self.code = self.code[1:]
                self.token += 'E'
                if self.code[0] != '+' and self.code[0] != '-':
                    self.say_error('syntax_err')
                self.token += self.code[0]
                self.code = self.code[1:]
                i = 0
                while i < len(self.code) and self.code[i].isdigit():
                    i += 1
                self.token += self.code[0:i]
                self.code = self.code[i:]
            self.token_type = 'digit'
            return
        # Обработка функкций и переменных
        if self.code[0].isalpha():
            self.token = self.code[0]
            self.code = self.code[1:]
            i = 0
            while i < len(self.code) and (self.code[i].isalpha() or self.code[i] == '_' or self.code[i].isdigit()):
                i += 1
            self.token += self.code[0:i] if (i < len(self.code)) else self.code[0:]
            self.code = self.code[i:]
            if self.token in function:
                self.token_type = 'function'
            elif self.token in boolean:
                self.token_type = 'delimiter'
            elif self.token in self.variables:
                self.token_type = 'variable'
            else:
                self.say_error('undef_err')
            return

    def run(self):
        return self.result.value()

    def say_error(self, err):
        self.error = err
        raise TFEMException(self.error)

    def compile(self):
        self.result = self.get_exp(self.result)
        if self.token_type == 'delimiter' and self.token != 'end':
            if self.token == '(' or self.token == ')':
                self.say_error('brackets_err')
            else:
                self.say_error('syntax_err')


# code = '               (        -    sin(pi * 0.5)^2)*cos(pi/2) + 2.5           '
# code = 'atan2(pi,2)<3.14159'
# code = 'sin(3.14)>=3.14159'
# code = '2*sin(pi/6)*cos(pi/6)'

# code = '((sin(x^2) + cos(y^2) <= 0) or x >= 0 and x < y)'
# parser = TParser()
# parser.add_variable('x', 0)
# parser.add_variable('y', 0)
# parser.add_variable('R', 2)
# parser.set_code(code)
# if parser.error == '':
#   print(parser.run())
