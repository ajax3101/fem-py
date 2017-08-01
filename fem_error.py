#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                   Реализация обработки ошибок
###################################################################


class TFEMException(Exception):
    def __init__(self, e):
        self.error = e

    def print_error(self):
        err_msg = '\nError: '
        if self.error == 'brackets_err':
            err_msg += 'unbalanced or unexpected bracket'
        elif self.error == 'syntax_err':
            err_msg += 'syntax error'
        elif self.error == 'undef_err':
            err_msg += 'undefined variable or function'
        elif self.error == 'redefinition_err':
            err_msg += 'redefinition variable'
        elif self.error == 'incorrect_fe_err':
            err_msg += 'incorrect finite element'
        elif self.error == 'read_file_err':
            err_msg += 'read file error'
        elif self.error == 'unknown_fe_err':
            err_msg += 'unknown finite element type'
        elif self.error == 'solve_method_err':
            err_msg += 'not specified method for solving linear systems'
        elif self.error == 'problem_type_err':
            err_msg += 'unknown problem type (static or dynamic)'
        elif self.error == 'elasticity_err':
            err_msg += 'not specified the parameters of elasticity'
        elif self.error == 'time_err':
            err_msg += 'incorrectly specified time parameters'
        else:
            err_msg += self.error
        print('\033[1;31m%s\033[1;m' % err_msg)
