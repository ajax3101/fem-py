#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fem_defs import DIR_X, DIR_Y, DIR_Z, INIT_U, INIT_V, INIT_U_T, INIT_V_T, INIT_U_T_T, INIT_V_T_T
from fem_object import TObject


def cylinder2():
    obj = TObject()
    # e = [6.5E+10]
    # m = [0.3]
    # E = 6.5E+10
    # nu = 0.3
    E = 5.9E+6
    nu = 0.489
    P = 0.4
    e1 = 9 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) * (
        (E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P))) / (
             3 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) + (
                 E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P)))
    m1 = (3 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) - 2 * (
        (E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P)))) / (
             6 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) + 2 * (
                 (E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P))))
    e = [float(e1)]
    m = [float(m1)]

    print(type(e))

    if obj.set_mesh('mesh/cyl2.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_variable('eps', 1.0E-6)
        obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y | DIR_Z)
        obj.add_boundary_condition('0', 'x=2', DIR_X | DIR_Y | DIR_Z)
        obj.add_concentrated_load('-1.0e+4*cos(atan2(y,z))', 'abs(y^2 + z^2 - 0.5^2) <= eps', DIR_Z)
        obj.add_concentrated_load('-1.0e+4*sin(atan2(y,z))', 'abs(y^2 + z^2 - 0.5^2) <= eps', DIR_Y)
        obj.add_concentrated_load('1.0e+4*cos(atan2(y,z))', 'abs(y^2 + z^2 - 0.25^2) <= eps', DIR_Z)
        obj.add_concentrated_load('1.0e+4*sin(atan2(y,z))', 'abs(y^2 + z^2 - 0.25^2) <= eps', DIR_Y)
        if obj.calc():
            obj.print_result()
            obj.plot('U')
            # obj.plot('V')
            # obj.plot('W')
            obj.plot('Exx')
            #        obj.plot('Eyy')
            #        obj.plot('Ezz')
            #        obj.plot('Exy')
            #        obj.plot('Exz')
            #        obj.plot('Eyz')
            obj.plot('Sxx')
            #        obj.plot('Syy')
            #        obj.plot('Szz')
            #        obj.plot('Sxy')
            #        obj.plot('Sxz')
            #        obj.plot('Syz')


def cylinder3():
    obj = TObject()
    #e = [6.5E+10]
    #m = [0.3]
    E = 5.9E+6
    nu = 0.489
    P = 0.4
    e1 = 9 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) * (
        (E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P))) / (
             3 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) + (
                 E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P)))
    m1 = (3 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) - 2 * (
        (E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P)))) / (
             6 * ((E / (3 - 6 * nu)) - (P * (E / (3 - 6 * nu))) / (1 - ((1 + nu) / (3 - 3 * nu)) * (1 - P))) + 2 * (
                 (E / (2 + 2 * nu)) - (P * (E / (2 + 2 * nu))) / (1 - ((8 - 10 * nu) / (15 - 15 * nu)) * (1 - P))))
    e = [float(e1)]
    m = [float(m1)]

    print(type(e))

    if obj.set_mesh('mesh/shell2.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'z=0', DIR_X | DIR_Y | DIR_Z)
        obj.add_surface_load('-1.0E+6*cos(atan2(y,x))', 'abs(x^2 + y^2 - 0.005^2) < 1.0e-6', DIR_X)
        obj.add_surface_load('-1.0E+6*sin(atan2(y,x))', 'abs(x^2 + y^2 - 0.005^2) < 1.0e-6', DIR_Y)
        if obj.calc():
            #obj.print_result()
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.plot('U')
            #obj.plot('V')
            #obj.plot('W')
            obj.plot('Exx')
            #        obj.plot('Eyy')
            #        obj.plot('Ezz')
            #        obj.plot('Exy')
            #        obj.plot('Exz')
            #        obj.plot('Eyz')
            obj.plot('Sxx')
            #        obj.plot('Syy')
            #        obj.plot('Szz')
            #        obj.plot('Sxy')
            #        obj.plot('Sxz')
            #        obj.plot('Syz')


if __name__ == "__main__":
    # beam()
    # head3d()
    # body1d()
    # cube()
    # console()
    # tank3()
    # cylinder()
    # cylinder2()
    cylinder3()
    # quad()
    # console4()
    # cube_test()
    # console_dynamic()

"""
1. Добавить загрузку названий функций в объект
3. OpenGL
"""
