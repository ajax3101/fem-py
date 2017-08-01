#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fem_defs import DIR_X, DIR_Y, DIR_Z, INIT_U, INIT_V, INIT_U_T, INIT_V_T, INIT_U_T_T, INIT_V_T_T
from fem_object import TObject


def body1d():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/body1d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_X)
        #    obj.add_volume_load('-1.0E+5', '', DIR_X)
        obj.add_concentrated_load('-1.0E+5', 'x=1', DIR_X)
        if obj.calc():
            obj.print_result()
            obj.plot('U')
            obj.plot('Exx')
            obj.plot('Sxx')


def cube():
    obj = TObject()
    e = [203200]
    m = [0.27]
    if obj.set_mesh('mesh/cube.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'z=0', DIR_X | DIR_Y | DIR_Z)
        #    obj.add_volume_load('-1000', '', DIR_Z)
        obj.add_surface_load('-1000', 'z=1', DIR_Z)
        #    obj.add_concentrated_load('-1000', 'z=1', DIR_Z)
        if obj.calc():
            obj.print_result()
            obj.plot('U')
            obj.plot('V')
            obj.plot('W')
            obj.plot('Exx')
            obj.plot('Eyy')
            obj.plot('Ezz')
            obj.plot('Exy')
            obj.plot('Exz')
            obj.plot('Eyz')
            obj.plot('Sxx')
            obj.plot('Syy')
            obj.plot('Szz')
            obj.plot('Sxy')
            obj.plot('Sxz')
            obj.plot('Syz')


def cube_test():
    obj = TObject()
    e = [203200]
    m = [0.27]
    if obj.set_mesh('mesh/cube_test.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Y)
        obj.add_surface_load('-1000', 'y=1', DIR_Y)
        if obj.calc():
            obj.print_result()


def beam():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/beam.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('iterative')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Y | DIR_Z)
        obj.add_volume_load('-1.0E+5', '', DIR_Y)
        if obj.calc():
            obj.print_result()
            obj.plot('U')
            obj.plot('V')
            obj.plot('W')
            obj.plot('Exx')
            obj.plot('Eyy')
            obj.plot('Ezz')
            obj.plot('Exy')
            obj.plot('Exz')
            obj.plot('Eyz')
            obj.plot('Sxx')
            obj.plot('Syy')
            obj.plot('Szz')
            obj.plot('Sxy')
            obj.plot('Sxz')
            obj.plot('Syz')


def console():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/console.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y)
        obj.add_concentrated_load('-1.0E+6', 'x=10', DIR_Y)
        if obj.calc():
            obj.print_result()
            obj.plot('U')
            obj.plot('V')
            obj.plot('Exx')
            obj.plot('Eyy')
            obj.plot('Exy')
            obj.plot('Sxx')
            obj.plot('Syy')
            obj.plot('Sxy')


def console4():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/console4.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y)
        #    obj.add_concentrated_load('-1.0E+5', 'x=10', DIR_Y)
        obj.add_volume_load('-1.0E+5', '', DIR_Y)
        if obj.calc():
            obj.print_result()
            obj.plot('U')
            obj.plot('V')
            obj.plot('Exx')


def quad():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/quad.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Y)
        obj.add_concentrated_load('-1.0E+5', 'y=1', DIR_Y)
        if obj.calc():
            obj.print_result()


def cylinder():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/cyl.trpa'):
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
            obj.plot('V')
            obj.plot('W')
            obj.plot('Exx')
            obj.plot('Eyy')
            obj.plot('Ezz')
            obj.plot('Exy')
            obj.plot('Exz')
            obj.plot('Eyz')
            obj.plot('Sxx')
            obj.plot('Syy')
            obj.plot('Szz')
            obj.plot('Sxy')
            obj.plot('Sxz')
            obj.plot('Syz')



def tank3():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/tank3.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_elasticity(e, m)
        obj.add_variable('eps', 1.0E-6)
        obj.add_variable('min', 0.0015)
        obj.set_width(10)
        obj.set_precision(5)
        obj.add_boundary_condition('0', 'y=-0.598 and abs(x^2+z^2-1.6635^2)<=eps', DIR_X | DIR_Y | DIR_Z)
        obj.add_boundary_condition('0', 'x=0', DIR_X)
        obj.add_boundary_condition('0', 'z=0', DIR_Z)
        obj.add_surface_load('1.0e+4*cos(atan2(z,x))',
                             '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037-min)^2) <= eps)', DIR_X)
        obj.add_surface_load('1.0e+4*sin(atan2(z,x))',
                             '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037-min)^2) <= eps)', DIR_Z)
        obj.add_surface_load('1.0e+4*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)', DIR_X)
        obj.add_surface_load('1.0e+4*cos(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)', DIR_Y)
        obj.add_surface_load('1.0e+4*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037-min)^2) <= eps)', DIR_Z)
        obj.add_surface_load('1.0e+4*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)', DIR_X)
        obj.add_surface_load('1.0e+4*cos(atan2((x^2+z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)', DIR_Y)
        obj.add_surface_load('1.0e+4*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,y))',
                             '(y > 0) and (abs(x^2 + y^2 + z^2 - (1.037-min)^2) <= eps)', DIR_Z)
        obj.add_surface_load('-5.0e+3*cos(atan2(z,x))',
                             '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037)^2) <= eps)', DIR_X)
        obj.add_surface_load('-5.0e+3*sin(atan2(z,x))',
                             '(y <= 0 and y>=-0.2690) and (abs(x^2 + z^2 - (1.037)^2) <= eps)', DIR_Z)
        obj.add_surface_load('-5.0e+3*cos(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)', DIR_X)
        obj.add_surface_load('-5.0e+3*cos(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)', DIR_Y)
        obj.add_surface_load('-5.0e+3*sin(atan2(z,x))*sin(atan2((x^2+z^2)^0.5,(y + 0.2690)))',
                             '(y < -0.2690) and (abs(x^2 + z^2 + (y + 0.2690)^2 - (1.037)^2) <= eps)', DIR_Z)
        obj.add_surface_load('-5.e+3', '(y=-1.7235) and (x^2+z^2 - 0.34205^2 <= eps)', DIR_Y)
        obj.add_surface_load('-5.e+3', '(y=-1.944) and (x^2+z^2 - 0.657857^2 <= eps and x^2+z^2 - 0.562143^2 >= eps)',
                             DIR_Y)
        obj.add_surface_load('5.0e+3*cos(atan2(z,x))', 'abs(y+0.6431) <= eps and abs(x^2 + z^2 - 1.6389^2) <= eps',
                             DIR_X)
        obj.add_surface_load('5.0e+3*sin(atan2(z,x))', 'abs(y+0.6431) <= eps and abs(x^2 + z^2 - 1.6389^2) <= eps',
                             DIR_Z)
        obj.add_surface_load(
            '5.0e+3*x*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+(z*(1.0644108554^2))^2)^0.5)',
            '(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps', DIR_X)
        obj.add_surface_load(
            '5.0e+3*(y+1.1013629509)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+(z*(1.0644108554^2))^2)^0.5)',
            '(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps', DIR_Y)
        obj.add_surface_load(
            '5.0e+3*z*(1.0644108554^2)/(((x*(1.0644108554^2))^2+(y+1.1013629509)^2+(z*(1.0644108554^2))^2)^0.5)',
            '(y>-0.6431 and y <-0.0234) and abs(y-((x^2+z^2)^0.5)*(-1.0644108554)-1.1013629509)<=eps', DIR_Z)
        obj.add_surface_load(
            '-5.0e+3*x*(1.0018498686^2)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+(y-1.3808172524)^2)^0.5)',
            '(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps', DIR_X)
        obj.add_surface_load(
            '5.0e+3*(y-1.3808172524)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+(y-1.3808172524)^2)^0.5)',
            '(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps', DIR_Y)
        obj.add_surface_load(
            '-5.0e+3*z*(1.0018498686^2)/(((x*(1.0018498686^2))^2+(z*(1.0018498686^2))^2+(y-1.3808172524)^2)^0.5)',
            '(y>-1.944 and y <-1.7235) and abs(y - ((x^2+z^2)^0.5)*(-1.0018498686)+1.3808172524)<=eps', DIR_Z)
        obj.add_surface_load(
            '5.0e+3*x*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+(3*z*(1.3260378897^2))^2)^0.5)',
            '(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps', DIR_X)
        obj.add_surface_load(
            '5.0e+3*(y-2.8163434974)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+(3*z*(1.3260378897^2))^2)^0.5)',
            '(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps', DIR_Y)
        obj.add_surface_load(
            '5.0e+3*z*(1.3260378897^2)/(((3*x*(1.3260378897^2))^2+(y-2.8163434974)^2+(3*z*(1.3260378897^2))^2)^0.5)',
            '(y>-1.944 and y < -0.6431) and abs(y-((x^2+z^2)^0.5)*(1.3260378897)+2.8163434974)<=eps', DIR_Z)
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.plot('U')
            obj.plot('V')
            obj.plot('W')


def head3d():
    obj = TObject()
    e = [1000]
    m = [0.3]
    if obj.set_mesh('mesh/head3d.trpa'):
        obj.set_problem_type('static')
        obj.set_solve_method('direct')
        obj.set_width(10)
        obj.set_precision(5)
        #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_variable('eps', 1.0E-6)
        obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Z)
        obj.add_boundary_condition('0', 'y=991.3', DIR_X | DIR_Y | DIR_Z)
        obj.add_surface_load('-1*cos(atan2(z,x))', 'abs(x^2 + z^2 - 210^2) <=0.001', DIR_X)
        obj.add_surface_load('-1*sin(atan2(z,x))', 'abs(x^2 + z^2 - 210^2) <= 0.001', DIR_Z)
        if obj.calc():
            obj.print_result()
            obj.plot('U')


def console_dynamic():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    if obj.set_mesh('mesh/console.trpa'):
        obj.set_problem_type('dynamic')
        obj.set_solve_method('direct')
        obj.set_damping(1.0E+3)
        obj.set_time(0, 1.0, 0.25)
        obj.set_width(10)
        obj.set_precision(5)
        #    obj.set_solve_method('iterative')
        obj.set_elasticity(e, m)
        obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y)
        obj.add_concentrated_load('-1.0E+5*cos(t)', 'x=10', DIR_X)
        obj.add_initial_condition('0', INIT_U)
        obj.add_initial_condition('0', INIT_V)
        obj.add_initial_condition('0', INIT_U_T)
        obj.add_initial_condition('0', INIT_V_T)
        obj.add_initial_condition('0', INIT_U_T_T)
        obj.add_initial_condition('0', INIT_V_T_T)
        if obj.calc():
            obj.print_result('mesh/' + obj.object_name() + '.res')
            obj.plot('U', 0.0)
            obj.plot('U', 0.25)
            obj.plot('U', 0.5)
            obj.plot('U', 0.75)
            obj.plot('U', 1.0)


if __name__ == "__main__":
    # beam()
    # head3d()
    # body1d()
    # cube()
    # console()
    # tank3()
    # cylinder()
    # cylinder2()
    # quad()
    # console4()
    # cube_test()
    # console_dynamic()

"""
1. Добавить загрузку названий функций в объект
3. OpenGL
"""
