#!/usr/bin/python
import numpy as nu
import scipy.integrate as intgr


def rectangle_method(func, a, b, n):
    h = (b - a)/n
    xs = nu.arange(a + h/2, b, h)
    return h*nu.sum(map(func, xs))


def trapezoidal_rule(func, a, b, n):
    h = (b - a)/n
    xs = nu.arange(a + h, b, h)
    return h*((func(a) + func(b))/2 + nu.sum(map(func, xs)))


def simpson_rule(func, a, b, n):
    return (2*rectangle_method(func, a, b, n) + trapezoidal_rule(func, a, b, n))/3


def runge(j, j2, k,):
    return (2**k*j2 - j)/(2**k - 1)


def max_der2(a, b):
    #return max(map(lambda x: (nu.abs(2*nu.sin(x + 0.05)/(nu.cos(2*(x + 0.05)) - 1.0))), nu.arange(a, b, 0.0001)))
    return 20.008


def printing(func, a, b, n):

    def method_print(method, k):
        one_precision = method(func, a, b, n)
        double_precision = method(func, a, b, 2*n)

        print one_precision, "  ", double_precision, "  ",\
            runge(one_precision, double_precision, k), \
            " ", nu.abs(one_precision - intgr.romberg(func, a, b)),

    print "METHOD OF INTEGRATION       N =", n,\
        "           N =", 2*n, "            RUNGE             ERROR                 A"
    print "============================================================================================================"
    print "RECTANGLE METHOD    ",
    method_print(rectangle_method, 2)
    print "   ", (b - a)**3 * max_der2(a, b)/(24 * n**2), "\n"
    print "TRAPEZOIDAL METHOD  ",
    method_print(trapezoidal_rule, 2)
    print "  ", ((b - a)**3 * 20.008)/(12 * n**2), "\n"
    print "SIMPSON'S RULE      ",
    method_print(simpson_rule, 4)
    print "  ", (b - a)**5 * max_der2(a, b)/(2880 * n**4), "\n"


start = 0.0
finish = 0.4
f = lambda x: 1/(0.25 + nu.sin(0.05 + x))

printing(f, start, finish, 8)