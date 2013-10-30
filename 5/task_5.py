#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import operator


def lagrange_maker(ys, xs):
    def w_tmp(x, j):
        def remove_element(xss, k):
            return xss[:k] + xss[(k+1):]

        return reduce(operator.mul, map(lambda z: (x - z)/(xs[j] - z), remove_element(xs, j)), 1.0)

    return lambda x: sum([ys[i]*w_tmp(x, i) for i in range(0, len(xs))])


def w(xs, x):
    return reduce(operator.mul, map(lambda z: (x - z), xs), 1.0)


def max_value(func, a, b):
    return max(map(func, np.arange(a, b, (b - a)/20.0)))


def error_up_limit(xs, derivation, x, a, b):
    return np.abs(w(xs, x) / reduce(operator.mul, range(2, len(xs) + 1)) * max_value(derivation, a, b))


def cos_error(xs, x):
    der = lambda t: (np.cos(t) if len(xs) % 2 else np.sin(t)) * (-1)**((len(xs) + 2) / 2)

    return error_up_limit(xs, der, x, -3*np.pi/2, 4*np.pi/3)


def sin9_error(xs, x):
    der = lambda t: (np.cos(9*t) if (len(xs) + 1) % 2 else np.sin(9*t)) * (-1)**((len(xs) + 1) / 2) * 9**(len(xs) + 1)
    return error_up_limit(xs, der, x, -3*np.pi/2, 4*np.pi/3)


def sin9(x):
    return np.sin(9*x)


def full_lagrange(xs, func, error_func):
    ys = map(func, xs)
    lagrange = lagrange_maker(ys, xs)

    a = min(xs)
    b = max(xs)
    xi = np.arange(a, b, (b - a)/20.0)
    for i in range(len(xi)):
        print_str = map(lambda x: str.zfill(str(x), 9)[:9],
                        [xi[i], func(xi[i]), lagrange(xi[i]), np.abs(func(xi[i]) - lagrange(xi[i])), error_func(xs, xi[i])])
        print 'X = {0}   F(x) = {1}   L(x) = {2}  R(x) = {3}  Rn(x) = {4}'\
            .format(print_str[0], print_str[1], print_str[2], print_str[3], print_str[4])
    say = raw_input("\n\nPress any key to plot all curves\n")
    s = np.arange(-2*np.pi, 2*np.pi, 0.01)
    plt.plot(s, map(lagrange, s))
    plt.plot(xs, ys, 'ro')
    plt.plot(s, map(func, s))
    plt.show()


x_table = [np.pi*-3/2, -0.5*np.pi, np.pi, np.pi*4/3]

start = min(x_table)
finish = max(x_table)

mid_table = map(lambda x: x/2.0, x_table)
right_table = map(lambda x: finish - x, mid_table)
left_table = map(lambda x: start - x, mid_table)
double_table = [np.pi*3/2, -1*np.pi/5, np.pi, np.pi*4/3, np.pi*1/7, np.pi*4/9, np.pi*6/13, np.pi*16/45, np.pi*31/21]


def run_all(func, error_func):

    full_lagrange(x_table, func, error_func)

    raw_input("\nMiddle-based table")
    full_lagrange(mid_table, func, error_func)

    raw_input("\nLagrange Polynomial for left side table")
    full_lagrange(left_table, func, error_func)

    raw_input("\nLagrange Polynomial for right side table")
    full_lagrange(right_table, func, error_func)

    raw_input("\nLagrange Polynomial for double sized table")
    full_lagrange(double_table, func, error_func)


#lagrange = lagrange_maker(map(sin9, x_table), x_table)
raw_input("\n---------------------------------------\n Lagrange Polynomial for  f(x) = cos(x)\n")
run_all(np.cos, cos_error)

raw_input("\n---------------------------------------\n Lagrange Polynomial for  f(x) = sin(9x)\n")
run_all(sin9, sin9_error)
