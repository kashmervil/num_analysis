#!/usr/bin/python
import numpy as np
import operator
import scipy.integrate as int

alpha = -0.37
a = 0.0
b = 0.5


def solve_system2():
    a1 = -(b - a)**(2+alpha)/(2+alpha)
    b1 = (b - a)**(1 + alpha)/(1 + alpha)
    c1 = (b - a)**(3 + alpha)/(3 + alpha)
    a2 = -(b - a)**(3+alpha)/(3+alpha)
    b2 = (b - a)**(2 + alpha)/(2 + alpha)
    c2 = (b - a)**(4 + alpha)/(4 + alpha)

    solution = np.linalg.solve(np.array([[a1, b1], [a2, b2]]), np.array([-c1, -c2]))
    return np.roots([1., -solution[0], solution[1]])


def solve_system4():

    a1 = -(b - a)**(alpha + 4)/(alpha + 4)
    b1 = (b - a)**(alpha + 3)/(alpha + 3)
    c1 = -(b - a)**(alpha + 2)/(alpha + 2)
    d1 = (b - a)**(alpha + 1)/(alpha + 1)
    e1 = (b - a)**(alpha + 5)/(alpha + 5)

    a2 = -e1
    b2 = -a1
    c2 = -b1
    d2 = -c1
    e2 = (b - a)**(alpha + 6)/(alpha + 6)

    a3 = -e2
    b3 = -a2
    c3 = -b2
    d3 = -c2
    e3 = (b - a)**(alpha + 7)/(alpha + 7)

    a4 = -e3
    b4 = -a3
    c4 = -b3
    d4 = -c3
    e4 = (b - a)**(alpha + 8)/(alpha + 8)

    s = np.linalg.solve(np.array([[a1, b1, c1, d1], [a2, b2, c2, d2],
                                  [a3, b3, c3, d3], [a4, b4, c4, d4]]), np.array([-e1, -e2, -e3, -e4]))
    return np.roots([1., -s[0], s[1], -s[2], s[3]])


def w(xs):
        return lambda x, j: reduce(operator.mul, map(lambda z: (x - z)/(xs[j] - z), np.delete(xs, j)), 1.0)


def quadrature(xs, f):
    l_func = w(xs)

    def a_k(k):
        return int.quad(lambda x: x**alpha*l_func(x, k), a, b)[0]
    for i in range(len(xs)):
        if i == 0:
            print len(xs),
        else:
            print " ",
        print '\n      x%d= %.6f    A%d= %.8f' % ((i+1), xs[i], (i+1), a_k(i)),

    return sum(map(lambda x: a_k(x) * f(xs[x]), range(len(xs))))

print "======================================================================="
print "       Homework 7 \n " \
      "  Computation of a univariate definite integral "
print "   function is f(x) = x^(%.2f)*cos(x)" % alpha
print "N         X              A               Quadrature     Error"
truesol = int.quad(lambda x: x**alpha * np.cos(x), a, b)[0]
sol2 = quadrature(solve_system2(), lambda x: np.cos(x))
print "   %.8f   " % sol2, np.abs(sol2 - truesol)
sol4 = quadrature(solve_system4(), lambda x: np.cos(x))
print "   %.8f   " % sol4, np.abs(sol4 - truesol)




