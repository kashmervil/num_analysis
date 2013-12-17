import numpy as np
import matplotlib.pyplot as plot
lowBound = 0
upBound = 1
a = 0.5 #0.25 1.5
k = 1. #0.25 2.
h = 0.1


def euler_method(func, step):
    xs = np.arange(lowBound, upBound, step)
    ys = np.zeros(len(xs))
    for i in range(len(xs) - 1):
        ys[i + 1] = ys[i] + step * func(xs[i], ys[i])
    return xs, ys


def helper(xs):
    tmp = []
    for i in range(len(xs) - 1):
        tmp.append(xs[i + 1] - xs[i])
    return tmp


def pp(xs):
    if len(xs[-1]) == 1:
        return xs
    else:
        xs.append(helper(xs[-1]))
        return pp(xs)


def runge_kutta_method(func, step):
    xs = np.arange(lowBound, upBound, step)
    ys = np.zeros(len(xs))
    for i in range(len(xs) - 1):
        k1 = step * func(xs[i], ys[i])
        k2 = step * func(xs[i] + step/2, ys[i] + k1/2)
        k3 = step * func(xs[i] + step/2, ys[i] + k2/2)
        k4 = step * func(xs[i] + step, k3)
        ys[i + 1] = ys[i] + 1./6 * (k1 + 2*k2 + 2*k3 + k4)
    return xs, ys


def adams_method(func, step):
    xs, ys = runge_kutta_method(func, step)
    etta = map(lambda x, y: step*func(x, y), xs, ys)
    table = pp([etta])
    for i in range(5, len(xs)):
        ys[i] = ys[i - 1] + table[0][i - 1] + 0.5*table[1][i - 2]
        ys[i] += (5.0/12)*table[2][i - 3] + (3.0/8)*table[3][i - 4] + (251.0/720)*table[4][i - 5]
    return xs, ys, table


def error(step):

    ff = lambda x, y: (a - x + y**2)/(1 + k*y + x**2)
    fx = lambda x, y: (a-2*x-a*x**2+a*k*y-2*k*x*y + y**2-x**2*y**2+k*y**3)/(1+x**2+k*y)**2
    fy = lambda x, y: -(y*(1+2*a*x-x**2+k*y+2*x*y**2))/(1+x**2+k*y)**2
    xs, ys, _ = adams_method(ff, step)
    m1 = max(map(lambda x, y: np.abs(fx(x, y)), xs, ys))
    m2 = max(map(lambda x, y: np.abs(fx(x, y)), xs, ys))
    m3 = max(map(lambda x, y: np.abs(fy(x, y)), xs, ys))
    m4 = m2 + m1*m3
    er = map(lambda x, y: np.abs(x - y), euler_method(ff, step)[1], ys)
    est = map(lambda x: (m4/(2*m3))*step*np.exp(m3*x), xs)
    return er, est


def euler_print_error(step):
    print "\nEuler method for h =%.3f" % step
    print "Step        Error       Estimation      "
    er = error(step)
    for i in range(len(np.arange(lowBound, upBound, step))):
        print "  %d         %.4f       %.4f" % (i, er[0][i], er[1][i])


def euler_plot(step):
    xs, ys = euler_method(f, step)
    plot.plot(xs, ys, "g")
    xs, ys, _ = adams_method(f, h)
    plot.plot(xs, ys, "b")
    plot.show()


def print_table(xs):
    for j in range(len(xs[0])):
        for w in range(len(xs) - j):
            if w < 6:
                print ' ', '{0:.6f}'.format(xs[w][j]),
        print ''

########################################################

f = lambda x, y: (a - x + y**2)/(1 + k*y + x**2)

print "======================================================================================"
print "                    Homework 8 Cauchy problem "
euler_print_error(h*2)
euler_plot(h*2)
raw_input()
euler_print_error(h)
euler_plot(h)
raw_input()
euler_print_error(h/2)
euler_plot(h/2)

runge = runge_kutta_method(f, h)
adams = adams_method(f, h)
print "\n\nAdams          Runge "
for i in range(len(adams[1])):
    print "%.8f " % adams[1][i], runge[1][i]
print "\nAdams table"
print_table(adams[2])