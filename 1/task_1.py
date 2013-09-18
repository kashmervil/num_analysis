#! /usr/bin/env python

import sys
import math

NUM_EPS = 1.0e-10
NUM_SPLITER = 1.0e-1

def polynom_neg(xs):
    r = []
    for i in range(len(xs)):
        if i % 2:
            r.append(-xs[i])
        else:
            r.append(xs[i])
    return r

def polynom_reverse(x):
    r = []
    n = len(x) - 1
    w = cmp(x[n],0)
    for i in range(len(x)):
        r.append(w*x[n-i])
    return r
     
def polynom_up_limit(x):
    first = filter(lambda x: x < 0,x)
    return 1 + (math.fabs(min(x))*1.0/x[0])**(1.0/x.index(first[0]))



def polynom_limits(x):
    negat= polynom_neg(x)
    negrev = polynom_reverse(x)
    rev = polynom_neg(negrev)
    return [-polynom_up_limit(negat),-1.0/polynom_up_limit(rev),
           1.0/polynom_up_limit(negrev), polynom_up_limit(x)]

def complex_limits(_):
    return [-math.pi,-math.pi/2.0,-0.2,0.2]

def polynom_value(p,x0):
    r = 0
    for i in range(len(p)):
        r = r*x0 + p[i]
    return r

def complex_value(_,x):
        return x*x + 4*math.sin(x)

def polynom_derivation(p,x):
    n = len(p) - 1
    d = 0
    for i in range(n):
        d = d*x + p[i]*(n - i)
    return d

def complex_derivation(_,x):
    return 2*x + 4*math.cos(x)



def sgn(x):
    return math.copysign(1,x) 

def splitter(x,get_value,start,end):
    counter = start
    change = []
    last = 0
    current = get_value(x,start)
    while counter < end:
        counter += NUM_SPLITER
        current, last = get_value(x,counter),current
	if sgn(last)*sgn(current) <= NUM_EPS:
            change.append(counter-NUM_SPLITER)
    return change

def root_segments(x,get_value,limits):
    lim = limits(x)
    print 'Roots bounds',(lim[0],lim[1]),(lim[2],lim[3]),'\n' 
    return splitter(x,get_value,lim[0],lim[1]) + splitter(x,get_value,lim[2],lim[3]) 


def newton(x,get_value,limits,derivation):
    print "Newton's method for finding approximations to the roots/n"
    print 'of ',x 
    segments = root_segments(x,get_value,limits)
    for i in range(len(segments)):
        root = segments[i]
        counter = 0
        while (math.fabs(get_value(x,root)) > NUM_EPS):
            counter += 1
            print ' STEP ',counter ,'  X=',root,'  F(X)=',get_value(x,root)
            root -=get_value(x,root)/derivation(x,root)
        print '\n' 

def polynomial_newton(x): return newton(x,polynom_value,polynom_limits,polynom_derivation)

def complex_newton(): return newton('funtion f(x) = x^2 +4*sin(x)',complex_value,complex_limits,complex_derivation)
y = [1,0.75,-12.25,3]
x = [2048,0,-6144,0,6912,0,-3584,0,840,0,-72,0,1]

print """
*****  Numerical Analysis, Homework 1 
*****  Tasks 1-3, Newton's method for finding approximations to the roots of functions
*****  1) 2048*x^12 - 6144*x^10 + 6912*x^8 - 3584*x^6 + 840*x^4 -72*x^2 + 1
*****  2) x^3 + 0.75*x^2 + 3
*****  3) x^2 + 4*sin(x) 
"""
        
polynomial_newton(y)

NUM_SPLITER = 1.0e-2

polynomial_newton(x)

complex_newton()




