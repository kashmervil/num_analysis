#! /usr/bin/env python

#usage: pass parameters a and k to the command line

import math
import sys

def F(k,x,y): return math.exp(k*x + y) - x*y -1.4

def G(a,x,y): return (x/a)**2.0 + y**2.0 - 4

def dFx(k,x,y): return k*math.exp(k*x + y) - y

def dFy(k,x,y): return math.exp(k*x + y) - x

def dGx(a,x): return 2.0*x/(a**2.0)

def dGy(y): return 2.0*y

def dx(a,k,x,y): return -F(k,x,y)*dGy(y) + G(a,x,y)*dFy(k,x,y)

def dy(a,k,x,y): return -dFx(k,x,y)*G(a,x,y) + dGx(a,x)*F(k,x,y)

def d(a,k,x,y): return dFx(k,x,y)*dGy(y) - dGx(a,x)*dFy(k,x,y)

def X(a,k,x,y): return dx(a,k,x,y)/d(a,k,x,y)

def Y(a,k,x,y): return dy(a,k,x,y)/d(a,k,x,y)

NUM_EPS = 1.0e-9
NUM_SPLITTER = 1000.0

try:   
    a = float(sys.argv[1])
    k = float(sys.argv[2])
except:
    print "    Wrong command line parameters"
    raise

def sgn(x):
    return math.copysign(1,x)

def splitter():
    change = []
    xs = range(0,int(NUM_SPLITTER*a)*2)
    pairs  = map(lambda x: (x/NUM_SPLITTER, 4 - (x/(a*NUM_SPLITTER))**2),xs)
    pairs = filter(lambda x: ( 0 <= x[1]),pairs)
    pairs = map(lambda x: (x[0],math.sqrt(x[1])),pairs)
    pairs2 = map(lambda x: (-x[0],x[1]),pairs)
    pairs2.reverse()
    pairs = pairs2 + pairs
    pairs2 = map(lambda x: (x[0],-x[1]),pairs)
    pairs2.reverse()
    pairs = pairs + pairs2

    last = pairs[0]
    for current in pairs:
        if sgn(F(k,last[0],last[1]))*sgn(F(k,current[0],current[1])) < 0:
            change.append((last,current))
        last = current
    return change


def newton():
    segments = splitter()
    segments = map(lambda x: [x[0][0],x[0][1]],segments)
    for i in range(len(segments)):
        root = segments[i]
        counter = 0
        while (math.fabs(F(k,root[0],root[1])) > NUM_EPS):
            counter += 1
            print ' STEP ',counter ,'  X=',root[0],
            '  Y=',root[1],'  F(X,Y)=',F(k,root[0],root[1])
           
            root[0] += X(a,k,root[0],root[1])
            root[1] += Y(a,k,root[0],root[1])
        print '\n'

print """

*****    Numerical Analysis, Homework 2 
*****    Task 1, Newton's method for finding approximations to the roots of nonlinear systems
"""
print ''' Alternating-sign segments ''',splitter(), '\n'

newton()
