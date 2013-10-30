
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


def get_up_bound_index(xs, x0):
    tmp = filter(lambda x: x <= x0, xs)
    return len(tmp)


def get_down_bound_index(xs, x0):
    tmp = filter(lambda x: x <= x0, xs)
    return len(tmp) - 1


def beg_of_table(xs, x0):
    index = get_down_bound_index(xs[0], x0)
    tmp = 0.0
    t = (x0 - xs[0][index])*10.0
    k = 1.0
    for i in range(len(xs) - 1 - index):
        tmp += k * xs[i + 1][index]
        k *= (t - i)/(i + 1.0)
    return tmp


def end_of_table(xs, x0):
    index = get_up_bound_index(xs[0], x0)
    tmp = 0.0
    t = (x0 - xs[0][index])*10.0
    k = 1.0
    for i in range(index + 1):
        tmp += k * xs[i + 1][index - i]
        k *= (t + i)/(i + 1.0)
    return tmp


def mid_of_table(xs, x0):
    index = get_down_bound_index(xs[0], x0)
    tmp = 0.0
    t = (x0 - xs[0][index])*10.0
    k = 1.0
    for i in range(index + 1):
        tmp += k * xs[i + 1][index - i/2]
        k *= (t + (-1)**i * i)/(i + 1.0)
    return tmp


def reverse_interpolation(xs, y0):
    index = get_up_bound_index(xs[1], y0)
    t_cur = -0.5

    def foo(t):
        return (y0 + t * xs[2][index - 1] - end_of_table(xs, t/10.0 + xs[0][index]))/xs[2][index - 1]
    

    counter = 0
    while abs(foo(t_cur) - t_cur) > 1e-7:
        counter += 1
        print counter 
        t_cur = foo(t_cur)
    return t_cur/10.0 + table[0][index]

def derivation(xs,x0):
    index = get_up_bound_index(xs[1], y0)
    t_cur = -0.5

    def foo(t):
        return (y0 + t * xs[2][index - 1] - end_of_table(xs, t/10.0 + xs[0][index]))/xs[2][index - 1]
    return (foo(x0 + 1e-3) - foo(x0))/1e-3


def all_der(xs,x0):
    index = get_up_bound_index(xs[0], x0)
    mi = ma = derivation(xs,x0)
    t = -1.0
    while (t < 0.0):
        x0 = t/10.0 + xs[0][index]
        t += 1e-3
        cur = derivation(xs, x0)
        if ma < cur:
            ma = cur
        
        if mi > cur:
            mi = cur
    return (min, max)
       
table = [[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
         [1.623250, 1.664792, 1.701977, 1.734832, 1.763404, 1.787764, 1.808002, 1.824230, 1.836580, 1.845201]]
table = pp(table)

for w in range(len(table)):
        if w > 5:
            table[w] = map(lambda x: 0.0, table[w])

for j in range(len(table[0])):
    for w in range(len(table) - j):
        if w < 6:
            print ' ', '{0:.6f}'.format(table[w][j]),
    print ''


x1 = 0.172544
x2 = 0.815445
x3 = 0.269765
y = 1.739372
#y = 1.787764
print "Beginning of table"

print 'X3 = ', x3, '    T = ', (x3 - table[0][2]) * 10, 'P(X3)= ', beg_of_table(table, x3)
print '\nEnd of table'
print 'X2 = ', x3, '    T = ', (x3 - table[0][3]) * 10, 'P(X3)= ', end_of_table(table, x3)
print '\nMiddle of table'
print 'X3 = ', x3, '    T = ', (x3 - table[0][2]) * 10, 'P(X3)= ', mid_of_table(table, x3)

x4 = reverse_interpolation(table, y)
print 'P(X*) = ', y, 'T = ', x4 - table[0][get_up_bound_index(table[0], x4)], 'X* = ', x4

print "\n\n\n"
print "Final Check with End of Table algoritm", end_of_table(table, x4)
