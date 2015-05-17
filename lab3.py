__author__ = 'vlad'
from math import *
from functools import reduce
import sys
import numpy
sys.setrecursionlimit(150)

# Var 6

f_6 = lambda x: x * cosh(x) - x**3 + x**2 * log(x)
d2f_6 = lambda x: -6*x + 2*log(x) + 2*sinh(x) + x*cosh(x) + 3

LEFT, RIGHT = 3, 5
EPS = 10**-9


def inner_x(start, stop, n_steps):
    step = (stop - start)/n_steps
    start+=step
    while start < stop:
        yield start
        start += step


def trapeze(n, a=LEFT, b=RIGHT, y=f_6, d2y=d2f_6):
    h = (b - a)/n
    x = inner_x(a, b, n)
    zeta = (b-a)/2
    r = -((b-a)/12) * h**2 * d2y(zeta)
    return h * ((y(a) + y(b))/2 + fsum(map(y, x)))


def trapeze2(e, a=LEFT, b=RIGHT, y=f_6, d2y=d2f_6):

    n = ceil(sqrt(((b-a)**3 * fabs(d2y(b)))/(12*e)))
    h = (b - a)/n
    x = inner_x(a, b, n)
    print(n)
    return h * ((y(a) + y(b))/2 + fsum(map(y, x)))


# todo: understand what is what
def romberg(e,f=f_6, a=LEFT, b=RIGHT, n=10):
    """Estimate the integral of f(x) from a to b using Romberg Integration.

    USAGE:
        r = romberg( f, a, b, n )

    INPUT:
        f       - function to integrate,
        [a, b]  - the interval of integration,
        n       - number of levels of recursion

    OUTPUT:
        numpy float array - Romberg integration array; most accurate
                            answer should be at bottom of right-most column,

    NOTES:
        Based on an algorithm in "Numerical Mathematics and Computing"
        4th Edition, by Cheney and Kincaid, Brooks-Cole, 1999.

    """

    r = numpy.array([[0] * (n+1)] * (n+1), float)
    h = b - a
    # hc = []
    # hc.append(h)

    r[0,0] = 0.5 * h * (f(a) + f(b))
    len = n

    powerOf2 = 1
    for i in range(1, n + 1):

        # Compute the halved stepsize and use this to sum the function at
        # all the new points (in between the points already computed)
        h = 0.5 * h

        sum = 0.0
        powerOf2 = 2 * powerOf2
        for k in range(1, powerOf2, 2):
            sum = sum + f(a + k * h)

        # Compute the composite trapezoid rule for the next level of
        # subdivision.  Use Richardson extrapolation to refine these values
        # into a more accurate form.

        r[i,0] = 0.5 * r[i-1,0] + sum * h

        powerOf4 = 1
        for j in range( 1, i + 1 ):
            powerOf4 = 4 * powerOf4
            r[i,j] = r[i,j-1] + ( r[i,j-1] - r[i-1,j-1] ) / ( powerOf4 - 1 )

        if fabs(r[i,i] - r[i,i-1]) <= e:
            len = i+1
            break

    # ra = r[0:len, 0:len]
    return r[0:len, 0:len], r[len-1,len-1]


if __name__ == '__main__':
    print(trapeze2(10**-9))
    # print(trapeze(1000000))
    a,v = romberg(10**-9)
    print(v)
    print(a)
