__author__ = 'vlad'
from math import *


# Var 14

# Global
E = 10**-5

dy = lambda x,y: (cos(x) - y**2)/y

NODES = 21
X0, Y0 = 1, 1
B = 2.0
H = (B - X0) / NODES

def average_node(f, x, y):
    for i in range(NODES):
        # print(i)
        k1 = f(x,y)
        k2 = f(x + H/2, y + k1*H/2)
        y += H*k2
        x += H
        yield x,y

def kutta_merson(f, x, y, h=H):
    for i in range(NODES):
        k1 = f(x,y)
        h = H #????????
        while 1:
            k2 = f(x + h/3, y + k1*h/3)
            k3 = f(x + h/3, y + k1*h/6 + k2*h/6)
            k4 = f(x + h/2, y + k1*h/8 + k2*3*h/8)
            k5 = f(x + h, y + k1*h/2 - k3*3*h/2 + 2*h*k4)
            y1t = y + h/2 * (k1 - 3*k3 + 4*k4)
            y1 = y + h/6 * (k1 + 4*k4 + k5)
            R = 0.2 * fabs(y1 - y1t)
            if R <= E:
                x += H #???????
                y = y1
                break
            else:
                h /= 2
        yield x,y


# methods Ser tochki, kutta mersona

if __name__ == "__main__":
    print(H)
    for e1,e2 in zip(average_node(dy,X0,Y0),kutta_merson(dy,X0,Y0)): print(e1,e2)



