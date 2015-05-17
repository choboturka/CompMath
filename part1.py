__author__ = 'vlad'
from methods import *
from numpy.linalg import eigvals
from functools import reduce

if __name__ == '__main__':

    ex = np.mat("4 2 1; 2 5 3; 1 3 6")

    a,v = LA.eig(M3)
    print("Correct eigenvalues \n",a)

    a, v = yacobi(M3)
    print("My yacobi \n",a)

    a = lu(M3)
    print("My LU \n",a)

    a = qr(M3)
    print("My QR \n",a)
