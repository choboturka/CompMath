__author__ = 'vlad'
import numpy as np
import numpy.linalg as LA
from math import sqrt, fabs, copysign, log10
from functools import reduce

# region CONSTANTS

E = 10**-3
M2 = np.mat(np.loadtxt('table2.txt'))
M3 = np.mat(np.loadtxt('table3.txt'))

# endregion

# region auxiliary functions
def max_nondiag(mat):
    n = mat.shape[0]
    l = mat.A1.tolist()
    nl = l.copy()
    del nl[0::n+1]
    ans = [(pos//n,pos % n) for pos,x in enumerate(l) if pos not in range(0,n**2+1,n+1) and x == max(nl, key=fabs)]
    return ans[0]
# endregion


#2 3 1 0 1 0


# PARTIAL
#todo 1. exponential (alg 2)

def exp_high(A, aitken=False, second=False):
    # init

    delta = E/10

    n = A.shape[0]
    y = np.mat(np.ones(n)).T   # y(0)

    eig_val = 10               # lambda(0)
    z = y / LA.norm(y)
    p = 3
    r = 1
    B = LA.matrix_power(A, 2**r)

    for k in range(1,10000):
        y = B * z


#todo 2. eitken (alg 2)

#3. reverse iterations (alg 3)
def rev_iter(A):
    #init
    eig_val = 10                        # lambda(0)
    y = np.mat(np.ones(A.shape[0])).T  # y(0) = [1,1,1,...].T
    # y += 5
    s = np.dot(y.A1, y.A1)
    p = sqrt(s)
    z = y / p
    ''':type : numpy.matrixlib.defmatrix.matrix'''

    #iterations

    for k in range(1,100):
        y = LA.solve(A, z)
        s = np.dot(y.A1, y.A1)
        t = np.dot(y.A1, z.A1)
        p = sqrt(s)

        # eig_val = s/t
        # nz = y/p
        # ''':type : numpy.matrixlib.defmatrix.matrix'''
        # stop = LA.norm(nz) - LA.norm(z) < E
        # if stop <= E:
        #     z = nz
        #     break
        # z = nz

        new_eig_val = s / t
        z = y / p
        if abs(1/new_eig_val - 1/eig_val) <= E:
            break
        eig_val = new_eig_val
        # print("{0} iteration, stop = {1}".format(k,stop))
        # print(1/eig_val, z.A1, 'iter {0}'.format(k))
    return 1/eig_val, z.A1

# FULL

#todo 4. dot product (alg 1)

def dot_product(M):
    pass
    


#5. yacobi (alg 0) todo fix formulas
def yacobi(A):
    B = A
    n = A.shape[0]
    CT = np.mat(np.eye(n,n))
    for k in range(1,1000):
        # NB = B

        # print('iter {0}'.format(k))

        i,j = max_nondiag(B)
        p = 2*B[i,j]
        q = B[i,i] - B[j,j]
        d = sqrt(p**2 + q**2)
        r = fabs(q)/(2*d)
        c = sqrt(0.5 + r)

        if p!=0 and q !=0 and log10(fabs(q)) - log10(fabs(p)) >=2:
            s = copysign(fabs(p)/(2*d*c),p*q)
        else:
            s = copysign(sqrt(0.5 - r), p*q)

        T = np.mat(np.eye(n,n))
        T[i,i] = T[j,j] = c
        T[i,j] , T[j,i] = -s, s
        CT = CT * T

        # NB[i,i] = c**2 * B[i,i] + s**2 * B[j,j] + 2*c*s*B[i,j]
        # NB[j,j] = s**2 * B[i,i] + c**2 * B[j,j] - 2*c*s*B[i,j]
        # NB[i,j] = NB[j,i] = 0
        # for l in range(0,n):
        #     if l != i and l != j:
        #         NB[i,l] = NB[l,i] = c*B[l,i] + s*B[l,j]
        #         NB[j,l] = NB[l,j] = -s*B[l,i] + c*B[l,j]
        B = T.T * B * T

        def stop():
            M = B.copy()
            np.fill_diagonal(M,0)
            if LA.norm(M, 'fro') <= E:
                return True
            else:
                return False

        if stop():
            break

    return np.diag(B), CT*-1


#6. lu (1 diagonal in u)
def lu(A):
    n=A.shape[0]
    for k in range(1,100):
        L = np.mat(np.zeros((n,n)))
        U = np.mat(np.eye(n))
        for index in np.ndindex(n,n):
            i,j = index[0], index[1]
            if i >= j:
                L[i,j] = A[i,j] - sum([L[i,k]*U[k,j] for k in range(j)])
            else:
                U[i,j] = (A[i,j] - sum([L[i,k]*U[k,j] for k in range(i)]))/L[i,i]

        A = U * L

        if LA.norm(np.triu(A,1)) < E:
            break

    return np.diag(A)


#7.  qr (alg 0)
def qr(A):
    n = A.shape[0]
    for k in range(1,1000):
        # print(k)
        q, r = LA.qr(A, mode='complete')
        NA = r * q

        stop = []
        for i in range(0,n):
            if NA[i,i] - A[i,i] < E:
                stop.append(True)
            else:
                stop.append(False)
        A = NA

        if False not in stop:
            break

    return np.diag(A)

