__author__ = 'vlad'
#region lab3

def frange(start, stop, step):
    while start < stop:
        yield start
        start += step

def intervals(start, stop, n_steps):
    step = (stop - start)/n_steps
    while start <= stop:
        yield start
        start += step

def romberg_brute(e=0, a=LEFT, b=RIGHT, y=f_6):
    h = b-a
    It = lambda h: h/2 * (y(a)+ y(b))
    I = It
    # i = 1
    # while True:
    for i in range(1,10):
        Ihi_1 = It(h)
        h /=2
        Ihi = It(h)
        for k in range(1,i):

            Rhi = (Ihi - Ihi_1)/(2**(2*k)-1)

            Ihi = Ihi + Rhi

            # def Rk_1(hi):
            #     return I(h) - I(h/2)/(2**(2*k) - 1)
            #
            # def Ik(hi):
            #     return Rk_1(hi) - I(hi)
            #
            # I = Ik

            # if fabs(R(h)) <= e:
            #     return I(h)
        # i += 1

    return I(h)

def romberg_er(y, a=LEFT, b=RIGHT):
    def R(n,m):
        if n == m == 0:
            return 0.5 * (b - a) * (y(a) + y(b))
        elif m == 0:
            hn = ((b - a)/(2**n))
            return 0.5 * R(n - 1, 0) + hn * sum((y(a+(2*k - 1) * hn) for k in range(1, 2**(n-1))))
        else:
            return R(n, m-1) + (1/(4**m - 1)) * (R(n, m-1) - R(n-1, m-1))
    n = 21#ceil(log2((b-a)) + 14)
    print(n)
    return R(n, 1)

def romberg2(func, a=LEFT, b=RIGHT):
    def I(p, h):
        Bs = [0.5, 1/3, 3/8, 2/45, 5/228]
        cofs = [ [1,1],
                 [1,4,1],
                 [1,3,3,1],
                 [7,32,12,32,7],
                 [19,75,50,50,75,19]
        ]
        k = p - 1
        return Bs[k] * h * sum( a * func(x) for a, x in zip(a[k], frange(a,b,h)) )

    hs = {0:  b - a}
    tabl = {(0,0): 0.5 * (b - a) * (func(a) + func(b))}
    i=1
    while True:
        hs[i] = hs[i-1]/2
        xpairs = zip(frange(a,b,hs[i]), frange(a+hs[i],b,hs[i]))
        tabl[(i,0)] = 0.5 * hs[i] * sum( func(x0) + func(x1) for x0, x1 in xpairs )
        for k in range(1, i+1):
            R = 1/(2**(2*k)) * (tabl[(i,k-1)] - tabl[(i-1,k-1)])
            tabl[(i,k)] = tabl[(i,k-1)] + R
            # print('r =', R)
            if abs(R) < EPSILON:
                print('n =', (b-a)/hs[i])
                return tabl[(i,k)]
        i+=1



#endregion