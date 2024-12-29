import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as spi
import time as ti
import scipy.special as sps
mob = pd.read_csv('mobius.csv')
zer = pd.read_csv('zeroes.csv')
dfm = pd.DataFrame(mob)
dfz = pd.DataFrame(zer)
#print(df.iloc[2]) col value = iloc[col - 2]
q = 0.05
start = 2.0
end = 10.0
t = np.arange(start, end, q)
c = 20
eSumTerms = 10
eLowerBound = -1000
#def T(x):
    #return np.real(sum(((m(n-1)/n)*spi.quad(lambda x: (np.e**(0.5 + z(n-1))/(0.5 + z(n-1))), -np.inf + z(n-1)*1j*np.log(x)/n, (0.5 + z(n-1)*1j)*np.log(x)/n) for n in (1, 3))))
def T(x):
    return Li(x) - sumg(x)
def sumf(x):
    output = 0
    for n in range(1, c+1):
        output += f(x,z(n))[0]
    return 2*np.sqrt(x)*output/(np.log(x))
def f(x,b):
    return (1/np.sqrt(0.25 + b**2))*np.cos(b*np.log(x) - np.arctan(2*b))
def sumg(x):
    output2 = 0
    for n in range(1, c+1):
        output2 += g(x,z(n))
    return 2*output2
def g(x,b):
    output3 = 0.0
    for k in range(1, 10):
        output3 += (mobius(k)/k)*(np.real(U(x,b,k)[0]))
    return output3
def z(n):
    return dfz.iloc[n-1]
def m(n):
    return dfm.iloc[n-1]
def integ(x,b,n):
    return spi.quad(lambda z: (np.e**(z + (b*1j*np.log(x)/n))/(z + (b*1j*np.log(x)/n))), -np.inf, np.log(x)/(2*n), complex_func=True)
def integ2(z):
    output = 0
    for n in range(0, eSumTerms):
        output += sps.factorial(n)/(np.pow(z,(n+1)))
    return np.e**(z)*output
def evalInteg(x,b,n):
    uBound = (0.5 + b*1j)*np.log(x)/n 
    lBound = -np.inf + b*1j*np.log(x)/n
    return integ2(uBound)
def Li(x):
    return spi.quad(lambda k: (1/np.log(k)), 2, x)[0]
def mobius(n):
    output = 0+0j
    arr = coprime(n)
    arr = np.pow(np.e, 2j*np.pi*arr/n)
    return np.round(np.real(np.sum(arr)))
def coprime(n):
    output = np.arange(1,n+1,1)
    return output[np.gcd(output,n) == 1]
Y = np.vectorize(T)
U = np.vectorize(integ)
V = np.vectorize(mobius)
#print(z(2)*1j*np.log(10)/2)
#print(spi.quad(lambda x: (np.e**(0.5 + z(2))/(0.5 + z(2))), z(2)*1j*np.log(10)/2, (0.5 + z(2)*1j)*np.log(10)/2))
#print("n: " + str(n) for n in range(1,c+1))
#print("Joi: " + str(T(4)))
#print(np.real(U([9,10],0.5,14.14,2)[0][0]))
start = ti.time_ns()
s = Y(t)
print(ti.time_ns() - start)
plt.plot(t,s)
plt.ylabel('some numbers')
plt.grid(True)
plt.show()