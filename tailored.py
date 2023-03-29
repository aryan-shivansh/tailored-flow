# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#Input parameters
T1 = 300    #K
y1 = 1.4  #gamma1
m1 = 28.9647e-3  #MW1
T4 = 300    #K
y4 = 1.667   #gamma4
m4 = 4.003e-3     #MW4
R = 8.314
delta_p = 27 # bars, delta_p = p4 - p1, pressure diff at which diaphragm raptures 

def f(M):
    #a = p52
    #b = T32
    # c = T21
    # d =T34
    # e = T32 calculated by T21 and T34, from first two eq of 3 
    a = ( (3*y1 - 1)*(M**2) - 2*(y1 - 1) )/( (y1 - 1)*(M**2) + 2)
    b = (m4*(y4 -1)*(1+ ((y4+1)/(y4-1))*a ))/(m1*(y1 -1)*(1+ ((y1+1)/(y1-1))*a ))
    c = ( ((y1-1)/(y1+1))*(M**2 - 1) + M**2 )*( ((y1-1)/(y1+1))*(M**2 - 1) + 1 )/(M**2)
    
    a1 = (y1*R*T1/m1)**0.5
    # print(a1)
    a4 = (y4*R*T4/m4)**0.5
    # print(a4)
    d = ( 1 - ((y4-1)/(y1+1))*(a1/a4)*(M - 1/M))**2
    e = d*T4/(c*T1)
    return(e-b)

# print(f(2))
def M_tailored(X):
    [T1,y1,m1,T4,y4,m4,R] = X
    A = 0.9
    B = 20
    j = 0
    while j<100:
        M = (A+B)/2
        #print(M)
        if f(M)>0:
            A = M
        else:
            B = M
        j = j+1
    return(M)
# n = 100
# temp = np.zeros(n)
# M_ = np.linspace(1.0005,10,n)
# for j in range(n):
#     temp[j] = f(M_[j])

# plt.plot(M_,temp)
X = [T1,y1,m1,T4,y4,m4,R]
# print(M_tailored(X))

M = M_tailored(X)
print("\nIncident shock mach no = ", M )


# T1 = 300
# T4 = 300
gamma1 = y1
# m1 = 39.948e-3
gamma4 = y4
Ms = M
# m4 = 16.043e-3
#From normal shock relations we get
a1 = np.sqrt(gamma1*R*T1/m1) #speed of sound for driven
a4 = np.sqrt(gamma4*R*T4/m4) # driver
c1 = 2*gamma1/(gamma1+1)
c2 = ((gamma4-1)/(gamma1+1))*(a1/a4)
c3 = 2*gamma4/(gamma4-1)
r = (1 + c1*(Ms**2 - 1))/(1 - c2*(Ms**2-1)/Ms)**c3 

p1 = delta_p/(r-1)
p4 = r*p1
# print(p4*1e-5)
# print(p1*1e-5)
print("ratio, p4/p1 = ", r)
print("p1 = ", p1, " bars")
print("p4 = ", p4, " bars")
