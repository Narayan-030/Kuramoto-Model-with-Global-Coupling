import numpy as np
import matplotlib.pyplot as plt

#Write a code to integrate a second order ODE using RK4 method:

N = 100      #Number of oscillators
Initial_condition = np.random.uniform(0,2*np.pi,N).tolist()
Natural_freq = np.random.standard_cauchy(N).tolist()
h = 0.05
Number_of_iter = 30000
t = np.linspace(0,Number_of_iter*h,Number_of_iter).tolist()

Coupling_constant = np.linspace(0,5,100).tolist()


A = np.zeros((len(t),N))
A[0] = Initial_condition

def Kuramoto(y, t,w,K, r,psi):
    return  w + K*r*np.sin(psi - y)
#dummy = 0

Order_parameter = []
for K in Coupling_constant:
    r_temp = []
    for i in range(len(t)-1):
        r = (np.sqrt((np.sum(np.cos(A[i])))**2+(np.sum(np.sin(A[i])))**2))/N
        psi = np.arctan((np.sum(np.sin(A[i])))/(np.sum(np.cos(A[i]))))
        r_temp.append(r)
        for j in range(N):#for jth oscillator at ith time
            yi = A[i][j]
            T = t[i]
            wj = Natural_freq[j]                  
            k1 = Kuramoto(yi,T,wj,K,r,psi)
            k2 = Kuramoto(yi + k1 * h / 2., T + h / 2.,wj,K,r,psi)
            k3 = Kuramoto(yi + k2 * h / 2., T + h / 2.,wj,K,r,psi)
            k4 = Kuramoto(yi + k3 * h, T + h,wj,K,r,psi)
            yi1 = yi + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)
            A[i+1][j] = yi1
    Order_parameter.append(np.average(r_temp[-5000:]))
    #dummy = dummy+1
    #print(dummy*100/len(Coupling_constant))
