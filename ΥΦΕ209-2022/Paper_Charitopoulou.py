# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:57:27 2022

@author: des_p
"""
# Paper Quantum Information
# Charitopoulou Despoina 4405

import numpy as np
import warnings
from scipy.integrate import solve_ivp
from scipy.integrate import simpson
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Lane-Emden equation

def L_E(t, z, gamma):
    theta, x = z
    return [x,-(2/t)*x-(theta**(1/(gamma-1)))]

# Parameters

g = np.linspace(1.2, 1.8, 100)

tmin = 0.0001
tmax = 20

S_ce = np.zeros(len(g)) # Configurational entropy 
S3_ce = np.zeros((3,len(g))) # Sa^3 for each gamma value and  0.95π/R, 1.00π/R, 1.05π/R
M_mass = np.zeros(len(g)) # Mass 

con = np.array([1,0])
g_le = np.array([1.2,1.4,1.7])
kappa = np.array([0.95,1.00,1.05])
i = 10000   

for j in g_le: 
    
    sol = solve_ivp(L_E, [tmin, tmax], con ,t_eval=np.linspace(tmin,tmax,i),method='RK45',args=(j,),atol=1e-8,rtol=1e-8)
    k_min = np.pi/sol.t[-1]
    k = np.linspace(k_min,100*k_min,i)                                                                        
    h_kmin = (simpson(((sol.y[0])**(1/(j-1)))*np.sin(k_min*sol.t)*sol.t,sol.t,dx=0.001)*(1/k_min))**2                  
    h_k = np.zeros(i) 
    
    for i1,j1 in enumerate(k):
            
        h_k[i1] = (simpson(((sol.y[0])**(1/(j-1)))*np.sin(j1*sol.t)*sol.t,sol.t,dx=0.001)*(1/j1))**2 
            
    f1 = h_k/h_kmin  
   
    plt.figure(1)
    plt.plot(k/np.sqrt(j/(j-1)),f1,label="γ="+str(j)) # Fig.1
    

plt.xlim([0,1.5])
plt.ylim([0,1.1])
plt.ylabel(r"$f(|k|)$")
plt.xlabel(r"$k/\sqrt{4\pi G/K}$")
plt.grid()
plt.legend()

for i2,j2 in enumerate(g):
    
    sol = solve_ivp(L_E, [tmin, tmax], con ,t_eval=np.linspace(tmin,tmax,i),method='RK45',args=(j2,),atol=1e-8,rtol=1e-8)
    k_min = (np.pi/sol.t[-1])/kappa
    
    for i3,j3 in enumerate(k_min):
        
        k = np.linspace(j3,100*j3,i)
        h_kmin = (simpson(((sol.y[0])**(1/(j2-1)))*np.sin(j3*sol.t)*sol.t,sol.t,dx=0.001)*(1/j3))**2 
        h_k = np.zeros(i)
        
        for i4,j4 in enumerate(k):
            
            h_k[i4] = (simpson(((sol.y[0])**(1/(j2-1)))*np.sin(j4*sol.t)*sol.t,sol.t,dx=0.001)*(1/j4))**2
            
        f2 = h_k/h_kmin
        
        M_mass[i2] = ((4*np.pi*(j2/(j2-1))**(3/2))*simpson(((sol.y[0])**(1/(j2-1)))*(sol.t)**2,sol.t,dx=0.001))/200
        
        if i3==1:
            S_ce[i2] = -4*np.pi*((j2/(j2-1)))**(-3/2)*simpson(f2*np.log(f2)*(k**2),k,dx=0.001)
            
        S3_ce[i3][i2] = -4*np.pi*simpson(f2*np.log(f2)*(k**2),k,dx=0.001)
        
        
# Fig.2        
plt.figure(2)
plt.plot(g,S_ce,color='r',label=r'$S\rho_0^{-1}/\left( \left( \frac{K}{4\pi G}\right)^{-\frac{3}{2}} \rho_c^{2-\frac{3}{2}\gamma} \right)$')
plt.plot(g,M_mass,'--',color='g',label=r'$M/\left( 200\left( \frac{K}{4\pi G}\right)^{\frac{3}{2}} \rho_c^{\frac{3}{2}\gamma -2} \right)$')
plt.xlim([1.25,1.7])
plt.ylim([0.4,1.3])
plt.legend()
plt.xlabel("γ")
plt.grid()   

# Fig.5
plt.figure(3)
plt.plot(g,S3_ce[0],'--',color='b',label='π/(0.95R)')
plt.plot(g,S3_ce[1],color='black',label='π/(1.00R)')
plt.plot(g,S3_ce[2],'-.',color='r',label='π/(1.05R)')
plt.axvline(4/3,color='grey',ls='--')
plt.text(1.35, 4.64, "4/3", rotation=0)
plt.axvline(5/3,color='grey',ls='--')
plt.text(1.63, 4.64, "5/3", rotation=0)
plt.xlim([1.25,1.75])
plt.ylim([4.6,5.6])
plt.xlabel("γ")
plt.ylabel(r"$S\alpha^{3}$")
plt.legend()
plt.grid()
            
            
plt.show()