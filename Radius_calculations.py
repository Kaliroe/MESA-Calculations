#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:44:28 2019

@author: kaliroe
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from Radius_functions import *


q = 4

L1 = find_L1(q)[0]
L2 = find_L2(q)
L3 = find_L3(q)
L = find_outer_limit(q)[0]
###############################################################################
############################## Contour Map ####################################
###############################################################################


x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
x_cm = (1/(1+q))
Z = Phi_z(q,X,Y,0)

fig, ax = plt.subplots()#(figsize=(15,15))

ax.tick_params(axis='both', labelsize=15)

plt.contour(X, Y, Z, 250, cmap = 'gray')#cmap='Wistia')# cmap='viridis')

#sctr = plt.contour(X, Y, Z, 250, cmap = 'Pastel1')

#plt.colorbar(sctr, ax=ax, format='%d')
#make this just say high potential energy to low potenital energy.

CS = plt.contour(X, Y, Z, 100, levels = (Phi_z(q,L1,0,0),Phi_z(q,L2,0,0), Phi_z(q,L3,0,0)),cmap ="Dark2")#cmap ="tab10")# cmap ='Set1')#cmap='Wistia')# "Oranges")

fmt = {}
strs = ['L1', 'L2', 'L3']
for l, s in zip(CS.levels, strs):
    fmt[l] = s
    
plt.clabel(CS,CS.levels[:],inline=1, fmt=fmt, fontsize=18)
#plt.title("Equipotential surfaces of the Roche Potential")
plt.xlabel("x (units of separation)",fontsize=20)
plt.ylabel("y (units of separation)",fontsize=20)
plt.tight_layout()
plt.savefig("contour.pdf")

###############################################################################
########################### Radius Calculations ###############################
###############################################################################

df = pd.read_csv('Outer_radius_dataframe.csv', dtype=np.float64, sep=',', header=0)

q_array = df.values[:,0]
radius_array = df.values[:,1]
rl_radius_array = df.values[:,2]
rl_alt_array = df.values[:,3]


rl_error = np.abs(rl_radius_array-rl_alt_array)/rl_radius_array

#plt.figure()
#plt.plot(q_array,rl_radius_array, ".",label = "eggleton")
#plt.plot(q_array,rl_alt_array, ".",label = "data")
#plt.xlabel("q value")
#plt.ylabel("Roche Lobe radius")
#plt.title("level = 8")
#plt.xscale("log")
#plt.legend()

#plt.figure()
#plt.plot(q_array,rl_error, ".")
#plt.xlabel("q value")
#plt.ylabel("Roche Lobe Error")
#plt.title("level = 8")
#plt.xscale("log")




def cofunc(x,c,d,e):
    a = 0.8149/2
    b = 0.8149/np.pi
    return  a + b*np.arctan(c*np.log10(x)+d + e*np.power(np.log10(x),3))

#def cofunc(x,a,b,c,d):
 #   return  a + b*np.arctan(c*np.log10(x)+d)

coef1, other1 = curve_fit(cofunc,q_array, radius_array)

y1 = cofunc(q_array,coef1[0],coef1[1],coef1[2])#,coef1[3])

#plt.plot(q_array,y1, ".",label = "fit")
#plt.xscale("log")
#plt.legend()
#plt.title("Outer Lagrangian Radius")
#plt.xlabel("q value")
#plt.ylabel("Radius in units of separation")

error = np.abs(radius_array-y1)/radius_array

#plt.figure()
#plt.plot(q_array,error,".")
#plt.xscale("log")
#plt.title("Error")
#plt.xlabel("q value")
#plt.ylabel("Error")

r_rl = radius_array/rl_alt_array

n = len(q_array) 
mean = sum(np.log10(q_array)*r_rl)/n 
sigma = np.sqrt(sum((r_rl*np.log10(q_array)-mean)**2)/n)/2

a = 1/(np.sqrt(2*np.pi)*sigma)

def gaus(x,a,b,c):
    n = len(x)                          #the number of data
    x = np.log10(x)
    data = x*r_rl
    mean = sum(data)/n
    print(mean)
    #mean = .25
    sigma = np.sqrt(sum(np.power((data-mean),2))/n)#/2
    print(sigma)
    #sigma = 1.38
    return a*np.exp(-np.power((c*(data-mean)),2)/(2*sigma**2)) + b

popt,pcov = curve_fit(gaus,q_array,r_rl)

y_gaus = gaus(q_array,popt[0],popt[1],popt[2])

g_error = np.abs(r_rl-y_gaus)/r_rl

#plt.figure()
#plt.title("radius / roche lobe radius error")
#plt.xlabel("q value")
#plt.ylabel("error")
#plt.plot(np.log10(q_array),g_error)

def func(x):
    a = 2.7412
    b = 0.4426
    c = 21.4669
    d = -0.478
    e = 0.8894
    f = 7.1311
    g = 12.1788
    sigma = c/(g+np.exp(d*x))
    arg = np.power((x-b)/sigma,2)
    func = a/(1+arg)/(f+np.exp(e*x))+1
    return func

y_pablo = func(np.log10(q_array))

p_error = np.abs(r_rl-y_pablo)/r_rl
"""
plt.figure()
plt.title("radius / roche lobe radius error")
plt.xlabel("q value")
plt.ylabel("error")
plt.plot(np.log10(q_array),p_error)
"""
fig, ax = plt.subplots()#(figsize=(15,15))

ax.tick_params(axis='both', labelsize=20)
#plt.figure()
#plt.xscale("log")
#plt.title("radius / roche lobe radius", fontsize=20)
plt.xlabel("log q", fontsize=20)
plt.ylabel(r'$R / R_{rl}$', fontsize=20)
plt.plot(np.log10(q_array),r_rl, ".", ms=8,label = "data")
plt.plot(np.log10(q_array),y_pablo, lw=2, label = "fit")
plt.legend(fontsize=18)
plt.tight_layout()
plt.savefig("rl.pdf")

###############################################################################
####################### Optically thin calculations ###########################
###############################################################################

def Jackson(q,rv):
    Phi = -1/(1+q)*(1 + (1/2)*(1/(1+q))) - 1/(rv*(1+q))*(q + \
             (1/3)*(q+1)*np.power(rv,3) + (4/45)*(q*np.power((1+1/q),2) + \
              9/q + 3*(1 + 1/q))*np.power(rv,6))
    return Phi

Roche_array = []
Outer_array = []
Difference = []

for i in range(0,len(q_array)):
    
    q = q_array[i]
    
    L1 = find_L1(q)[0]

    if q <= 1:
        L = find_L2(q)
    else:
        L = find_L3(q)
        
    Phi_L1 = Phi_z(q,L1,0,0)
    Phi_outer = Phi_z(q,L,0,0)
    
    Outer_Jackson = Jackson(q_array[i],radius_array[i])
    Roche_Jackson = Jackson(q_array[i],rl_alt_array[i])
    
    Roche = Roche_Jackson/Phi_L1
    Outer = Outer_Jackson/Phi_outer
    
    diff = (Outer_Jackson - Roche_Jackson)/(Phi_outer - Phi_L1)
    
    Roche_array.append(Roche)
    Outer_array.append(Outer)
    Difference.append(diff)
    

fig, ax = plt.subplots()#(figsize=(15,15))
ax.tick_params(axis='both', labelsize=15)
#plt.figure()
#plt.xscale("log")
#plt.title("Phi_jackson/Phi_real", fontsize=20)
plt.xlabel("log q", fontsize=20)
plt.ylabel(r'$\Phi_{jackson}/\Phi_{real}$', fontsize=20)
plt.plot(np.log10(q_array),Roche_array, ".", label = "Roche")
plt.plot(np.log10(q_array),Outer_array, ".", label = "Outer")
plt.legend(fontsize=18)
plt.tight_layout()
plt.savefig("JackvsReal.pdf")

fig, ax = plt.subplots()#(figsize=(15,15))
ax.tick_params(axis='both', labelsize=15)
#plt.figure()
#plt.xscale("log")
#plt.title("Phi jackson Phi real Difference", fontsize=20)
plt.xlabel("log q", fontsize=20)
plt.ylim(0.9,1.2)
plt.xlim(-2,2)
plt.ylabel(r'$\Delta \Phi_{jackson} \Phi_{real}$', fontsize=20)
plt.plot(np.log10(q_array),Difference, ".")
plt.tight_layout()
plt.savefig("JackDifference.pdf")







