#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:13:40 2019

@author: kaliroe
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
import random as rd
from scipy.optimize import curve_fit
import time

G = 6.67428e-8 #cm^3⋅g−1⋅s−2
Rsun = 6.9598e10

###############################################################################
############################# Phi and Derviatives #############################
###############################################################################

def dPhi(q,x,y):
    x_cm = 1/(1+q)
    r = np.sqrt(x**2+y**2)
    dPhi_dx = q*x/((1+q)*r**3) -(1-x)/((1+q)*(((1-x)**2)+y**2)**(3/2)) -(x-x_cm)
    return dPhi_dx

def Phi(q,x,y):
    r = np.sqrt(x**2+y**2)
    x_cm = 1/(1+q)
    phi = -q/(r*(1+q)) -1/((q+1)*np.sqrt((1-x)**2+y**2)) -(1/2)*((x-x_cm)**2+y**2)
    return phi

def Phi_z(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    x_cm = 1/(1+q)
    phi = -q/(r*(1+q)) -1/((q+1)*np.sqrt((1-x)**2+y**2+z**2)) -(1/2)*((x-x_cm)**2+y**2)
    return phi

def dPhi_y(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    dPhi_dy = q*y/((1+q)*r**3) + y/((1+q)*(((1-x)**2)+y**2+z**2)**(3/2)) - y
    return dPhi_dy

def dPhi_x(q,x,y,z):
    x_cm = 1/(1+q)
    r = np.sqrt(x**2+y**2+z**2)
    dPhi_dx = q*x/((1+q)*r**3) -(1-x)/((1+q)*(((1-x)**2)+y**2+z**2)**(3/2)) - (x-x_cm)
    return dPhi_dx

def ddPhi_y(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    ddPhi_dy = q/((1+q)*r**3) - 3*q*y**2/((1+q)*r**5) + \
    1/((1+q)*(((1-x)**2)+y**2+z**2)**(3/2)) - \
    3*y**2/((1+q)*(((1-x)**2)+y**2+z**2)**(5/2)) - 1
    return ddPhi_dy

def dPhi_z(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    dPhi_dz = q*z/((1+q)*r**3) + z/((1+q)*(((1-x)**2)+y**2+z**2)**(3/2))
    return dPhi_dz

def ddPhi_z(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    ddPhi_dz = q/((1+q)*r**3) - 3*q*z**2/((1+q)*r**5) + \
    1/((1+q)*(((1-x)**2)+y**2+z**2)**(3/2)) - \
    3*z**2/((1+q)*(((1-x)**2)+y**2+z**2)**(5/2))
    return ddPhi_dz

def ddddPhi_y(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    ddddPhi_dy = -(1/(1+q))*(q*(3/r**5 - 15*y**2/r**7 + 6/r**5 \
    - 30*y**2/r**7 - 45*y**2/r**7 + 105*y**4/r**9) + 3/((((1-x)**2)+y**2+z**2)**(5/2))\
    - 15*y**2/((((1-x)**2)+y**2+z**2)**(7/2)) + 6/((((1-x)**2)+y**2+z**2)**(5/2))\
    - 30*y**2/((((1-x)**2)+y**2+z**2)**(7/2)) - 45*y**2/((((1-x)**2)+y**2+z**2)**(7/2))\
    + 105*y**4/((((1-x)**2)+y**2+z**2)**(9/2)))
    return ddddPhi_dy

def ddddPhi_z(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    ddddPhi_dz = -(1/(1+q))*(q*(3/r**5 - 15*z**2/r**7 + 6/r**5 \
    - 30*z**2/r**7 - 45*z**2/r**7 + 105*z**4/r**9) + 3/((((1-x)**2)+y**2+z**2)**(5/2))\
    - 15*z**2/((((1-x)**2)+y**2+z**2)**(7/2)) + 6/((((1-x)**2)+y**2+z**2)**(5/2))\
    - 30*z**2/((((1-x)**2)+y**2+z**2)**(7/2)) - 45*z**2/((((1-x)**2)+y**2+z**2)**(7/2))\
    + 105*z**4/((((1-x)**2)+y**2+z**2)**(9/2)))
    return ddddPhi_dz


def dPhi_dydz(q,x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    dydz = (1/(1+q))*(q*(-3/r**5 + 15*z**2/r**7 + 15*y**2/r**7 \
    - 75*z**2*y**2/r**7) - 3/((((1-x)**2)+y**2+z**2)**(5/2))\
    + 15*z**2/((((1-x)**2)+y**2+z**2)**(7/2)) + 15*y**2/((((1-x)**2)+y**2+z**2)**(7/2))\
    - 105*z**2*y**2/((((1-x)**2)+y**2+z**2)**(9/2)))
    return dydz

###############################################################################
############################### Critial Points ################################
###############################################################################

def find_L1(q):
    
    x_cm = (1/(1+q))
    
    upper_bound = 1
    lower_bound = 0
    
    limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
    
    tolerance = 1e-6
    
    while limit > tolerance:
        x = (lower_bound+upper_bound)/2
        
        dPhi_new = dPhi(q,x,0)
        
        if dPhi_new > 0:
            lower_bound = x
        elif dPhi_new < 0:
            upper_bound = x
        else:
            #print("no converge L1, q = ", q)
            break
        
        limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
        
    L1 = (upper_bound + lower_bound)/2
    
    dist = L1 - x_cm
    
    return L1,dist,lower_bound,upper_bound

def find_L2(q):
            
    if q < 1:
        upper_bound = 0
        lower_bound = -1
        limit = np.abs(upper_bound-lower_bound)/np.abs(lower_bound)
           
    if q >= 1:    
        upper_bound = 2
        lower_bound = 1
        limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
    
    
    tolerance = 1e-6
    
    while limit > tolerance:
        x = (lower_bound+upper_bound)/2
        
        dPhi_new = dPhi(q,x,0)
        
        if dPhi_new > 0:
            lower_bound = x
        elif dPhi_new < 0:
            upper_bound = x
        else:
            print("no converge L2, q = ", q)
            break
        
        if q < 1:
            limit = np.abs(upper_bound-lower_bound)/np.abs(lower_bound)
        else:
            limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
        
    L2 = (upper_bound + lower_bound)/2
    
    return L2


def find_L3(q):
                
    if q >= 1:
        upper_bound = 0
        lower_bound = -1
        limit = np.abs(upper_bound-lower_bound)/np.abs(lower_bound)
           
    if q < 1:    
        upper_bound = 2
        lower_bound = 1
        limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
    
    limit = np.abs(upper_bound-lower_bound)
    
    tolerance = 1e-6
    
    while limit > tolerance:
        x = (lower_bound+upper_bound)/2
        
        dPhi_new = dPhi(q,x,0)
        
        if dPhi_new > 0:
            lower_bound = x
        elif dPhi_new < 0:
            upper_bound = x
        else:
            print("no converge L3, q = ", q)
            break
        
        if q >= 1:
            limit = np.abs(upper_bound-lower_bound)/np.abs(lower_bound)
        else:
            limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
        
    L3 = (upper_bound + lower_bound)/2
    
    return L3


def find_outer_limit(q):
    
    upper_bound = 2
    lower_bound = 0
    
    limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
    L1 = find_L1(q)[0]
    
    tolerance = 1e-8
    
    while limit > tolerance:
        x = (lower_bound+upper_bound)/2
        
        dPhi_new = dPhi_y(q,L1,x,0)
        
        if dPhi_new > 0:
            lower_bound = x
        elif dPhi_new < 0:
            upper_bound = x
        else:
            print("no converge L, q = ", q)
            break
        
        limit = np.abs(upper_bound-lower_bound)/np.abs(upper_bound)
        
    L = (upper_bound + lower_bound)/2
    
    return L,upper_bound,lower_bound



