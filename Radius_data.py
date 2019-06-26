#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:16:46 2019

@author: kaliroe
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
import time
from Radius_functions import *


def eggleton_r(q):
    eggleton_r = (0.49*q**(2/3))/(0.6*q**(2/3)+np.log(1+q**(1/3)))
    return eggleton_r

#"adaptive" stratified montecarlo
def volb(q,xmin, xmax, ymin, ymax, zmin, zmax, level, max_level,R_rl=False):
    #print("entering level", level, xmin, xmax, ymin, ymax)
    num_inside = 0
    L2 = find_L2(q)
    L3 = find_L3(q)

    if q <= 1:
        L = L2
    else:
        L = L3
        
    if R_rl:
        L = L1
    
    p = Phi_z(q,L,0,0)
        
    full_vol = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
    coords = np.array([[xmin,ymin,zmin],[xmin,ymax,zmin],[xmax,ymin,zmin],\
                       [xmax,ymax,zmin],[xmin,ymin,zmax],[xmin,ymax,zmax],[xmax,ymin,zmax],\
                       [xmax,ymax,zmax],[0.5*(xmax+xmin),0.5*(ymax+ymin),0.5*(zmax+zmin)]])
    for coord2 in coords:
        
        x_act = coord2[0]
        y_act = coord2[1]
        z_act = coord2[2]
                
        x_vector = np.array([x_act,y_act,z_act])
        
        x_g = dPhi_x(q,x_act,y_act,z_act)
        y_g = dPhi_y(q,x_act,y_act,z_act)
        z_g = dPhi_z(q,x_act,y_act,z_act)
        
        gradient = (x_g, y_g, z_g)
        
        product = np.dot(x_vector,gradient)
        
        value = Phi_z(q,x_act,y_act,z_act)
        if value <= p and product >= 0:
        #if value <= p and product > 0:
            num_inside += 1
            plot2_x.append(coord2[0])
            plot2_y.append(coord2[1])
            plot2_z.append(coord2[2])
    if num_inside == 9:
        #print("full level",level)
        return (full_vol,9)
    elif level == max_level:
        #print("bottom level",level)
        #do a random sample
        x = xmin + (xmax-xmin)*np.random.rand()
        y = ymin + (ymax-ymin)*np.random.rand()
        z = zmin + (zmax-zmin)*np.random.rand()
        
        x_vector = np.array([x,y,z])
        
        x_g = dPhi_x(q,x,y,z)
        y_g = dPhi_y(q,x,y,z)
        z_g = dPhi_z(q,x,y,z)
        
        gradient = (x_g, y_g, z_g)
        
        product = np.dot(x_vector,gradient)
        
        if Phi_z(q,x,y,z) <= p and product >= 0:
        #if Phi_z(q,x,y,z) <= p and product > 0:
            #print("return full")
            #num_inside += 1
            return (full_vol,10)
        else:
            #print("return none",x,y)
            return(0,9)
    elif num_inside == 0:
        #print("empty level",level)
        return(0,8)
    else:
        #print("need extra level", num_inside, level)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)
        if R_rl:
            area1 = volb(q,xmin,xmid,ymin,ymid,zmin,zmid,level+1,max_level,R_rl=True)
            area2 = volb(q,xmid,xmax,ymin,ymid,zmin,zmid,level+1,max_level,R_rl=True)
            area3 = volb(q,xmin,xmid,ymid,ymax,zmin,zmid,level+1,max_level,R_rl=True)
            area4 = volb(q,xmid,xmax,ymid,ymax,zmin,zmid,level+1,max_level,R_rl=True)
            area5 = volb(q,xmin,xmid,ymin,ymid,zmid,zmax,level+1,max_level,R_rl=True)
            area6 = volb(q,xmid,xmax,ymin,ymid,zmid,zmax,level+1,max_level,R_rl=True)
            area7 = volb(q,xmin,xmid,ymid,ymax,zmid,zmax,level+1,max_level,R_rl=True)
            area8 = volb(q,xmid,xmax,ymid,ymax,zmid,zmax,level+1,max_level,R_rl=True)
        else:
            area1 = volb(q,xmin,xmid,ymin,ymid,zmin,zmid,level+1,max_level)
            area2 = volb(q,xmid,xmax,ymin,ymid,zmin,zmid,level+1,max_level)
            area3 = volb(q,xmin,xmid,ymid,ymax,zmin,zmid,level+1,max_level)
            area4 = volb(q,xmid,xmax,ymid,ymax,zmin,zmid,level+1,max_level)
            area5 = volb(q,xmin,xmid,ymin,ymid,zmid,zmax,level+1,max_level)
            area6 = volb(q,xmid,xmax,ymin,ymid,zmid,zmax,level+1,max_level)
            area7 = volb(q,xmin,xmid,ymid,ymax,zmid,zmax,level+1,max_level)
            area8 = volb(q,xmid,xmax,ymid,ymax,zmid,zmax,level+1,max_level)
        A = area1[0] + area2[0] + area3[0] + area4[0] + area5[0] + area6[0] +\
            area7[0] + area8[0]
        full_pts = area1[-1] + area2[-1] + area3[-1] + area4[-1] + area5[-1] + \
            area6[-1] + area7[-1] + area8[-1]
            
        R_v1 = np.abs(np.cbrt(3*A/(4*np.pi)))
        #print(R_v1)
        #print("returning area")
        return(A,R_v1,plot2_x,plot2_y,plot2_z,full_pts)
    
start_time = time.time()

radius_array = []
rl_radius_array = []
rl_alt_array = []
q_array = np.logspace(-10,10,20)

for q in q_array:
    L1 = find_L1(q)[0]
    L2 = find_L2(q)
    L3 = find_L3(q)
    
    plot2_x = []
    plot2_y = []
    plot2_z = []

    if q < 1:
        L = L2
    else:
        L = L3
    print("q = %f"%q)
    
    dist = np.abs(L1-L)/2
    A,r, pltx, plty, pltz, full = volb(q,L1,L,-dist, dist, -dist, dist, 1, 7)
    rl = eggleton_r(q)
    rl_alt = volb(q,L1,L,-dist, dist, -dist, dist, 1, 7, R_rl=True)[1]
    rl_radius_array.append(rl)
    radius_array.append(r)
    rl_alt_array.append(rl_alt)
    
    #plt.figure()
    #ax = plt.axes(projection='3d')

    #ax.plot3D(plot2_x, plot2_y, plot2_z, ".")

    
radius_array = np.array(radius_array)
rl_radius_array = np.array(rl_radius_array)
rl_alt_array = np.array(rl_alt_array)

dictionary = {"q value": q_array, "L radius value": radius_array, "Roche Lobe": rl_radius_array,"Alt Roche Lobe": rl_alt_array}

#part1 = pd.DataFrame(dictionary)

#export_csv1 = part1.to_csv (r'Outer_radius_dataframe.csv', index = None, header=True)


end = time.time()
print("total minutes:", (end - start_time)/60)