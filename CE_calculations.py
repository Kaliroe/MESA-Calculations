#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 14:36:39 2019

@author: Kaliroe Pappas, kmpappa2@illinois.edu

"""

import numpy as np
import pandas as pd
import mesa_data as md
import matplotlib.pyplot as plt
from constants import *
"""
from matplotlib import rcParams
#plt.rcParams('text', usetex=True)
rcParams['font.family'] = 'stixgeneral'
rcParams['font.size'] = 14
"""

"""The following is specific to my personal computer"""

#Run this script on terminal in the Downloads directory


data1 = pd.read_csv("mesa-dev/massive_bins_alpha01.dat", delim_whitespace=True, \
        header=1, names=['log10(M_1i)(Msun)', 'qratio(M_2i/M_1i)', 'log10(P_i)(days)',\
        'result', 'CE_flag', 'log_L_1', 'log_T_1', 'M_1f(Msun)', 'M_2f(Msun)', 'Porb_f(d)',\
        'tmerge(Gyr)', 'He_core_1(Msun)', 'C_core_1(Msun)', 'runtime(min)'])

data_0100 = data1[867:918]
data_0150 = data1[816:867]
data_0200 = data1[765:816]
data_0250 = data1[714:765]
data_0300 = data1[663:714]
data_0350 = data1[612:663]
data_0400 = data1[561:612]
data_0450 = data1[510:561]
data_0500 = data1[459:510]
data_0550 = data1[408:459]
data_0600 = data1[357:408]
data_0650 = data1[306:357]
data_0700 = data1[255:306]
data_0750 = data1[204:255]
data_0800 = data1[153:204]
data_0850 = data1[102:153]
data_0900 = data1[51:102]
data_0950 = data1[0:51]


mesa_data1 = pd.read_csv("mesa-dev/massive_bins_1.500.dat", delim_whitespace=True, \
        header=1, names=['log10(M_1i)(Msun)', 'qratio(M_2i/M_1i)',\
        'log10(P_i)(days)', 'result', 'CE_flag', 'log_L_1', 'log_T_1',\
        'M_1f(Msun)', 'M_2f(Msun)', 'Porb_f(d)', 'tmerge(Gyr)', 'He_core_1(Msun)',\
        'C_core_1(Msun)', 'runtime(min)'])

mesa_0300 = mesa_data1[0:46]
mesa_0350 = mesa_data1[46:]
#mesa_0350 = mesa_data1[46:92]
#mesa_0400 = mesa_data1[92:]

history_error = '0.100_3.370_error/LOGS1/history.data'
history_in = '0.100_3.390_no_interaction/LOGS1/history.data'
history_smt = '0.100_3.380_stable_MT/LOGS1/history.data'
history_merger = '0.100_3.280_CE_merger/LOGS1/history.data'
history_ejection_m = '0.100_3.320_CE_ejection_m/LOGS1/history.data'
history_ejection = '0.100_3.350_CE_ejection/LOGS1/history.data'


"""General functions for reading history.data and binary_history.data files"""


def read(file):
        data = pd.read_csv(file, delim_whitespace=True, \
        header=1, names=['log10(M_1i)(Msun)', 'qratio(M_2i/M_1i)',\
        'log10(P_i)(days)', 'result', 'CE_flag', 'log_L_1', 'log_T_1',\
        'M_1f(Msun)', 'M_2f(Msun)', 'Porb_f(d)', 'tmerge(Gyr)', 'He_core_1(Msun)',\
        'C_core_1(Msun)', 'runtime(min)'])
        return data


    
def separation_ratio(string, save = False, name = "000"):
    history = md.mesa_data(string)
    model = history.get("model_number")
    a = history.get("Separation_Ratio")
    #CE = history.get("CE_flag")
    rl1 = history.get("rl_relative_overflow_1")
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Model Number")
    ax1.set_ylabel("Separation Ratio")
    ax1.plot(model, a, label="Binary Separation", color="c")
    ax1.tick_params(axis='y', labelcolor="c")
    plt.legend()
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    #ax2.set_ylabel("Common Envelope Flag")  # we already handled the x-label with ax1
    #ax2.plot(model, CE, label="CE flag", color="m")
    ax2.tick_params(axis='y', labelcolor="m")
    ax2.set_ylabel("Relative Overflow")
    ax2.plot(model, rl1, label="rl 1", color="m")
    plt.legend()
    plt.title("Separation Ratio " + name)
    
    fig.tight_layout()
    if save:
        plt.savefig( "Separation_Ratio_" + name + ".png")
    plt.show()


def b_separation(string, save = False, name = "000"):
    history = md.mesa_data(string)
    model = history.get("model_number")
    a = history.get("binary_separation")
    #CE = history.get("CE_flag")
    rl1 = history.get("rl_relative_overflow_1")
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Model Number")
    ax1.set_ylabel("Binary Separation")
    ax1.plot(model, a, label="Binary Separation", color="b")
    ax1.tick_params(axis='y', labelcolor="c")
    plt.legend()
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    #ax2.set_ylabel("Common Envelope Flag")  # we already handled the x-label with ax1
    #ax2.plot(model, CE, label="CE flag", color="m")
    ax2.tick_params(axis='y', labelcolor="m")
    ax2.set_ylabel("Relative Overflow")
    ax2.plot(model, rl1, label="rl 1", color="m")
    plt.title("Binary Separation " + name)
    
    fig.tight_layout()
    if save:
        plt.savefig( "binary_separation_" + name + ".png")
    plt.show()


def rl_overflow(string, save = False, name = "000"):
    history = md.mesa_data(string)
    #model = history.get("model_number")
    model = history.get("star_1_mass")
    rl1 = history.get("rl_relative_overflow_1")
    Radius_flag = history.get("Radius_Lagrangian_flag")
    L3_flag = history.get("Outer_Lagrangian_flag")
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Donor Mass")
    ax1.set_ylabel("Relative Roche Lobe overflow")
    ax1.plot(model, rl1, label="Roche Lobe overflow", color="c")
    ax1.tick_params(axis='y', labelcolor="c")

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("Mass Transfer Flag")  # we already handled the x-label with ax1
    ax2.plot(model, L3_flag, label="Outer Mass Transfer Flag", color="m")
    ax2.plot(model, Radius_flag, label="Radius Outer Mass Transfer Flag")
    ax2.tick_params(axis='y', labelcolor="m")
    plt.legend()

    plt.title("Mass Transfer Rate")
    if save:
        plt.savefig( "rl_relative_overflow_" + name + ".png")
    #plt.show()

def mt_types(string, save = False, name = "000"):
    
    #take abs and put on log scale
    
    history = md.mesa_data(string)
    #model = history.get("model_number")
    model = history.get("star_1_mass")
    rl1 = history.get("rl_relative_overflow_1")
    CE = history.get("CE_flag")
    Outer_thick = history.get("Outer_MT_thick")
    Outer_thick = np.abs(Outer_thick)
    L1_thick = history.get("L1_MT_thick")
    L1_thick = np.abs(L1_thick)
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Donor Mass")
    ax1.set_ylabel("Relative Roche Lobe overflow")
    ax1.plot(model, rl1, label="Roche Lobe overflow", color="c")
    ax1.plot(model, CE, label="CE flag", color = "darkorange")
    ax1.tick_params(axis='y', labelcolor="c")
    plt.legend(loc = 'lower left')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("Mass Transfer")  # we already handled the x-label with ax1
    ax2.plot(model, Outer_thick, label="Outer MT thick", color="m")
    ax2.plot(model, L1_thick, label="L1 MT thick")
    ax2.tick_params(axis='y', labelcolor="m")
    ax2.set_yscale("log")
    plt.legend()

    plt.title("Mass Transfer Rate")
    if save:
        plt.savefig( "MT_rates_" + name + ".png")
    #plt.show()

def mass_comp(string):
    history = md.mesa_data(string)
    model = history.get("model_number")
    m1 = history.get("star_1_mass")
    m2 = history.get("star_2_mass")
    plt.plot(model,m1,label="Mass 1")
    plt.plot(model,m2,label="Mass 2")
    plt.xlabel("Model Number")
    plt.ylabel("Mass in Msun")
    plt.title("Mass Comparison")
    plt.legend()
    
def mass_ce(string):
    history = md.mesa_data(string)
    model = history.get("model_number")
    m1 = history.get("star_1_mass")
    m2 = history.get("star_2_mass")
    CE = history.get("CE_flag")
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("model number")
    ax1.set_ylabel("Mass in Msun")
    ax1.plot(model, m1, label="Mass 1", color="b")
    ax1.plot(model, m2, label="Mass 2", color="g")
    ax1.tick_params(axis='y', labelcolor="c")
    plt.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("Common Envelope Flag")  # we already handled the x-label with ax1
    ax2.plot(model, CE, label="CE flag", color="m")
    ax2.tick_params(axis='y', labelcolor="m")
    
    plt.title("Mass to Common Envelope Comparison")

    fig.tight_layout()

def mass_rate_comp(string):
    history = md.mesa_data(string)
    model = history.get("model_number")
    m1 = history.get("lg_mstar_dot_1")
    m2 = history.get("lg_mstar_dot_2")
    plt.plot(model,m1,label="Mass 1")
    plt.plot(model,m2,label="Mass 2")
    plt.xlabel("Model Number")
    plt.ylabel("Log Rate of Mass Loss")
    plt.title("Mass Loss")
    plt.legend()
    
def mass_trans(string, string2= "name"):
    history = md.mesa_data(string)
    #model = history.get("model_number")
    model = history.get("star_1_mass")
    m1 = history.get("lg_mtransfer_rate")
    plt.plot(model, m1, label = string2)
    plt.xlabel("Donor Mass")
    plt.ylabel("Log Rate of Mass Transfer")
    plt.title("Mass Transfer Rate")
    
def mass_types(string):
    history = md.mesa_data(string)
    #model = history.get("model_number")
    model = history.get("star_1_mass")
    L1_thin = history.get("L1_MT_thin")
    L1_thin = - L1_thin + 1e-8
    L1_thick = history.get("L1_MT_thick")
    L1_thick = - L1_thick + 1e-8
    Outer_thick = history.get("Outer_MT_thick")
    Outer_thick = - Outer_thick + 1e-8
    Outer_thin = history.get("Outer_MT_thin")
    Outer_thin = - Outer_thin + 1e-8
    plt.plot(model, np.log10(L1_thin), linewidth = 2,label = "Inner MT thin")
    plt.plot(model, np.log10(L1_thick), linewidth = 2, label = "Inner MT thick")
    plt.plot(model, np.log10(Outer_thick), linewidth = 2, label = "Outer MT thick")
    plt.plot(model, np.log10(Outer_thin), linewidth = 2, label = "Outer MT thin")
    #plt.gca().invert_yaxis()
    #plt.yscale("log")
    plt.rcParams['xtick.labelsize']=20
    plt.rcParams['ytick.labelsize']=20
    plt.ylim(-6,-2)
    plt.xlabel(r'Donor Mass $(M_{\odot})$',fontsize=20)
    plt.ylabel(r'Log Rate of Mass Transfer $(M_{\odot}/year)$',fontsize=20)
    #plt.title("Types of Mass Transfer Rate",fontsize=20)
    
def mass_trans_CE(string, save = False, name = "000"):
    history = md.mesa_data(string)
    model = history.get("model_number")
    mtrans = history.get("lg_mtransfer_rate")
    CE = history.get("CE_flag")
    fig, ax1 = plt.subplots(figsize=(9,5))
    ax1.set_xlabel("model number")
    ax1.set_ylabel("Log mass transfer rate")
    ax1.plot(model, mtrans, ".",label="mass transfer rate", color="c")
    ax1.tick_params(axis='y', labelcolor="c")
    ax1.set_ylim(-10,0)
    plt.legend()
    plt.title("mass_transfer_rate_" + name)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("Common Envelope Flag")  # we already handled the x-label with ax1
    ax2.plot(model, CE, label="CE flag", color="m")
    ax2.tick_params(axis='y', labelcolor="m")
    plt.legend()
    if save:
        plt.savefig( "mass_transfer_rate_" + name + ".png")
    plt.show()

def radius_types(string):
    history = md.mesa_data(string)
    model = history.get("star_1_mass")
    R_outer = history.get("Outer_radius")
    R_outer = R_outer/Rsun
    R_Roche = history.get("rl_1")
    R_star = history.get("log_R")
    R_star = 10**R_star
    plt.plot(model, np.log10(R_outer), label = "Outer Radius")
    plt.plot(model, np.log10(R_Roche), label = "Roche lobe Radius")
    plt.plot(model, np.log10(R_star), label = "Star Radius")
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18
    plt.xlabel("Donor Mass",fontsize=18)
    plt.ylabel("Binary Radii",fontsize=18)
    plt.title("Different Radii",fontsize=20)
    
def star_radius(string):
    history = md.mesa_data(string)
    model = history.get("star_1_mass")
    R_star = history.get("log_R")
    R_star = 10**R_star
    plt.plot(model, np.log10(R_star), label = "Star Radius no outer MT")
    
def separation(string1,string2):
    history = md.mesa_data(string1)
    history_no = md.mesa_data(string2)
    model1 = history.get("star_1_mass")
    model2 = history_no.get("star_1_mass")
    sep = history.get("binary_separation")
    old_sep = history_no.get("binary_separation")
    plt.plot(model1, np.log10(sep), label = "Separation with outer Mass Loss")
    plt.plot(model2, np.log10(old_sep), label = "Separation without outer Mass Loss")
    plt.hlines(0)
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18
    plt.xlabel("Donor Mass (Rsun)",fontsize=18)
    plt.ylabel("log Binary separation (Rsun)",fontsize=18)
    plt.title("Binary Separation",fontsize=20)
    
    
    
def rel_radius_types(string1,string2):
    history1 = md.mesa_data(string1)
    history2 = md.mesa_data(string2)
    model1 = history1.get("star_1_mass")
    model2 = history2.get("star_1_mass")
    R_outer1 = history1.get("outer_relative_overflow")
    R_Roche1 = history1.get("rl_relative_overflow_1")
    R_outer2 = history2.get("outer_relative_overflow")
    R_Roche2 = history2.get("rl_relative_overflow_1")
    plt.plot(model1, R_outer1, "-.", color = 'orange',label = "Outer Lagrangian w/outer MT")
    plt.plot(model1, (R_Roche1), ":", color = 'orange', label = "Roche Lobe w/outer MT")
    plt.plot(model2, (R_outer2), "-.", color = 'b', label = "Outer Lagrangian w/out outer MT")
    plt.plot(model2, (R_Roche2), ":", color = 'b', label = "Roche Lobe w/out outer MT")
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18
    plt.xlabel("Donor Mass (Rsun)",fontsize=18)
    plt.ylabel("Relative Overflow",fontsize=18)
    plt.title("Relative Overflow",fontsize=20)
    
    
def duplex(string1,string2):
    history1 = md.mesa_data(string1)
    history2 = md.mesa_data(string2)
    model1 = history1.get("star_1_mass")
    model2 = history2.get("star_1_mass")
    sep = history1.get("binary_separation")
    old_sep = history2.get("binary_separation")
    R_outer1 = history1.get("outer_relative_overflow")
    R_Roche1 = history1.get("rl_relative_overflow_1")
    R_outer2 = history2.get("outer_relative_overflow")
    R_Roche2 = history2.get("rl_relative_overflow_1")
    
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0})
    ax1.plot(model1, R_outer1, "-.", color = 'orange',linewidth=2)#,label = "Outer Lagrangian w/outer MT")
    ax1.plot(model1, (R_Roche1), color = 'orange',linewidth=2)#, label = "Roche Lobe w/outer MT")
    ax1.plot(model2, (R_outer2), "-.", color = 'b',linewidth=2)#, label = "Outer Lagrangian w/out outer MT")
    ax1.plot(model2, (R_Roche2),color = 'b',linewidth=2)#, label = "Roche Lobe w/out outer MT")
    ax2.plot(model1, np.log10(sep), ':',color = 'orange',linewidth=2)#, label = "Separation with outer Mass Loss")
    ax2.plot(model2, np.log10(old_sep), ':',color='b',linewidth=2)# label = "Separation without outer Mass Loss")
    ax1.hlines(0, 13.5,30.5,'grey','dashed')
    
    
    legend_elements = [plt.Line2D([0], [0], color = 'orange',lw = 6, label = "Including outer MT"),
                       plt.Line2D([0], [0], color = 'b', lw=6,label = "Not including outer MT"),
                       plt.Line2D([0], [0], color = 'k',label = "Roche Lobe relative overflow"),
                       plt.Line2D([0], [0], linestyle='-.',color = 'k',label = "Outer Lagrangian relative overflow"),
                       plt.Line2D([0], [0], linestyle=':',color = 'k', label = "Separation")]

    #plt.legend(handles=legend_elements,loc='upper left')
    plt.legend(handles=legend_elements,loc='upper center',bbox_to_anchor=(0.5, 2.3),ncol=3,fontsize=16)
    
    
    #legend_elements2 = [plt.Line2D([0], [0], color = 'orange',lw = 6, label = "Separation with outer Mass Loss"),
    #                   plt.Line2D([0], [0], color = 'b', lw=6,label = "Separation without outer Mass Loss")]
    
    #ax2.legend(handles=legend_elements2)
    
    
    ax1.set_ylim(-.3,.5)
    #ax2.set_ylim(0,3.5)
    
    #ax1.set_ylabel(r'$1/q$')
    
    plt.rcParams['xtick.labelsize']=18
    plt.rcParams['ytick.labelsize']=18
    
    plt.xlim(13.5,30.5)
    plt.xlabel(r'Donor mass $(M_{\odot})$',fontsize=20)
    ax1.set_ylabel("Relative Overflow",fontsize=20)
    ax2.set_ylabel(r'log Binary separation $(R_{\odot})$',fontsize=20)
        
    plt.show()


    
    
    
    
    
    
    
    
    
    
    
