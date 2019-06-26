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
    m1 = history.get("L1_MT_thin_flag")
    m2 = history.get("L1_MT_thin")
    #m2 = np.abs(m2)
    m2 = -m2
    #m2 = np.log10(m2)
    plt.plot(model, m1, label = "MT thin flag")
    plt.plot(model, m2, label = "MT thin rate")
    plt.xlabel("Donor Mass")
    plt.ylabel("Log Rate of Mass Transfer")
    plt.title("Thin Mass Transfer Rate")
    
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

def plot_shit():
    merger = md.mesa_data('0.100_3.280_CE_merger/binary_history.data')
    ejection = md.mesa_data('0.100_3.350_CE_ejection/binary_history.data')
    error = md.mesa_data('0.100_3.370_error/binary_history.data')
    stable = md.mesa_data('0.100_3.380_stable_MT/binary_history.data')
    no_interaction = md.mesa_data('0.100_3.390_no_interaction/binary_history.data')
    
    
    plt.plot(merger.get("model_number"),merger.get("binary_separation"),label="Merger")
    plt.plot(ejection.get("model_number"),ejection.get("binary_separation"),label="Ejection")
    plt.plot(error.get("model_number"),error.get("binary_separation"),label="Error")
    plt.plot(stable.get("model_number"),stable.get("binary_separation"),label="Stable")
    plt.plot(no_interaction.get("model_number"),no_interaction.get("binary_separation"),label="No Interaction")
    plt.legend()
    plt.savefig("separation.png")
    plt.show()
    
    plt.plot(merger.get("model_number"),merger.get("rl_relative_overflow_1"),label="Merger")
    plt.plot(ejection.get("model_number"),ejection.get("rl_relative_overflow_1"),label="Ejection")
    plt.plot(error.get("model_number"),error.get("rl_relative_overflow_1"),label="Error")
    plt.plot(stable.get("model_number"),stable.get("rl_relative_overflow_1"),label="Stable")
    plt.plot(no_interaction.get("model_number"),no_interaction.get("rl_relative_overflow_1"),label="No Interaction")
    plt.legend()
    plt.savefig("rl_overflow.png")
    plt.show()
    
    
    plt.plot(merger.get("model_number"),merger.get("lg_mtransfer_rate"),label="Merger")
    plt.plot(ejection.get("model_number"),ejection.get("lg_mtransfer_rate"),label="Ejection")
    plt.plot(error.get("model_number"),error.get("lg_mtransfer_rate"),label="Error")
    plt.plot(stable.get("model_number"),stable.get("lg_mtransfer_rate"),label="Stable")
    plt.plot(no_interaction.get("model_number"),no_interaction.get("lg_mtransfer_rate"),label="No Interaction")
    plt.legend()
    plt.savefig("mass_transfer.png")
    plt.show()


def log_dt():
    history_ni = md.mesa_data('0.100_3.390_no_interaction/LOGS1/history.data')
    history_smt = md.mesa_data('0.100_3.380_stable_MT/LOGS1/history.data')
    history_merger = md.mesa_data('0.100_3.280_CE_merger/LOGS1/history.data')
    history_error = md.mesa_data('0.100_3.370_error/LOGS1/history.data')
    history_ejection = md.mesa_data('0.100_3.350_CE_ejection/LOGS1/history.data')
    history_ejection_m = md.mesa_data('0.100_3.320_CE_ejection_m/LOGS1/history.data')
    
    log_dt_error = history_error.get("log_dt")
    model_error = history_error.get("model_number")
    
    log_dt_merger = history_merger.get("log_dt")
    model_merger = history_merger.get("model_number")
    
    log_dt_smt = history_smt.get("log_dt")
    model_smt = history_smt.get("model_number")
    
    log_dt_ni = history_ni.get("log_dt")
    model_ni = history_ni.get("model_number")
    
    log_dt_ejection = history_ejection.get("log_dt")
    model_ejection = history_ejection.get("model_number")
    
    log_dt_ejection_m = history_ejection_m.get("log_dt")
    model_ejection_m = history_ejection_m.get("model_number")
    
    plt.plot(model_smt, log_dt_smt, label = "stable mass transfer")
    
    plt.plot(model_ejection, log_dt_ejection, label = "ejection")
    
    plt.plot(model_ejection_m, log_dt_ejection_m, label = "ejection merger")
    
    plt.plot(model_ni, log_dt_ni, label = "no interaction")
    
    plt.plot(model_error, log_dt_error, label = "error")
    
    plt.plot(model_merger, log_dt_merger, label = "merger")
    
    plt.legend()
    
    plt.xlabel("Model Number")
    plt.ylabel("Log(dt)")
    
    plt.savefig("log_dt.png")
    
    plt.show()

    
def dt_retries(string, save = False , name = "000"):
    history = md.mesa_data(string)
    log_dt = history.get("log_dt")
    model = history.get("model_number")
    num_re = history.get("num_retries")
    #CE = history.get("CE_flag")

    fig, ax1 = plt.subplots(figsize=(9,5))
    ax1.set_xlabel("Model Number")
    ax1.set_ylabel("Time Step (log years)")
    ax1.plot(model, log_dt, ".", label = "log dt", color = "c")
    ax1.tick_params(axis='y', labelcolor="c")
    plt.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("Number of Retries")  # we already handled the x-label with ax1
    ax2.plot(model, num_re, ".", label="number of retries", color="m")
    ax2.tick_params(axis='y', labelcolor="m")
    plt.title("Time Step and Retries " + name)
    plt.legend()
    fig.tight_layout()
    
    if save:
        plt.savefig( "Time_Step_and_Retries_" + name + ".png")
    plt.show()


def ce_comp():
    
    history_ni = md.mesa_data('0.100_3.390_no_interaction/LOGS1/history.data')
    history_smt = md.mesa_data('0.100_3.380_stable_MT/LOGS1/history.data')
    history_merger = md.mesa_data('0.100_3.280_CE_merger/LOGS1/history.data')
    history_error = md.mesa_data('0.100_3.370_error/LOGS1/history.data')
    history_ejection = md.mesa_data('0.100_3.350_CE_ejection/LOGS1/history.data')
    history_ejection_m = md.mesa_data('0.100_3.320_CE_ejection_m/LOGS1/history.data')
    
    mass_rate_error = history_error.get("lg_mtransfer_rate")
    model_error = history_error.get("model_number")
    CE_error = history_error.get("CE_flag")
    
    mass_rate_merger = history_merger.get("lg_mtransfer_rate")
    model_merger = history_merger.get("model_number")
    CE_merger = history_merger.get("CE_flag")
    
    mass_rate_smt = history_smt.get("lg_mtransfer_rate")
    model_smt = history_smt.get("model_number")
    CE_smt = history_smt.get("CE_flag")
    
    mass_rate_ni = history_ni.get("lg_mtransfer_rate")
    model_ni = history_ni.get("model_number")
    CE_ni = history_ni.get("CE_flag")
    
    mass_rate_ejection = history_ejection.get("lg_mtransfer_rate")
    model_ejection = history_ejection.get("model_number")
    CE_ejection = history_ejection.get("CE_flag")
    
    mass_rate_ejection_m = history_ejection_m.get("lg_mtransfer_rate")
    model_ejection_m = history_ejection_m.get("model_number")
    CE_ejection_m = history_ejection_m.get("CE_flag")
    
    fig, ax1 = plt.subplots(figsize=(13,6))
    ax1.set_xlabel("Model Number")
    ax1.set_ylabel("Log Rate of Mass Loss")

    
    ax1.plot(model_smt, mass_rate_smt,".", label = "stable MT (transfer)")
    
    ax1.plot(model_ejection, mass_rate_ejection,".", label = "ejection (transfer)")

    ax1.plot(model_ejection_m, mass_rate_ejection_m,".", label = "ejection merger (transfer)")
    
    ax1.plot(model_ni, mass_rate_ni,".", label = "no interaction (transfer)")
    
    ax1.plot(model_error, mass_rate_error,".", label = "error (transfer)")
    
    ax1.plot(model_merger, mass_rate_merger,".", label = "merger (transfer)")
    
    ax1.set_ylim(-10,0)

    ax1.tick_params(axis='y', labelcolor="c")
    plt.legend()
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("CE flag")  # we already handled the x-label with ax1
    
    ax2.plot(model_smt, CE_smt, label = "stable MT (CE)")
    
    ax2.plot(model_ejection, CE_ejection, label = "ejection (CE)")
    
    ax2.plot(model_ejection_m, CE_ejection_m, label = "ejection merger (CE)")
    
    ax2.plot(model_ni, CE_ni, label = "no interaction (CE)")
    
    ax2.plot(model_error, CE_error, label = "error (CE)")
    
    ax2.plot(model_merger, CE_merger, label = "merger (CE)")
    
    ax2.tick_params(axis='y', labelcolor="m")
    
    
    plt.title("Mass Transfer Rate and CE Flag")
    plt.legend()
    fig.tight_layout()
    plt.savefig("Mass_Transfer_Rate_and_CE_Flag.png")
    plt.show()

    
    
    
