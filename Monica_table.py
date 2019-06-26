"""Monica's code for creating the tables"""


import matplotlib.lines as mlines
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import ascii

from matplotlib import rcParams
rcParams['font.family'] = 'stixgeneral'
rcParams['font.size'] = 14

def make_axis(qunique, punique): #def make_axis(qunique, p):
    '''
    requires: all unique xvalues, all unique y values
    returns: x, Xgrid, y, Ygrid
    '''
    #q = np.arange(qunique[0], qunique[-1], 0.05)
    #p = np.arange(punique[0], punique[-1] + 0.01, 0.01)
    
    q_step = qunique[1] - qunique[0]
    p_step = punique[1] - punique[0]
    
    q = np.arange(qunique[0], qunique[-1] + q_step, q_step) #for the small table
    #q = np.arange(qunique[0], qunique[-1] - q_step, q_step) #for the full table
    p = np.arange(punique[0], punique[-1], p_step)
    
    # x axis
    #m2_masses = np.linspace(5, 30, 17)
    shift = ((q[1]-q[0]))/2
    new_q = q-shift

    last = q[-1]+shift
    #new_m = list(m2_masses)
    new_q = list(new_q)
    new_q.append(last)
    
    # y axis
    #T_log = np.log10(np.logspace(Psmall,Plarge,p_num))
    shift = ((p[1]-p[0]))/2
    new_p = p-shift
    last = p[-1]+shift
    new_p = list(new_p)
    new_p.append(last)    

    # make mesh
    QQ, PP = np.meshgrid(new_q, new_p)
    return q, QQ, p, PP

def get_unique_values(data):
    '''
    requires: data
    returns: unique values for P, q, and results
    '''
    
    p = data['log10(P_i)(days)']
    punique = np.unique(p)
    
    q = data['qratio(M_2i/M_1i)']
    qunique = np.unique(q)

    res = data['result']
    resunique = np.unique(res)    
    
    return punique, qunique, resunique

def get_outer_values(data):
    '''
    requires: data
    returns: unique values for P, q, and results
    '''
    
    p = data['log10(P_i)(days)']
    punique = np.unique(p)
    
    q = data['qratio(M_2i/M_1i)']
    qunique = np.unique(q)

    res = data['Outer_flag']
    resunique = np.unique(res)    
    
    return punique, qunique, resunique

def get_stable_values(data):
    '''
    requires: data
    returns: unique values for P, q, and results
    '''
    
    p = data['log10(P_i)(days)']
    punique = np.unique(p)
    
    q = data['qratio(M_2i/M_1i)']
    qunique = np.unique(q)

    res = data['CE_flag']
    resunique = np.unique(res)    
    
    return punique, qunique, resunique

## edit here 

my_data = '/Users/kaliroe/Downloads/mesa-dev/massive_bins_1.500.dat'
alpha1 = '/Users/kaliroe/Downloads/mesa-dev/massive_bins_alpha01.dat'
alpha5 = '/Users/kaliroe/Downloads/mesa-dev/massive_bins_alpha05.dat'
full_table = '/Users/kaliroe/Downloads/mesa-dev/full_table_massive_bins_1.500.dat'
cutoff = '/Users/kaliroe/Downloads/mesa-dev/cutoff_table_massive_bins_1.500.dat'

# to combine two different data sets
#data1 = pd.read_csv('./1.500__0.1_0.15.dat', comment='#', delim_whitespace=True)
#data2 = pd.read_csv('./1.500__0.3.dat', comment='#', delim_whitespace=True)
#data = data1.append(data2, ignore_index=True)

def make_graph(file, show = True, save = False, name = 'massive_bins_graph'):
    data = ascii.read(file)
    #data = pd.read_csv(file, comment='#', delim_whitespace=True)
    
    punique, qunique, resunique = get_unique_values(data)
    
    qvals, QQ, pvals, PP = make_axis(qunique, punique)
    
    # determine colors
    color_dictionary = {'CE_ejection':   '#EC6666',
                      'CE_ejection_m': '#90A245',
                      'CE_merger':     '#F5C258',
                      'error':         '#1668E8',
                      'no_interaction':'#473335',
                      'stable_MT_to_wide_binary': '#98C0CB',
                      'stable_MT_to_merger': 'thistle',
                      'None':           'white',
                      'No_Run':          'white' } #changed
    
    final_colors = []
    res_dictionary = {}
    for i in range(0,len(resunique)+1):
        
        if i < (len(resunique)):
            result = resunique[i]
            res_dictionary[result] = i         
            
        # always add a result = "None" corresponding to grids without data
        else:
            result = 'None'
            res_dictionary[result] = i
        final_colors.append( color_dictionary[result] )
        
    # create color boundaries to be used in plot    
    color_bounds = []
    vals = np.arange(0,len(res_dictionary))
    shift = ((vals[1]-vals[0]))/2
    new_vals = vals-shift
    last = vals[-1]+shift
    color_bounds = list(new_vals)
    color_bounds.append(last)    
    
    print((res_dictionary))
    
    # here is where things can break
    # np.arange does weird things with the precision of floats. 
    # Something this cell works perfectly and other times the values/shape of 
    # pvals is wrong
    
    # if the code crashes because the shapes of  X, Y, and Z in the plot are not 
    # compatible, this cell might be wrong
    
    # In that case, make sure that pvals has a min = punique[0] and max = punique[-1]
    # with increments of 0.01
    
    print( (punique[-1]))
    #pvals = np.arange(start = punique[0], stop = punique[-1] + 0.01, step=0.01 )
    p_step = punique[1] - punique[0]
    pvals = np.arange(start = punique[0], stop = punique[-1], step=p_step )
    print(pvals)
    
    evol_all = []
    for i, mass_ratio in enumerate(qvals):
        
        m_i_idx = np.where(data['qratio(M_2i/M_1i)'] ==   np.around(mass_ratio,2))[0]
        p_for_mi = data['log10(P_i)(days)'][m_i_idx]
        result_for_mi = data['result'][m_i_idx]
        
        # check if the index array is empty
        # if empty, result at all P for that mass_ratio = 100
        if m_i_idx.size == 0:
            evol = np.ones_like(pvals)*res_dictionary['None']
            
        else: 
            evol = np.ones_like(pvals)
            for j, period in enumerate(pvals):
                # np.arange includes too many decimals (e.g. 0.1 = 0.99999998) and np.where does not work
                # np.where will result indicies starting from 0
                # result_for_mi and p_for_mi have can have indicies starting with values != 0 because p_for_mi
                #          is a pandas DataFrame
                # I change the pandas DataFrame to arrays to that their inidicies will start from 0
    
                # periods and results in data
                p_for_mi_array      = np.array(p_for_mi)
                result_for_mi_array = np.array(result_for_mi)
                
                # grab index for the result of mass_ratio_i and period_j
                result_ind = np.where(p_for_mi_array ==  np.around(period,2))[0] 
            
                if result_ind.size ==0:
                    evol[j] = res_dictionary['None']
                else:
                    result_string = result_for_mi_array[result_ind][0]
                    result_ID = res_dictionary[result_string]
                    evol[j] = result_ID
                    
        evol_all.append(evol)
        
    # QQ and PP should have the same shape
    # np.array(evol_all).T should have oness element in its columns
    
    
    
    print(np.shape(QQ))
    print(np.shape(PP))
    print(np.shape(np.array(evol_all).T))
    
    cmap = colors.ListedColormap(final_colors)
    norm = colors.BoundaryNorm(color_bounds, cmap.N)
    
    fig, ax = plt.subplots(figsize=(10,11))
    h = plt.pcolormesh(QQ, PP, np.array(evol_all).T, cmap=cmap, norm=norm, edgecolors = 'white')
    
    # edit this to make the x axis tickmarks match the actual values of q 
    #ax.set_xticks(np.linspace(.3,.4,3))
    
    
    ax.set_xlabel(r'$q$', fontsize=18)
    ax.set_ylabel(r'$\log_{10}P_{\rm orb,i} \ [\rm{days}]$', fontsize=18)
    
    legend_list = []
    for i, result in enumerate(res_dictionary.items()):
        color_ = color_dictionary[result[0]]
        print(color_)
        legend_list.append( mlines.Line2D([], [], color=color_, linestyle = 'None', marker='s', markersize=11, label=result[0] )  )
    
    lgd = ax.legend( handles = legend_list ,bbox_to_anchor=(0.5, 1.05, 0.0, 0),loc = 'center', fontsize = 15., frameon = False, handletextpad=1.0, handlelength = 1.3, ncol = 3, numpoints = 1 )
    
    if save:
        plt.savefig(name + '.png')
    if show:
        plt.show()


#Graph that only shows outer mass transfer

def make_outer_graph(file, show = True, save = False, name = 'massive_bins_graph'):
    data = ascii.read(file)
    
    punique, qunique, resunique = get_outer_values(data)
    
    qvals, QQ, pvals, PP = make_axis(qunique, punique)
    
    # determine colors
    color_dictionary = {'Outer_MT':   'gold',
                      'No_outer_MT': 'tab:cyan',
                      'No_data':     'white'} #changed
    
    final_colors = []
    res_dictionary = {}
    for i in range(0,len(resunique)+1):
        
        if i < (len(resunique)):
            result = resunique[i]
            res_dictionary[result] = i         
            
        # always add a result = "None" corresponding to grids without data
        else:
            result = 'None'
            res_dictionary[result] = i
        
        if result == 0.0:
            color = 'tab:cyan'
        elif result == 1.0:
            color = 'gold'
        else:
            color = 'white'
        final_colors.append( color )
        
    # create color boundaries to be used in plot    
    color_bounds = []
    vals = np.arange(0,len(res_dictionary))
    shift = ((vals[1]-vals[0]))/2
    new_vals = vals-shift
    last = vals[-1]+shift
    color_bounds = list(new_vals)
    color_bounds.append(last)    
    
    print((res_dictionary))
    
    print( (punique[-1]))
    #pvals = np.arange(start = punique[0], stop = punique[-1] + 0.01, step=0.01 )
    p_step = punique[1] - punique[0]
    pvals = np.arange(start = punique[0], stop = punique[-1], step=p_step )
    print(pvals)
    
    evol_all = []
    for i, mass_ratio in enumerate(qvals):
        
        m_i_idx = np.where(data['qratio(M_2i/M_1i)'] ==   np.around(mass_ratio,2))[0]
        p_for_mi = data['log10(P_i)(days)'][m_i_idx]
        result_for_mi = data['Outer_flag'][m_i_idx]
        
        # check if the index array is empty
        # if empty, result at all P for that mass_ratio = 100
        if m_i_idx.size == 0:
            evol = np.ones_like(pvals)*res_dictionary['None']
            
        else: 
            evol = np.ones_like(pvals)
            for j, period in enumerate(pvals):
                # np.arange includes too many decimals (e.g. 0.1 = 0.99999998) and np.where does not work
                # np.where will result indicies starting from 0
                # result_for_mi and p_for_mi have can have indicies starting with values != 0 because p_for_mi
                #          is a pandas DataFrame
                # I change the pandas DataFrame to arrays to that their inidicies will start from 0
    
                # periods and results in data
                p_for_mi_array      = np.array(p_for_mi)
                result_for_mi_array = np.array(result_for_mi)
                
                # grab index for the result of mass_ratio_i and period_j
                result_ind = np.where(p_for_mi_array ==  np.around(period,2))[0] 
            
                if result_ind.size ==0:
                    evol[j] = res_dictionary['None']
                else:
                    result_string = result_for_mi_array[result_ind][0]
                    result_ID = res_dictionary[result_string]
                    evol[j] = result_ID
                    
        evol_all.append(evol)
        
    # QQ and PP should have the same shape
    # np.array(evol_all).T should have oness element in its columns
    
    
    
    print(np.shape(QQ))
    print(np.shape(PP))
    print(np.shape(np.array(evol_all).T))
    
    cmap = colors.ListedColormap(final_colors)
    norm = colors.BoundaryNorm(color_bounds, cmap.N)
    
    fig, ax = plt.subplots(figsize=(10,11))
    h = plt.pcolormesh(QQ, PP, np.array(evol_all).T, cmap=cmap, norm=norm, edgecolors = 'white')
    
    # edit this to make the x axis tickmarks match the actual values of q 
    #ax.set_xticks(np.linspace(.3,.4,3))
    
    
    ax.set_xlabel(r'$q$', fontsize=18)
    ax.set_ylabel(r'$\log_{10}P_{\rm orb,i} \ [\rm{days}]$', fontsize=18)
    
    legend_list = []
    for i, result in enumerate(color_dictionary.items()):
        color_ = color_dictionary[result[0]]
        legend_list.append( mlines.Line2D([], [], color=color_, linestyle = 'None', marker='s', markersize=11, label=result[0] )  )
    
    lgd = ax.legend( handles = legend_list ,bbox_to_anchor=(0.5, 1.05, 0.0, 0),loc = 'center', fontsize = 15., frameon = False, handletextpad=1.0, handlelength = 1.3, ncol = 3, numpoints = 1 )
    
    if save:
        plt.savefig(name + '.png')
    if show:
        plt.show()
        
#Graph that only shows outer mass transfer

def make_stability_graph(file, show = True, save = False, name = 'massive_bins_graph'):
    data = ascii.read(file)
    
    punique, qunique, resunique = get_stable_values(data)
    
    qvals, QQ, pvals, PP = make_axis(qunique, punique)
    
    # determine colors
    color_dictionary = {'Stable':   'lemonchiffon',
                      'Unstable': 'coral',#'cornflowerblue',
                      'No_data':     'white'} #changed
    
    final_colors = []
    res_dictionary = {}
    for i in range(0,len(resunique)+1):
        
        if i < (len(resunique)):
            result = resunique[i]
            res_dictionary[result] = i         
            
        # always add a result = "None" corresponding to grids without data
        else:
            result = 'None'
            res_dictionary[result] = i
        
        if result == 0:
            color = 'lemonchiffon'
        elif result == 1:
            color = 'coral'#'cornflowerblue'
        else:
            color = 'white'
        final_colors.append( color )
        
    # create color boundaries to be used in plot    
    color_bounds = []
    vals = np.arange(0,len(res_dictionary))
    shift = ((vals[1]-vals[0]))/2
    new_vals = vals-shift
    last = vals[-1]+shift
    color_bounds = list(new_vals)
    color_bounds.append(last)    
    
    print((res_dictionary))
    
    print( (punique[-1]))
    #pvals = np.arange(start = punique[0], stop = punique[-1] + 0.01, step=0.01 )
    p_step = punique[1] - punique[0]
    pvals = np.arange(start = punique[0], stop = punique[-1], step=p_step )
    print(pvals)
    
    evol_all = []
    for i, mass_ratio in enumerate(qvals):
        
        m_i_idx = np.where(data['qratio(M_2i/M_1i)'] ==   np.around(mass_ratio,2))[0]
        p_for_mi = data['log10(P_i)(days)'][m_i_idx]
        result_for_mi = data['CE_flag'][m_i_idx]
        
        # check if the index array is empty
        # if empty, result at all P for that mass_ratio = 100
        if m_i_idx.size == 0:
            evol = np.ones_like(pvals)*res_dictionary['None']
            
        else: 
            evol = np.ones_like(pvals)
            for j, period in enumerate(pvals):
                p_for_mi_array      = np.array(p_for_mi)
                result_for_mi_array = np.array(result_for_mi)
                
                # grab index for the result of mass_ratio_i and period_j
                result_ind = np.where(p_for_mi_array ==  np.around(period,2))[0] 
            
                if result_ind.size ==0:
                    evol[j] = res_dictionary['None']
                else:
                    result_string = result_for_mi_array[result_ind][0]
                    result_ID = res_dictionary[result_string]
                    evol[j] = result_ID
                    
        evol_all.append(evol)
        
    
    cmap = colors.ListedColormap(final_colors)
    norm = colors.BoundaryNorm(color_bounds, cmap.N)
    
    fig, ax = plt.subplots(figsize=(10,11))
    h = plt.pcolormesh(QQ, PP, np.array(evol_all).T, cmap=cmap, norm=norm, edgecolors = 'white')

    
    ax.set_xlabel(r'$q$', fontsize=18)
    ax.set_ylabel(r'$\log_{10}P_{\rm orb,i} \ [\rm{days}]$', fontsize=18)
    
    legend_list = []
    for i, result in enumerate(color_dictionary.items()):
        color_ = color_dictionary[result[0]]
        legend_list.append( mlines.Line2D([], [], color=color_, linestyle = 'None', marker='s', markersize=11, label=result[0] )  )
    
    lgd = ax.legend( handles = legend_list ,bbox_to_anchor=(0.5, 1.05, 0.0, 0),loc = 'center', fontsize = 15., frameon = False, handletextpad=1.0, handlelength = 1.3, ncol = 3, numpoints = 1 )
    
    if save:
        plt.savefig(name + '.png')
    if show:
        plt.show()
