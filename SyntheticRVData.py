#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:11:10 2019

@author: ZhexingLi
"""

# Script that produces synthetic RV data at any epoch for any time baseline, with
# some randomness in time of observation and rv scatter to mimic real data.

# import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import root
import random

# Input options to toggle

# First element as first planet, second element as second planet etc.
P = [29.8]               # orbital period of the planet(s) (unit of days), 29.8
Tp = [2455000]          # time of periastron, in unit of days (JD), in the past (Tp < ref1)
e = [0.0]                 # eccentricity
K = [1]               # semi amplitude, in unit of m/s, 0.0024
w = [0]               # argument of periastron of the STAR, in unit of radians
V_0 = 0.0                     #systematic velocity, in unit of m/s, set here universally as 0 m/s

rms = 1               # rms from original fitting , 0.00279
unc = 1               # rv data uncertainty

time_base = 2453000       # number to subtract for display purpose

# type of time stamp, either 'complete', i.e. provided by user input, or 'reference', i.e.
# only start and end time stamps provided by user
time_stamp = 'complete'

if time_stamp == 'complete':
    # path to the dataset containing the needed timestamp
    path = '/Users/zhexingli/Desktop/UCR/01-RESEARCH/PROJECTS/AAT/Targets/08-HD134606/HD134606AATbinned.txt'
    data = pd.read_csv(path,header=None,skiprows=1,delim_whitespace=True,usecols=[0],names=['time'])
    stamps = list(data['time'])
    num = len(stamps)
elif time_stamp == 'reference':
    # reference time 1, all new data starts after this
    # reference time 2, all new data ends before this. If set to 0, then end time determined by outer period
    ref1 = 2452984.7100400003
    ref2 = 2459483.121358  
    stamps = [ref1,ref2]
    num = 66               # number of data points
else:
    print("Please provide a valid time stamp type.")
# path to store created data
savepath = '/Users/zhexingli/Desktop/UCR/01-RESEARCH/PROJECTS/AAT/Targets/08-HD134606/AATActivity/'
name = 'HD134606'


#######################################################################################

# compute the ending time if ref2 given is 0
if stamps[1] == 0:
    start = stamps[0]
    end = start + np.max(P)
    stamps = [start,end]
else:
    pass
    
def Synth_Data (error,scatter,number):
    '''
    Function that calculates required number of synthetic data points at the required
    time specified above. Also randomizes and scatters the calculated RV points.
    
    Input: 1) error: Error estimate for each data point
           2) scatter: RMS from the original fitting as gaussian sigma to simulate scatter
           3) number: Number of data points for each run
    Output: One text file that stores time of data, rv, and err in the same directory
            as other data files and setup file.
    '''
    
    startT, endT = stamps[0],stamps[-1]
    
    length = int(endT - startT)
    
    # plot out individual and total rv curves
    datapoint(stamps,int(length*2))
    
    if number > 2 and time_stamp == 'reference':
    # get individual data points
    # add extra number of points so some could be randomly excluded for the purpose
    # of adding more randomness in the time of each data point
        add = int(number*1)
        number = number + add
        vel, time_random, vel_ind = datapoint(stamps,number,data = 'Yes')
    
        '''
        # individual planetary components (noiseless)
        pl_b = list(vel_ind[0])
        pl_c = list(vel_ind[1])
        pl_d = list(vel_ind[2])
        '''
        vel = list(vel)
        time_random = list(time_random)
        
        # randomly exclude number of data points that are added earlier to create more
        # randomness in the time
        point = np.arange(1,len(vel)-1,1)     # number indices to be randomly selected except the first and the last data points
        exclude = random.sample(list(point),add)      # list containing the indices of the data points to be excluded
        for item in exclude:
            ind = item
            vel[ind] = 'NaN'
            time_random[ind] = 'NaN'
            #pl_b[ind] = 'NaN'          # individual planetary components are noiseless
            #pl_c[ind] = 'NaN'
            #pl_d[ind] = 'NaN'
    
        vel = [x for x in vel if x != 'NaN']
        time_random = [x for x in time_random if x != 'NaN']
        #pl_b = [x for x in pl_b if x != 'NaN']
        #pl_c = [x for x in pl_c if x != 'NaN']
        #pl_d = [x for x in pl_d if x != 'NaN']
        vel = np.array(vel)
        time_random = np.array(time_random)
    else:
        vel, time_random, vel_ind = datapoint(stamps,number,data = 'Yes')
        
    # Randonmize calculated total velocity data with the error provided using a Gaussian
    # filter.
    vel_random = []
    for i in range (len(vel)):
        new_vel = random.gauss(vel[i],scatter)
        vel_random.append(new_vel)
        
    vel_random = np.array(vel_random)

    #pl_b = np.array(pl_b)      # again, individual component is noiseless
    #pl_c = np.array(pl_c)
    #pl_d = np.array(pl_d)
    
    plt.errorbar(time_random - time_base,vel_random,yerr=unc,marker = 'o', mfc = 'k',\
                 mec = 'cornflowerblue',markersize=5,ls='none',label='Total Randomized')
    #plt.legend(loc='upper left',prop={'size':5},fontsize='large')
    plt.legend(loc='upper left',fontsize='small')
    plt.xlabel('Time - {base} (JD)'.format(base=time_base))
    plt.ylabel('RV (m/s)')
    plt.show()
    
    df = pd.DataFrame()
    df['time'] = time_random
    df['rv'] = vel_random
    df['err'] = error
    '''
    df1 = pd.DataFrame()
    df1['time'] = time_random
    df1['rv'] = pl_b
    df1['err'] = error
    
    df2 = pd.DataFrame()
    df2['time'] = time_random
    df2['rv'] = pl_c
    df2['err'] = error
    
    df3 = pd.DataFrame()
    df3['time'] = time_random
    df3['rv'] = pl_d
    df3['err'] = error
    '''
    path0 = savepath + name + '_SynthRV.txt'
    df.to_csv(path0,sep='\t',index=False, header=True)
    '''
    path0 = path + 'rvdata_noiseless_b.txt'
    df1.to_csv(path0,sep='\t',index=False, header=True)
    path0 = path + 'rvdata_noiseless_c.txt'
    df2.to_csv(path0,sep='\t',index=False, header=True)
    path0 = path + 'rvdata_noiseless_d.txt'
    df3.to_csv(path0,sep='\t',index=False, header=True)
    '''
    figpath = savepath + name + '_SynthRV.png'
    #plt.savefig(figpath,dpi=500)


#######################################################################################

def datapoint(stamps, rv_num, data = None):
    '''
    Function that produce theoretical RV data points for individual planet and the
    entire system at the specified time stamps.
    
    Input: 1) stamps: full time stamps provided; or just the start and end points
           2) rv_num: number of rv data points to be produced
    Output: Non-randomized non-scattered RV data points at time stamps provided for
            the whole system as well as for the individual component
    '''
    
    # Determine the actual times to create data points based on the starting time and 
    # number of data points. Actual time depends on the one full period of the outer
    # most planet. Or just take time stamps provided.
    
    newtime = []              # Synthetic data points time in JD

    if time_stamp == 'reference' or data is None:
        start, end = stamps[0], stamps[-1]
        if rv_num > 1:
            time_step = (end - start)/(rv_num-1)
            for i in range (rv_num):
                t = start + time_step*i
                newtime.append(t)
        else:
            t = start       # number = 1 case, put data at the start
            newtime.append(t)
    elif time_stamp == 'complete':
        newtime = stamps
    else:
        print("Invalid type of time stamps.")
        
    # If just plotting individual rv curves, don't randomize time, otherwise,
    # slightly randomize the time of observation to avoid false periodicity (aliases)
    # to be picked up by periodogram
    if time_stamp == 'complete' or data is None:
        t_random = np.array(newtime)
    else:
        t_random = []
        for i in range (len(newtime)):
            new_time = random.gauss(newtime[i],3)       # sigma of 3 day
            t_random.append(new_time)
        t_random.sort()
        t_random = np.array(t_random)
    
    # Determine the total velocity values for all planets
    rv = 0         # total combined rv
    rv_ind = []    # individual planetary rv signals stored in lists
    for p in range (len(P)):
        # Determine orbital phase for each planet
        phi = []
        for t in t_random:
            ratio = (t-Tp[p])/P[p]
            if ratio > 1.0:  # when t is more than one period bigger than tp
                phase = ratio - int(ratio)
                phi.append(phase)
            elif 0 <= ratio <= 1:  # when t is within one period bigger than tp
                phi.append(ratio)
            elif -1 <= ratio < 0:  # when t is within one period less than tp
                phase = 1 + ratio
                phi.append(phase)
            else:   # when t is more than one period less than tp
                phase = 1 + (ratio - int(ratio))
                phi.append(phase)
                
        phi = np.array(phi)

        # Calculate the corresponding mean anomaly according to orbital phase
        m = 2*(np.pi)*phi/1.0              # mean anomaly

        # Calculate the eccentric anomaly according to mean anomaly
        eps = []
        for i in range(len(phi)):
            def func(x):
                return m[i] - x + e[p]*(np.sin(x))
            sol = root(func,1)
            eps.append(sol.x[0])
        eps = np.array(eps)

        # Calculate the true anomaly according to eccentric anomaly
        nu = []
        for i in range(len(phi)):
            v = np.arccos(((np.cos(eps[i])) - e[p])/(1 - e[p]*(np.cos(eps[i]))))
            if eps[i] <= np.pi:
                nu.append(v)
            else:
                v = np.pi + abs(np.pi - v)
                nu.append(v)
        nu = np.array(nu)
        
        # Calculate the velocity according to true anomaly
        V = V_0 + K[p]*(np.cos(w[p]+nu)+e[p]*(np.cos(w[p])))
        rv_ind.append(V)
        if data is None:
            if p == 0:
                plt.plot(t_random - time_base,V,'m-',alpha=0.7,label='Planet b')
            elif p == 1:
                plt.plot(t_random - time_base,V,'r-',alpha=0.7,label='Planet c')
            elif p == 2:
                plt.plot(t_random - time_base,V,'g-',alpha=0.7,label='Planet d')
            else:
                pass
            rv = rv + V
        else:
            rv = rv + V
    
    if data is None:
        plt.plot(t_random - time_base,rv,'k-',label='Total')
    else:
        pass
    
    rv = np.array(rv)
    
    # store unrandomized rv curve data
    df = pd.DataFrame()
    df['time'] = t_random
    df['rv'] = rv
    
    path0 = savepath + name + '_SynthRV.txt'
    #df.to_csv(path0,sep='\t',index=False, header=True)
    
    return rv, t_random, rv_ind

Synth_Data(unc,rms,num)
