# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:15:09 2023

@author: Cristiano
"""

import pandas as pd
from pathlib import Path
import csv
import time
import os
import random
from os import listdir
from os import makedirs
from os.path import isfile, join, exists
import numpy as np
import tkinter as tk
from tkinter import filedialog
from rhofun import rhofun
from fr_rothC import fr_rothC
# Vsevolod: for some reasons it didn't work on my Ubuntu 22.04. It is not used further so I commented it out
#from NS import NS 
import matplotlib.pyplot as plt
import scipy as sp
import scipy.linalg
import time # Vsevolod: added to measure execution time of the algorithms

from sklearn.metrics import mean_absolute_error,r2_score,mean_squared_error
import numpy as np

def getfile():
    import tkinter as Tkinter, tkinter.filedialog as tkFileDialog
    root = Tkinter.Tk()
    root.after(100, root.focus_force)
    root.after(200,root.withdraw)    
    file_path = tkFileDialog.askopenfilename(parent=root,title='Pick a file')    
    return file_path 

print('Select scenario Excel file')

scenario="polonia"

filename=getfile()
tab = pd.read_excel(filename, header=None)

# Read scenario
outname =tab[1][0]
compare = 'OUTPUT_%s/compare_%s.xlsx'%(outname,outname)
comp_table = pd.read_excel(compare, header=0)
years=comp_table.values[:,0]
for i in range(len(years)):
	years[i]=float(years[i])+float(i%12)/12.0
    
if scenario=="ucraina":
#use this for ucraina

    to_compare=np.array(comp_table.values[:,2], dtype='float')
    plt.plot(years, comp_table.values[:,3], 'gold', label='q=0.97 corr=1')
    plt.plot(years, comp_table.values[:,4], 'navy', label='q=0.97 corr=0')
    plt.plot(years, comp_table.values[:,5], 'orangered', label='q=1 model 1')
    plt.plot(years, to_compare, 'D',color='green', label='Measured values')
    #�plt.plot(years,soc_mean)
    plt.xlabel('Year')
    plt.ylabel('SOC')
    plt.grid(axis = 'y')
    #plt.ylim(20, 80)
    plt.legend(loc="lower left")
    
    measure=np.where(~np.isnan(to_compare))[0]
    obs=to_compare[measure]
    for i in range(3,6): 
        predict=pd.to_numeric(comp_table.values[:,i][measure])
        r2=r2_score(obs, predict)
        rmse=np.sqrt(mean_squared_error(obs, predict))
        RMSE='rmse_'+comp_table.columns[i]
        R2='r2_'+comp_table.columns[i]
        comp_table[RMSE]=rmse
        comp_table[R2]=r2
    comp_table.to_excel(compare,index=False)
    
    
elif scenario=="polonia":
    #use below for plot polonia
    years=years[years<2016]
    ind=np.where(years<2016)[0]
    
    to_compare=np.array(comp_table.values[:,2][ind], dtype='float')
    plt.plot(years, comp_table.values[:,3][ind], 'gold', label='q=0.95 corr=1')
    plt.plot(years, comp_table.values[:,4][ind], 'navy', label='q=0.95 corr=0')
    plt.plot(years, comp_table.values[:,5][ind], 'orangered', label='q=1')
    plt.plot(years, to_compare[ind], 'D',color='green', label='Measured values')
    
    plt.xlabel('Year')
    plt.ylabel('SOC')
    plt.grid(axis = 'y')
    #plt.ylim(np.min(to_compare[~np.isnan(to_compare)])-2, 60)
    #plt.ylim(20, 80)
    plt.legend(loc="lower left")
    
    #error and efficiency polonia
    measure=np.where(~np.isnan(to_compare))[0]
    obs=to_compare[measure]
    for i in (3,4,5): 
        predict=pd.to_numeric(comp_table.values[:,i][ind][measure])
        r2=r2_score(obs, predict)
        rmse=np.sqrt(mean_squared_error(obs, predict))
        RMSE='rmse_'+comp_table.columns[i]
        R2='r2_'+comp_table.columns[i]
        comp_table[RMSE]=rmse
        comp_table[R2]=r2
    comp_table.to_excel(compare,index=False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    