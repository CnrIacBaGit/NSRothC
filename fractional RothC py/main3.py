# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:39:02 2022

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


def getfile():
    import tkinter as Tkinter, tkinter.filedialog as tkFileDialog
    root = Tkinter.Tk()
    root.after(100, root.focus_force)
    root.after(200,root.withdraw)    
    file_path = tkFileDialog.askopenfilename(parent=root,title='Pick a file')    
    return file_path 

print('Select scenario Excel file')



filename=getfile()
tab = pd.read_excel(filename, header=None)

# Read scenario
outname =tab[1][0]
dataname = tab[1][1]
parname = tab[1][2]
# Vsevolod: added part of code to read fractional order model parameters
# fractional model parameters
use_fr_model=0
fr_q=1
fr_tau=1.0/(12.0*15.0) # time step in years (default - ~2 days)
fr_corr_type=0
# try:
#     if tab[1][3]!="":
#         use_fr_model=int(tab[1][3])
#     if tab[1][4]!="":
#         fr_q=float(tab[1][4])
#     if tab[1][5]!="":
#         fr_tau=float(tab[1][5])/12.0 # time step in a fraction of a month
#     if tab[1][6]!="":
#         fr_corr_type=int(tab[1][6])
# except:
#     pass



use_fr_model=int(input("Model= "))
if use_fr_model==1:
    use_fr_model=int(use_fr_model)
    fr_q=float(input("q= "))
    fr_corr_type=float(input("corr= "))
else:
    use_fr_model=int(use_fr_model)
    fr_q=1

    


    


if use_fr_model==1:
    print("fractional-order model. Q="+str(fr_q)+", tau="+str(fr_tau)+", corr_type="+str(fr_corr_type)) 
if use_fr_model==-1:
    print("standard RothC scheme")
if use_fr_model==0:
    print("non standard RothC scheme")
# Vsevolod: end of added part of code

print("current dir is: %s" % (os.getcwd()))
# Vsevolod: changed to current dir. 
#	Also all '\' slashes were changed to '/' and file/dir names made case sensitive to comply with unix rules
my_path = Path('./')

folder_out='OUTPUT_'+outname
path=my_path / folder_out

if not path.exists():
    path.mkdir()

# Read files
datatab   =  pd.read_excel(dataname, header=None)
partab = pd.read_excel(parname, header=None)

# Global parameters
k       = partab[2][0:4].values
r       = partab[2][4]
eta     = partab[2][5]
clay    = partab[2][6]
d       = partab[2][7]
meantemp=partab[2][8]

# Initial SOC content
C0      = partab[2][9:13].values
# C0[0]=C0[3]*0.004919
# C0[1]=C0[3]*0.143926
# C0[2]=C0[3]*0.021407
# C0[3]=C0[3]*0.829748
IOM     = partab[2][13];

# Read data
years   = datatab[0].values[1:].astype(int)
months  = datatab[1].values[1:].astype(str)
temp    = datatab[2].values[1:].astype(float)
rain    = datatab[3].values[1:].astype(float)
pet     = datatab[4].values[1:].astype(float)
g       = datatab[5].values[1:].astype(float)
f       = datatab[6].values[1:].astype(float)
cover   = datatab[7].values[1:].astype(float)
try:
    to_compare = datatab[9].values[1:].astype(float)
except:
    pass


unique_y=years#cris add

# Settings
gamma= r/(r+1)
avecg   = np.array([gamma, 1-gamma, 0, 0]).reshape(1,-1)
avecf   = np.array([eta, eta, 0, 1-2*eta]).reshape(1,-1)

b         = np.dot(g.reshape(g.size,1),avecg) + np.dot(f.reshape(f.size,1),avecf)
#b  = np.dot(f.reshape(f.size,1),avecf)
b         = b.T
meantemp=temp[0]

x = 1.67 * ( 1.85 + 1.6* np.exp(-0.0786*clay) )
biohumfac = 1/(x+1)
alpha= 0.46*biohumfac
beta = 0.54*biohumfac

# Time parameters
t0=1
T =len(years)
tspan   = np.array([t0,T]).reshape(1,-1)
dt=1 #month

# Rate modifying factors
[ka,kb,kc,acc] = rhofun(temp,meantemp, rain,pet,cover,clay,d);
#kb=np.ones(len(kb)).reshape(-1,1)#aggiunto da me
rho=ka*kb*kc

rho_table = pd.DataFrame({'years': years, 'months': months, 'ka':ka.reshape(-1,), 'kb':kb.reshape(-1,), 'kc':kc.reshape(-1,), 'acc TSMD': acc.reshape(-1,)})
#rho_table.to_excel('rho_table.xlsx', index=False)
file_out='OUTPUT_%s/RHO_%s.xlsx'%(outname,outname)
rho_table.to_excel(file_out, index=False)

# Run
#-------------------------------------
Lam=np.array([(0, 0, 0, 0),(0, 0, 0, 0), (alpha, alpha, alpha, alpha), (beta, beta, beta, beta)])
N=len(rho)-1
# Vsevolod: branch to three algorithms - non-standard scheme (use_fr_model==0), Crank-Nicholson for fractional-order model (use_fr_model==1)
#           and standard scheme (use_fr_model==-1) 
t1=time.time()
if use_fr_model==0: 
    A=-np.dot((np.eye(4)-Lam),np.diag(k))
    M=np.eye(4)-Lam
    Minv=np.linalg.inv(M)
    At=-(np.dot(np.dot(M,np.diag(k)),Minv))
    Atinv=-np.dot(np.dot(M,np.diag(1./k)),Minv)
    
    Cout=[C0]
    tout=[t0] 
    #F=np.zeros([4,4,N])
    for n in range(N):
        if rho[n]==0:
            C0=C0 + b[:,n]
        else:
            F= np.dot(Atinv, ((sp.linalg.expm(dt* rho[n]*At) -np.eye(4)) / rho[n]))
            C0=C0 + np.dot(F,(rho[n]*np.dot(A,C0)+b[:,n]) )
        Cout.append(C0)
    #t0=tspan[0][0] 
    #tend=tspan[0][1] 
    
    
    # while  t0 <= tend-N*dt:  
    #     for n in range(N):
    #         C0=C0 + np.dot(F[:,:,n],(rho[n]*np.dot(A,C0)+b[:,n]) )
    #         t0=t0+dt
    #         Cout.append(C0)
    #         tout.append(t0)
	
# call external fr_rothC solver
if use_fr_model==1:
    [tout,Cout]=fr_rothC(alpha,beta,12*k,rho,12*b,C0,years[0],T,fr_q,fr_tau,fr_corr_type)

# standard RothC scheme
# c(n+1)=(Lam+(I-Lam)exp(-dt*rho(n)*diag(k))c(n)+dt*b(n)
if use_fr_model==-1:
    Cout=[C0]
    tout=[t0] 
    for n in range(N):
        C0=np.dot(Lam+ np.dot(np.eye(4)-Lam, sp.linalg.expm(-dt* rho[n]*np.diag(k))),C0)+dt*b[:,n]
        Cout.append(C0)
t2=time.time()
print("execution time - "+str(t2-t1)+"s")

tout=np.stack(tout)
Cout=np.stack(Cout)
#NT=len(tout)
#np=(NT-1)/N

#----------------------------
#[tout  ,Cout]  = NS(tspan,dt,C0,alpha,beta,k,rho,b,0);  
soc  = np.sum(Cout,1) +IOM;

#cris add
averagey=[]
for j in np.unique(unique_y):
    averagey.append(soc[np.where(years==j)[0]].mean())
averagey=np.stack(averagey)
table_mean_soc=pd.DataFrame({'years':np.unique(unique_y), 'mean_soc':averagey})
#end cris add


# Save table
out_table = pd.DataFrame({'years': years, 'months': months, 'DPM':Cout[:,0],'RPM': Cout[:,1],'BIO':Cout[:,2], 'HUM':Cout[:,3],'IOM':IOM*np.ones(T), 'SOCa': soc})

file_out='OUTPUT_%s/SOC_%s.xlsx'%(outname,outname)
out_table.to_excel(file_out,index=False)

#begin cris----> add part to write dataframe with different model and parameter
soc_type='SOC_%s_%s_%s'%(use_fr_model,fr_q,fr_corr_type)

compare = 'OUTPUT_%s/compare_%s.xlsx'%(outname,outname)
if exists(compare):
    comp_table = pd.read_excel(compare, header=0)
    comp_table[soc_type]=soc
    comp_table[soc_type].astype('float')
    comp_table.to_excel(compare,index=False)
else:
    
    comp_table = pd.DataFrame({'years': years, 'months': months, 'to_compare':to_compare})
    comp_table[soc_type]=soc
    comp_table[soc_type].astype('float')
    comp_table.to_excel(compare,index=False)

#end cris -----------------------------------------

### Plot SOC

years=years.astype(float)
for i in range(len(years)):
	years[i]=float(years[i])+float(i%12)/12.0
plt.plot(years,soc)
try:
    plt.plot(years,to_compare,'bo')
except:
    pass
plt.xlabel('Year')
plt.ylabel('SOC')
plt.title(('%s'%outname))
#plt.ylim(120, 160)
plt.show()

#begin cris plot with mean
# soc_series=pd.Series(soc)
# soc_mean=soc_series.rolling(6, min_periods=1).mean()
# soc_mean=soc_mean.values
# comp_table[soc_type+'mean']=soc_mean
# comp_table.to_excel(compare,index=False)
# # maxx=np.where(soc_mean==np.max(soc_mean))[0]
# # soc_mean[maxx]=np.mean(soc_mean[0:11])
# plt.plot(years,soc_mean)
# try:
#     plt.plot(years,to_compare,'bo')
# except:
#     pass
# plt.xlabel('Year')
# plt.ylabel('SOC')
# plt.title(('%s'%outname))

# plt.show()
#\end cris

#plt.savefig( ('OUTPUT_%s/SOC_%s.png'%(outname,outname)) )

print(('Simulation ended with success. Results have been saved in the folder %s' %(outname)))

#################cris######################
#compare polonia
# years=years[years<2016]
# ind=np.where(years<2016)[0]
# measure=np.where(~np.isnan(to_compare[ind]))[0]
# obs=to_compare[measure]

# predict=soc[ind][measure]
# r2=r2_score(obs, predict)
# rmse=np.sqrt(mean_squared_error(obs, predict))
# print('RMSE:', rmse)
# print('r2:', r2)


#compare ucraina

measure=np.where(~np.isnan(to_compare))[0]
obs=to_compare[measure]

predict=soc[measure]
r2=r2_score(obs, predict)
rmse=np.sqrt(mean_squared_error(obs, predict))
print('RMSE:', rmse)
print('r2:', r2)

#compare with mean

measure=np.where(~np.isnan(to_compare))[0]
obs=to_compare[measure]

aa=(datatab[0].values[1:][~np.isnan(to_compare)]).astype(float)
predict=[]
for h in aa:
    predict.append(table_mean_soc.mean_soc[np.where(table_mean_soc.years==h)[0]])
    

predict=np.stack(predict)
r2=r2_score(obs, predict)
rmse=np.sqrt(mean_squared_error(obs, predict))
print('RMSE:', rmse)
print('r2:', r2)
############end cris#####################


