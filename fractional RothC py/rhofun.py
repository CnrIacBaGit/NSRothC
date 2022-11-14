# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 07:03:33 2022

@author: Cristiano
"""
import numpy as np

def rhofun(temp, meantemp, rain,pet,cover,clay,d):

    n=len(temp)
    
    # the rate modifying factor for temperature
    a = (47.91/ ( 1+ np.exp( 106.06 / (temp-meantemp+106.6/np.log(46.91)) ) )).reshape(-1,1)
    
    # the rate modifying factor for moisture
    Mcover = -(20 + 1.3*clay-0.01*clay**2) * d/23
    Mbare= Mcover/1.8
    vec=rain-pet
    accTSMD=np.zeros(n).reshape(-1,1)
    i=0
    while vec[i]>0:
        i=i+1
        
    for j in range(i,n):
        if cover[j]==0:
            accTSMD[j]= min( max( accTSMD[j-1]+vec[j] , Mcover ), 0)
        else:
            accTSMD[j]= min( max( accTSMD[j-1]+vec[j] , Mbare ), 0)
        
    
    
    b=np.ones(n).reshape(-1,1)
    indc =np.intersect1d(np.where((cover==0))[0], np.where((accTSMD<=0.444*Mcover))[0])
    indb = np.intersect1d(np.where(cover==1)[0], np.where( accTSMD<=Mbare)[0]) #0.444 ???
    b[indc] = 0.2+ (1-0.2)* (Mcover-accTSMD[indc])/(Mcover-0.444*Mcover)
    b[indb] = 0.2+ (1-0.2)* (Mbare-accTSMD[indb])/(Mbare-0.444*Mbare)
    
    c=np.ones(n).reshape(-1,1)
    c[np.where(cover==1)]=0.6

    return [a,b,c,accTSMD]