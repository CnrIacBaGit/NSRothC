# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:23:45 2022

@author: Cristiano
"""


def NS(tspan,dt,C0,alpha,beta,k,rho,b,fig):
    import numpy as np
    import scipy as sp

    Lam=np.array([(0, 0, 0, 0),(0, 0, 0, 0), (alpha, alpha, alpha, alpha), (beta, beta, beta, beta)])
    A=-np.dot((np.eye(4)-Lam),np.diag(k))
    M=np.eye(4)-Lam
    Minv=np.linalg.inv(M)
    At=-(np.dot(np.dot(M,np.diag(k)),Minv))
    Atinv=-np.dot(np.dot(M,np.diag(1./k)),Minv)
    N=len(rho)-1
    
    
    F=np.zeros([4,4,N])
    for n in range(N):
        if rho[n]==0:
            F[:,:,n] = np.eye(4)
        else:
            F[:,:,n]= np.dot(Atinv, ((sp.linalg.expm(dt* rho[n]*At) -np.eye(4)) / rho[n]))
        
    t0=tspan[0][0] 
    tend=tspan[0][1] 
    Cout=[C0]
    tout=[t0] 
    
    while  t0 <= tend-N*dt:  
        for n in range(N):
            C0=C0 + np.dot(F[:,:,n],(rho[n]*np.dot(A,C0)+b[:,n]) )
            t0=t0+dt
            Cout.append(C0)
            tout.append(t0)
       
       
    tout=np.stack(tout)
    Cout=np.stack(Cout)
    
    NT=len(tout)
    np=(NT-1)/N
    
    '''
    
    if fig==1   
        #fig, axs = plt.subplots(2)
        figure()
        subplot(2,2,1)
        plot(tout,Cout(:,1),'g','LineWidth',2)
        title('DPM')
        xlim(tspan)
        xticks(tout(1:N:end))
        xticklabels(split(num2str(0:np)))
        xlabel('Year')
        subplot(2,2,2)
        plot(tout,Cout(:,2),'g','LineWidth',2)
        title('RPM')
        xlim(tspan)
        xticks(tout(1:N:end))
        xticklabels(split(num2str(0:np)))
        xlabel('Year')
        subplot(2,2,3)
        plot(tout,Cout(:,3),'g','LineWidth',2)
        title('BIO')
        xlim(tspan)
        xticks(tout(1:N:end))
        xticklabels(split(num2str(0:np)))
        xlabel('Year')
        subplot(2,2,4)
        plot(tout,Cout(:,4),'g','LineWidth',2)
        title('HUM')
        xlim(tspan)
        xticks(tout(1:N:end))
        xticklabels(split(num2str(0:np)))
        xlabel('Year')
    '''
     
    return [tout,Cout]