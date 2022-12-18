#!/usr/bin/python
import math
import functools
from functools import partial

# correction factor (left-hand size - corr(t)*Da C
def corr_factor(t,q,corr_type):
	if (corr_type==1):
		return pow(t,q-1)
	return 1
# out=mult*rho(t)*A*in
def mmult(t,_in,out,k,alpha,beta,rho,corr_type,q,t0):
	month=int(t*12.0)
	if month>=len(rho):
		month=len(rho)-1
	out[0]=-k[0]*_in[0]
	out[1]=-k[1]*_in[1]
	out[2]=alpha*(k[0]*_in[0]+k[1]*_in[1]+k[3]*_in[3])+(alpha-1.0)*k[2]*_in[2]
	out[3]=beta*(k[0]*_in[0]+k[1]*_in[1]+k[2]*_in[2])+(beta-1.0)*k[3]*_in[3]
	for i in range(4):
		out[i]=out[i]*(rho[month][0]/corr_factor(t+t0,q,corr_type))
# Crank-Nicholson finite difference with L1 Caputo derivative approximation
# for step j+1, in - c(t_j), old_c - prev solutions
ds={}
def d(j,k,q):
	global ds
# read from dictionary if present
	if ds.get(j-k)!=None:
		return ds[j-k]
	if (j==k):
		 return 1.0
	val=pow((float)(j-k+1),1.0-q)-pow((float)(j-k),1.0-q)
# save in dictionary
	ds[j-k]=val
	return val
def scheme2_one_step(t,tau,j,q,_in,out,idx,alpha,beta,rho,k,b,fr_corr_type,t0,old_c):
	month=int(t*12.0)
	month_tau=int((t+tau)*12.0)
	if q!=1.0:
		mult=1.0/(math.gamma(2.0-q)*pow(tau,q))
	else:
		mult=1.0/tau
	# vb1 - for right-hand side, vb2,vb3 - auxiliary, lM - matrix in the left-hand size
	vb1=[0,0,0,0]
	vb2=[0,0,0,0]
	vb3=[0,0,0,0]
	lM=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
	# right-hand side vector
	# b
	if month_tau>=len(b[0,:]):
		month_tau=len(b[0,:])-1
	if month>=len(b[0,:]):
		month=len(b[0,:])-1
	for i in range(4):
		vb1[i]=b[i,month]/corr_factor(t+t0,q,fr_corr_type)
		vb2[i]=b[i,month_tau]/corr_factor(t+t0,q,fr_corr_type)
		vb1[i]=0.5*(vb1[i]+vb2[i])
	# Mc(t_j)
	mmult(t,_in,vb2,k,alpha,beta,rho,fr_corr_type,q,t0)
	for i in range(4):
		vb1[i]=vb1[i]+0.5*vb2[i]
	# c_j part from left-hand size
	for i in range(4):
		vb1[i]=vb1[i]+mult*d(j,j,q)*_in[i]
	# non-local part from left-hand side
	if q!=1.0:
		for i in range(4):
# simple
			for kk in range(j):
				vb1[i]=vb1[i]-(mult*d(j,kk,q)*(old_c[4*(kk+1)+i]-old_c[4*kk+i]))
# through map/reduce
#			if j!=0:
#				vb1[i]=vb1[i]-mult*functools.reduce(lambda a,b:a+b,map(lambda k:d(j,k,q)*(old_c[4*(k+1)+i]-old_c[4*k+i]),idx[:j]))
	# left-hand side matrix
	# M(t_j+1)
	for kk in range(4):
		vb2[0]=vb2[1]=vb2[2]=vb2[3]=0
		vb2[kk]=1.0
		mmult(t+tau,vb2,vb3,k,alpha,beta,rho,fr_corr_type,q,t0)
		for i in range(4):
			lM[i][kk]=-0.5*vb3[i]
	# diagonal part
	for i in range(4):
		lM[i][i]=lM[i][i]+mult*d(j,j,q)
	# solve lM * out = vb1 (simple Gaussian elimination)
	for h in range(4):
		for i in range(4):
			if i>=h+1:
				f = lM[i][h] / lM[h][h]
				lM[i][h]= 0
				for j in range(4):
					if j>=h+1:
						lM[i][j]= lM[i][j]-lM[h][j] * f
				vb1[i]=vb1[i]-vb1[h]*f
	out[3]=vb1[3]/lM[3][3]
	out[2]=(vb1[2]-out[3]*lM[2][3])/lM[2][2]
	out[1]=(vb1[1]-out[3]*lM[1][3]-out[2]*lM[1][2])/lM[1][1]
	out[0]=(vb1[0]-out[3]*lM[0][3]-out[2]*lM[0][2]-out[1]*lM[0][1])/lM[0][0]
# ////////////////// main ///////////////////
def fr_rothC(alpha,beta,k,rho,b,C0,t0,tend,fr_q,fr_tau,fr_corr_type):
	old_c=[]
	idx=[]
	_in=[0,0,0,0]
	out=[0,0,0,0]
	tout=[]
	cout=[]
	month=0
	t=0 # start from 0, not t0!
	n_steps=int((((tend)/12.0)+1.0)/fr_tau)
	for i in range(4*(n_steps+1)):
		old_c.append(0.0);
	for i in range(n_steps+1):
		idx.append(i)
	for i in range(4):
		_in[i]=old_c[i]=C0[i]
	tout.append(t)
	cout.append([C0[0],C0[1],C0[2],C0[3]])
	for i in range(n_steps):
		scheme2_one_step(t,fr_tau,i,fr_q,_in,out,idx,alpha,beta,rho,k,b,fr_corr_type,t0,old_c)
		for _k in range(4):
			_in[_k]=old_c[4*(i+1)+_k]=out[_k]
		m=int(t*12)
		if m!=month:
			tout.append(t)
			cout.append([out[0],out[1],out[2],out[3]])
			month=m
#			print(out)
			if len(cout)==tend:
				break
		t=t+fr_tau
	return [tout,cout]
