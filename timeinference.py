import sys

import _numpypy as np
from math import *
#from numarsy import *
#import numarsy.linear_algebs as la
import copy

#import sndom
from math import sqrt,exp
import random

def maximum(vec):
	m1=0				
	for i in range(len(vec)):
		m1=max(max(vec[i]),m1)	
	return m1

##############################################################################
#	this progsm geneste a bipartite network, using networkX
##############################################################################


Number=5
K=Number
L=Number
S=Number
#G=5
T=365 #dies sep OCt 2010
p=118

num=sys.argv[1]	#links observed

#time={}
#for t in range(T):
#	time[t]=[]
fh=open('/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Format_Tonyi/Format_Toñi/N115-Enginyeria_Electronica,_Electrica_i_Automatica-2010_118/linkstrain0.dat','r')
#N107-Quimica_Fisica_i_Inorganica-2010_116
#N117-Enginyeria_Mecanica-2010_104
#fout=open('FilMEC10.dat','w')
igot=fh.readlines()
links=[]
for line in igot:
	about = line.strip().split(' ')
	#if int(about[1])!=int(about[2]):
	#	#if (int(about[1]),int(about[2])) not in time[int(about[0])]:
	links.append((int(about[0]),int(about[1]),int(about[2])))
			#time[int(about[0])].append((int(about[1]),int(about[2])))
			#time[int(about[0])].append((int(about[2]),int(about[1])))
			#fout.write('%s' % line)
#fout.close()

fh.close()

fh=open('/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Format_Tonyi/Format_Toñi/N115-Enginyeria_Electronica,_Electrica_i_Automatica-2010_118/nolinkstrain0.dat','r')
#fout=open('FilMEC10.dat','w')
igot=fh.readlines()
nolinks=[]
for line in igot:
	about = line.strip().split(' ')
	#if int(about[1])!=int(about[2]):
	#	#if (int(about[1]),int(about[2])) not in time[int(about[0])]:
	nolinks.append((int(about[0]),int(about[1]),int(about[2])))
			#time[int(about[0])].append((int(about[1]),int(about[2])))
			#time[int(about[0])].append((int(about[2]),int(about[1])))
			#fout.write('%s' % line)
#fout.close()

fh.close()

'''
#fout=open('nolinksFilMEC10.dat','w')
c=0
nc=0
nolinks=[]
for t in range(T):
	for i in range(p):
		for j in range(p-(i+1)):
			k=i+(j+1)
			b=(i,k)
			#b2=(k,i)
			if b not in time[t]:
				nc=nc+1
				nolinks.append((t,i,k))
				#fout.write('%s %s %s\n' % (t,i,k))
			else:
				c=c+1
print 'links',c,'nolinks', nc, c+nc, 'links list', len(links), len(list(set(links)))
#fout.close()
'''
############################################

sampling=1	

for c in range(sampling):
	print c
	theta=[]
	ntheta=[]
	for i in range(p):
		vec=[]
		for k in range(K):
			vec.append(random.random())
		theta.append(vec)
		ntheta.append([0.0]*K)
	#theta=np.arsy(theta)
	tau=[]
	ntau=[]
	for t in range(T):
		vec=[]
		for s in range(S):
			vec.append(random.random())
		tau.append(vec)
		ntau.append([0.0]*S)
	#tau=np.arsy(tau)
	pr=[]
	npr=[]
	nonpr=[]
	for k in range(K):
		pr.append([])
		npr.append([])
		nonpr.append([])
		for l in range(L):
			pr[k].append([])
			npr[k].append([])
			nonpr[k].append([])
	for k in range(K):
	    for z in range(L-(k)):
		l=k+z
		if k == l:
			vec=[]
			for s in range(S):
				vec.append(random.random())
		    	pr[k][l] = vec
			npr[k][l]=[0.]*S
			nonpr[k][l]=[0.]*S
		else:
			vec=[]
			for s in range(S):
				vec.append(random.random())
			pr[k][l]=vec
			pr[l][k]=vec
			npr[k][l]=[0.]*S
			npr[l][k]=[0.]*S
			nonpr[k][l]=[0.]*S
			nonpr[l][k]=[0.]*S
	
	#Normalizations:
	for i in range(p):
		D=0.
		for k in range(K):
			D=D+theta[i][k]
		for k in range(K):
			theta[i][k]=theta[i][k]/D

	for t in range(T):
		D=0.
		for s in range(S):
			D=D+tau[t][s]
		for s in range(S):
			tau[t][s]=tau[t][s]/D

#########################################################################################
	Runs=500
	Lant=-100000000000
	cnt=0
	count=0
	ver=0
	#fout=open('likeli'+str(Runs)+'.dat','w')
	for g in range(Runs):
		#print g
		for e in links:
			t=e[0]
			n1=e[1]
			n2=e[2]
			#s=linksr[n]
			D=0.	
			for s in range(S):
				for l in range(L):
					for k in range(K):
						D=D+theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s]
			for s in range(S):
				for l in range(L):
					for k in range(K):
						a=(theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s])/D
						ntheta[n1][k]=ntheta[n1][k]+a
						ntheta[n2][l]=ntheta[n2][l]+a
						ntau[t][s]=ntau[t][s]+a
						npr[k][l][s]=npr[k][l][s]+a
		for e in nolinks:
			t=e[0]
			n1=e[1]
			n2=e[2]
			D=0.	
			for s in range(S):
				for l in range(L):
					for k in range(K):
						D=D+theta[n1][k]*theta[n2][l]*tau[t][s]*(1-pr[k][l][s])
			for s in range(S): 
				for l in range(L):
					for k in range(K):
						a=(theta[n1][k]*theta[n2][l]*tau[t][s]*(1-pr[k][l][s]))/D
						ntheta[n1][k]=ntheta[n1][k]+a
						ntheta[n2][l]=ntheta[n2][l]+a
						ntau[t][s]=ntau[t][s]+a
						nonpr[k][l][s]=nonpr[k][l][s]+a
		#print 'nolink', nonpr[0][0][0], npr[0][0][0]
		#Normalizations:
		err=0.
		for i in range(p):
			D=0.
			for k in range(K):
				D=D+ntheta[i][k]
			for k in range(K):
				ntheta[i][k]=ntheta[i][k]/D

		for t in range(T):
			D=0.
			for s in range(S):
				D=D+ntau[t][s]
			for s in range(S):
				ntau[t][s]=ntau[t][s]/D
		for s in range(S):
			for l in range(L):
				for k in range(K):
					npr[k][l][s]=npr[k][l][s]/(nonpr[k][l][s]+npr[k][l][s])

		theta=copy.copy(ntheta)
		tau=copy.copy(ntau)	
		for k in range(K):
			for l in range(L):
				pr[k][l]=npr[k][l]
		for i in range(p):
			ntheta[i]=[0.]*K
		for t in range(T):
			ntau[t]=[0.]*S
		for k in range(K):
			for l in range(L):
				npr[k][l]=[0.]*S
				nonpr[k][l]=[0.]*S
		#print theta[0], pr[0][0]
	Like=0.
	for e in links:
		t=e[0]
		n1=e[1]
		n2=e[2]
		D=0.
		for s in range(S):
			for l in range(L):
				for k in range(K):
					D=D+theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s]
		Like=Like+log(D)
	for e in nolinks:
		t=e[0]
		n1=e[1]
		n2=e[2]		
		D=0.
		for s in range(S):
			for l in range(L):
				for k in range(K):
					D=D+theta[n1][k]*theta[n2][l]*tau[t][s]*(1-pr[k][l][s])
		Like=Like+log(D)
	print Like
	#fout.write('%s %s\n' % (g,Like))
	#fout.close()
	fout2=open('/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Format_Tonyi/Format_Toñi/N115-Enginyeria_Electronica,_Electrica_i_Automatica-2010_118/Test0/'+str(K)+'_'+str(S)+'nodesresults'+str(num)+'.dat','w')
	fout2.write('%s\n' % Like)
	for i in range(p):
		for kk in range(K):
			fout2.write('%s ' % (theta[i][kk]))
		fout2.write('\n')
		#theta[i][6],theta[i][7],theta[i][8],theta[i][9]))
	for t in range(T):
		for ss in range(S):
			fout2.write('%s ' % (tau[t][ss]))
		fout2.write('\n')
		#,tau[t][2],tau[t][3],tau[t][4]))
		#eta[j][5], eta[j][6],eta[j][7],eta[j][8],eta[j][9]))
	#for s in range(S):
	for k in range(K):
		for l in range(K):
			for s in range(S):
				fout2.write('%s ' % (pr[k][l][s]))
			fout2.write('\n')
			#,pr[k][l][2],pr[k][l][3],pr[k][l][4]))
	fout2.close()

		
'''
#k num of pasmeters:
k=p*K+m*L
BIC=-2*Lold+k*log(len(B.edges()))
print 'best likelihood', Lold, 'num topics', -2*Lold, k*log(len(B.edges())),k, len(B.edges())
print 'BIC?', BIC
#print 'theta', theta
'''
