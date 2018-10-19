import sys

import _numpypy as np
from math import *
#from numarsy import *
#import numarsy.linear_algebs as la
import copy

#import sndom
from math import sqrt,exp
import random
import time

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
p=104

num=sys.argv[1]	#links observed

#time={}
#for t in range(T):
#	time[t]=[]
fh=open('/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Format_Tonyi/Format_Toñi/N117-Enginyeria_Mecanica-2010_104/linkstrain0.dat','r')
#fh=open('/export/home/shared/Temp/train_links.dat','r')
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

fh=open('/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Format_Tonyi/Format_Toñi/N117-Enginyeria_Mecanica-2010_104/nolinkstrain0.dat','r')
#fh=open('/export/home/shared/Temp/train_nolinks.dat','r')
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
#fout3=open('likelinks2000.dat','w')
############################################
chances=1
for c in range(chances):
	print c
	theta={}
	ntheta={}
	for i in range(p):
		for j in range(p-(i+1)):
			z=i+j+1
			vec=[]
			for k in range(K):
				vec.append(random.random())
			theta[(i,z)]=vec
			theta[(z,i)]=vec
			ntheta[(i,z)]=[0.0]*K
			ntheta[(z,i)]=[0.0]*K
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
	for i in range(K):
		b=[]
		for j in range(S):
			vec.append(random.random())
		pr.append(vec)
		npr.append([0.0]*S)
		nonpr.append([0.0]*S)
	'''
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
	'''
	#Normalizations:
	for i in range(p):
		for j in range(p-(i+1)):
			z=i+j+1
			D=0.
			for k in range(K):
				D=D+theta[(i,z)][k]
			for k in range(K):
				theta[(i,z)][k]=theta[(i,z)][k]/D
				theta[(z,i)][k]=theta[(i,z)][k]

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
	#time0=time.time()
	for g in range(Runs):
		#print g
		for e in links:
			t=e[0]
			n=(e[1],e[2])
			D=0.	
			for s in range(S):
				for k in range(K):
					D=D+theta[n][k]*tau[t][s]*pr[k][s]
			for s in range(S):
				for k in range(K):
					a=(theta[n][k]*tau[t][s]*pr[k][s])/D
					ntheta[n][k]=ntheta[n][k]+a
					ntau[t][s]=ntau[t][s]+a
					npr[k][s]=npr[k][s]+a
			
		for e in nolinks:
			t=e[0]
			n=(e[1],e[2])
			D=0.	
			for s in range(S):
				for k in range(K):
					D=D+theta[n][k]*tau[t][s]*(1-pr[k][s])
			for s in range(S): 
				for k in range(K):
					a=(theta[n][k]*tau[t][s]*(1-pr[k][s]))/D
					ntheta[n][k]=ntheta[n][k]+a
					ntau[t][s]=ntau[t][s]+a
					nonpr[k][s]=nonpr[k][s]+a
			
		#print 'nolink', nonpr[0][0][0], npr[0][0][0]
		#Normalizations:
		#for z in range(L-(k)):
		#l=k+z
		err=0.
		for i in range(p):
			for j in range(p-(i+1)):
				z=i+j+1
				D=0.
				for k in range(K):
					D=D+ntheta[(i,z)][k]+ntheta[(z,i)][k]
				for k in range(K):
					ntheta[(i,z)][k]=(ntheta[(i,z)][k]+ntheta[(z,i)][k])/(D+0.00000000001)
					ntheta[(z,i)][k]=ntheta[(i,z)][k]
		#print 'theta',sum(ntheta[(1,2)]), sum(ntheta[(2,1)]), ntheta[(1,2)], ntheta[(2,1)]
		for t in range(T):
			D=0.
			for s in range(S):
				D=D+ntau[t][s]
			for s in range(S):
				ntau[t][s]=ntau[t][s]/(D+0.0000000000001)
		
		for s in range(S):
			for k in range(K):
				npr[k][s]=npr[k][s]/(nonpr[k][s]+npr[k][s])
		
		#print npr[0][0]
		theta=copy.copy(ntheta)
		tau=copy.copy(ntau)	
		pr=copy.copy(npr)
		for i in range(p):
			for j in range(p-(i+1)):
				z=i+j+1
				ntheta[(i,z)]=[0.]*K
				ntheta[(z,i)]=[0.]*K
		for t in range(T):
			ntau[t]=[0.]*S
		for k in range(K):
			npr[k]=[0.]*S
			nonpr[k]=[0.]*S
		#print theta[0], pr[0][0]
	#time1=time.time()
	#print time1-time0, (time1-time0)/20.
	
	Like=0.
		
	for e in links:
		t=e[0]
		n=(e[1],e[2])
		D=0.
		for s in range(S):
				for k in range(K):
					D=D+theta[n][k]*tau[t][s]*pr[k][s]
		Like=Like+log(D)
	
	for e in nolinks:
		t=e[0]
		n=(e[1],e[2])
		D=0.
		for s in range(S):
			for k in range(K):
				D=D+theta[n][k]*tau[t][s]*(1-pr[k][s])
		Like=Like+log(D)
	print Like

	#fout3.write('%s %s\n' % (g,Like))
	
	#fout3.close()
	#print g, Like
	#fout2=open('/export/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Rogerresults.dat','w')
	fout2=open('/home/shared/Projects/Agodoy/Bipartite/RealTimeinference/Paper/Format_Tonyi/Format_Toñi/N117-Enginyeria_Mecanica-2010_104/Test0/'+str(K)+'_'+str(S)+'linksresults'+str(num)+'.dat','w')
	fout2.write('%s\n' % Like)
	for i in range(p):
		for j in range(p-(i+1)):
			z=i+j+1
			for kk in range(K):
				fout2.write('%s ' % (theta[(i,z)][kk]))
			fout2.write('\n')
			#,theta[(i,z)][2],theta[(i,z)][3],theta[(i,z)][4]))
			#theta[i][6],theta[i][7],theta[i][8],theta[i][9]))
	for t in range(T):
		for ss in range(S):
			fout2.write('%s ' % (tau[t][ss]))
		fout2.write('\n')
		#,tau[t][2],tau[t][3],tau[t][4]))
		#eta[j][5], eta[j][6],eta[j][7],eta[j][8],eta[j][9]))
	#for s in range(S):
	for k in range(K):
			for s in range(S):
				fout2.write('%s ' % (pr[k][s]))
			fout2.write('\n')
			#,pr[k][2],pr[k][3],pr[k][4]))
	fout2.close()
