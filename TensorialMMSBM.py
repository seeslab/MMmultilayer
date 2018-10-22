import sys

import _numpypy as np
from math import *
#from numarsy import *
#import numarsy.linear_algebs as la
import copy


#import sndom
from math import sqrt,exp
import random


train=sys.argv[1] #train
test=sys.argv[2] #test
p=int(sys.argv[3]) #nodes drugs69
T=int(sys.argv[4]) #layers drugs85
K=int(sys.argv[5]) #groups of nodes
L=K
S=int(sys.argv[6]) #groups of layers
R=int(sys.argv[7]) #different labels
sampling=int(sys.argv[8]) #different initializations
printp=int(sys.argv[9]) #0/1 0noprint 1 print params
############################################

fh=open(train,'r')
igot=fh.readlines()
trainn=[]

for line in igot:
	about = line.strip().split(' ')
	trainn.append((int(about[0]),int(about[1]),int(about[2]),int(about[3])))
fh.close()

fh2=open(test,'r')
igot2=fh2.readlines()
testt={}
praij={}
for line in igot2:
	about = line.strip().split(' ')
	testt[(int(about[0]),int(about[1]),int(about[2]))]=int(about[3])
	praij[(int(about[0]),int(about[1]),int(about[2]))]=[0.]*R
fh2.close()
###################################################
for w in range(sampling):
	theta=[]
	ntheta=[]
	for i in range(p):
		vec=[]
		for k in range(K):
			vec.append(random.random())
		theta.append(vec)
		ntheta.append([0.0]*K)
	tau=[]
	ntau=[]
	for t in range(T):
		vec=[]
		for s in range(S):
			vec.append(random.random())
		tau.append(vec)
		ntau.append([0.0]*S)
	pr=[]
	npr=[]
	for k in range(K):
		pr.append([])
		npr.append([])
		for l in range(L):
			pr[k].append([])
			npr[k].append([])
	for k in range(K):
		for z in range(L-(k)):
			l=k+z
			if k == l:
				vec=[]
				b=[]
				for s in range(S):
					a=[]
					for r in range(R):
						a.append(random.random())
					vec.append(a)
					b.append([0.]*R)
				pr[k][l] = vec
				npr[k][l]=b
			else:
				vec=[]
				b=[]
				for s in range(S):
					a=[]
					for r in range(R):
						a.append(random.random())
					vec.append(a)
					b.append([0.]*R)
				pr[k][l]=vec
				pr[l][k]=vec
				npr[k][l]=b
				npr[l][k]=b

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

	for k in range(K):
		for l in range(L):
			for a in range(S):
				D=0.
				for r in range(R):
					D=D+pr[k][l][s][r]
				for r in range(R):
					pr[k][l][s][r]=pr[k][l][s][r]/D

	#########################################################################################
	Runs=1000
	#Al=1.
	#While Al>0.00000001:
	for g in range(Runs):
		for e in trainn:
			t=int(e[0])
			n1=int(e[1])
			n2=int(e[2])
			ra=int(e[3])
			D=0.	
			for s in range(S):
				for l in range(L):
					for k in range(K):
						D=D+theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s][ra]
			for s in range(S):
				for l in range(L):
					for k in range(K):
						a=(theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s][ra])/D
						ntheta[n1][k]=ntheta[n1][k]+a
						ntheta[n2][l]=ntheta[n2][l]+a
						ntau[t][s]=ntau[t][s]+a
						npr[k][l][s][ra]=npr[k][l][s][ra]+a

		#Normalizations:
		err=0.
		for i in range(p):
			D=0.
			for k in range(K):
				D=D+ntheta[i][k]
			for k in range(K):
				ntheta[i][k]=ntheta[i][k]/(D+0.00000000001)

		for t in range(T):
			D=0.
			for s in range(S):
				D=D+ntau[t][s]
			for s in range(S):
				ntau[t][s]=ntau[t][s]/(D+0.0000000000001)

		for k in range(K):
			for l in range(L):
				for s in range(S):
					D=0.
					for r in range(R):
						D=D+npr[k][l][s][r]
					for r in range(R):
						npr[k][l][s][r]=npr[k][l][s][r]/(D+0.000000001)

		theta=copy.copy(ntheta)
		tau=copy.copy(ntau)
		for k in range(K):
			for l in range(L):
				for s in range(S):
					pr[k][l][s]=npr[k][l][s]
		for i in range(p):
			ntheta[i]=[0.]*K
		for t in range(T):
			ntau[t]=[0.]*S
		for k in range(K):
			for l in range(L):
				for s in range(S):
					npr[k][l][s]=[0.]*R

	Like=0.
	for e in trainn:
		t=int(e[0])
		n1=int(e[1])
		n2=int(e[2])
		ra=int(e[3])
		D=0.
		for s in range(S):
			for l in range(L):
				for k in range(K):
					D=D+theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s][ra]
		Like=Like+log(D)
		
	#scores:
	for e in praij.keys():
		t=int(e[0])
		n1=int(e[1])
		n2=int(e[2])
		for rr in range(R):
			pra=0.
			for s in range(S):
				for k in range(K):
					for l in range(L):
						pra=pra+theta[n1][k]*theta[n2][l]*tau[t][s]*pr[k][l][s][rr]
			praij[(t,n1,n2)][rr]=praij[(t,n1,n2)][rr]+pra/sampling
	
	if printp!=0: 
		fout2=open('TMMSBMparamsK'+str(K)+'L'+str(S)+'_'+str(w)+'.dat','w')
		fout2.write('%s\n' % Like)
		for i in range(p):
			for kk in range(K):
					fout2.write('%s ' % (theta[(i)][kk]))
			fout2.write('\n')
		for t in range(T):
			for ss in range(S):
					fout2.write('%s ' % (tau[t][ss]))
			fout2.write('\n')
		#for s in range(S):
		for k in range(K):
			for l in range(L):
				for s in range(S):
					for rr in range(R):
						fout2.write('%s ' % (pr[k][l][s][rr]))
					fout2.write('\n')
		fout2.close()
####################
## Test scores
###################
fout=open('TMMSBMscoresK'+str(K)+'L'+str(S)+'.dat','w')
for e in praij.keys():
	fout.write('%s %s %s %s ' % (e[0],e[1],e[2],testt[e]))
	for rr in range(R):
		fout.write('%s ' % praij[e][rr])
	fout.write('\n')
fout.close()



