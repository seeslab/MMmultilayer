import sys

import _numpypy as np
from math import *
#from numarsy import *
#import numarsy.linear_algebs as la
import copy


#import sndom
from math import sqrt,exp
import random
#from scipy.stats.stats import pearsonr

def maximum(vec):
	m1=0				
	for i in range(len(vec)):
		m1=max(max(vec[i]),m1)	
	return m1

def pearson_1d(x, y):
    """ Pearson product-moment correlation of vectors `x` and `y`

    Parameters
    ----------
    x : array shape (N,)
        One-dimensional array to correlate with `y`
    y : array shape (N,)
        One dimensional array to correlate with `x`

    Returns
    -------
    r_xy : scalar
        Pearson product-moment correlation of vectors `x` and `y`.
    """
    # Mean-center x -> mc_x
    # Mean-center y -> mc_y
    # a : Get sum of products of mc_x, mc_y
    # b : Get sum of products of mc_x on mc_x
    # c : Get sum of products of mc_y on mc_y
    # return a / (sqrt(b) * sqrt(c))
    mc_x=[]
    mc_y=[]
    mx=sum(x)/float(len(x))
    my=sum(y)/float(len(y))
    for i in range(len(x)):
	mc_x.append(x[i]-mx)
	mc_y.append(y[i]-my)
    #mc_x = x - np.mean(x)
    #mc_y = y - np.mean(y)
    a=0.
    b=0.
    c=0.
    for i in range(len(x)):
	a=a+mc_x[i]*mc_y[i]
	b=b+mc_x[i]*mc_x[i]
	c=c+mc_y[i]*mc_y[i]
    #a = mc_x.dot(mc_y)
    #b = mc_x.dot(mc_x)
    #c = mc_y.dot(mc_y)
    return a / (sqrt(b) * sqrt(c))

##############################################################################
#	this progsm geneste a bipartite network, using networkX
##############################################################################


Number=20
K=Number #drugs
L=Number #other drugs
S=Number

T=85 #cell lines
p=69 #drugs

num=sys.argv[1]	#links observed

############################################
for CV in range(5):
	fh=open('/export/home/shared/Projects/Agodoy/Dream/Paper/train'+str(CV)+'.dat','r')
	igot=fh.readlines()
	inter=[]
	#label=[-6.38738637,4.5191740,6800]
	#label=[-1.5,17,6800]
	#label=[0,6800]
	syner={}
	label=[-20.,20,6800]
	#label=[-60,-20.,20.,60,6800.]
	#label=[20.,6800.]
	#label=[-20.,6800.]
	for line in igot:
		about = line.strip().split(' ')
		#inter.append((int(about[0]),int(about[1]),int(about[2]),float(about[3])))
		b=float(about[3])
		syner[(int(about[0]),int(about[1]),int(about[2]))]=b
		g=0
		while b>label[g]:
			#print g,b, label[g],'max', max(b, label[g])
			g=g+1
		inter.append((int(about[0]),int(about[1]),int(about[2]),g))

	fh.close()
	N=len(label)
	R=len(label)

	###################################################sampling!!!!!!
	sampling=1
	#fout2=open('sampling27.6.dat','w')
	for w in range(sampling):
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
		theta[(7,7)]=[0.2]*K
		ntheta[(7,7)]=[0.]*K
		#theta=np.arsy(theta)
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
			for s in range(S):
				pr[k].append([])
				npr[k].append([])
		for k in range(K):
			for s in range(S):
				a=[]
				for r in range(R):
					a.append(random.random())
				pr[k][s] = a
				npr[k][s]=[0.]*R
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
				tau[t][s]=tau[t][s]/(D+0.0000000001)
		for k in range(K):
			for s in range(S):
				D=0.
				for r in range(R):
					D=D+pr[k][s][r]
				for r in range(R):
					pr[k][s][r]=pr[k][s][r]/D


	#########################################################################################
		Runs=1000
		cnt=0
		count=0
		ver=0
		#fout=open('likeli'+str(Runs)+'.dat','w')
		for g in range(Runs):
			#print g
			for e in inter:
				t=int(e[0])
				ra=int(e[3])
				n=(int(e[1]),int(e[2]))
				D=0.	
				for s in range(S):
					for k in range(K):
						D=D+theta[n][k]*tau[t][s]*pr[k][s][ra]
				for s in range(S):
					for k in range(K):
						a=(theta[n][k]*tau[t][s]*pr[k][s][ra])/D
						ntheta[n][k]=ntheta[n][k]+a
						ntau[t][s]=ntau[t][s]+a
						npr[k][s][ra]=npr[k][s][ra]+a
			#normalization:
			for i in range(p):
				for j in range(p-(i+1)):
					z=i+j+1
					D=0.
					for k in range(K):
						D=D+ntheta[(i,z)][k]+ntheta[(z,i)][k]
					for k in range(K):
						ntheta[(i,z)][k]=(ntheta[(i,z)][k]+ntheta[(z,i)][k])/(D+0.00000000001)
						ntheta[(z,i)][k]=ntheta[(i,z)][k]

			for t in range(T):
				D=0.
				for s in range(S):
					D=D+ntau[t][s]
				for s in range(S):
					ntau[t][s]=ntau[t][s]/(D+0.0000000000001)

			for s in range(S):
				for k in range(K):
					D=0.
					for r in range(R):
						D=D+npr[k][s][r]
					for r in range(R):
						npr[k][s][r]=npr[k][s][r]/(D+0.0000000000000000001)
		
			theta=copy.copy(ntheta)
			tau=copy.copy(ntau)
			for k in range(K):
				for s in range(S):
					for r in range(R):
						pr[k][s][r]=npr[k][s][r]	

			for i in range(p):
				for j in range(p-(i+1)):
					z=i+j+1
					ntheta[(i,z)]=[0.]*K
					ntheta[(z,i)]=[0.]*K
			for t in range(T):
				ntau[t]=[0.]*S

			for k in range(K):
				for s in range(S):
					npr[k][s]=[0.]*R

		Like=0.
		for e in inter:
			t=int(e[0])
			n=(int(e[1]),int(e[2]))
			ra=int(e[3])
			D=0.
			for s in range(S):
				for k in range(K):
					D=D+theta[n][k]*tau[t][s]*pr[k][s][ra]
			Like=Like+log(D)
		print Like

		#print g, Like
		#fout2=open('/export/home/shared/Projects/Agodoy/Dream/Paper/AntaOthers_Links/Test'+str(CV)+'/'+str(K)+'_'+str(S)+'results'+str(num)+'.dat','w')
		fout2=open('/export/home/shared/Projects/Agodoy/Dream/Paper/Links_taugroups_3/Test'+str(CV)+'/'+str(K)+'_'+str(S)+'results'+str(num)+'.dat','w')
		#fout2=open('/export/home/shared/Projects/Agodoy/Dream/Paper/Links_taufixa_3/results'+str(w)+'.dat','w')
		#fout2=open('prova.dat','w')
		fout2.write('%s\n' % Like)
		for i in range(p):
			for j in range(p-(i+1)):
				z=i+j+1
				#fout2.write('%s %s %s %s %s\n' % (theta[(i,z)][0], theta[(i,z)][1],theta[(i,z)][2],theta[(i,z)][3],theta[(i,z)][4]))
				for kk in range(K):
					fout2.write('%s ' % (theta[(i,z)][kk]))
				fout2.write('\n')
			#,theta[i][5],theta[i][6],theta[i][7],theta[i][8]))
			#,theta[i][7],theta[i][8],theta[i][9]))
		for t in range(T):
				#fout2.write('%s %s %s %s %s\n' % (tau[t][0], tau[t][1],tau[t][2],tau[t][3],tau[t][4]))
				for ss in range(S):
					fout2.write('%s ' % (tau[t][ss]))
				fout2.write('\n')
		for k in range(K):
			for s in range(S):
				#for t in range(T):
				fout2.write('%s %s %s\n' % (pr[k][s][0], pr[k][s][1],pr[k][s][2]))
					#,pr[k][s][2],pr[k][s][3],pr[k][s][4]))
					#,pr[k][l][s][3],pr[k][l][s][4],
					#pr[k][l][s][5]))
					#,pr[k][l][s][7],pr[k][l][s][8],pr[k][l][s][9]))
		fout2.close()

