import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import csv
import shutil 
import numpy as numpy
import time
import numpy.linalg as la
from mpl_toolkits.mplot3d import *
import scipy as sp
from scipy.interpolate import interp1d
from matplotlib.pyplot import * 
from pylab import *
import urllib
import cmath
import math as math 
from random import *
from scipy import optimize 
clf #clears plot
import itertools as itertools
import shlex
import pickle 
import os 
from scipy.integrate import quad
import scipy.special as sps
from pylab import figure, show, rand
from matplotlib.patches import Ellipse
import os.path
from heapq import nsmallest

import Zlm2 as Zlm2
Zlm2=reload(Zlm2)

import Functs as Functs
Functs=reload(Functs)


####Inputs  
twopi=2.0*pi
hc=197.3269
I=(0+1j)
eps=0.00000000000000001

mpi=0.03928
mK=0.08344 
aniso=3.4534 
Loas=32
L=Loas*aniso

Ecm=arange(0.07,0.24,.0001)
 
"""
********************
jack knife code
********************
"""
def jack(X):
	nn=len(X)
	JACK=zeros(nn)
	for kk in range(nn):
		JACK[kk]=(sum(X)-X[kk])/float(nn-1.)
	return JACK
 

###standard deviation and convariance matrix
	
def jack_sig(X):
	nn=float(len(X))
	ave=mean(X)
	return pow((((nn-1.)/nn)*sum((X-ave)**2.)),.5)

def inflate(X):
	NN=len(X)
	amp=float(sqrt(NN-1.0))
	mx=mean(X)
	return mx+((X-mx)*amp)


def deflate(X):
	NN=len(X)
	amp=1.0/float(sqrt(NN-1.0))
	mx=mean(X)
	return mx+((X-mx)*amp)


def loadjack(files):
	with open(files, 'r') as fin:
  		data = fin.read().splitlines(True)


	NG=int(data[0].split()[0])
	out=ones(NG)
	for ng in range(NG):
		out[ng]=data[ng+1].split()[1]
	return inflate(jack(out))



#########################################
####################Inputs###############
#########################################


def Couts(part_type,irrep,dname,Loas):
	"part_type= pipi or KKbar"
	filename="couts/"+part_type+"_d"+dname+"_"+irrep+"_L."+str(int(Loas))+".txt"
	cout=loadtxt(filename)
	
	return cout
	 
	 

"""
we need smart interpolating functions for the 
phase shift, QCs and LL factors...
"""
def S(Z):
	return (Z+I)/(Z-I)
	
def Z(S):
	return I*(S+1.0)/(S-1.0)
	
def clever_Z_interp(Eout,Vinp):
	S_VinpR,S_VinpI=real(S(Vinp)),imag(S(Vinp))
	fS_VinpR_tmp=interp1d(Ecm,S_VinpR)
	fS_VinpI_tmp=interp1d(Ecm,S_VinpI)
 

	return real(Z(fS_VinpR_tmp(Eout)+I*fS_VinpI_tmp(Eout)))					
 

def fLuscher_pcotdelta(Eout,part_type,irrep,dname,Loas): 
	"""change total energy to CM energy"""
	EEcm=pow(Eout,2.0)
	
	if part_type=="pipi":
		qq=(EEcm/4.0)-pow(mpi,2.0)
		qout=sqrt(qq)	
		
	if part_type=="KKbar":
		qq=(EEcm/4.0)-pow(mK,2.0)
		qout=sqrt(qq)

	"""cout=(c00+alpha20*c20+alpha22*c22)"""
	cout=Couts(part_type,irrep,dname,Loas)

	"""interpolating the function"""
	fcout=clever_Z_interp(Eout,cout)
	
 	"cot*delta"
 	cotd=(4.0*pi*fcout)/qout
	
 	"delta"
 	delta=arctan(1.0/cotd)*180.0/pi
 	if mean(delta)<0 and part_type=="pipi":
 		delta=delta+180.0

	return delta
 
 
filename='../spec_final_onlypipi.list'

with open(filename, 'r') as fin:
  	states = fin.read().splitlines(True)


for i0 in range(len(states)):
	
	dname,irrep,L_label,jack_energy_file=states[i0].split()
	
	Etot=loadjack(jack_energy_file)
	ps=fLuscher_pcotdelta(Etot,"pipi",irrep,dname,Loas)
	
	d=array([float(dname[0]),float(dname[1]),float(dname[2])])
	P=twopi*d/L
	Eout=sqrt(pow(Etot,2.0)-dot(P,P))
	
	x,dx=mean(Eout/mpi),std(Eout/mpi)
	y,dy=mean(ps),std(ps)
 	errorbar(x,y,yerr=dy,xerr=dx,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
	

 
filename='../spec_final_onlykkbar.list'

with open(filename, 'r') as fin:
  	states = fin.read().splitlines(True)


for i0 in range(len(states)):
	
	dname,irrep,L_label,jack_energy_file=states[i0].split()
	
	Etot=loadjack(jack_energy_file)
	ps=fLuscher_pcotdelta(Etot,"KKbar",irrep,dname,Loas)
	
	d=array([float(dname[0]),float(dname[1]),float(dname[2])])
	P=twopi*d/L
	Eout=sqrt(pow(Etot,2.0)-dot(P,P))
	
	x,dx=mean(Eout/mpi),std(Eout/mpi)
	y,dy=mean(ps),std(ps)
 	errorbar(x,y,yerr=dy,xerr=dx,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
	

#############################################################################
####################### Interpolating functions #################################
#############################################################################
 

for i in range(len(Einps)):
	
	print ""
	print "**************************"
	print ""
	print "Einp"
	print "Einps[i]=",mean(Einps[i])
	print "P",P
	print "L",L 
	PP=pow(P,2.0)
	Einps[i]=jack(Einps[i])
	EE=pow(Einps[i],2.0)

	Einps[i]=inflate(sqrt(EE-PP))


"""
************************************************************************************************************************************************************************************************************
************************************************************************************************************************************************************************************************************
************************************************************************************************************************************************************************************************************
"""



"""
************************************************************************************************************************************************************************************************************
************************************************************************************************************************************************************************************************************
************************************************************************************************************************************************************************************************************
"""



 
#####################################################
####################Breit Wigner Phase###############
#####################################################
"""
 # |                                parameter name |         value +/- error         |       ||      0     1
-----------------------------------------------------------------------------------------------  ----- -----
 0 |                                         J1_mR |       0.15047 +/- 0.0002684     |       ||   1.00 -0.06
 1 |                                         J1_gR |        5.7848 +/- 0.10622       |       ||         1.00
-----------------------------------------------------------------------------------------------  ----- -----

"""
mR=0.15047
grho=5.7848
dmR=0.0002684
dgrho=0.10622
corr=-0.06

mpi=mpis[0]


psq=(pow(mR,2.0)-pow(2.0*mpi,2.0))/4.0
pcm=sqrt(psq)
qrho=pcm	
Gammarho=pow(grho,2.0)*pow(pcm,3.0)/(6.0*pi*pow(mR,2.0))

print "Gammarho=",Gammarho
print "Gammarho/mr=",Gammarho/mR

Covf=array([[pow(dmR,2.0),dmR*dgrho*corr],[dmR*dgrho*corr,pow(dgrho,2.0)]])
Mean=[mR,grho]
"""
corr=1.00 -0.22
	-0.22  1.00
"""

NG=603
mRens,grhoens=np.random.multivariate_normal(Mean,Covf,NG).T
runfirsttime="n"
if runfirsttime=="y":
	savetxt('jack_masses/mRho_ens.txt',mRens)
	savetxt('jack_masses/gRho_ens.txt',grhoens)
	end
else:	
	mRens=loadtxt('jack_masses/mRho_ens.txt')
	grhoens=loadtxt('jack_masses/gRho_ens.txt')

"""
here is the Breit wigner phase shift
"""

def pitorho(Ecm): 
	return Functs.gen_pitorho(Ecm,mpi,type,Params,ens) 

Params=[grhoens,mRens]
type="BW"
modellabel=r'$\rm{Breit}$'+'-'+r'$\rm{Wigner}$'
ens='n'
	
Ecm=arange(.132,0.29,.0001)

EEcm=pow(Ecm,2.0)
qq=(EEcm/4.0)-pow(mpi,2.0)
qout=sqrt(qq)

"""
here we make phaseshift (ps) and qcot(ps) 
as a function of the energy and as a function of the LECs of the phase shift
the latter will be labeled as "ens" meaning ensemble
"""

psens =zeros((NG,len(Ecm)),dtype=float)		
ps =zeros(len(Ecm),dtype=float)		
qcotpsens =zeros((NG,len(Ecm)),dtype=float)		
psmax =zeros(len(Ecm))
psmin =zeros(len(Ecm))
qcotpsmax =zeros(len(Ecm))
qcotpsmin =zeros(len(Ecm))
qcotps =zeros(len(Ecm))

subplot(211)
for n0 in range(len(Ecm)):
	psens[:,n0],qcotpsens[:,n0],ps[n0],psmax[n0],psmin[n0],qcotps[n0],qcotpsmax[n0],qcotpsmin[n0]=Functs.phase_shift(Ecm[n0],mpi,type,Params,ens) 

units="lattice"
if units=="phys":
	fill_between(Ecm/at,psmax*180.0/pi, psmin*180.0/pi,facecolor='orange',alpha='.5', interpolate=True)
if units=="lattice":
	fill_between(Ecm/mpi,psmax*180.0/pi, psmin*180.0/pi,facecolor='orange',alpha='.5', interpolate=True)



print ""
print "**************************"
print 'irrep =',irrep


narrow="n"

for j in range(len(Einp)):

	print ""
	print "**************************"
	print "Einp=",Einp[j]
	print 'fQC(Einp)=',fQC(Einp[j])
	print "phase[%] at Einp=",fps(Einp[j])*180.0/pi,fLuscher_delta(Einp[j])*180.0/pi
	print 
	Eout=optimize.fsolve(fQC,Einp[j])
	while Eout >0.23 or Eout <0.138:
		Ej=numpy.random.normal(Einp[j], dEinp[j], size=1)
		while Ej >0.23 or Ej <0.138:
			Ej=numpy.random.normal(Einp[j], dEinp[j], size=1)
		Eout=optimize.fsolve(fQC,Ej)
 

	print 'fQC(Eout)=',fQC(Eout)
	print "Ecm=",Eout
	print "phase[%] at Ecm=",fps(Eout)*180.0/pi ,fLuscher_delta(Eout)*180.0/pi
	print "assymetric error (sign ambiguity?)"
	#print "LL factor [Lfree]=",abs(fRinvf(Eout)),'+',(abs(fRinvfmax(Eout)-fRinvf(Eout))),'+',abs(fRinvf(Eout)-fRinvfmin(Eout))
	#print "LL factor=",abs(fRinvf(Eout))*fRinvfree(Eout),'+',(abs(fRinvfmax(Eout)-fRinvf(Eout)))*fRinvfree(Eout),'+',abs(fRinvf(Eout)-fRinvfmin(Eout))*fRinvfree(Eout)



	"""
	here we plot the phase shift
	"""
	subplot(211)
	y=fps(Eout)*180.0/pi 
	"""model prediction of spectrum"""
	if units=="phys":
		errorbar(Eout/at,y,yerr=0,xerr=0,markersize=8,fmt='o',color='k', mfc='white',mec='k',elinewidth=2, capsize=4, mew=1.4)	
	if units=="lattice":
		errorbar(Eout/mpi,y,yerr=0,xerr=0/mpi,markersize=8,fmt='o',color='k',mfc='white',mec='k', elinewidth=2, capsize=4, mew=1.4)	

	
	Ej=numpy.random.normal(Einp[j], dEinp[j], size=5000)
	print shape(Einps[j])

	ym=fLuscher_delta(Einps[j])*180.0/pi
	ym,dy=mean(ym),std(ym)
	
	if ym<0:
		ym=ym+180.0
		
	print ""
	print "**********************************"
	print "**********************************"
	print "shape shift at lattice energy using clever interpolation"
	print ym,dy
	print "**********************************"
	print "**********************************"
	print ""
	
	"""lattice determination of spectrum"""
	if units=="phys":
		errorbar(Einps[j]/at,ym,yerr=dy,xerr=dEinp[j],markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
	if units=="lattice":
		errorbar(mean(Einps[j]/mpi),ym,yerr=dy,xerr=std(Einps[j]/mpi),markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
	hack=0
	if hack==1:
		ym=Functs.Luscher_ps_no_interp_vec(Einps[j],mpi,irrep,d,L)*180.0/pi
		ym,dy=mean(ym),std(ym)
		if ym<0:
			ym=ym+180.0
		print ""
		print "**********************************"
		print "**********************************"
		print "No interpolation!!!!!!"
		print "phase shift=",ym,dy
		print "**********************************"
		print "**********************************"
		print ""
		"""lattice determination of spectrum"""
		if units=="lattice":
			errorbar(mean(Einps[j]/mpi),ym,yerr=dy,xerr=std(Einps[j]/mpi),markersize=8,fmt='d',color='g',mfc='white',mec='g', elinewidth=2, capsize=4, mew=1.4)	
	
	plt.ylabel(r'$\delta_1/{\circ}$',size=25)
	
 
	"""
	here we plot the ratio of the LL factors to the free one
	"""
	subplot(212)
	 
	ym_no_inter=fRinv(Einps[j])
	
	
	ylim([0,5.0])


	#########################
	####Saving LL factors####
	#########################
	
	print ""
	print "**************************"
	print "LL factor in absolute units"

	LLs[:,3*j]=deflate(fRinv(Einps[j]))
	LLs[:,3*j+1]=mean(LLs[:,3*j])*ones(len(Einps[j]))
	LLs[:,3*j+2]=(fRinv(Eout))*ones(len(Einps[j]))

	print mean(LLs[:,3*j]), jack_sig(LLs[:,3*j])

	

	########################################
	####Saving pi_to_rho proportionality####
	########################################

	ens="n"
	deltaps=fps(Einps[j])
	out=pitorho(Einps[j]) 
	
	RhotoPi[:,3*j]=abs(deflate(out)*LLs[:,3*j])
	RhotoPi[:,3*j+1]=abs(mean(RhotoPi[:,3*j])*ones(len(Einps[j])))
	
	out0=pitorho(Eout) 
	RhotoPi[:,3*j+2]=abs(out0*ones(len(Einps[j]))*LLs[:,3*j+2])

	out=pitorho(Einps[j]) 

	RhotoPi2[:,j]=abs(deflate(ym_no_inter*out))

	"""
	RhotoPinarrow[:,3*j]=abs(LLs[:,3*j]/fRinvrho(Einps[j]))
	RhotoPinarrow[:,3*j+1]=abs(mean(RhotoPinarrow[:,3*j])*ones(len(Einps[j])))
	RhotoPinarrow[:,3*j+2]=abs((1.0/fRinvrho(Eout))*ones(len(Einps[j]))*LLs[:,3*j+2])
	"""
	

	"""
	Here I determine the LL and RhotoPi as ensembles 
	to propagate the systematic error due to the fit of the phase shift
	"""	

	for ng in range(NG):
		Rinvens=sqrt(fcos2ps(Ecm)*abs(16.0*pi*fEftot(Ecm)*Einitial*derRens[ng,:]))/fNif(Ecm)
		fRinvsys=interp1d(Ecm,Rinvens)
		LLssys[ng,2*j]=fRinvsys(mean(Einps[j]))
	ens="y"

	LLssys[:,2*j+1]=abs(pitorho(mean(Einps[j]))*LLssys[:,2*j])
	
	print "**************************"
	print "proportionality between pi pi and rho"	
	print mean(RhotoPi[:,3*j]), jack_sig(RhotoPi[:,3*j])
	print mean(RhotoPi2[:,j]), jack_sig(RhotoPi2[:,j])
	

	"""
	here we plot the ratio of the LL factors to the rho value
	"""
		
	ens="n"
	
	if narrow=="y":
		ym,dy=mean(RhotoPinarrow[:,3*j]), jack_sig(RhotoPinarrow[:,3*j])
		if units=="phys":
			errorbar(mean(Einps[j])/at,ym,yerr=dy,xerr=std(Einps[j]),markersize=8,fmt='d',color='r', elinewidth=2, capsize=4, mew=1.4)	
		if units=="lattice":
			errorbar(mean(Einps[j]),ym,yerr=dy,xerr=std(Einps[j]),markersize=8,fmt='d',color='r', elinewidth=2, capsize=4, mew=1.4)	
	
		ym=abs(fRinv(Eout)/fRinvrho(Eout))
 
		if units=="phys":
			errorbar(Eout/at,ym,yerr=0,xerr=0,markersize=8,fmt='o',color='k', elinewidth=2, capsize=4, mew=1.4)	
		if units=="lattice":
			errorbar(Eout,ym,yerr=0,xerr=0,markersize=8,fmt='o',color='k', elinewidth=2, capsize=4, mew=1.4)	
	else:
		ym=abs(fRinv(Eout)*pitorho(Eout))

		if units=="phys":
			errorbar(Eout/at,ym,yerr=dy,xerr=dEout,markersize=8,fmt='o',color='k',mfc='white',mec='k', elinewidth=2, capsize=4, mew=1.4)	
		if units=="lattice":
			if j==0:
				modeldata=errorbar(Eout/mpi,ym,yerr=0,xerr=0,markersize=8,fmt='o',color='k', mfc='white',mec='k',elinewidth=2, capsize=4, mew=1.4)	
			else:
				errorbar(Eout/mpi,ym,yerr=0,xerr=0,markersize=8,fmt='o',color='k', mfc='white',mec='k',elinewidth=2, capsize=4, mew=1.4)	
		
		
		"""here we plot the LL factor without interpolating"""
		ym,dystat=mean(abs(ym_no_inter*pitorho(Einps[j]))),std(ym_no_inter*pitorho(Einps[j]))

		dysys=std(LLssys[:,2*j+1])
		dyt=sqrt(pow(dystat,2.0)+pow(dysys,2.0))
		if units=="phys":
			errorbar(mean(Einps[j])/at,ym,yerr=dy,xerr=std(Einps[j])/at,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
		if units=="lattice":
			if j==0:
				latticedata=errorbar(mean(Einps[j])/mpi,ym,yerr=dystat,xerr=std(Einps[j])/mpi,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
			else:
				errorbar(mean(Einps[j])/mpi,ym,yerr=dyt,xerr=std(Einps[j])/mpi,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=1, capsize=4, mew=1.4)	
				errorbar(mean(Einps[j])/mpi,ym,yerr=dystat,xerr=std(Einps[j])/mpi,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	

 	if hack==1:
	 	ens='n'
	 	ym_no_inter_2=Functs.Rinv_no_inter(Einps[j],mpi,irrep,d,ps,qcotps,Ecm,ens,L)
		
		print ""
		print "*****************************************"
		print "*****************************************"
		print "with interpolation:" ,ym,dystat
		print "without interpolation:" ,mean(abs(ym_no_inter_2*pitorho(Einps[j]))),std(ym_no_inter_2*pitorho(Einps[j]))
		print "*****************************************"
		print "*****************************************"
		ym,dystat=mean(abs(ym_no_inter_2*pitorho(Einps[j]))),std(ym_no_inter_2*pitorho(Einps[j]))
		dysys=std(LLssys[:,2*j+1])
		dyt=sqrt(pow(dystat,2.0)+pow(dysys,2.0))
		if j==0:
			latticedata=errorbar(mean(Einps[j])/mpi,ym,yerr=dystat,xerr=std(Einps[j])/mpi,markersize=8,fmt='d',color='r',mfc='white',mec='r', elinewidth=2, capsize=4, mew=1.4)	
		else:
			errorbar(mean(Einps[j])/mpi,ym,yerr=dyt,xerr=std(Einps[j])/mpi,markersize=8,fmt='d',color='g',mfc='white',mec='g', elinewidth=1, capsize=4, mew=1.4)	
			errorbar(mean(Einps[j])/mpi,ym,yerr=dystat,xerr=std(Einps[j])/mpi,markersize=8,fmt='d',color='g',mfc='white',mec='g', elinewidth=2, capsize=4, mew=1.4)	
		

		print "******************************************************************************************************************************"
		print "******************************************************************************************************************************"


#########################
####Saving LL factors####
#########################
print "Saving LL factors in ",final_dir_LL
dirjk=dirjk+"/LLfactors/"
savetxt(dirjk+final_dir_LL,LLs)
savetxt(dirjk+final_dir_rhotopi,RhotoPi)
savetxt(dirjk+final_dir_rhotopi2,RhotoPi2)
savetxt(dirjk+final_dir_LL_sys,LLssys)

##########################################################################################
##########################################################################################
"""
Here we plot the phase shift
"""
subplot(211) 
plt.yticks(arange(0,180,40),size=12)
if units=="phys":
	fill_between(Ecm[2:len(Ecm)-2]/at,psmax[2:len(Ecm)-2]*180.0/pi, psmin[2:len(Ecm)-2]*180.0/pi,facecolor='orange',alpha='.5', interpolate=True)
	plt.text(0.18/at, 80, irrepname,size=20, ha="center", va="center",bbox = dict(boxstyle="round",ec=(0., 0, 0), fc=(1., 1, 1)))
	plt.xticks(arange(780,1120,40),size=12)
	xlim([2.0*mpi/at,.2/at]) 
if units=="lattice":
	fill_between(Ecm[2:len(Ecm)-2]/mpi,psmax[2:len(Ecm)-2]*180.0/pi, psmin[2:len(Ecm)-2]*180.0/pi,facecolor='orange',alpha='.5', interpolate=True)
	

	plt.xticks(arange(2.0,3,.2),size=12)
	xlim([.14/mpi,.198/mpi]) 

ylim([0,180]) 

##########################################################################################
##########################################################################################
"""
Here we plot the LL factor bands
"""
subplot(212)
if narrow=="y":
	if units=="phys":
		fill_between(Ecm[2:len(Ecm)-2]/at,Rinvmax[2:len(Ecm)-2]/abs(fRinvrho(Ecm[2:len(Ecm)-2])), Rinvmin[2:len(Ecm)-2]/abs(fRinvrho(Ecm[2:len(Ecm)-2])),facecolor='g',alpha='.5', interpolate=True)
	if units=="lattice":
		fill_between(Ecm[2:len(Ecm)-2],Rinvmax[2:len(Ecm)-2]/abs(fRinvrho(Ecm[2:len(Ecm)-2])), Rinvmin[2:len(Ecm)-2]/abs(fRinvrho(Ecm[2:len(Ecm)-2])),facecolor='g',alpha='.5', interpolate=True)
else:
	if units=="phys":
		fill_between(Ecm[2:len(Ecm)-2]/at,fRinvmax(Ecm[2:len(Ecm)-2])*abs(pitorho(Ecm[2:len(Ecm)-2])), fRinvmin(Ecm[2:len(Ecm)-2])*abs(pitorho(Ecm[2:len(Ecm)-2])),facecolor='b',alpha='.4', interpolate=True)
	if units=="lattice":
		fill_between(Ecm[2:len(Ecm)-2]/mpi,fRinvmax(Ecm[2:len(Ecm)-2])*abs(pitorho(Ecm[2:len(Ecm)-2])), fRinvmin(Ecm[2:len(Ecm)-2])*abs(pitorho(Ecm[2:len(Ecm)-2])),facecolor='b',alpha='.4', interpolate=True)
if units=="phys":
	fill_between(Ecm[2:len(Ecm)-2]/at,fRinvfmax(Ecm[2:len(Ecm)-2]), fRinvfmin(Ecm[2:len(Ecm)-2]),facecolor='orange',alpha='.5', interpolate=True)
plt.xticks(arange(780,1120,40),size=12)
xlim([2.0*mpi/at,.2/at]) 
if units=="lattice":
	#fill_between(Ecm[2:len(Ecm)-2]/mpi,Rinvfmax[2:len(Ecm)-2], Rinvfmin[2:len(Ecm)-2],facecolor='orange',alpha='.5', interpolate=True)
	
	plt.xticks(arange(2.0,3,.2),size=12)
	xlim([.14/mpi,.198/mpi]) 
	plt.xlabel(r'$\rm{E}^*_{\pi\pi}/m_\pi$',size=25)


subplot(211)
extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
latticelabel=r'$\rm{lattice}$'+' ' +r'$\rm{QCD}$'
legend([extra, latticedata, modeldata], (irrepname,latticelabel, modellabel),numpoints=1, fancybox=True, shadow=True, loc=4)


#plt.text(2.5, 80, irrepname,size=20, ha="center", va="center",bbox = dict(boxstyle="round",ec=(0., 0, 0), fc=(1., 1, 1)))

subplot(212)
plt.ylabel(r'$\sqrt{\mathcal{K}_{\pi\pi}/\mathcal{R}}$',size=25)

plt.yticks(arange(0,6,1),size=12)
ylim([0,3.5]) 


	
gcf().subplots_adjust(left=0.20)
gcf().subplots_adjust(bottom=0.20)
if narrow=="y":
	plt.savefig('figures/'+irrep+"_narrow",dpi=150)
else:
	plt.savefig('figures/'+irrep,dpi=150)