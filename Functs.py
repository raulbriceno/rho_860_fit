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
import h5py as h5py
import os.path
from heapq import nsmallest

import Zlm2 as Zlm2
Zlm2=reload(Zlm2)

####Inputs  
twopi=2.0*pi
hc=197.3269
I=(0+1j)
eps=0.00000000000000001

at=0.2951/1672.0

xi=1.0/2.0
 

"""
***********************
	code that uses all the functions defined below
***********************
""" 

def phase_shift(Ecm,mpi,type,Params,ens):
	
	if type=="BW":
		grhoens,mRens=Params
		return ps_BW(Ecm,mpi,grhoens,mRens,ens)

	if type=="BW_HQ":
		grhoens,mRens,Rangeens=Params
		return ps_BWHQ(Ecm,mpi,grhoens,mRens,Rangeens,ens)

	if type=="BW_Gauss":
		grhoens,mRens,Rangeens=Params
		return ps_BWGauss(Ecm,mpi,grhoens,mRens,Rangeens,ens)

	if type=="Kmat":
		grhoens,mRens=Params
		return ps_Kmat(Ecm,mpi,grhoens,mRens,ens)
	
	
	
def gen_pitorho(Ecm,mpi,type,Params,ens):
	
	if type=="BW":
		grhoens,mRens=Params
		return pitorho_BW(Ecm,mpi,grhoens,mRens,ens)

	if type=="BW_HQ":
		grhoens,mRens,Rangeens=Params
		return pitorho_BWHQ(Ecm,mpi,grhoens,mRens,Rangeens,ens)

	if type=="BW_Gauss":
		grhoens,mRens,Rangeens=Params
		return pitorho_BWGauss(Ecm,mpi,grhoens,mRens,Rangeens,ens)

	if type=="Kmat":
		grhoens,mRens=Params
		return pitorho_Kmat(Ecm,mpi,grhoens,mRens,ens)
		
		
"""
##############################################
##############################################
############################################## 
"""


"""
***********************
	BW Wigner Codes
***********************
""" 

def ps_BW(Ecm,mpi,grhoens,mRens,ens):
	gsq=grhoens*grhoens
	mRsq=mRens*mRens
	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	pcm=sqrt(psq)
	
	#Gamma_1
	G1=gsq*pow(pcm,3.0)/(6.0*pi*Esq)
	#numerator	
	num=Ecm*G1
	#denumerator
	dnum=mRsq-Esq
	#tan(phi)
	tanphi=num/dnum

	#Gamma_1
	G1=gsq*psq/(6.0*pi*Esq)
	#numerator	
	num=Ecm*G1
	#denumerator
	dnum=mRsq-Esq
	#qcot(phi)
	qcotphi=dnum/num
	eps=.000001
	qcotphi=qcotphi+eps
	phi=arctan(tanphi)
	if shape(phi) != ():
		for j in range(len(phi)):
			
			if phi[j]<0:
				phi[j]=phi[j]+pi
	else:
		if phi<0:
			phi=phi+pi
	return phi,qcotphi, mean(phi),mean(phi)+std(phi),mean(phi)-std(phi),mean(qcotphi),mean(qcotphi)+std(qcotphi),mean(qcotphi)-std(qcotphi)


def pitorho_BW(Ecm,mpi,grhoens,mRens,ens):
	if ens=="n":
		gsq=mean(grhoens*grhoens)
		mRsq=mean(mRens*mRens)
	if ens=="y":
		gsq=grhoens*grhoens
		mRsq=mRens*mRens
	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	pcm=sqrt(psq)
	
	#Gamma_1
	G1=gsq*pow(pcm,3.0)/(6.0*pi*Esq)

	rho=sqrt(8.0*pi*Ecm/(pcm*xi))

	tand=Ecm*G1/(mRsq-Esq)
	delta=arctan(tand)
	num=sqrt(Ecm*G1)*sin(delta)

	dnum=Ecm*G1
	amp=rho*num/dnum
	return 1.0/amp
	



"""
********************************
	BW with Gaussian Barrier
********************************
""" 

def ps_BWGauss(Ecm,mpi,grhoens,mRens,Rangeens,ens):
	gsq=grhoens*grhoens
	mRsq=mRens*mRens
	R2=Rangeens*Rangeens
	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	psR=(mRsq-pow(2.0*mpi,2.0))/4.0
	pcm=sqrt(psq)
	
	#Gamma_1
	"""
	the only difference with the regular BW is the scaling factor
	"""
	amp_Gaus=exp(-psq*R2)/exp(-psR*R2)
	Gamma=gsq*pow(pcm,3.0)/(6.0*pi*Esq)
	
	#numerator	
	
	num=Ecm*Gamma*amp_Gaus
	#denumerator
	dnum=mRsq-Esq
	#tan(phi)
	tanphi=num/dnum

	#Gamma_1

	G1_over_pcm=gsq*psq/(6.0*pi*Esq)
	#numerator	
	num=Ecm*G1_over_pcm*amp_Gaus
	#denumerator
	dnum=mRsq-Esq
	#qcot(phi)
	qcotphi=dnum/num
	eps=.000001
	qcotphi=qcotphi+eps
	phi=arctan(tanphi)
	if shape(phi) != ():
		for j in range(len(phi)):
			
			if phi[j]<0:
				phi[j]=phi[j]+pi
	else:
		if phi<0:
			phi=phi+pi
	return phi,qcotphi, mean(phi),mean(phi)+std(phi),mean(phi)-std(phi),mean(qcotphi),mean(qcotphi)+std(qcotphi),mean(qcotphi)-std(qcotphi)



def pitorho_BWGauss(Ecm,mpi,grhoens,mRens,Rangeens,ens):

	if ens=="n":
		gsq=mean(grhoens*grhoens)
		mRsq=mean(mRens*mRens)
		R2=mean(Rangeens*Rangeens)
	if ens=="y":
		gsq=grhoens*grhoens
		mRsq=mRens*mRens
		R2=Rangeens*Rangeens

	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	psR=(mRsq-pow(2.0*mpi,2.0))/4.0

	pcm=sqrt(psq)
	
	#Gamma_1
	"""
	the only difference with the regular BW is the scaling factor
	"""
	amp_Gaus=exp(-psq*R2)/exp(-psR*R2)
	
	#numerator	
	
	G1=gsq*pow(pcm,3.0)/(6.0*pi*Esq)
	G1=G1*amp_Gaus
	
	rho=sqrt(8.0*pi*Ecm/(pcm*xi))

	tand=Ecm*G1/(mRsq-Esq)
	delta=arctan(tand)
	num=sqrt(Ecm*G1)*sin(delta)

	dnum=Ecm*G1
	amp=rho*num/dnum
	return 1.0/amp
	


"""
********************************
	BW with HQ Barrier
********************************
""" 


def ps_BWHQ(Ecm,mpi,grhoens,mRens,Rangeens,ens):
	gsq=grhoens*grhoens
	mRsq=mRens*mRens
	R2=Rangeens*Rangeens
	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	psR=(mRsq-pow(2.0*mpi,2.0))/4.0
	pcm=sqrt(psq)
	
	#Gamma_1
	"""
	the only difference with the regular BW is the scaling factor
	"""
	amp_HQ=(1.0+psR*R2)/(1.0+psq*R2)
	Gamma=gsq*pow(pcm,3.0)/(6.0*pi*Esq)
	
	#numerator	
	
	num=Ecm*Gamma*amp_HQ
	#denumerator
	dnum=mRsq-Esq
	#tan(phi)
	tanphi=num/dnum

	#Gamma_1

	G1_over_pcm=gsq*psq/(6.0*pi*Esq)
	#numerator	
	num=Ecm*G1_over_pcm*amp_HQ
	#denumerator
	dnum=mRsq-Esq
	#qcot(phi)
	qcotphi=dnum/num
	eps=.000001
	qcotphi=qcotphi+eps
	phi=arctan(tanphi)
	if shape(phi) != ():
		for j in range(len(phi)):
			
			if phi[j]<0:
				phi[j]=phi[j]+pi
	else:
		if phi<0:
			phi=phi+pi
	return phi,qcotphi, mean(phi),mean(phi)+std(phi),mean(phi)-std(phi),mean(qcotphi),mean(qcotphi)+std(qcotphi),mean(qcotphi)-std(qcotphi)


def pitorho_BWHQ(Ecm,mpi,grhoens,mRens,Rangeens,ens):
	if ens=="n":
		gsq=mean(grhoens*grhoens)
		mRsq=mean(mRens*mRens)
		R2=mean(Rangeens*Rangeens)

	if ens=="y":
		gsq=grhoens*grhoens
		mRsq=mRens*mRens
		R2=Rangeens*Rangeens

	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	psR=(mRsq-pow(2.0*mpi,2.0))/4.0

	pcm=sqrt(psq)
	
	#Gamma_1
	"""
	the only difference with the regular BW is the scaling factor
	"""
	amp_HQ=(1.0+psR*R2)/(1.0+psq*R2)
	
	#numerator	
	
	G1=gsq*pow(pcm,3.0)/(6.0*pi*Esq)
	G1=G1*amp_HQ
	
 	rho=sqrt(8.0*pi*Ecm/(pcm*xi))

	tand=Ecm*G1/(mRsq-Esq)
	delta=arctan(tand)
	num=sqrt(Ecm*G1)*sin(delta)

	dnum=Ecm*G1
	amp=rho*num/dnum
	return 1.0/amp



"""
********************************
	Single Channel K-matrix
********************************
""" 
	
def ps_Kmat(Ecm,mpi,grhoens,mRens,ens):

	gsq=grhoens*grhoens
	mRsq=mRens*mRens
	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	pcm=sqrt(psq)
	
	#Gamma_1
	G1=8.0*gsq*pow(pcm,3.0)/Esq
	#rho_phase space
	rhops=2.0*pcm/sqrt(Esq)
	xips=pow(rhops,2.0)
	logterm=gsq*pow(2.0*pcm,2.0)*(rhops/pi)*(log((rhops+xips)/(rhops-xips)))
	#numerator	
	num=Ecm*G1
	#denumerator
	dnum=mRsq-Esq+logterm
	#tan(phi)
	tanphi=num/dnum

	#Gamma_1
	G1=8.0*gsq*psq/Esq

	#numerator	
	num=Ecm*G1
	#denumerator
	dnum=mRsq-Esq+logterm
	#qcot(phi)
	qcotphi=dnum/num
	eps=.000001
	qcotphi=qcotphi+eps
	phi=arctan(tanphi)
	if shape(phi) != ():
		for j in range(len(phi)):
			
			if phi[j]<0:
				phi[j]=phi[j]+pi
	else:
		if phi<0:
			phi=phi+pi
	return phi,qcotphi, mean(phi),mean(phi)+std(phi),mean(phi)-std(phi),mean(qcotphi),mean(qcotphi)+std(qcotphi),mean(qcotphi)-std(qcotphi)


def pitorho_Kmat(Ecm,mpi,grhoens,mRens,ens):
	if ens=="n":
		gsq=mean(grhoens*grhoens)
		mRsq=mean(mRens*mRens)
	if ens=="y":
		gsq=grhoens*grhoens
		mRsq=mRens*mRens
	Esq=pow(Ecm,2.0) 
	psq=(Esq-pow(2.0*mpi,2.0))/4.0
	pcm=sqrt(psq)
	
	#Gamma_1
	
	G1=8.0*gsq*pow(pcm,3.0)/Esq
	#rho_phase space
	rhops=2.0*pcm/sqrt(Esq)
	xips=pow(rhops,2.0)
	logterm=gsq*pow(2.0*pcm,2.0)*(rhops/pi)*(log((rhops+xips)/(rhops-xips)))
	
	"""
	num=sqrt(Ecm*G1)
	dnum=mRsq-Esq-I*Ecm*G1
	amp=rho*num/dnum
	"""
	
	
	"""
	alternatively, I unambiguiously separate the real piece and throw away 
	the overall phase that exactly cancels with the scattering phase in the LL-factor
	"""
	rho=sqrt(8.0*pi*Ecm/(pcm*xi))

	tand=Ecm*G1/(mRsq-Esq+logterm)
	delta=arctan(tand)
	num=sqrt(Ecm*G1)*sin(delta)

	dnum=Ecm*G1
	amp=rho*num/dnum
	return 1.0/amp







"""
#################################################
#################################################
#################################################
"""


def no_interp_QC(Ecm0,mpi,irrep,d,qcotps,Eps,ens,L):

	"""
	####################
	first, we define an interpolating function for the q*cot*(detal), which have been evaluated in the energy range Eps
	"""
	
	fqcotps=interp1d(Eps,qcotps)
	
	"""
	####################
	"""
	
	gamma0=Zlm2.gamma(Ecm0,d,L)
	alphapipi=1.0/2.0
	
	EEcm=pow(Ecm0,2.0)
	qq=(EEcm-pow(2.0*mpi,2.0))/4.0
	qout=sqrt(qq)
	x=qq*pow(L/(twopi),2.0)
	##################################
	(l,m)=(0,0)
	
	Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
	amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
	c00=amp*Z0
	if irrep !='T1mm':
		##################################
		(l,m)=(2,0)
		Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
		amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
		c20=amp*Z0
		##################################
		(l,m)=(2,2)
		Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
		amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
		c22=amp*Z0
		################################## 
	else:
		c20,c22=0,0
	
 
	alpha0=alphas(irrep)
	cout=c00+((c20*alpha0[0]+c22*alpha0[1])/qq)
	qcotphi=-4.0*pi*cout
		
	QC=(1.0/qcotphi)+(1.0/fqcotps(Ecm0))
	return QC


def Luscher_ps_no_interp_vec(Ecms,mpi,irrep,d,L):
	
	if len(Ecms)>0:
		out=[]
		for i0 in range(len(Ecms)):
			Ecm0=Ecms[i0]
			out.append(Luscher_ps_no_interp(Ecm0,mpi,irrep,d,L))
		return array(out)
		
	else:
		return Luscher_ps_no_interp(Ecms,mpi,irrep,d,L)
		
def Luscher_ps_no_interp(Ecm0,mpi,irrep,d,L):

 
	"""
	####################
	"""
	
	gamma0=Zlm2.gamma(Ecm0,d,L)
	alphapipi=1.0/2.0
	
	EEcm=pow(Ecm0,2.0)
	qq=(EEcm-pow(2.0*mpi,2.0))/4.0
	qout=sqrt(qq)
	x=qq*pow(L/(twopi),2.0)
	##################################
	(l,m)=(0,0)
	
	Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
	amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
	c00=amp*Z0
	if irrep !='T1mm':
		##################################
		(l,m)=(2,0)
		Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
		amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
		c20=amp*Z0
		##################################
		(l,m)=(2,2)
		Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
		amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
		c22=amp*Z0
		################################## 
	else:
		c20,c22=0,0
	
 
	alpha0=alphas(irrep)
	cout=c00+((c20*alpha0[0]+c22*alpha0[1])/qq)
	qcotphi=-4.0*pi*cout
	
	cotd=4.0*pi*cout/qout	

	tand=1.0/cotd
	Luscher_delta=arctan(tand)
 	return Luscher_delta




def Zboost_E_no_interp(Ecm0,mpi,irrep,d,L,l,m):

 
	"""
	####################
	"""
	
	gamma0=Zlm2.gamma(Ecm0,d,L)
	alphapipi=1.0/2.0
	
	EEcm=pow(Ecm0,2.0)
	qq=(EEcm-pow(2.0*mpi,2.0))/4.0
	qout=sqrt(qq)
	x=qq*pow(L/(twopi),2.0)
	##################################
	amp=sqrt(4.0*pi)*pow(twopi/L,l-2)/(gamma0*pow(L,3.0))
	Z0=Zlm2.Zboost(x,d,l,m,alphapipi,gamma0)
	return  amp*Z0
	
 



def Rinv_no_inter(Ecminp,mpi,irrep,d,ps,qcotps,Eps,ens,L):

	"""
	####################
	first, we define an interpolating function for the phase shift, which have been evaluated in the energy range Eps
	"""
	fps=interp1d(Eps,ps)
	
	"""
	####################
	"""
	
	
	dEcm=0.00001
	out=ones(len(Ecminp))
	
	for n0 in range(len(Ecminp)):
		Ecm0=Ecminp[n0]
		derR=(no_interp_QC(Ecm0+dEcm,mpi,irrep,d,qcotps,Eps,ens,L)-no_interp_QC(Ecm0,mpi,irrep,d,qcotps,Eps,ens,L))/dEcm
		cos2ps=pow(cos(fps(Ecm0)),2.0)
		P=d*twopi/(L)
		PP=dot(P,P)
		Eftot=sqrt(pow(Ecm0,2.0)+PP)
		Einitial=mpi
		
		Ni=2.0*Einitial
		Nf=2.0*Eftot
		xi=1.0/2.0
		Nif=sqrt(Ni*Nf*xi)
		
		
		out[n0]=sqrt(cos2ps*abs(16.0*pi*Eftot*Einitial*derR))/Nif
		Ecm0

	return out 



"""
#################################################
#################################################
#################################################
"""

"""
Interpolating function for the raw LL factor
"""
def persona_interp(EcmRange,Vinp,Eout,ens):
	"""
	here I define a brute force interpolation which just takes the average of nearing neighbors
	EcmRange = range of input energies where Vinp is defined
	Vinp = is the function to interpolate
	Eout = is the energy where it is going to be evaluated 
	"""
	EcmL=list(EcmRange)

	if shape(Eout)==():
		"meaning, if Eout is a number"
		Ecm0s=nsmallest(2, EcmL, key=lambda x: abs(x-Eout))
		"Ecm0s is the two energies in EcmL that are closest to Eout"
		out=0
		ne0=EcmL.index(Ecm0s[0])
		ne1=EcmL.index(Ecm0s[1])
		return (Vinp[ne0]+Vinp[ne1])/2.0

	else:
		"is Eout is a list, we just loop over its components and repeat what we did above"
		out=zeros(len(Eout))
		for i0 in range(len(Eout)):	
			E0=Eout[i0]
			Ecm0s=nsmallest(2, EcmL, key=lambda x: abs(x-E0))
			ne0=EcmL.index(Ecm0s[0])
			ne1=EcmL.index(Ecm0s[1])
			out[i0]=(Vinp[ne0]+Vinp[ne1])/2.0
		return out
		


def alphas(irrep):
	if irrep=='T1mm':
		alpha0=[0,0]

	if irrep=='D4A1d001':
		alpha0=[2.0/sqrt(5.0),0]

	if irrep=='D4E2d001':
		alpha0=[-1.0/sqrt(5.0),0]

	if irrep=='D2A1d011':
		alpha0=[-1.0/sqrt(5.0),-I*sqrt(6.0/5.0)]

	if irrep=='D2B1d011':
 		alpha0=[-1.0/sqrt(5.0),I*sqrt(6.0/5.0)]	

	if irrep=='D2B2d011':
		alpha0=[2.0/sqrt(5.0),0]
	
	if irrep=='D3A1d111':
		alpha0=[0,-2.0*I*sqrt(6.0/sqrt(5.0))]
	
	if irrep=='D3E2d111':
		alpha0=[0,I*sqrt(6.0/sqrt(5.0))]
	
	return alpha0	