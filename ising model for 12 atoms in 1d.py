# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 18:53:44 2021

@author: Pc
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import sparse

def kroner (marray):
  i=1 
  result = marray[0].copy()
  while i<len (marray):
    result =np.kron(result,marray[i])
    i+=1
  return result

def math_ave(operator,psi):
  psi_dagger = np.conj(psi)
  result=np.matmul(operator,psi)
  result=np.matmul(psi_dagger, result)
  return result

sigma_x = (np.array([[0,1],[1,0]]))*0.5
sigma_y = (np.array([[0,-1j],[1j,0]]))*0.5
sigma_z = (np.array([[1,0],[0,-1]]))*0.5
#sigma_x= sparse (s_x)
identity=np.eye(2)

def sigma12_x (i):
  n=12
  s_x=kroner([np.eye(2**(i)),sigma_x,np.eye(2**(n-1-i))])
  return s_x

def sigma12_y(j):
    n=12
    s_y=kroner([np.eye(2**(j-1)),sigma_y,np.eye(2**(n-j))])
    return s_y

def sigma12_z(k):
    n=12
    s_z=kroner([np.eye(2**(k)),sigma_z,np.eye(2**(n-1-k))])
    return s_z

i=1
H1=np.matmul(sigma12_z(0),sigma12_z(1))+np.matmul(sigma12_z(11),sigma12_z(0))
F= sigma12_x(0) + sigma12_x(11)
while i<11:
    H1=H1+np.matmul(sigma12_z(i),sigma12_z(i+1))
    F = F+sigma12_x(i)
    i+=1
K=H1
L=F

def SzSz(l):
    h = l*0.1
    H = -K-h*L
    eigval,eigvec=np.linalg.eigh(H)
    szsz = np.zeros ((12))
    for i in range (0,12):
        szsz[i] =math_ave(np.matmul(sigma12_z(0),sigma12_z(i)),eigvec[:,0])- math_ave(sigma12_z(0),eigvec[:,0])*math_ave(sigma12_z(i),eigvec[:,0])
        
    return szsz

x=np.arange (0,12)
plt.plot(x,SzSz(1))
plt.plot(x,SzSz(2))
plt.plot(x,SzSz(3))
plt.plot(x,SzSz(4))
plt.plot(x,SzSz(5))
plt.plot(x,SzSz(6))
plt.plot(x,SzSz(7))
plt.plot(x,SzSz(8))
plt.plot(x,SzSz(9))
plt.plot(x,SzSz(10))




plt.show()