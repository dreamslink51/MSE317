#!/usr/bin/env python
import numpy as np
from numpy import linalg as LA
import argparse

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to make calculate mobility from snapshot_KMC.x.')
parser.add_argument('-dfile', action="store", default="distance.dat", dest='d_file',
                    help='Name of the distance file. (Default=distance.dat)')
parser.add_argument('-tfile', action="store", default="tcheck.dat", dest="t_file",
                    help='Name of the tcheck file. (Default=tcheck.dat)')
parser.add_argument('-temp', action="store", default=283, dest="temp",
                    help='Temperatuer. (Default=283K)')
parser.add_argument('-nsample', action="store", default=10, dest="step",
                    help='Number of KMC runs. (Default=10)')

prm = parser.parse_args()


# load displacement
dist = np.loadtxt(prm.d_file)
# load time
tck = np.loadtxt(prm.t_file)


step=np.int(prm.step) # number of KMC runs
T=np.float(prm.temp) # Temperature
n=np.int(np.floor(dist.shape[0]/step)) # number of MD snapshots

# calculated mobility eigen values for n MD snapshot
# unlike matlab, numpy outputs eigenval in an 1d array.
matrices = np.zeros([n,3],dtype=np.float)
matrices_sorted = np.zeros([n,3],dtype=np.float)
# calculated mobility eiegen vectors for n MD snapshot
matrices_V = np.zeros([3*n,3],dtype=np.float)

# averaged mobility eigen values for n MD snapshot
avedir = np.zeros(3)
# averaged mobility eiegen vectors for n MD snapshot
aveeig = np.zeros(3)

for i in range(n):
    sta = i*step
    fini = (i+1)*step
    # print dist[sta:fini,:].T
    # MSD/(2*time*k*Temperature)
    mobmat = np.cov(dist[sta:fini,:].T)*1e-8**2 \
             /(2*np.mean(tck[sta:fini])*8.617343e-5*T)

    # Diag mobmat
    # note that the order of eigval is randomized
    D, V = LA.eig(mobmat)

    # fill in the mobility matrix
    matrices[i,:] = D
    # sort eigenvals
    D.sort()
    matrices_sorted[i,:] = D
    matrices_V[i*3:(i+1)*3,:] = V.T

    # add all to ave (not sure why?)
    avedir=avedir+V
    aveeig=aveeig+D


# print eigenvalue and eigenvector for each MD snapshot
# note that the eigenvec may change order.
for i in range(n):
    print "MD step:",i
    print 'Mobility[1]: %15.12f ; Eigenvec: %15.8f, %15.8f, %15.8f'\
        %(matrices[i,0],matrices_V[i*3,  0],matrices_V[i*3,  1],matrices_V[i*3,  2])
    print 'Mobility[2]: %15.12f ; Eigenvec: %15.8f, %15.8f, %15.8f'\
        %(matrices[i,1],matrices_V[i*3+1,0],matrices_V[i*3+1,1],matrices_V[i*3+1,2])
    print 'Mobility[3]: %15.12f ; Eigenvec: %15.8f, %15.8f, %15.8f'\
            %(matrices[i,2],matrices_V[i*3+2,0],matrices_V[i*3+2,1],matrices_V[i*3+2,2])

# calculate mean value
mx, my, mz = np.mean(matrices_sorted[:,:],axis=0)
print '\nSorted (small to large) Mobilities:'
print 'Mobility[1]: %15.12f'%(mx)
print 'Mobility[2]: %15.12f'%(my)
print 'Mobility[3]: %15.12f'%(mz)
