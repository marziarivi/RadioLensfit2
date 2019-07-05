"""
  shear.py
  Author: Marzia Rivi (2018)

  arguments: -nf number of files with name ellipticities<n>.txt

  Compute shear as a weighted mean of the galaxies ellipticity
  use bootstrap to compute standard deviation
  Use data simulations where shapes are generated to be zero on average.
  Opposite ellipticities are consecutives. Therefore for each bad data,
  the opposite is removed too from the shear computation.

"""

import sys
import argparse
import math
import numpy as np
#from pylab import *
import astropy.stats as astro

parser = argparse.ArgumentParser(description='bootstrap_std')
parser.add_argument('-nf',dest='nfiles', type=int, default=1, help='Number of measurement files')

args = parser.parse_args(sys.argv[1:])

e_max = 0.804   # ellipticity cutoff
e_0 = 0.0256    # circularity parameter
a = 0.2539      # dispersion
Nfactor = 2.43180252985281 # A=1/0.4112176 normalization factor

def prior_ellipticity(ee1,ee2):
    emod = np.sqrt(np.multiply(ee1,ee1)+np.multiply(ee2,ee2))
    p_e = Nfactor*np.divide(np.multiply(emod,(1.-np.exp((emod-e_max)/a))),np.multiply(1.+emod,np.sqrt(np.multiply(emod,emod)+e_0*e_0)))

# read ellipticities
me1=[]
me2=[]
err1=[]
err2=[]
SNR=[]
w=[]
k=0
remove = 0
nfiles=args.nfiles
bad = 0
num=0
ngal = 10000

try:
    while num<nfiles:
        name = "ellipticities%d.txt"%(num)
        print name
        data=np.loadtxt(name, skiprows=1, delimiter='|',usecols=(1,2,3,4,5,6))
        i=0
        while i<len(data) and k<ngal:
            error1 = data[i,1]
            error2 = data[i,3]
            var = data[i,4]
            SNRvalue = data[i,5]
            
            if var > 1e-5 and error1>1e-3 and error2>1e-3 and SNRvalue >= 10:
                if remove == 1:
                    me1[k] = data[i,0]
                    err1[k] = error2
                    me2[k] = data[i,2]
                    err1[k] = error1
                    w[k] = var*e_max*e_max/(e_max*e_max-2*var)
                    SNR[k] = SNRvalue
                else:
                    me1.append(data[i,0])
                    err1.append(error2)
                    me2.append(data[i,2])
                    err1.append(error1)
                    w.append(var*e_max*e_max/(e_max*e_max-2*var))
                    SNR.append(SNRvalue) # 8
                k = k+1
                remove = 0
            else:  # remove the opposite too
                if i%2 == 0:  # even line: opposite is the previous one
                    remove = 1  # remove previous element by overlapping it with the following
                    k = k-1
                else:         # odd line: opposite is the next one
                    i = i+1     # skip next data line
                print "bad measure:",data[i,:]
                bad = bad + 1
            i=i+1
        num = num + 1
except:
    print 'ERROR!'

# if last line is bad measure
if remove == 1:
    k = k-2
    me1 = np.delete(me1,k)
    me2 = np.delete(me2,k)
    SNR = np.delete(SNR,k)
    w = np.delete(w,k)

print "mean SNR ",np.mean(SNR), "median SNR ",np.median(SNR)
print "ngal: ", len(me1)," bad: ",bad
print "min SNR",np.min(SNR)

print 'measured mean: ',np.mean(me1), np.mean(me2)

# compute shape noise
pe = prior_ellipticity(me1,me2)
mean1 = np.average(me1,weights=pe)
mean2 = np.average(me2,weights=pe)
mean12 = np.average(np.multiply(me1,me2),weights=pe)
var12 = mean12 - mean1*mean2
mean11 = np.average(np.multiply(me1,me1),weights=pe)
mean22 = np.average(np.multiply(me2,me2),weights=pe)
var1 = mean11 - mean1*mean1
var2 = mean22 - mean2*mean2

sigma_shape = np.sqrt(var1*var2-var12*var12)
print "shape noise: ",sigma_shape

# compute shear
norient = 2
ngal = len(me1)/norient
print norient*ngal
w = np.divide(1.,w+sigma_shape)
data = np.array(me1).reshape(ngal,norient)
weight = np.array(w).reshape(ngal,norient)  #np.ones((ngal,norient))
alldata=np.concatenate((data,weight),axis=1)
data = np.array(me2).reshape(ngal,norient)
weight = np.array(w).reshape(ngal,norient)
# concatenate all and reshape as 4-dim array of dimensions (ngal,ncomponents,ntypes,norient), where types = measure, weight
alldata=np.concatenate((alldata,data,weight),axis=1).reshape(ngal,2,2,norient)

bootnum = 1000
bootsamples = astro.bootstrap(alldata,bootnum)

e1_shears=np.zeros(bootnum)
e2_shears=np.zeros(bootnum)
for i in range(bootnum):
    e1_shears[i] = np.average(bootsamples[i,:,0,0,:],weights=bootsamples[i,:,0,1,:])
    e2_shears[i] = np.average(bootsamples[i,:,1,0,:],weights=bootsamples[i,:,1,1,:])

print 'e1_shear: mean = ',np.mean(e1_shears),' std = ',np.std(e1_shears)
print 'e2_shear: mean = ',np.mean(e2_shears),' std = ',np.std(e2_shears)



