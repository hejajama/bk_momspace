#!/usr/bin/python
#-*- coding: UTF-8 -*-
# BK Equation solver
# Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011

# FT given datafile from k space to r space
# Data is given in form
# ktsqr N(ktsqr)
# N(r) = r^2/(2\pi) \int d^2k exp(-i k.r) N(k)
# As N(k)=N(ktsqr),
# N(r) = r^2 \int dk k J_0(kr) N(k)

# Filename is given as a first argument

import sys
sys.path.append("./src/")
from matplotlibhelper import *
from scipy.integrate import simps
from scipy.special import j0

kt=[]
ktn=[]

rpoints = 100   # Evaluate for rpoitns different dipole sizes
minr = 0.05
maxr = 10.0
r=[]
rn=[]
dr = (maxr-minr)/rpoints
for i in range(rpoints):
    r.append(minr+dr*i)
        
fig = figure()
p1=fig.add_subplot(111)

files = []
for i in range(len(sys.argv)):
    if i==0:
        continue
    kt=[]
    ktn=[]
    readfile(sys.argv[i], kt, ktn, 1)
    rn=[]


    for i in range(len(kt)):
        kt[i]=sqrt(kt[i])

    for i in range(rpoints):
        tmpktn=[]
        for j in range(len(ktn)):
            tmpktn.append(r[i]*r[i]*kt[j]*j0(kt[j]*r[i])*ktn[j])
        res = simps(tmpktn, kt)
        rn.append(res)

    p1.semilogx(r, rn)

pp = PdfPages("out.pdf")
savefig(pp, format="pdf")
pp.close()
