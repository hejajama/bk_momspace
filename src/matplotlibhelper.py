#!/usr/bin/python
# -*- coding: UTF-8 -*- 

# Reads data from a file in format
# x y

import sys
import pdb
sys.path.append('/home/hejajama/lib/python');

from matplotlib.pyplot import *
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')

# 
dashes = ['--', #    : dashed line
          '-', #     : solid line
          '-.', #   : dash-dot line
          ':', #    : dotted line
           '-']
colors = [ "black", "blue", "red", '0.35', '0.75', '0.5', '0.1', '0.2', '0.6', '0.8']

textsize = 16

# Possibility to use , as a desimal separator
import locale
locale.setlocale(locale.LC_ALL,"fi_FI")
fnx = lambda x : locale.format("%.1f", x)
fnx2 = lambda x : locale.format("%.2f",x)  
fnx3 = lambda x : locale.format("%.3f",x) 

# Read data and multiply all y values by factor m
def readfile(file, xlist, ylist,m=1.0,err=[]): 
	f = open(file,"r")
	lines=f.readlines()
	n=len(lines)

	for i in range(n):
		s=lines[i].split()
		if (len(s)==2 and s[0][0]!="#"):
			xlist.append(float(s[0]))
			ylist.append(float(s[1])*m)
		if (len(s)==3 and s[0]!="#"):
			xlist.append(float(s[0]))
			ylist.append(float(s[1])*m)
			err.append(float(s[2])*m)
	f.close()

# Read data from two first columns, ignore others
def readfile_xy(file, xlist, ylist, ycol=1):
    f = open(file,"r")
    lines=f.readlines()
    n=len(lines)

    for i in range(n):
        #print i #str(len(lines[i])) + str(i)
        s=lines[i].split()
        if (len(lines[i])>1 and s[0][0]!="#"):
            #pdb.set_trace()
            xlist.append(float(s[0]))
            ylist.append(float(s[ycol]))
    f.close()

# Read parameters from the file
# Parameters are strings/numbers appearing on a line formatted as
# "### parameter"
def read_parameters(file):
    f = open(file, "r")
    lines = f.readlines()
    n = len(lines)
    params=[]
    for i in range(n):
        if (len(lines[i])<3):
            continue
        s = lines[i].split()
        if (s[0]=="###"):
            params.append(s[1])
    f.close()
    return params

# Read data and multiply all y values by factor m
# syntax in file: x y ystaterr posyerr negyerr
def readfile_errorarray(file, xlist, ylist,m=1.0,errorarray=[]): 
	f = open(file,"r")
	lines=f.readlines()
	n=len(lines)
	tmpposerr=[]
	tmpnegerr=[]
	
	for i in range(n):
		s=lines[i].split()
		if (len(s)==5 and s[0]!="#"):
			xlist.append(float(s[0])*m)
			ylist.append(float(s[1])*m)
			staterr = float(s[2])
			possysterr = float(s[3])
			negyerr = float(s[4])
			tmpposerr.append(sqrt(staterr**2+possysterr**2)*m)
			tmpnegerr.append(sqrt(staterr**2+negyerr**2)*m)
	
	f.close()
	errorarray.append([])
	for j in tmpnegerr:
		errorarray[0].append(j)
	errorarray.append([])
	for j in tmpposerr:
		errorarray[1].append(j)

# Scale list by a constant
def scale_list(list, const):
    for i in range(len(list)):
        list[i] = list[i]*const
