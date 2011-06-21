#!/usr/bin/python
#-*- coding: UTF-8 -*-
import os
import sys
import math
import pdb
sys.path.append("/home/hejajama/lib/")
sys.path.append("../src/")

from matplotlibhelper import *

#yvals = [1,2,4,6,8,10,20,30,40]
yvals = [2,6,10,30]

# [text id style, width]
modes = [
    ["", "plain", dashes[1], 1],
    ["RC", "rc", dashes[0], 2] ]

ic = "invpower"
icfile = ic + "/loglogder/ic"

minx=1e-4
maxx=1e10
miny=-1.1
maxy=0.1

scale=0  # 1: scale with Q_s, 0: don't


def sqrt_data(x):
    for i in range(len(x)):
        x[i]=math.sqrt(x[i])

fig = figure()
p1=fig.add_subplot(111)
if (scale==1):
    xlabel(r"$k_T/Q_s$")
else:
    xlabel(r"$k_T$")
ylabel(r"d ln $N(k^2)/$ d ln $k^2$")

# Plot initial condition
xdata=[]
ydata=[]
readfile(icfile, xdata, ydata)
sqrt_data(xdata)
p1.semilogx(xdata, ydata, label="IC", linestyle=dashes[1], linewidth=3)


for y in yvals:
    for mode in modes:
        data = ic + "/loglogder/" + mode[1] + "/loglogder_y" + str(y)
        xdata=[]
        ydata=[]
        readfile(data, xdata, ydata)
        sqrt_data(xdata)

        
        if (scale==1):
            params = read_parameters(data)
            satscale = float(params[0])
            print "satscale at y=" + str(y) + " for " + mode[0] + " is " + str(satscale) \
                + " minkt/q=" + str(xdata[0]/satscale) + " maxkt/q=" + str(xdata[len(xdata)-1]/satscale)
            
            # Scale x by saturation scale
            scale_list(xdata, 1.0/satscale)
        
        lbl = "y=" + str(y)
        if mode[0]!="":
            lbl = lbl + ", " + mode[0]
        p1.semilogx(xdata, ydata, label=lbl, linestyle = mode[2], linewidth = mode[3])

fig.suptitle(r"Initial condition " + ic)
axis([minx,maxx,miny,maxy])
leg=legend(prop=dict(size=textsize),labelspacing=0.001)
leg.draw_frame(False)
file = "loglogder.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()
