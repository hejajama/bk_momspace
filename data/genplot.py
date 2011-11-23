#!/usr/bin/python
#-*- coding: UTF-8 -*-
import os
import sys
import math
sys.path.append("/home/hejajama/lib/")
sys.path.append("../src/")

from matplotlibhelper import *

#yvals = [1,2,4,6,8,10,20,30,40]
yvals = [0,2,10,40,50]

# [text id style, width]
modes = [
    #["invpower", "plain", dashes[1], 1],
    #["invpower rc", "rc", dashes[1], 2],
    ["ftipsat", "../../ftipsat/plot/plain", dashes[2], 1],
    #["ftipsat rc", "../../ftipsat/plot/rc", dashes[2], 2]
    #["KC", "kc", dashes[0], 2]#, 
    #["parent", "rc", dashes[3], 2]
]

ic = "ftipsat"

icfile = ic + "/plot/ic"

minx=1e-4
maxx=1e20
miny=1e-15
maxy=100


def sqrt_data(x):
    for i in range(len(x)):
        x[i]=math.sqrt(x[i])

fig = figure()
p1=fig.add_subplot(111)
xlabel(r"$k_T$")
ylabel(r"$N(k_T)$")

# Plot initial condition
xdata=[]
ydata=[]
readfile_xy(icfile, xdata, ydata)
sqrt_data(xdata)
#p1.semilogx(xdata, ydata, label="IC", linestyle=dashes[1], linewidth=3)

color=0
for y in yvals:
    for mode in modes:
        data = ic + "/plot/" + mode[1] + "/plot_y" + str(y)
        xdata=[]
        ydata=[]
        print data
        readfile_xy(data, xdata, ydata)
        sqrt_data(xdata)
        lbl = "y=" + str(y)
        if mode[0]!="":
            lbl = lbl + ", " + mode[0]
        p1.loglog(xdata, ydata, label=lbl, linestyle = mode[2], linewidth = mode[3]) #, color=colors[color])
        color=color+1

fig.suptitle(r"Initial condition " + ic)
axis([minx,maxx,miny,maxy])
leg=legend(prop=dict(size=textsize),labelspacing=0.001, loc=1)
leg.draw_frame(False)
file = "plot.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()
