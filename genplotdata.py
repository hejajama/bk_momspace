#!/usr/bin/python
import os

yvals = [2,4,6,8,10,20,30,40,50]
postfix="rc"
threads=1

#ic: parameter for ./bk , dirname
#ic = ["INVPOWER", "invpower"]
ic = ["FTIPSAT", "ftipsat"]

method = "BRUTEFORCE"
datafile = "output_ftipsat_bruteforce_maxy50"

for y in yvals:
    cmd = "OMP_NUM_THREADS=" + str(threads) + " ./bk -ic " + ic[0] \
        + " -method " + method + " -mode SINGLE_PLOT -y " + \
        str(y) + " -data " + datafile
    if postfix!="":
        cmd = cmd + "_" + postfix
    cmd = cmd + ".dat"

    dir = postfix
    if postfix=="":
        dir = "plain"
    output = "data/" + ic[1] + "/plot/" + dir + "/plot_y" + str(y)

    cmd = cmd + " > " + output

    print cmd
    os.system(cmd)
