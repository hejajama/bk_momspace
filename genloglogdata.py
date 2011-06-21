#!/usr/bin/python
import os

yvals = [1,2,4,6,8,10,20,30]
postfix="rc"
threads=1
#ic: parameter for ./bk , dirname
ic = ["INVPOWER", "invpower"]
#ic = ["FTIPSAT", "ftipsat"]

method = "BRUTEFORCE"
datafile = "output_" + ic[1] + "_bruteforce_maxy40"

for y in yvals:
    cmd = "OMP_NUM_THREADS=" + str(threads) + " ./bk -ic " + ic[0] \
        + " -method " + method + " -mode LOGLOG_DERIVATIVE -y " + \
        str(y) + " -data " + datafile
    if postfix!="":
        cmd = cmd + "_" + postfix
    cmd = cmd + ".dat"

    dir = postfix
    if postfix=="":
        dir = "plain"
    output = "data/" + ic[1] + "/loglogder/" + dir + "/loglogder_y" + str(y)

    cmd = cmd + " > " + output

    print cmd
    os.system(cmd)
