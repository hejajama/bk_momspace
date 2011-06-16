#!/usr/bin/python
import os

yvals = [1,2,4,6,8,10,20,30,40]
postfix=""
threads=1
ic = "INVPOWER"
method = "BRUTEFORCE"
datafile = "output_invpower_bruteforce_maxy40"

for y in yvals:
    cmd = "OMP_NUM_THREADS=" + str(threads) + " ./bk -ic " + ic \
        + " -method " + method + " -mode LOGLOG_DERIVATIVE -y " + \
        str(y) + " -data output_invpower_bruteforce_maxy40"
    if postfix!="":
        cmd = cmd + "_" + postfix
    cmd = cmd + ".dat"

    dir = postfix
    if postfix=="":
        dir = "plain"
    output = "data/loglogder/" + dir + "/loglogder_y" + str(y)

    cmd = cmd + " > " + output

    print cmd
    os.system(cmd)
