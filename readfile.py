# this python code is used for read data file
#and use the data to plot a 3D picture of particles at that moment
#intend to use for lammps dump file "dumpfile.txt"

import re
import matplotlib.pyplot as plt
from scipy import *

##x=list()
##y=list()
DD=list()

log= open("MSD_fixprint.txt")
for line in log:
    line=line.rstrip()
    if re.search("[0-9]+[.][0-9]+",line):
        num=line.split()
        ##x.append(float(num[0]))
        ##y.append(float(num[1]))
        DD.append(float(num[0]))

X=zeros(len(DD))
for i in range(len(DD)):
    X[i]=i*100

plt.loglog(DD,X)
plt.show()
