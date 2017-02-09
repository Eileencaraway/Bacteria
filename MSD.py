import glob
from scipy import *
import re
import matplotlib.pyplot as plt

x=list()
y=list()
count=0
file = open('MSD.txt')
for line in file:
    line= line.rstrip()
    count+=1
    if count >1:
        num= line.split()
        x.append(float(num[0]))
        y.append(float(num[1]))



## read from multiple files
#filename = glob.glob(('dumpfile.*.txt'))
#filename = sort(filename)
## list is extendable while array cannot
#x=list()
#y=list()
#print(filename)
## f refer to each file
#for f in filename:
#    file= open(f)
    ##print(f)
#    count=0  ## count for the line , so this file only suitable for the same pattern data
#    for line in file:
#        line= line.rstrip()
#        count+=1
#        if count==10:
#            num= line.split()
#            x.append(float(num[0]))
#            y.append(float(num[1]))

## change list to array
x=array(x)
print(x)
y=array(y)
print(y)


plt.loglog(x,y)
plt.show()
