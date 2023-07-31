#!/usr/bin/python

import numpy as np
from operator import itemgetter

filename="predstr.out-nicocr"

#data=np.genfromtxt("predstr.out-nicocr")
data=[]
with open(filename,'r') as in_file:
  for line in in_file:
    ll=line.split()
    ll[3]=float(ll[3])
    data.append(ll)
data_sorted=sorted(data,key=itemgetter(3))
for line in data_sorted:
  for element in line:
    print element,"	",
  print '\n'
#data_array=np.asarray(data)
lb,ub=data_sorted[0][3],data_sorted[len(data_sorted)-1][3] #np.min(data_array[:,3]),np.max(data_array[:,3])
N=60
ddos=(ub-lb)/(N-1)
dos=np.linspace(0,0,N)
index=0
for line in data_sorted:  
  if line[3]>lb+index*ddos:
    index +=1
    print dos[index-1]
  dos[index] +=1

print dos.transpose()

