
import numpy as np

begin,end=-0.1013,-0.0695
bin_num=30
bin_size=(end-begin)/(bin_num-1)
E_arr=np.linspace(begin,end,30)
dos_arr=np.genfromtxt("dos.dat")
for i in range(0,len(E_arr)):
  print(str(E_arr[i])+"	"+str(dos_arr[i]))
