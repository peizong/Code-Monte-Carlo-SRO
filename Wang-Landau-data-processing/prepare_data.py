
import numpy as np

def prepare_data(outcar,begin,end):
  bin_num=30
  E_arr=np.linspace(begin,end,bin_num)
  dos,sro=[],[]
  col_file_head=2
  count=0
  with open(outcar,'r') as in_file: #sys.argv[1]
    for line in in_file:
      if count<bin_num+col_file_head and count > col_file_head-1:
        le=line.split()
        dos.append(le[0])
      if count>bin_num*2+col_file_head+2-1:
        le=line.split(",")
        for i in range(0,len(le)-1):
          le[i]=round(float(le[i]),4)
        sro.append(le[0:len(le)-1])
      count += 1
  for i in range(0,len(E_arr)):
    E_arr[i]=round(E_arr[i],4)
    line_i=str(E_arr[i])+"	"+str(dos[i])
    for j in [0,1,2,4,5,8]: #range(0,len(sro[i])-1):
      line_i +="	"+str(sro[i][j])
    line_i +="	"+str(sro[i][len(sro[i])-1])
    print(line_i)
if __name__ == "__main__":
  begin,end=-0.1213,-0.0795
  outcar="log.wlce"
  prepare_data(outcar,begin,end)
