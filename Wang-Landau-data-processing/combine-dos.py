
#!/usr/bin/python
import numpy as np
from scipy.interpolate import splev,splrep,interp1d #, spline
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats

class ppdos:
  """post processing for one dos file """
  def __init__(self, N_bin, N_sec,sec_ID0,  N_and, dos_for_N_atoms,skip_header_file):
    self.N_bin=N_bin
    self.N_sec=N_sec
    self.sec_ID0=sec_ID0
    self.N_and=N_and
    self.N=dos_for_N_atoms
    self.E_range=[-0.1213,-0.0595]
    self.data=[]
    self.whole_dos=[]
    self.dos_spline=[]
    self.skip_header=skip_header_file
  def read_data(self):
    N_bin=self.N_bin
    N_whole=self.N_bin*self.N_sec-self.N_and*(self.N_sec-1)
    dE=(self.E_range[1]-self.E_range[0])/(N_whole-1)
    for fi in range(self.sec_ID0,self.sec_ID0+self.N_sec):
      self.filename="dos_"+str(fi)
      # find local energy range for each dos section
      El=self.E_range[0]+(self.N_bin-self.N_and)*fi*dE
      Eu=El+dE*(N_bin-1)
      E_range=[El,Eu]
      dos=np.array([np.zeros(N_bin),np.zeros(N_bin),np.zeros(N_bin),np.zeros(N_bin),\
           np.zeros(N_bin),np.zeros(N_bin),np.zeros(N_bin),np.zeros(N_bin)])
      for i in range(0,N_bin):
        dos[0][i]=E_range[0]+i*dE 
      with open(self.filename,'r') as dosfile:
       #speical line locators:
        dos_locator=3
        sro_locator=65
        i=0
        for line in dosfile:
          if i>dos_locator-2 and i<dos_locator+(self.N_bin-1):
            ll=line.split()
            dos[1][i-(dos_locator-1)]=ll[0]
            #set the starting point of the first dos section as zero
            #if fi==self.sec_ID0: 
            #  dos[1][i-(dos_locator-1)]-=min(dos[1])
          if i>sro_locator-2:
            ll=line.split(",")
            #dos[2][i-(sro_locator-1)],dos[3][i-(sro_locator-1)],dos[4][i-(sro_locator-1)],\
            #dos[5][i-(sro_locator-1)]=ll[0],ll[1],ll[2],ll[3]
            for k in range(0,6):
              k_selected=[0,1,2,4,5,8] #choose CoCo,CoCr,CoNi,CrCr,CrNi,NiNi paris
              kp=k_selected[k]
              dos[k+2][i-(sro_locator-1)]=ll[kp]
          i+=1
      dat=np.asarray(dos).transpose()
      dat[:,0]=dat[:,0] #*self.N
      self.data.append(dat)
    #set the starting point of the first dos section as zero
    
    for i in range(0,len(dat[0,:])):
      self.whole_dos.append(np.zeros(N_whole)) #,np.zeros(N_whole)]
      self.whole_dos[i][0:self.N_bin]=self.data[0][:,i]
  def line_regression(self):
    for i in range(0, self.N_sec-1):
      x,y=self.data[i][(self.N_bin-self.N_and):self.N_bin,1], self.data[i+1][0:self.N_and,1]
      solution=stats.linregress(x,y)
      new_dos=(self.data[i+1][:,1]-solution[1])/solution[0]
      new_dos[0:self.N_and]=(new_dos[0:self.N_and]+ \
                            self.data[i][(self.N_bin-self.N_and):self.N_bin,1])/2.0
      self.data[i+1][:,1]=new_dos
      offset_i=(self.N_bin-self.N_and)*(i+1)
      self.whole_dos[0][offset_i:(offset_i+self.N_bin)]=self.data[i+1][0:self.N_bin,0]
      self.whole_dos[1][offset_i:(offset_i+self.N_bin)]=self.data[i+1][0:self.N_bin,1]
   # slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)    
   # print self.whole_dos #solution[0],solution[1]  #slop, intercept, r_value
  def combine_sro(self):
    for i in range(0,self.N_sec-1):
      for j in range(2,len(self.whole_dos)):
        x,y=self.data[i][(self.N_bin-self.N_and):self.N_bin,j], self.data[i+1][0:self.N_and,j]
        new_sro=(x+y)/2.0
        self.data[i+1][0:self.N_and,j]=new_sro
        offset_i=(self.N_bin-self.N_and)*(i+1)  
        self.whole_dos[j][offset_i:(offset_i+self.N_bin)]=self.data[i+1][0:self.N_bin,j]
  def print_combined_data(self):
    print("E","	","dos","	","Co-Co","	","Co-Cr","	","Co-Ni","	","Cr-Cr","	","Cr-Ni","	","Ni-Ni")
    for i in range(0,len(self.whole_dos[0])):
      line_i=''
      for j in range(0,len(self.whole_dos)):
        if j<len(self.whole_dos)-1:
          line_i +=str('{:.4f}'.format(self.whole_dos[j][i]))+"	"
        else: line_i +=str('{:.4f}'.format(self.whole_dos[j][i])) #+"\n"
      print(line_i)
        #if j==len(self.whole_dos)-1: print "\n"
  def cal_dos_spline(self):
    x_new =np.linspace(np.min(self.whole_dos[0]), np.max(self.whole_dos[0]), num=0.2*len(self.whole_dos[0]), endpoint=True)
    #mid=splrep(self.whole_dos[0],self.whole_dos[1])
    #y_new=splev(x_new,mid)
    #y_new=spline(self.whole_dos[0],self.whole_dos[1],x_new)
    tck = splrep(self.whole_dos[0],self.whole_dos[1], s=0)
    y_new=splev(x_new,tck) #,der=0)
    self.dos_spline= np.transpose(np.array([x_new,y_new]))
  def plot(self):
    plt.subplot(121)
    plt.plot(self.whole_dos[0],self.whole_dos[1])
    plt.plot(self.dos_spline[:,0],self.dos_spline[:,1])
    plt.subplot(122)
    for i in range(2,len(self.whole_dos)):
      plt.plot(self.whole_dos[0],self.whole_dos[i])
    plt.show()
#-----------ternary-------------------
caseNiCoCr=ppdos(30,6,0,6,192,1)
caseNiCoCr.read_data()
caseNiCoCr.line_regression()
caseNiCoCr.cal_dos_spline()
caseNiCoCr.combine_sro()
caseNiCoCr.print_combined_data()
caseNiCoCr.plot()
#------------end----------------------



