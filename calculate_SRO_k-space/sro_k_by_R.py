#!/nfs/apps/Compilers/Python/Anaconda/2.7/bin/python
##!/usr/bin/python

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import copy as copy
import time
import operator
from collections import OrderedDict

class sro:
  def __init__(self,data_file,N,NR): #,period):
    self.in_file=data_file
    self.labels=[] #labels
    self.c=[] 
    self.labels=[]
    self.N=N
    self.NR=NR
    self.period=[]
    #self.solute_pair=solute_pair
    self.coord=[]
    self.atoms_pos=[]
    self.atoms_pos_spins=[]
    self.neighbor_list=[]
    self.sro_real=[]
    self.alpha_k_arr=[]
  def __del__(self):
    print("class deleted")
  def read_data(self):
    count=0
    with open(self.in_file,'r') as in_file: #sys.argv[1]
      for line in in_file:
        le=line.split()
        if count<5 and count>1:
          le[0],le[1],le[2]=float(le[0]),float(le[1]),float(le[2])
          le[0],le[1],le[2]=round(le[0],4),round(le[1],4),round(le[2],4)
          self.coord.append(le)
        if count==5: 
          for i in range(0,len(le)):
            le[i]=float(le[i])
            self.c.append(le[i])
          self.c=1.0*np.asarray(self.c)/np.sum(self.c)
        if count==6: coord_type=le[0]
        if count>6:
          if le !=list():
            self.labels.append(le[3])
            pos=np.array([float(le[0]),float(le[1]),float(le[2])])
            if coord_type=='Direct':
              pos=pos.dot(self.coord)
            #print pos
            le[0],le[1],le[2]=pos[0],pos[1],pos[2]
            le[0:3]=self.move_in_supercell(le[0:3])
            self.atoms_pos.append(le)
        count +=1
    self.labels=list(OrderedDict.fromkeys(self.labels))
    print(self.labels)
    self.period=[self.coord[0][0],self.coord[1][1],self.coord[2][2]]
    #self.atoms_pos_spins=copy.copy(self.atoms_pos)
  def manipulate_sign(self,pos,ith):
    if ith==0: return [pos[0],pos[1],pos[2]]
    elif ith==1: return [-pos[0],pos[1],pos[2]]
    elif ith==2: return [pos[0],-pos[1],pos[2]]
    elif ith==3: return [pos[0],pos[1],-pos[2]]
    elif ith==4: return [-pos[0],-pos[1],pos[2]]
    elif ith==5: return [-pos[0],pos[1],-pos[2]]
    elif ith==6: return [pos[0],-pos[1],-pos[2]]
    elif ith==7: return [-pos[0],-pos[1],-pos[2]]
    else: print("Error! Index out of range!")
  def find_neighbors(self):
    neighbor_list=[]
    for i in self.atoms_pos:
      if np.linalg.norm(i[0:3])==0:
        R0=i
        break
    for i in self.atoms_pos:
      vec_length = np.linalg.norm(i[0:3])
      vec_length = round(vec_length,4)
      if vec_length !=0 and vec_length <= np.min(self.period):
        for j in range(0,8):
          neighbor_list.append([vec_length,self.manipulate_sign(i[0:3],j)])
    neighbor_list.sort(key=operator.itemgetter(0))
    reg_list=[]
    element=[]
    dist=neighbor_list[0][0]
    for i in neighbor_list:
      if dist==i[0]:
       element.append(i[1])
      else:
        reg_list.append([dist,element])
        dist=i[0]
        element=[]
        element.append(i[1])
    for i in range(0,len(reg_list)):      
      id_to_be_removed=[]
      for j in range(0,len(reg_list[i][1])-1):
         for k in range(j+1,len(reg_list[i][1])):
           if np.linalg.norm(np.asarray(reg_list[i][1][j])-np.asarray(reg_list[i][1][k]))==0:
             id_to_be_removed.append(k)
      id_to_be_removed=list(set(id_to_be_removed))
      id_to_be_removed.sort(reverse=True)
      for j in id_to_be_removed:
        reg_list[i][1].pop(j)
      print(i,"	",reg_list[i][0],"	",len(reg_list[i][1]))
    self.neighbor_list=reg_list
  def Kronecker(self,a,b):
    if a==b: return 1
    else: return 0
  def move_in_supercell(self,neighbor_pos):
    pos=np.zeros(3)
    for i in range(0,3):          
      pos[i]=round(neighbor_pos[i],4)
      while(pos[i]-self.coord[i][i] >= 0): pos[i] -= self.coord[i][i]
      while(pos[i] <0 ): pos[i] += self.coord[i][i]
    return pos
  def find_spin(self,NNN_pos):
    spins=np.zeros(len(self.labels))
    for atom in self.atoms_pos:
      atom_pos=np.array(atom[0:3]) #([atom[0],atom[1],atom[2]]) # latt*np.array([i[0],i[1],i[2]])
      dist=np.linalg.norm(atom_pos-NNN_pos)
      if dist<1e-4:
        #for i in range(0,len(self.labels)):
        spins=atom[3] #self.spin_of(atom)
        break
      #else: print "no match"
    #print "check here:", atom,atom_pos, NNN_pos, spins
    if np.sum(spins)==0: print(atom, NNN_pos, "No Match of atomic positions!")
    return spins
  def spin_of(self,atom_pos):
    spins=np.zeros(len(self.labels))
    for i in range(0,len(self.labels)):
      if atom_pos[3]==self.labels[i]: spins[i]=1
    #print atom_pos[3], self.labels, spins
    if np.sum(spins)==0: 
      print("label wrong, here is the details:", atom_pos[3], self.labels)
      return None
    else:
      return spins
  def cal_spin2(self):
    count=np.zeros(self.NR)
    for i in range(0,len(self.atoms_pos)):
      self.atoms_pos[i][3]=self.spin_of(self.atoms_pos[i])
    sro_real=np.zeros([len(self.labels)**2,self.NR])
    #c_loc=np.zeros([len(self.labels),self.NR])
    for i in range(0,len(self.atoms_pos)): 
      spin_i3=self.atoms_pos[i][3]
      for k in range(0,self.NR):
        for j in range(0,len(self.neighbor_list[k][1])):
          neighbor_pos=np.asarray(self.atoms_pos[i][0:3]) + np.asarray(self.neighbor_list[k][1][j])
          neighbor_pos=self.move_in_supercell(neighbor_pos)
          spin_j3=self.find_spin(neighbor_pos)
          for k1 in range(0,len(spin_i3)): 
            for k2 in range(0,len(spin_j3)):
              spin_ik1=spin_i3[k1] 
              spin_jk2=spin_j3[k2] 
              spin_ij=float(spin_ik1)*float(spin_jk2)
              sro_real[k1*len(spin_i3)+k2][k] += 1.0*spin_ij/len(self.neighbor_list[k][1])/len(self.atoms_pos)
    for k in range(0,self.NR):
      for k1 in range(0,len(self.labels)):
        for k2 in range(0,len(self.labels)):
          kp=k1*len(self.labels)+k2
          Warren_Cowley_factor=self.c[k1]*(self.Kronecker(k1,k2)-self.c[k2])
          sro_real[kp][k]=(sro_real[kp][k]-self.c[k1]*self.c[k2])/Warren_Cowley_factor 
    self.sro_real=sro_real
    print(self.sro_real)
  def cal_spin(self):
    count=np.zeros(self.NR) #0
    #spin_i_avg=np.zeros(len(self.labels)**2)
    for i in range(0,len(self.atoms_pos)):
      self.atoms_pos[i][0:3]=self.move_in_supercell(self.atoms_pos[i][0:3])
      self.atoms_pos[i][3]=self.spin_of(self.atoms_pos[i])
      print(self.atoms_pos[i][0:3])
    sro_real=np.zeros([len(self.labels)**2,self.NR])
    for i in range(0,len(self.atoms_pos)-1):
      spin_i3=self.atoms_pos[i][3] #self.spin_of(self.atoms_pos[i])
      for j in range(i+1,len(self.atoms_pos)):
        dist_mu_nu=np.linalg.norm(np.asarray(self.atoms_pos[i][0:3])-np.asarray(self.atoms_pos[j][0:3]))
        dist_mu_nu=round(dist_mu_nu,4)
        if dist_mu_nu<self.neighbor_list[self.NR][0]:
          spin_j3=self.atoms_pos[j][3]
          for k in range(0,self.NR):
            #print dist_mu_nu, self.neighbor_list[k][0]
            if dist_mu_nu==self.neighbor_list[k][0]:
              Ri=k
              break
    #      print("Ri:",Ri)
          for k1 in range(0,len(spin_i3)): 
            for k2 in range(0,len(spin_j3)):
              spin_ik1=spin_i3[k1] 
              spin_jk2=spin_j3[k2] 
              spin_ij=float(spin_ik1)*float(spin_jk2)
             # if k1*len(spin_i3)+k2==2 and Ri==1: print spin_i3[k1], spin_j3[k2], spin_ij
              sro_real[k1*len(spin_i3)+k2][Ri] += (spin_ij-self.c[k1]*self.c[k2])/self.c[k1]/(self.Kronecker(k1,k2)-self.c[k2])       #spin_ij
          count[k]+=1
    self.sro_real=1.0*sro_real/count
    print("SRO_R:", self.sro_real)
  def sro_k(self,k,pair):
    index=pair[0]*len(self.labels)+pair[1]
    alpha_k=0 #complex(0,0)
    N_atom=0
    V=self.period[0]*self.period[1]*self.period[2]*1.0
    for nR in range(0,self.NR): #len(self.neighbor_list)):
      for nRi in range(0,len(self.neighbor_list[nR][1])):
        phase=np.pi*np.dot(k,self.neighbor_list[nR][1][nRi])
        #print "phase:", np.cos(phase)
        alpha_k += self.sro_real[index][nR]*np.cos(phase) #/V
    #alpha_k = alpha_k/(self.c[pair[0]]*(self.Kronecker(pair[0],pair[1])-self.c[pair[1]]))
    return alpha_k
  def sro(self,pair):
    alpha_k_list=[]
    for ix in range(0,self.N[0]):
      alpha_k_list.append([0])
      for iy in range(0,self.N[1]):
        k=np.array([4.0*ix/(self.N[0]-1),4.0*iy/(self.N[1]-1),0]) #/np.array([coord[0][0],coord[1][1],coord[2][2]])
        alpha_k_list[ix].append(self.sro_k(k,pair))
     # if iy==19 or iy==39: print sro_k(atoms_pos,k,c)
      alpha_k_list[ix].pop(0)
    return alpha_k_list
  def plot_sro(self,solute_pair,new_labels):
    alpha_k_list=self.sro(solute_pair)
    self.alpha_k_arr=np.asarray(alpha_k_list)
    x=np.linspace(0,2.0,self.N[0])
    X,Y=np.meshgrid(x,x)
    Z=np.asarray(self.alpha_k_arr)
   #level=np.linspace(np.min(Z),np.max(Z)/2.0,40)
   # print Z
    plt.contour(X,Y,Z,80,cmap=cm.coolwarm) #, vmin=0, vmax=60)
   # plt.annotate(new_labels, (8, 1), backgroundcolor='w')
    plt.annotate(new_labels, (.7, .1), xycoords='axes fraction', backgroundcolor='w')
    #plt.legend()
   #plt.scatter(X,Y,Z)
   # plt.plot(X[0,:],Z[0,:])
    plt.colorbar()
  def plot_sro2(self,alpha_k_list):
    self.alpha_k_arr=np.asarray(alpha_k_list)
    x=np.linspace(0,2.0,self.N[0])
    X,Y=np.meshgrid(x,x)
    Z=np.asarray(self.alpha_k_arr)
   #level=np.linspace(np.min(Z),np.max(Z)/2.0,40)
    plt.contour(X,Y,Z,40,cmap=cm.coolwarm) #, vmin=0, vmax=60)
   #plt.scatter(X,Y,Z)
   # plt.plot(X[0,:],Z[0,:])
    plt.colorbar()
  def plot_sro_line(self,Label,solute_pair,line_color):
    alpha_k_list=self.sro(solute_pair)
    self.alpha_k_arr=np.asarray(alpha_k_list)
    Z=np.asarray(self.alpha_k_arr) 
    lenZ2=len(Z[0,:])/2 #(len(Z[0,:]+1))/2
    x1=np.linspace(0,1,lenZ2)
    x2=np.linspace(1,2,lenZ2)
    x3=np.linspace(2,3,lenZ2)
    x4=np.linspace(2,3,lenZ2)
    y1=Z[0,0:lenZ2]
    y2=Z[0:lenZ2,lenZ2] 
    #y3=Z[0:(len(Z[0,:])/2+0),(len(Z[:,0])/2+0)]
    #y3=y3[::-1]
    y4=[]
    for i in range(0,lenZ2):
      y4.append(Z[i,i])
    y4=y4[::-1]
    x5=np.linspace(3,4,lenZ2)
    y5=[]
    for i in range(0,lenZ2):
      k=np.array([4.0*i/(self.N[0]-1),4.0*i/(self.N[1]-1),4.0*i/(self.N[1]-1)])
      y5.append(self.sro_k(k,solute_pair))
    #X,Y=np.meshgrid(x,x)
   #level=np.linspace(np.min(Z),np.max(Z)/2.0,40)
   # plt.contour(X,Y,Z,40,cmap=cm.coolwarm,label=Label) #, vmin=0, vmax=60)
   #plt.scatter(X,Y,Z)
   # plt.plot(X[0,:],Z[0,:])
    plt.plot(x1,y1,line_color, x2,y2,line_color,x4,y4,line_color,x5,y5,line_color)
   # plt.colorbar()

prex_r="random_"
prex_o="order_"
#write_sro_k(sys.argv[1],prex_o+"sro4.dat")
#write_sro_k(sys.argv[2],prex_r+"sro9.dat")

#sro_dat=read_sro("sro9.dat")
#---------------start calculations-------------------
def plt_sro_k(rfile,sro_source):
  NR=10 
  N=np.array([40,40])
  case=sro(rfile,N,NR) #,period)
  case.read_data()
  case.find_neighbors() #3 # number of period
  labels=case.labels
  case.cal_spin2()
  fig_id=0
  for i in range(0,len(labels)):
    for j in range(0,len(labels)):
      new_labels=labels[i]+"-"+labels[j]
      fig_id+=1
      if j>=i:
        plt.subplot(len(labels),len(labels),fig_id)
        solute_pair=np.array([i,j])
        case.plot_sro(solute_pair,new_labels)
  labels_all=str(labels[0])
  for i in range(1,len(labels)):
    labels_all += str(labels[i])
  plt.savefig(str(sys.argv[1])+'.png', dpi=300)
  plt.show()
  #--------plot-lines---------
  #line_color_o="-b"
  #line_color_r="-r"
  #plt.subplot(3410)
  #order.plot_sro_line("Ni-Co",solute_pair1,line_color_o)
  #random.plot_sro_line("Ni-Co",solute_pair1,line_color_r)
  #plt.subplot(3411)
  #order.plot_sro_line("Co-Cr",solute_pair2,line_color_o) 
  #random.plot_sro_line("Co-Cr",solute_pair2,line_color_r)
  #plt.subplot(3412)
  #order.plot_sro_line("Ni-Cr",solute_pair3,line_color_o) 
  #random.plot_sro_line("Ni-Cr",solute_pair3,line_color_r)

start=time.time()
plt_sro_k(sys.argv[1],sys.argv[2])
print("runtime in second:", time.time()-start)

