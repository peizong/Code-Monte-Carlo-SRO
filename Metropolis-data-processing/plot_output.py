#!/nfs/apps/Compilers/Python/Anaconda/2.7/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

from matplotlib import RcParams
latex_style_times = RcParams({'font.family': 'serif',
             'font.serif': ['Times'],
             'text.usetex': True,
             'axes.linewidth':2.0
             #'figsize':[6.4,4.8]
             })

plt.style.use(latex_style_times)
plt.rcParams.update({'font.size': 15})
plt.rcParams['lines.linewidth']=3.0
plt.rcParams["figure.figsize"]=[7,8] #[3.2,2.4]
plt.rcParams["figure.dpi"]=300
plt.tick_params(axis="both",direction="in",left=True,right=True,bottom=True,top=True)

#labels=np.array(["Al","Cr","Fe","Mn","Ti"])
def plot_lro(labels):
  #labels=np.array(["Ti","Fe","Al"]) #np.array(["Fe","Al","Ti"])
  data=np.genfromtxt("lro.dat",skip_header=0)
  for i in range(1,len(data[0])):
    plt.plot(data[:,0],data[:,i],label=labels[i-1])
  plt.legend()
  plt.xlabel("Temperature/K")
  plt.ylabel("Long range order parameter")
  plt.show()
def plot_sro(labels):
  #labels=np.array(["Ti","Fe","Al"]) #np.array(["Fe","Al","Ti"])
  mlabels=[]
  data=np.genfromtxt("sro.dat",skip_header=0)
  for i in range(0,len(labels)):
    idx=str(len(labels))+str(1)+str(i+1)
    plt.subplot(idx)
    for j in range(0,len(labels)):
      #mlabels.append(labels[i]+"-"+labels[j])
      new_label=labels[i]+"-"+labels[j]
      plt.plot(data[:,0],data[:,1+i*len(labels)+j],label=new_label) #mlabels[i*len(labels)+j])
    plt.legend()
    plt.xlabel("Temperature/K")
    if i==1: plt.ylabel("Renormalized SRO")
    plt.ylim(0,0.6)
  plt.savefig("metropolis666.png",dpi=300)
  #plt.show()
def plot_energy():
  data=np.genfromtxt("energies.dat-corrected",skip_header=0)
  x=np.linspace(1,len(data[:,1]),len(data[:,1]))
  plt.plot(x,data[:,1])
  plt.xlabel("MC steps")
  plt.ylabel("energies/eV")
  plt.savefig("energies-100K.png",dpi=300)
  plt.show()
labels=np.array(["Co","Cr","Ni"]) #np.array(["Ni","Co","Cr"])
#labels=np.array(["Ti","Fe","Al"])
#labels=np.array(["Al","Cr","Fe","Mn","Ti"])

plot_sro(labels)
#plot_lro(labels)
#plot_energy()
