
from dos import *
from matplotlib import RcParams
latex_style_times = RcParams({#'font.family': 'serif',
               #'font.serif': ['Times'],
               'text.usetex': True,
               'axes.linewidth':2.0
               #'figsize':[6.4,4.8]
               })

plt.style.use(latex_style_times)
plt.rcParams.update({'font.size': 15})
plt.rcParams["figure.figsize"]=[14.0,10.8]
#plt.rcParams["figure.figsize"]=[12.0,4.8]
plt.rcParams["figure.dpi"]=300
def plot_Cv():
  caseNiCoCr=ppdos(1E-5,192,'cdos.dat',1)
  T=np.linspace(200,2000,1000)
  caseNiCoCr.plot_specific_heat_can(T)
  #plt.show()
  plt.savefig('dos.png',dpi=300)
#plot_Cv()
def get_sros():
  case=ppdos(1E-5,240,'345/cdos.dat',1)
  T=np.linspace(100,3000,1000)
  case.print_sro_can(T)
get_sros()
def plot_sro():
  import shutil,os
  sizes=[108,192,240,480]
  file_dirs=['333','344','345','456']
  CE="3body"
  Tc_shift=0 #300
  for size,file_dir in zip(sizes,file_dirs):
    plot_sro3(CE,size,os.path.join(file_dir,'cdos.dat')) #,Tc_shift)
    #shutil.move("sro.png",file_dir)
    os.replace("sro.png",os.path.join(file_dir,file_dir+"_sro.png"))
plot_sro()

#plot binder cumulant
file_dir='.'
plot_test_supercell(file_dir)
#plot_binder_cumulant(file_dir)
