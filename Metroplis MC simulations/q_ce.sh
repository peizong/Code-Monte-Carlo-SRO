#! /usr/bin/env bash
#PBS -j oe
#PBS -N FeCo
#PBS -l walltime=24:00:00
#PBS -l nodes=6
#PBS -A MAT020

#-------------------------------------------

#--Change directory to the working directory.

cd $PBS_O_WORKDIR
# cd $MEMBERWORK 

#export OMP_NUM_THREADS=8
#module swap PrgEnv-intel PrgEnv-gnu
module load vasp5 #/5.4
oldword=XXX
#--Run the executable.
mmaps -m=3 &
date
time pollmach runstruct_vasp_magmomNiCoCr aprun -n 96 -N 16 vasp5
#time pollmach runstruct_vasp_magmom5 aprun -n 96 -N 16 -d 1 -j 2 vasp5
#time emc2 -gs=1 -mu0=1.5 -mu1=0.5 -dmu=0.04 -T0=300 -T1=5000 -dT=50 -k=8.617e-5 -dx=1e-3 -er=50
#time phb -T=900 -mu=5 -gs1=-1 -gs2=5 -dT=25 -dn -dx=1e-3 -er=50 -k=8.617e-5 -ltep=5e-3 -o=phm15.out 
#time phb -T=800 -mu=5 -gs1=-1 -gs2=5 -dT=25 -dx=1e-3 -er=50 -k=8.617e-5 -ltep=5e-3 -o=phm5d.out
#time ./cal_E0.sh
#for newword in 16
#do
#sed "s/$oldword/$newword/g" input/i_lsms_string > i_lsms
#time aprun -n1 -N1 -d8 -j1 /ccs/home/zp7/LSMS_3-eos-libxc/bin/lsms i_lsms
#cp k.out out/k.out_$newword
#done
date

