# Code-Monte-Carlo-SRO

This repository includes scripts/codes that we used to calculate temperature-dependent short-range-order parameters for complex concentrated alloys.

This is mainly a Python code tested on Python 3.6 with Mac OS 11.0.1. It is expected to run in different operation systems as long as Python3 works.

### Instruction
#### Calculate the contribution of the lattice distortion to the Hamiltonian
```shell
cd Lattice-Green-Function-NiCoCr
python LGF.py
```

#### Metroplis MC simulations
Once the Hamiltonian is ready, put the interaction parameters in the file Metroplis_MC_simulations/eci.out
Then submit the simulation as a job to a cluster of supercomputer.
```shell
cd Metroplis_MC_simulations
sbatch sb_mc.sh
```
Once the job finishes, change the format of restart.struct to POSCAR. Put the file in another folder calculate_SRO_k-space. The order parameters in the reciprocal space can be calculated by
```shell
cd calculate_SRO_k-space
sbatch sb_sro_k.sh
```
or submit it as a job to run it remotely in a cluster.
```shell
python sro_k_by_R.py restart.struct 0
```
The last parameter "0" means starting the calculation from scratch and then plotting a figure based on the result. If "1" is given, the code will search for the calculated short-range-order parameter in reciprocal and plot the figure based on the data. 
The calculation can take a while.

#### Wang-Landau simulations
Download the Wang-Landau Monte-Carlo code (see link below). Copy the interaction parameters (eci.out), structure file (init.struct) and other related files to a new folder or examples/example-ternary-NiCoCr/.
Then submit it as a job to a cluster.
```shell
cd examples/example-ternary-NiCoCr/
sbatch qsub.sh (This is for slurm manager; need to adjust this to specific job manager)
```
To accelerate the calculation, we divide the whole calculation into multiple ones and each one correspond to a specific energy range. This can be done with the following command, but the energy range (minimum of ordered and maximum of random states) need to be provided.
```shell
python gen_energy_array.py
```
One the multiple jobs are done, all important information is written to restart.log.wlce, including the density of state. Then copy the output files restart.log.wlce to Wang-Landau-data-processing.

First, use the code to generate a combined cdos.dat file for the density of state;
```shell
cd Wang-Landau-data-processing
python combine-dos. py
```
Then use dos.py to test the convergence and calculate the thermodynamic properties, such as specific heat Cv.
```shell
python dos.py
```
We can also plot the density of states (dos) with
```shell
python plot_dos.py
```
### Other related codes/tools
The drivers of the Metroplis Monte Carlo and Wang-Landau Monte Carlo can be found in other repositories. 

#### Metropolis
For Metropolis, please refer to https://github.com/peizong/Simulated-Annealing-Metropolis-MC-.

#### Wang-Landau

For Wang-Landau, please refer to https://github.com/peizong/WL-CE.
