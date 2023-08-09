### Code-Monte-Carlo-SRO

This repository includes scripts/codes that we used to calculate temperature-dependent short-range-order parameters for complex concentrated alloys.

This is mainly a Python code tested on Python 3.6 with Mac OS 11.0.1. It is expected to run in different operation systems as long as Python3 works.

#### Calculate the contribution of the lattice distortion to the Hamiltonian
'''shell
cd Lattice-Green-Function-NiCoCr
python LGF.py
'''

The drivers of the Metroplis Monte Carlo and Wang-Landau Monte Carlo can be found in other repositories. 

#### Metropolis
For Metropolis, please refer to https://github.com/peizong/Simulated-Annealing-Metropolis-MC-.

#### Wang-Landau

For Wang-Landau, please refer to https://github.com/peizong/WL-CE.
