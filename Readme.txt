Materials Modelling (19/20,20/21) MATE96009 / MATE97008  @ Imperial College London
Module: Research lecture - Kinetic Monte Carlo simulation - Charge transport in organics
Date: 15th March
GTAs: Mathias Golomb, Chengcheng Xiao, Ju Huang, Xinwei Wang

The full hand out is Exercises_python.pdf

## Usage
Get a Git repository by typing
git clone https://github.com/dreamslink51/MSE317.git

To go to KMC folder, you type
cd MSE317/KMC
The executable file is called snapshot_KMC.x and the input file is input.inp. 

To run snapshot_KMC.x you type
./snapshot_KMC.x input.inp

## Output analysis
Output files will be generated in the same folder including distance.dat and tcheck.dat. Those two files will be read by AnalyseDistance.py, a Python script, to compute mobilities.
./AnalyseDistance.py -nsample <your #KMC> -temp <your temerature>
nsample: number of KMC runs (Default: 10)
temp: temperature (K) (Default: 283)
For example  ./AnalyseDistance.py -nsample 10 -temp 283

### Extended Reading List
Yang et al, Intermolecular charge transfer parameters, electron–phonon couplings, and the validity of polaron hopping models in organic semiconducting crystals: rubrene, pentacene, and C60(2017)
(https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.7b00618)
