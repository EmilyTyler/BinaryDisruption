# BinaryDisruption

This code simulates the disruption of wide binary stars in the Galactic halo by compact objects. I wrote this alongside my PhD thesis, which provides the necessary context. My thesis can be found [here](http://eprints.nottingham.ac.uk/69110/).

The simulation code is in the `src` folder and the code for analysing the results is in the `scripts` folder.

## Python Scripts

To install the required python packages run:

`pip install -r requirements.txt`

### Analysis Script

The script that runs the analysis on the simulation output is `scripts/yoo_analysis.py`. This will read the observed binary data in `scripts/data` and the simulation output in `output` and perform a statistical analysis on it. This analysis will output the maximum allowed fraction of compact objects that could make up dark matter as a function of compact object mass.

The simulation takes a long time to run (up to two weeks for each perturber mass/dark matter density combination), already generated output files can be found [here](https://1drv.ms/f/s!AnBSDojWOifs3kIc1ps7gtET0jYq?e=f13UnM), just copy these into a folder called `output`.

To change the observed binaries on which to run the analysis, change the input file read in on line 56 of `scripts/yoo_analysis.py`.

### Calculating the Time-Averaged Dark Matter Density

The time-averaged dark matter density for binaries from the [AMR catalog](https://arxiv.org/abs/1406.5164) can be calculated and plotted by running this script: `scripts/average_dark_matter_density.py`.

## Simulation Code

It compiles uses g++, make sure you have a g++ compiler installed and then run `make` in the root directory. The executable will be built in the `bin` folder. Run `make clean` to delete the object files and the `bin` folder.

All of the input variables are set in the `main` function of `main.cpp`. The name of the output file is configured on line 95 of `main.cpp`.
