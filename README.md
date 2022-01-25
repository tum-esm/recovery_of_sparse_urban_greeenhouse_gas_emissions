# Recovery of sparse urban greenhouse gas emissions

[![DOI](https://zenodo.org/badge/450792204.svg)](https://zenodo.org/badge/latestdoi/450792204)

This is the code for the paper "Recovery of sparse urban greenhouse gas emissions"
and reproduces the results in the paper given the right input files.
The input files are not part of this repository.
The footprints used in the paper are found [here](https://www.dropbox.com/sh/xbjz02jnnihd4ha/AAAtx6RKs-YFZDV0HbtSYrata?dl=0).

The inventory files are not provided, for details please see the paper.

## Requirements to run the code

In order to run the code, Matlab 2021a and the Matlab libraries CVX and Gurobi 9.11 are needed.

CVX is available for free [here](http://cvxr.com/cvx/download/).
For [Gurobi](https://www.gurobi.com/products/gurobi-optimizer/) an academic version exists.


## Description of the Code

The code is located in the folder "source". All of the Matlab files in this folder are scripts which produce outputs similar to the figures in the paper.
The input data for each script can be modified in the "Input data" section of each script, where e.g. the emission inventory and amount of measurement stations can be chosen.

A detailed description for each script is found in the first lines of each script itself.

The footprints should be placed into a folder called "data/footprint/" and the inventory data into "data/emission_map/", however, this path can be changed in the input section of each script.

The code in subfolders are functions needed for the scripts and do not fulfill other purpose.

## How to cite

DOI of the github repository: [![DOI](https://zenodo.org/badge/450792204.svg)](https://zenodo.org/badge/latestdoi/450792204)

Please also cite our paper if you use this code.
