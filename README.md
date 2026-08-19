# MD_and_ML_chemistry_tools
Code for running and analyzing molecular dynamics simulations, including python analysis code, ML potential tools for MACE, bash scripts, and LAMMPS modifications.

## Python Analysis Scripts
The python scripts in this folder are primarily for post-simulation analysis of MD data. This includes code to calculate IR spectra from dipole or dipole derivative timeseries and analyze the locations and orientations of water molecules in clusters. Many of these files have been written to be self-contained, meaning they can be copied individually onto an HPC cluster and run without setting up a package or calling other python files. A list of package dependencies that covers all of the python code here can be found in requirements.txt.

## Bash Scripts
The bash scripts in this folder were written as tools for file management and processing and job management on HPC clusters. 

## LAMMPS Modifications
The files in this folder are modifications of the LAMMPS source code that I made to perform the functions needed for simulating charged hydrocarbon-water clusters. This includes a pair potential for soft spheres as described in https://doi.org/10.1021/acs.jpcb.3c07777.  
