# Cluster-Analysis
Perform cluster analysis on carbon-hydrogen GCMC simulations and export files for PIMD processing. 

Perform a cluster anaysis or generate and export all clusters in a lammpstrj file for analysis using PIMD. The program will generate a directory which contains:
- An .xyz file for the atoms within a bounding rectangle of the cluster for analysis within OVITO. The rectangle is extended by a distance which can be changed using the box_extension command. The hydrogen atoms of the cluster are additionally exported as a seperate atom type to make clear the contents of the cluster.
- A .lmp file to use as input for the PIMD calculation for each cluster.
- PIMD_Template.in, a lammps input file template which is used for each input file.
- run_cluster.sh, a bash script to queue a PIMD calculation on pawsey.
- runfiles.sh, a bash script to process and queue all PIMD calculations.
- script_settings.txt, a document which contains the settings used to generate the given file.

Note the .xyz and .lmp files are of the format cluster-id_number-of-atoms.lmp/xyz

## Usage

Run cluster_analysis.py to execute the script using a default file name and paramaters. Accepted arguments are:
