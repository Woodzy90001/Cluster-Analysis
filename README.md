# Cluster-Analysis
Perform cluster analysis on carbon-hydrogen GCMC simulations and export files for PIMD processing. 

This program will read a lammpstrj file and parse it allowing for a cluster anaysis to be performed or to generate and export all clusters for analysis using PIMD. The program will generate a directory which contains:
- An .xyz file for the atoms within a bounding rectangle of the cluster for analysis within OVITO. The rectangle is extended by a distance which can be changed using the box_extension command. The hydrogen atoms of the cluster are additionally exported as a seperate atom type to make clear the contents of the cluster.
- A .lmp file to use as input for the PIMD calculation for each cluster.
- PIMD_Template.in, a lammps input file template which is used for each input file.
- run_cluster.sh, a bash script to queue a PIMD calculation on pawsey.
- runfiles.sh, a bash script to process and queue all PIMD calculations.
- script_settings.txt, a document which contains the settings used to generate the given file.

Note the .xyz and .lmp files are of the format cluster-id_number-of-atoms.lmp/xyz. The program requires the github files cluster_analysis.py, run_cluster.sh abd PIMD-Template.in and the lammpstrj input file or frame file to be located within the same directory to run.

## Usage

Run cluster_analysis.py to execute the script using a default file name and paramaters. For arguments containing an = sign, ensure there is not whitespace around it. 

Accepted arguments are:

    INPUT-FILE-STRING.lammpstrj
If a filename is given followed by the suffix .lammpstrj, then the program will use that file if present in the folder.


    SAVED-FRAMES-STRING.fd
Load frame data: If a filename is given followed bt the suffix .fd, then the program will assume it is a saved frame and use it as the input frame if it is in the folder.


    --sf=FILE-NAME-STR
Save frames: If a filename is given the program will save the selected frame. 


    --cf-tf=TIMESTEP-INT
  Choose frame - Timestep number: Choose the frame of the lammpstrj to use by timestep number.

    --cf-fn=TIMESTEP-INT
  Choose frame - Frame number: Choose the frame of the lammpstrj to use by the frame count. eg. 4 would correspond to the 4th saved frame, which could be timestep 40000.

    --cs=CLUSTER-SIZE-FLOAT
  Cluster size: Set the size of the cluster. Default 4.5

    --cae=EXTENSION-SIZE-FLOAT
  Cluster area extension: Set how much to extend the bounding box of the cluster when saving atoms. Default 5

    --minc=MIN-CLUSTER-INT
  Minimum cluster: Sets the minimum cluster to use. Is used as a list indices with maxc. eg. --minc=3 --maxc=10 -> clusters[3:10]

    --maxc=MAX-CLUSTER-INT
  Maximum cluster: Sets the maximum cluster to use. Is used as a list indices with minc. eg. --minc=3 --maxc=10 -> clusters[3:10]
  
    --ca
Cluster Analysis: Perform a cluster analysis. This will calculate abd plot the cluster count, average coordination and average cluster size for a sweep of cluster radii.

    --mncs=CLUSTER-SIZE-FLOAT
  Minimum cluster size: Set the minimum cluster size for the cluster analysis scan. Default 2

    --mxcs=CLUSTER-SIZE-FLOAT
  Maximum cluster size: Set the maximum cluster size for the cluster analysis scan. Default 6

    --csnp=NUM-POINTS-INT
  Cluster scan number of points: Set how many points to scan for the cluster analysis between mncs and mxcs. Default 50

    --cssf
  Cluster scan save figure: Save the figures from the cluster scan.

    --st
  Show timings: Show timings for the lammpstrj parser.
