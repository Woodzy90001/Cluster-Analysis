'''

cluster_analysis.py

Parse a lammpstrj file into files for a PIMD calulcation. Also has support for
cluster analysis. 

'''
import os
import sys
import time
import pickle
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.spatial  import KDTree



PIMDTEMPLATENAME = "PIMD-Template.in"
RUNCLUSTERNAME = "run_cluster.sh"
DRIVERFILENAME = "runfiles.sh"









'''

Parse a lammpstrj into a dictionary containing the bounds, hydrogen and carbon
atoms

Parameters
----------
file_name : string
    Name of file to parse.
show_times : bool, optional
    Print timing info. The default is True.

Raises
------
Exception
    Raised if the atom count does match which means something has gone wrong.

Returns
-------
frame_holder : Dict{int : [hydrogen atom list, carbon atom list, bounds]}
        - hydrogen atom list : List<float, float, float, int> [x, y, z, ID]
        - carbon atom list   : List<float, float, float, int> [x, y, z, ID]
        - bounds             : List<float, float, float>      [x, y, z]
    Dictionary containing all the atoms for each frame, seperated into carbon
    and hydrogen and the bounds for each frame.

'''

def lammpstrj_parser(file_name, show_times=True):
    
    # Read file
    start_time = time.time()
    f = open(file_name, "r")
    data = f.readlines()
    f.close()
    
    if (show_times): print("{: <35} Time taken : {:.2f}s".format("File read.", time.time() - start_time))
    
    # Setup
    frame_holder = {}
    atoms = [[], []]
    frame_times = []
    timestep_next = False
    atom_num_next = False
    atoms_started = False
    bounds_counter = 0
    
    timestep = None
    atom_num = None
    bounds = None
    start_frame_time = time.time()
    for line in data:
        line_stripped = line.strip()
        
        # Process new timestep
        if (line_stripped == "ITEM: TIMESTEP"):
            frame_times.append(time.time() - start_frame_time)
            start_frame_time = time.time()
            timestep_next = True
            atoms_started = False
            
            # Don't do it for the first iteration
            if timestep is not None:
                # Save previous frame
                atoms.append(bounds)
                frame_holder[timestep] = atoms
            
                # Check number of atoms is correct
                if (atom_num != (len(atoms[0]) + len(atoms[1]))): raise Exception("Atom count doesnt match. \nExpected: {}\n\nHydrogen: {}\nCarbon: {}\n\nTotal: {}".format(atom_num, len(atoms[0]), len(atoms[1]), len(atoms[0]) + len(atoms[1]))) 
                
            # Make new atom array
            atoms = [[], []]
        
        # Process atom
        elif (atoms_started):
            atom_data = line_stripped.split(" ")
            ID = int(atom_data[0])
            atom_type = int(atom_data[1])
            x = float(atom_data[2])
            y = float(atom_data[3])
            z = float(atom_data[4])
            atoms[atom_type - 1].append([x, y, z, ID])
        
        # Process timestep, atom number or bounds
        elif (timestep_next) or (atom_num_next) or (bounds_counter > 0):
            if (timestep_next):
                timestep = int(line_stripped)
                timestep_next = False
            elif (atom_num_next):
                atom_num = int(line_stripped)
                atom_num_next = False
            elif (bounds_counter > 0):
                bounds[bounds_counter - 1] = float(line_stripped.split(" ")[1])
                bounds_counter += 1
                bounds_counter %= 4
            
        # Process headers
        else:
            if (line_stripped == "ITEM: NUMBER OF ATOMS"):
                atom_num_next = True
            elif (line_stripped == "ITEM: BOX BOUNDS pp pp pp"):
                bounds_counter = 1
                bounds = np.array([-1.1, -1.1, -1.1]) # Set dummy impossible value
            elif (line_stripped == "ITEM: ATOMS id type x y z"):
                atoms_started = True
                
    # Add the last iteration
    atoms.append(bounds)
    frame_holder[timestep] = atoms
    
    del frame_times[0]
    if (show_times): print("{: <35}           : {}".format("Total number of frames", len(frame_holder)))
    if (show_times): print("{: <35}           : {:.2f}s".format("Average frame time", sum(frame_times) / len(frame_times)))
    if (show_times): print("{: <35}           : {:.2f}s".format("Total time taken", time.time() - start_time))
    
    return frame_holder










'''

Generate and plot cluster statistics.

Parameters
----------
arg_information : ArgInformation
    Contains the frame data, min and max cluster size, if the plots are to be
    saved and cluster spacings.

Returns
-------
None.

'''

def AnalyseClusterSize(arg_information):
    
    # Extract params
    min_cluster_size = arg_information.min_cluster_size
    max_cluster_size = arg_information.max_cluster_size
    num_points = arg_information.num_points
    save_figures = arg_information.save_figures
    
    # Setup
    scan_sizes = np.linspace(min_cluster_size, max_cluster_size, num_points)
    cluster_count = []
    coodination = []
    avg_cluster_size = []
    count = 0
    cutoff = 0
    
    # get data for each cluster size
    print("Cluster scan - {: >3}% done.".format(0))
    for scan_radius in scan_sizes:
        arg_information.cluster_size = scan_radius
        hydrogen_atom_ids, hydrogen_atom_coords, coordination_numbers, clusters, clusters_ids = FindClusters(arg_information)
        cluster_count.append(len(clusters))
        avg_cluster_size.append(sum([len(i) for i in clusters]) / len(clusters))
        coodination.append(sum(coordination_numbers) / len(coordination_numbers))
        count += 1
                
        if (round(100 * count / num_points, -1) > cutoff): 
            cutoff = round(100 * count / num_points, -1)
            print("Cluster scan - {: >3}% done.".format(cutoff))

        #VerifyData(clusters, clusters_ids, data[0], atom_ids)

    figure_filename_prefix = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Plot data
    plt.figure()
    plt.plot(scan_sizes, cluster_count)
    plt.xlabel("Cutoff distance (A)")
    plt.ylabel("Cluster count (Number of clusters)")
    if save_figures: plt.savefig("{}_Cluster_Count.png".format(figure_filename_prefix), bbox_inches='tight')
    plt.show()


    plt.figure()
    plt.plot(scan_sizes, coodination)
    plt.xlabel("Cutoff distance (A)")
    plt.ylabel("Average coordination of atoms (atoms)")
    if save_figures: plt.savefig("{}_Coordination.png".format(figure_filename_prefix), bbox_inches='tight')
    plt.show()


    plt.figure()
    plt.plot(scan_sizes, avg_cluster_size)
    plt.xlabel("Cutoff distance (A)")
    plt.ylabel("Average cluster size (atoms)")
    if save_figures: plt.savefig("{}_Average_Cluster_Size.png".format(figure_filename_prefix), bbox_inches='tight')
    plt.show()










'''

Parse the command line arguments

Returns
-------
arg_information : ArgInformation
    A class containing all relevant variables for the cluster analysis.

'''

def ParseArgs():
    
    # Setup
    command_line_args = sys.argv
    
    print("Parsing input. ")
    
    file_name = "H2-1S_adsorb_-1.0.lammpstrj"
    frame_select = [None, -1]
    cluster_size = 4.5
    min_cluster_size = 2
    max_cluster_size = 6
    num_points = 50
    box_extension = 5
    min_cluster = None
    max_cluster = None
    cluster_analysis = False
    pickle_load_file = False
    show_timings = False
    save_data = False
    save_figures = False
    
    # iterate through all command line args
    if (len(command_line_args) > 1):
        for i in range(1, len(command_line_args)):
            current_arg = command_line_args[i]
            
            # Do cluster analysis
            if (current_arg == "--ca"): # Cluster Analysis
                cluster_analysis = True
                
            # Update file name
            elif (current_arg[-10:] == ".lammpstrj"):
                file_name = current_arg
                
            # Save the frame data 
            elif (current_arg[:4] == "--sf"): # Save Frames
                save_data = current_arg.split("=")
                if (len(save_data) > 1):
                    save_data = save_data[1]
                    if (save_data[-3:] != ".fd"): save_data += ".fd"
                else:
                    print("Error parsing frame filename. File data will not be saved. ")
                    save_data = False
            
            # Cut the clusters but stating a minimum cluster
            elif (current_arg[:6] == "--minc"): # Minimum cluster
                min_cluster = current_arg.split("=")
                try:
                    min_cluster = int(min_cluster[1])
                except:
                    print("Error parsing minimum cluster. 0 will be used. ")
                    min_cluster = None
                
            # Cut the clusters by stating a maximum cluster
            elif (current_arg[:6] == "--maxc"): # Maximum cluster
                max_cluster = current_arg.split("=")
                try:
                    max_cluster = int(max_cluster[1])
                except:
                    print("Error parsing maximum cluster. The last element will be used. ")
                    max_cluster = None
                if (type(min_cluster) == int) and (type(max_cluster) == int) and (max_cluster <= min_cluster): 
                    print("Error, maximum cluster must be strictly larger than the minimum cluster. Resetting to defaults. ")
                    min_cluster = None; max_cluster = None
                
            # Load frame data
            elif (current_arg[-3:] == ".fd"): 
                pickle_load_file = current_arg
                
            # Show timings
            elif (current_arg == "--st"): # Show Timings
                show_timings = True
            
            # Select a custom frame by timestep
            elif (current_arg[:7] == "--cf-ts"): # Choose Frame - Timestep Number
                save_data = current_arg.split("=")
                try:
                    frame_select[0] = int(save_data[1])
                except:
                    print("Error parsing frame number. The last frame will be used. ")
                    
            # Select a custom frame by frame number
            elif (current_arg[:7] == "--cf-fn"): # Choose Frame - Frame Number
                save_data = current_arg.split("=")
                try:
                    frame_select[1] = int(save_data[1])
                except:
                    print("Error parsing frame number. The last frame will be used. ")
            
            # Set cluster size
            elif (current_arg[:4] == "--cs"): # Cluster Size
                split_data = current_arg.split("=")
                try:
                    cluster_size = float(split_data[1])
                except:
                    print("Error parsing cluster size. 3.2A will be used. ")
                    
            # Set min cluster size for scan
            elif (current_arg[:6] == "--mncs"): # Minimum Cluster Size
                split_data = current_arg.split("=")
                try:
                    min_cluster_size = float(split_data[1])
                except:
                    print("Error parsing min cluster size. 2A will be used. ")
                    
            # Set min cluster size for scan
            elif (current_arg[:6] == "--mxcs"): # Maximum Cluster Size
                split_data = current_arg.split("=")
                try:
                    max_cluster_size = float(split_data[1])
                except:
                    print("Error parsing max cluster size. 6A will be used. ")
                    
            # Set min cluster size for scan
            elif (current_arg[:6] == "--csnp"): # Cluster Scan Number of Points
                split_data = current_arg.split("=")
                try:
                    num_points = int(split_data[1])
                except:
                    print("Error parsing cluster scan points. 50 will be used. ")
                    
            # Set min cluster size for scan
            elif (current_arg[:6] == "--cssf"): # Cluster Scan Save Figure
                save_figures = True
                
            # Set min cluster size for scan
            elif (current_arg[:6] == "--cae"): # Cluster Area Extension
                split_data = current_arg.split("=")
                try:
                    box_extension = float(split_data[1])
                except:
                    print("Error parsing cluster scan points. 50 will be used. ")
                
    
    print("Loading trj information. ")
    # Load pickled data if given, otherwise parse file
    if (pickle_load_file):
        frame_data = pickle.load(open(pickle_load_file, 'rb'))
    else:
        frame_data = lammpstrj_parser(file_name, show_times=show_timings)
        
        # Save the frame information if selected
        if (save_data): pickle.dump(frame_data, open(save_data, 'wb'))
    
    
    print("Getting frame. ")
    # Get all the timesteps in the dictionary
    timesteps = list(frame_data.keys())
    
    # Extract the correct frame
    if (frame_select[0] is None):
        try:
            frame_data = frame_data[timesteps[frame_select[1]]]
            frame_used = timesteps[frame_select[1]]
        except:
            print("Error getting frame, lest frame will be used. ")
            frame_data = frame_data[timesteps[-1]]
            frame_used = timesteps[-1]
    else:
        try:
            frame_data = frame_data[frame_select[0]]
            frame_used = frame_select[0]
        except:
            print("Error getting frame, lest frame will be used. ")
            frame_data = frame_data[timesteps[-1]]
            frame_used = timesteps[-1]
    
    arg_information = ArgInformation(cluster_analysis, show_timings, 
            cluster_size, frame_data, min_cluster_size, max_cluster_size, 
            num_points, save_figures, box_extension, frame_used, file_name,
            min_cluster, max_cluster)
    print("Input parsed. ")
    
    return arg_information










'''

Get the coorination of each hydrogen atom

Parameters
----------
hydrogen_atom_coords : np.array([np.array([float-x, floaty, float-z]), ..., np.array([float-x, floaty, float-z])])
    A numpy array containing all the hydrogen atoms.
radius : float
    Radius of cluster.
bounds : List<float, float, float>
    Bounds of the simulation box.

Returns
-------
coordination_numbers : np.array([int, ..., int])
    A numpy array containing the coordination number of all points.

'''

def GetCoordination(hydrogen_atom_coords, radius, bounds):
    
    # Make tree and get all indexes within the radius
    tree = KDTree(hydrogen_atom_coords, boxsize=bounds)
    indexes = tree.query_ball_point(hydrogen_atom_coords, radius)
    
    # Save the connected IDs, connected atoms and coordination number of each 
    # hydrogen atom
    coordination_numbers = []
    for i in range(len(indexes)):
        coordination_numbers.append(len(indexes[i]))
    
    coordination_numbers = np.asarray(coordination_numbers)
    
    return coordination_numbers










'''

Recursivley search all points within a radius around a point.

Parameters
----------
atom_copy : np.array([np.array([float-x, floaty, float-z]), ..., np.array([float-x, floaty, float-z])])
    Numpy list containing all points to process.
tree : kd-tree
    kd-tree containing all points.
radius : float
    Radius of cluster.
curr_cluster : List<np.array([float-x, floaty, float-z])>
    List of all elements found in the current cluster. Initially just the 
    original atom.
cluster_ids : List<int>
    List of all indexes of all elements found in the current cluster. 
    Initially just the original atoms index, 0.
atom_index : int
    Curr atom to processes index. Initially 0.

Returns
-------
cluster_ids : List<int>
    List of all the indexes of points found to be within the current cluster.
curr_cluster : List<np.array([float-x, floaty, float-z])>
    List of all elements found in the current cluster.

'''

def FindClusterRecursive(atom_copy, tree, radius, curr_cluster, cluster_ids, atom_index):
    
    # Get the indexes of all points within the radius of the current point
    indexes = tree.query_ball_point(atom_copy[atom_index], radius)
    # Remove the current point to have an array of any other points within
    # the radius of the current point
    indexes.remove(atom_index)
    # Remove any points already found to be in the current cluster
    indexes_cleaned = [x for x in indexes if x not in cluster_ids]
    # Add these new points to the list of all the clusters indexes
    curr_cluster.extend(atom_copy[indexes_cleaned])
    cluster_ids.extend(indexes_cleaned)
    
    # For each newly discovered point within the cluster, perform a search on
    # it
    for index in indexes_cleaned:
        cluster_ids, curr_cluster = FindClusterRecursive(atom_copy, tree, radius, curr_cluster, cluster_ids, index)
    
    return cluster_ids, curr_cluster










'''

Find all clusters for a given cluster radius

Parameters
----------
arg_information : ArgInformation
    A class containing all relevant variables for the cluster analysis.

Returns
-------
hydrogen_atom_ids : np.array([int, ..., int])
    A numpy array containing a map for all the hydrogen atom indexes. 
    Eg [4, 6, 1, 3] means the 1st hydrogen atom has an ID of 4 in the input 
    file.
hydrogen_atom_coords : np.array([np.array([float-x, floaty, float-z]), ..., np.array([float-x, floaty, float-z])])
    A numpy array containing all the hydrogen atoms.
hydrogen_coordination_number : np.array([int, ..., int])
    Coordination number of the hydrogen atoms.
clusters : List<List<np.array([float-x, floaty, float-z])>, ..., List<np.array([float-x, floaty, float-z])>>
    List of list of all points in a cluster.
clusters_ids : List<List<int>, ..., List<int>>
    List of list of indexes of points within the cluster. The indexes are from
    the reduced atom list whenever the cluster was parsed so all cluster IDs 
    begin from 0. This is only useful for counting cluster size and debugging.

'''

def FindClusters(arg_information):
    # Atom
    #   - in
    #       [x, y, z, ID]
    #   - out
    #       [ID], [x, y, z], [Coordinations], [Clusters - x, y, z]
    
    # Unpacking args
    frame = arg_information.frame_data
    radius = arg_information.cluster_size
    
    # Take hydrogen atoms and convert to numpy array
    hydrogen_atoms = np.asarray(frame[0])
    bounds = frame[2]                     # Bounds
    
    # Take only the xyz coordinates 
    hydrogen_atom_coords = hydrogen_atoms[:,:3]
    atom_copy = hydrogen_atom_coords.copy()

    # Get ids
    hydrogen_atom_ids = hydrogen_atoms[:,3:]

    # While there are still atoms which don't belong to a cluster, iterate.
    clusters = []
    clusters_ids = []
    while (len(atom_copy) > 0):
        
        # Take first atom of remaining atoms to process
        curr_cluster = [atom_copy[0]]
        # Build new tree of remaining atoms
        tree = KDTree(atom_copy, boxsize=bounds)
        # Recursivley process cluster
        cluster_ids, curr_cluster = FindClusterRecursive(atom_copy, tree, radius, curr_cluster, [0], 0)
        
        # Remove all atoms from copy list in this cluster
        atom_copy = np.delete(atom_copy, cluster_ids, axis=0)
        # For debugging puropses
        clusters_ids.append(cluster_ids)
        # Add cluster positions to cluster list
        clusters.append(curr_cluster)

    # get coordination stats about all points
    hydrogen_coordination_number = GetCoordination(hydrogen_atom_coords, radius, bounds)
    
    return hydrogen_atom_ids, hydrogen_atom_coords, hydrogen_coordination_number, clusters, clusters_ids










'''

Make a copy of the input file to make getting periodic atoms easier.
Inefficient on memory and time, but does take super long so it's acceptable

Parameters
----------
atoms : np.array([np.array([float, float, float]), ..., np.array([float, float, float])])
    Input array of all points to duplicate.
bounds : np.array([float, float, float]) - x-extent, y-extent, z-extent
    The extents of the box.

Returns
-------
all_atoms : np.array([np.array([float, float, float]), ..., np.array([float, float, float])])
    An array of all atoms and images of the atoms.

'''

def DuplicatePoints(atoms, bounds):
    
    all_atoms = []
    
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:
                
                offset_point = np.array([x, y, z]) * bounds
                for atom in atoms:
                    curr_position = atom + offset_point
                    all_atoms.append(curr_position)
                    
    all_atoms = np.asarray(all_atoms)
    
    return all_atoms










'''

Take a point and shift the xyz coordinates into images of the unit cell until
a position is found which minimises the distance to the previous point. 
This acts to ensure all points are close together and not split over the 
boundaries

Parameters
----------
last_point : np.array([float, float, float])
    The xyz coordinates of the last point.
curr_point : np.array([float, float, float])
    The xyz coordinates of the current point.
bounds : np.array([float, float, float])
    Extents of the box to shift the coordinate by.

Returns
-------
final_point : np.array([float, float, float])
    The shifted xyz coordinates of the current point.

'''

def DistanceMinimise(last_point, curr_point, bounds):
    
    # Initial setup
    best_distance = np.linalg.norm(last_point-curr_point)
    final_point = curr_point
    
    # For each periodic image
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:
                
                # Get offset point and distance
                offset_point = np.array([x, y, z]) * bounds
                curr_position = curr_point + offset_point
                curr_distance = np.linalg.norm(last_point-curr_position)
                
                # If distance is shorter than any previous coordinate, then 
                # save it 
                if (curr_distance < best_distance):
                    final_point = curr_position
                    best_distance = curr_distance

    return final_point










'''

Cull all the points outside the box of interest which is to be saved. 

Parameters
----------
atoms : np.array([np.array([float, float, float]), ..., np.array([float, float, float])])
    Input array of all points to refine.
analysis_mins : np.array([float, float, float])
    Minimum xyz coordinates.
analysis_maxs : np.array([float, float, float])
    Maximum xyz coordinates..

Returns
-------
final_atoms : List<np.array([float, float, float])>
    A list of atoms within the desired box.

'''

def GetBoxedAtoms(atoms, analysis_mins, analysis_maxs):
    
    return atoms[np.where(np.all(np.logical_and(atoms >= analysis_mins, atoms <= analysis_maxs), axis=1))]










'''

Save a copy of the driver files.

Returns
-------
PIMD_in : List<string>
    Contents of the PIMD input file.
run_cluster_sh : List<string>
    Contents of the LAMMPS sbatch script.

'''
    
def ReadTemplateFiles():

    f = open(PIMDTEMPLATENAME, "r")
    PIMD_in = f.readlines()
    f.close()
    
    f = open(RUNCLUSTERNAME, "r")
    run_cluster_sh = f.readlines()
    f.close()
    
    return PIMD_in, run_cluster_sh










'''

Write a copy of the driver files

Parameters
----------
PIMD_in : List<string>
    Contents of the PIMD input file.
run_cluster_sh : List<string>
    Contents of the LAMMPS sbatch script.
arg_information : ArgInformation
    A class containing all relevant variables for the cluster analysis.

Returns
-------
None.

'''
    
def WriteTemplateFiles(PIMD_in, run_cluster_sh, arg_information):

    
    f = open(PIMDTEMPLATENAME, "w")
    for line in PIMD_in: f.write(line)
    f.close()
    
    f = open(RUNCLUSTERNAME, "w")
    for line in run_cluster_sh: f.write(line)
    f.close()
    
    arg_information.WriteContents()










'''

Write the xyz file for the current cluster

Parameters
----------
curr_file_name : string
    The file name.
cluster_num : int
    Number of points in the current cluster.
total_atom_count : int
    Total hydrogen and carbon atoms to write in the current box.
bounds : np.array([float, float, float])
    The bounds of the current box.
mins : np.array([float, float, float])
    Zero point for the current box.
cluster : List<np.array([float-x, floaty, float-z])>
    List of all points in a cluster..
box_hydrogen_atoms : List<np.array([float, float, float])>
    List of all hydrogen atoms in the box.
box_carbon_atoms : List<np.array([float, float, float])>
    List of all carbon atoms in the box.

Returns
-------
None.

'''
    
def WriteXYZ(curr_file_name, cluster_num, total_atom_count, bounds, mins, cluster, box_hydrogen_atoms, box_carbon_atoms):
    
    # Make file and write header
    f = open("{}.xyz".format(curr_file_name), "w")
    f.write("{}\n{} {} {}\n".format(cluster_num + total_atom_count, bounds[0], bounds[1], bounds[2]))
    
    # Cluster points
    for atom in cluster:
        atom_shifted = atom - mins
        f.write("{} {} {} {}\n".format(0, atom_shifted[0], atom_shifted[1], atom_shifted[2]))
        
    # Hydrogen atoms
    for atom in box_hydrogen_atoms:
        atom_shifted = atom - mins
        f.write("{} {} {} {}\n".format(1, atom_shifted[0], atom_shifted[1], atom_shifted[2]))
        
    # Carbon atoms
    for atom in box_carbon_atoms:
        atom_shifted = atom - mins
        f.write("{} {} {} {}\n".format(2, atom_shifted[0], atom_shifted[1], atom_shifted[2]))
        
    f.close()
    
    
    
    
    
    
    
    
    
    
'''

Write the lmp file for the current cluster

Parameters
----------
curr_file_name : string
    The file name.
total_atom_count : int
    Total hydrogen and carbon atoms to write in the current box.
bounds : np.array([float, float, float])
    The bounds of the current box.
mins : np.array([float, float, float])
    Zero point for the current box.
box_hydrogen_atoms : List<np.array([float, float, float])>
    List of all hydrogen atoms in the box.
box_carbon_atoms : List<np.array([float, float, float])>
    List of all carbon atoms in the box.

Returns
-------
None.

'''
    
def WriteLMP(curr_file_name, total_atom_count, bounds, mins, box_hydrogen_atoms, box_carbon_atoms):
    
    # Make file and write header
    atom_number = 1
    f = open("{}.lmp".format(curr_file_name), "w")
    f.write("# LAMMPS data file written by cluster_analysis.py\n\n{} atoms\n\n0.0 {} xlo xhi\n0.0 {} ylo yhi\n0.0 {} zlo zhi\n\nAtoms  # full\n\n".format(total_atom_count, bounds[0], bounds[1], bounds[2]))

    # Hydrogen atoms
    for atom in box_hydrogen_atoms:
        atom_shifted = atom - mins
        f.write("{} 1 1 0.0 {} {} {}\n".format(atom_number, atom_shifted[0], atom_shifted[1], atom_shifted[2]))
        atom_number += 1
    
    # Carbon atoms
    for atom in box_carbon_atoms:
        atom_shifted = atom - mins
        f.write("{} 1 2 0.0 {} {} {}\n".format(atom_number, atom_shifted[0], atom_shifted[1], atom_shifted[2]))
        atom_number += 1
        
    f.close()












def WriteDriver(file_names):
    
    f = open(DRIVERFILENAME, "w")
    f.write('rm PIMD_*.in\n\n')
    input_file_formatted_string = 'input_files=('
    for file in file_names: input_file_formatted_string += '"{}" '.format(file)
    input_file_formatted_string = input_file_formatted_string.strip() + ')\n'
    f.write(input_file_formatted_string)
    f.write('for inp_file in "${input_files[@]}"\n')
    f.write('do\n')
    f.write("\tcluster_id=$(echo ${inp_file} | cut -d '_' -f 1)\n")
    f.write('\tcp PIMD-Template.in PIMD_${cluster_id}.in\n')
    f.write('\tsed -i -e  "s/INPUTFILE/${inp_file}/g" PIMD_${cluster_id}.in\n')
    f.write('\tsbatch {} '.format(RUNCLUSTERNAME))
    f.write('${inp_file}\n')
    f.write('done')
    f.close()







'''

Parse the arginformation and save the coodinates to be analysed. 

Parameters
----------
arg_information : ArgInformation
    A class containing all relevant variables for the cluster analysis.

Returns
-------
None.

'''
    
def WriteFiles(arg_information):

    # Save all clusters
    if not arg_information.cluster_analysis:
        
        # Get cluster information
        start_time = time.time()
        hydrogen_atom_ids, hydrogen_atom_coords, hydrogen_coordination_number, clusters, clusters_ids = FindClusters(arg_information)
        post_cluster = time.time()
        print("{: <35} Time taken : {:.2f}s".format("Clusters found.", post_cluster - start_time))
        
        # Sort the clusters and extract the carbon and hydrogen atoms
        clusters = sorted(clusters, key=lambda array : len(array), reverse=True)
        bounds = arg_information.frame_data[2]
        hydrogen_atoms = np.asarray(arg_information.frame_data[0])[:,:3]
        carbon_atoms = np.asarray(arg_information.frame_data[1])[:,:3]
        clusters_sorted_and_atoms_extracted = time.time()
        print("{: <35} Time taken : {:.2f}s".format("Clusters sorted.", clusters_sorted_and_atoms_extracted - post_cluster))
        
        # Make a periodic copy of the points for ease in producing a file over 
        # periodic boundaries
        all_hydrogen_atoms = DuplicatePoints(hydrogen_atoms, bounds)
        all_carbon_atoms = DuplicatePoints(carbon_atoms, bounds)
        points_duplicated = time.time()
        print("{: <35} Time taken : {:.2f}s".format("Points duplicated.", points_duplicated - clusters_sorted_and_atoms_extracted))
        
        # Save PIMD input and driver file
        PIMD_in, run_cluster_sh = ReadTemplateFiles()
        
        # Make save directory
        parent_dir = os.getcwd()
        save_dir = parent_dir + "/" + datetime.now().strftime("cluster_analysis_%Y-%m-%d_%H-%M-%S")
        os.mkdir(save_dir)
        os.chdir(save_dir)
        
        # Make settings file and copy PIMD script
        WriteTemplateFiles(PIMD_in, run_cluster_sh, arg_information)
        
        # Process each cluster
        file_num = 1
        file_names = []
        curr_progress = -1
        cluster_start = time.time()
        
        min_cluster_index = arg_information.min_cluster; min_cluster_index_offset = min_cluster_index if min_cluster_index is not None else 0
        file_num += min_cluster_index_offset
        max_cluster_index = arg_information.max_cluster
        if max_cluster_index is not None and (max_cluster_index >= len(clusters)): max_cluster_index = None; print("Max clusters specified larger than total number of clusters, using all clusters. ")
        tot_num_cluster = len(clusters[min_cluster_index:max_cluster_index])
        for cluster in clusters[min_cluster_index:max_cluster_index]:
            
            # Make the cluster not be split over a periodic boundary
            cluster = np.asarray(cluster)
            for i in range(1, len(cluster)):
                cluster[i] = DistanceMinimise(cluster[i - 1], cluster[i], bounds)
            
            # Get min and maxes of the cluster and bounds of the save box
            mins = np.min(cluster, axis=0)
            maxs = np.max(cluster, axis=0)
            analysis_extensions = np.ones(3) * arg_information.box_extension
            analysis_mins = mins - analysis_extensions
            analysis_maxs = maxs + analysis_extensions
            
            # Find points in the box
            box_hydrogen_atoms = GetBoxedAtoms(all_hydrogen_atoms, analysis_mins, analysis_maxs)
            box_carbon_atoms   = GetBoxedAtoms(all_carbon_atoms, analysis_mins, analysis_maxs)
            
            # Skip all clusters with no carbon interaction
            if (len(box_carbon_atoms) > 0):
                # Get bounds and min
                mins_hydrogen = np.min(box_hydrogen_atoms, axis=0)
                maxs_hydrogen = np.max(box_hydrogen_atoms, axis=0)
                mins_carbon = np.min(box_carbon_atoms, axis=0)
                maxs_carbon = np.max(box_carbon_atoms, axis=0)
                all_mins = np.array([mins, mins_hydrogen, mins_carbon])
                all_maxs = np.array([maxs, maxs_hydrogen, maxs_carbon])
                mins = np.min(all_mins, axis=0)
                maxs = np.min(all_maxs, axis=0)
                bounds = maxs-mins
                
                # File write setup
                total_atom_count = len(box_hydrogen_atoms) + len(box_carbon_atoms)
                curr_file_name = "{}_{}".format(file_num, total_atom_count)
                file_names.append(curr_file_name)
                
                # Write the files
                WriteXYZ(curr_file_name, len(cluster), total_atom_count, bounds, 
                         mins, cluster, box_hydrogen_atoms, box_carbon_atoms)
                WriteLMP(curr_file_name, total_atom_count, bounds, mins, 
                         box_hydrogen_atoms, box_carbon_atoms)
            
            # Give updates
            if (int(( file_num - min_cluster_index_offset ) / tot_num_cluster * 10) > curr_progress):
                curr_progress += 1
                
                print("{: <35} Time taken : {:.2f}s".format("{: >3}% done! Cluster {: <5} written.".format(curr_progress * 10, file_num), time.time() - cluster_start))
                cluster_start = time.time()
            
            file_num += 1
            
        # Copy driver
        WriteDriver(file_names)
        
    # Perform a cluster analysis
    else:
        AnalyseClusterSize(arg_information)










'''

Hold the command line arguments


'''

class ArgInformation:
    
    def __init__(self, cluster_analysis, show_timings, cluster_size, frame_data, min_cluster_size, max_cluster_size, num_points, save_figures, box_extension, frame_used, file_name, min_cluster, max_cluster):
        self.cluster_analysis = cluster_analysis        # Bool
        self.show_timings = show_timings                # Bool
        self.cluster_size = cluster_size                # Float
        self.frame_data = frame_data                    # np.array([,,])
        self.min_cluster_size = min_cluster_size        # Float
        self.max_cluster_size = max_cluster_size        # Float
        self.num_points = num_points                    # Int
        self.save_figures = save_figures                # Bool
        self.box_extension = box_extension              # Float
        self.frame_used = frame_used                    # Int
        self.file_name = file_name                      # String
        self.min_cluster = min_cluster                  # Int
        self.max_cluster = max_cluster                  # Int
        
    '''
    
    Save the paramaters used for the current file

    '''
    def WriteContents(self):

        f = open("script_settings.txt", "w")
        f.write(datetime.now().strftime("Execution date   : %Y-%m-%d_%H-%M-%S\n\n"))
        f.write("File used        : {}\n".format(self.file_name))
        f.write("Frame            : {}\n".format(self.frame_used))
        f.write("Box extension    : {}\n".format(self.box_extension))
        f.write("Min cluster      : {}\n".format("Default (0)" if self.min_cluster is None else self.min_cluster))
        f.write("Max cluster      : {}\n".format("Default (Max)" if self.max_cluster is None else self.max_cluster))
        f.write("Timings shown    : {}\n\n".format(self.show_timings))
        f.write("Cluster analysis : {}\n".format(self.cluster_analysis))
        f.write("Min cluster size : {}\n".format(self.min_cluster_size))
        f.write("Max cluster size : {}\n".format(self.max_cluster_size))
        f.write("#cluster points  : {}\n".format(self.num_points))
        f.write("Figures saved    : {}\n".format(self.save_figures))
        f.close()
        
        
if __name__ == "__main__":
    
    start_time = time.time()
    arg_information = ParseArgs()
    WriteFiles(arg_information)
    print("{: <35} Time taken : {:.2f}s".format("Done!", time.time() - start_time))
