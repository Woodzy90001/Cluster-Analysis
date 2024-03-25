rm PIMD_*.in
    
input_files=(1 0.1 0.01 0.001 0.0001 0.00001 0.000001)
for inp_file in "${input_files[@]}"
do
        cp PIMD-Template.in PIMD_${inp_file}.in
        sed -i -e  "s/TIMESTEP/${inp_file}/g" PIMD_${inp_file}.in
        sbatch run_timestep.sh ${inp_file}
done