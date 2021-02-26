#!/usr/bin/env bash

make 
([ $? -eq 0 ] && echo "Make success!") || (echo "Make failure!" && exit 1;) 


# deleting old data file
rm M1* M2* M3*

###################
num_time_steps=50
p_array=(16 36 49 64)
n_array=(256 1024 4096 16384 65536 262144 1048576)
#n_array=(256 1024 4096 16384 65536 262144)


for i in {1..5}; do
    for P in ${p_array[*]}; do
        python hostA.py	
        for N in ${n_array[*]};do
            printf "Iteration: %d Processes: %d - data points : %d\n" $i $P $N
            mpirun -np $P ./s1 $N $num_time_steps >> "M1_data_$P.txt" || exit  
            mpirun -np $P ./s2 $N $num_time_steps >> "M2_data_$P.txt" || exit 
            mpirun -np $P ./s3 $N $num_time_steps >> "M3_data_$P.txt" || exit 
        done
    done
done


# creating plots 

echo "Plotting Boxplot:\n"
python3 plot.py
echo "Created Plot.\n"
