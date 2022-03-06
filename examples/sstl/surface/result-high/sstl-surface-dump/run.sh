#!/bin/bash

../../../../../lammps/src/./lmp_mpi -in Inputscript.lammps 

for i in {1..78}
do
   pre=$((i-1))
   grid=$((i*1000))
   grid_pre=$(((i-1)*1000))

   echo $pre
   echo $grid
   echo $grid_pre
   echo 
   sed -i -- "s/a${grid_pre}.vti/a${grid}.vti/g" "Inputscript.lammps"
   sed -i -- "s/surface${pre}.txt/surface${i}.txt/g" "Inputscript.lammps"
   ../../../../../lammps/src/./lmp_mpi -in Inputscript.lammps 
done
