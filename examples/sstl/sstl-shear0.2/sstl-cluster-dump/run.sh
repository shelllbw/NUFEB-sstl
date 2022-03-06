#!/bin/bash
../../../../lammps/src/lmp_mpi -in Inputscript.lammps 
for i in {1..287}
do
   pre=$((i-1))
   echo $i
   pre_grid=$((pre*720))
   grid=$((i*720))
   
   sed -i -e "s/t${pre}.txt/t${i}.txt/g" ./Inputscript.lammps 
   sed -i -e "s/_a${pre_grid}.vti/_a${grid}.vti/g" ./Inputscript.lammps 
   ../../../../lammps/src/lmp_mpi -in Inputscript.lammps 
done
