#!/bin/bash


source=$1
aim=$2
#Initialization
gmx grompp -f ${source}minim.mdp -c ${source}init_box.gro -p ${source}topol_mod.top -o ${aim}minim.tpr
gmx mdrun -deffnm ${aim}minim -s ${aim}minim.tpr

#Equilibration
gmx grompp -f ${source}equil.mdp -c ${aim}minim.gro -p ${source}topol_mod.top -o ${aim}equil.tpr
gmx mdrun -v -deffnm ${aim}equil -s ${aim}equil.tpr

#Generate the desired files
printf '15 20 21 22 23 0' | gmx msd -f ${aim}equil.xtc -s ${aim}equil.tpr -o ${aim}density.xvg

echo 'Simulation finished.'
