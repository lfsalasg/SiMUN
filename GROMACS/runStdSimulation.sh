#!/bin/bash


source=$1
aim=$2
try=$3
#Initialization
gmx grompp -f ${source}minim.mdp -c ${source}init_box.gro -p ${source}topol.top -o ${aim}minim.tpr
gmx mdrun -deffnm ${aim}minim -s ${aim}minim.tpr

#Equilibration
gmx grompp -f ${source}equilibration/try${try}.mdp -c ${aim}minim.gro -p ${source}topol.top -o ${aim}equil.tpr
gmx mdrun -v -deffnm ${aim}equil -s ${aim}equil.tpr

#Production
gmx grompp -f ${source}final.mdp -c ${aim}equil.gro -p ${source}topol_mod.top -o ${aim}final.tpr
gmx mdrun -v -deffnm ${aim}final -s ${aim}final.tpr

#Generate the desired files
printf '3' | gmx msd -f ${aim}final.xtc -s ${aim}final.tpr -o ${aim}diff.xvg

echo 'Simulation finished.'

