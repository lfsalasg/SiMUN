#!/bin/bash


source=$1
aim=$2
#Initialization
gmx grompp -f ${source}minim.mdp -c ${source}init_box.gro -p ${source}topol_mod.top -o ${aim}minim.tpr
gmx mdrun -deffnm ${aim}minim -s ${aim}minim.tpr

#Equilibration
gmx grompp -f ${source}equil.mdp -c ${aim}minim.gro -p ${source}topol_mod.top -o ${aim}equil.tpr
gmx mdrun -v -deffnm ${aim}equil -s ${aim}equil.tpr
printf '19 20 21 23 0' | gmx energy -f ${aim}equil.edr -o ${aim}box.xvg

#Rescale the box
gmx analyze -f ${aim}box.xvg -w -b 800 > ${aim}box.log
n=0
while IFS= read -r line
do
    if [ $n -eq 9 ]
    then
        resultado=$line
    fi
    n=$((n+1))
done <${aim}"box.log"

read -a fields <<< $resultado
gmx editconf -f ${aim}equil.gro -density ${fields[1]} -o ${aim}scaled.gro

#Production
gmx grompp -f ${source}final.mdp -c ${aim}scaled.gro -p ${source}topol_mod.top -o ${aim}final.tpr
gmx mdrun -v -deffnm ${aim}final -s ${aim}final.tpr

echo 'Simulation finished.'

