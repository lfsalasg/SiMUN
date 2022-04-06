#!/bin/bash


initconf=$1
countVdW=$2
countCoulomb=$3
aim=$4
startVdW=$5
startCoulomb=$6

if [ $2 = 'VdW' ]; then 
	#Run equilibration 
	mpirun gmx_mpi grompp -f $target/equil_vdw.mdp -c $target/init_box.gro -p $target/topol.top -o $target/equil
	mpirun -np 2 gmx_mpi mdrun -deffnm $target/equil -s $target/equil.tpr -ntopm 4

	#Run production
	mpirun gmx_mpi grompp -f $target/final_vdw.mdp -c $target/equil.gro -p $target/topol.top -o $target/final
	mpirun -np 2 gmx_mpi mdrun -deffnm $target/final -s $target/final.tpr -ntopm 4

	#Run analysis
	printf "20 0" | mpirun gmx_mpi energy -f $final/final.edr -o $final/gibbs.xvg
	mpirun gmx_mpi analyze -f $target/gibbs.xvg > $final/gibbs.log
fi

if [ $2 = 'Coulomb' ]; then
	#Run equilibration 
	mpirun gmx_mpi grompp -f $target/equil_charge.mdp -c $target/init_box.gro -p $target/topol.top -o $target/equil
	mpirun -np 2 gmx_mpi mdrun -deffnm $target/equil -s $target/equil.tpr -ntopm 4

	#Run production
	mpirun gmx_mpi grompp -f $target/final_charge.mdp -c $target/equil.gro -p $target/topol.top -o $target/final
	mpirun -np 2 gmx_mpi mdrun -deffnm $target/final -s $target/final.tpr -ntopm 4

	#Run analysis
	printf "20 0" | mpirun gmx_mpi energy -f $final/final.edr -o $final/gibbs.xvg
	mpirun gmx_mpi analyze -f $target/gibbs.xvg > $final/gibbs.log
fi

