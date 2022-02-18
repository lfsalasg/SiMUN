#ABOUT

The current script runs a standard simulation following these steps:
1. A minimization step using the `minim.mdp` file.
2. A equilibration step using the `equil.mdp` file.
3. Scaling the box obtained during the equilibration stage to the average density. The script calculates the density starting at the 800ps. You may need to change this!
4. A production step using the `final.mdp` file. Be aware that this script was thought to produce on the NVT ensemble.

## How to use?

1. Prepare all your files (topologic, structure, and the mdp files) on a folder (this folder is the `<root_dir>`). You may need to edit the given mdp files.
2. Create a folder to save all the results (This is the `<target_dir>`). Results will have the same name as the mdp files.
3. Run the `startSimulation.sh` script as:
```
bash startSimulation.sh <root_dir> <target_dir>
```
### Required files in the root_dir
- minim.mdp: The minimization config file
- equil.mdp: The equilibration config file
- final.mdp: The production config file
- topol_mod.top: The file with all the topologic information
- init_box.gro: The initial box for the MD simulation

### Output files 
- edr, log, tpr, ttr, xtc and gro files of the minimization, equilibration and production steps.
- The average size of the box between 800 and the final point of the equilibration (`box.log`)

## Stuff to do

- Add `help` option
- Let choose NVT or NPT

## Contributors

Luis Salas lfsalasg@unal.edu.co
