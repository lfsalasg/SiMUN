#ABOUT

The current script runs a series of standard simulations following these steps:
1. A equilibration step using the `equil.mdp` file.
2. Scaling the box obtained during the equilibration stage to the average density. The script calculates the density starting at the 800ps. You may need to change this!
3. A production step using the `final.mdp` file. Be aware that this script was thought to produce on the NVT ensemble.
4. The calculation of the MSD from 0 to the last frame with default step (You can change this in the source code)

## How to use?
1. This scripts generates random initial configurations by giving a random initial temperature around the target temperature of the simulation. However you have to provide the **minimized** initial boxes.
1. Prepare all your files (topologic, structure, and the mdp files) on a folder (this folder is the `<root_dir>`). You may need to edit the given mdp files. All the initial boxes must be named as `init_box_<size>` and the size must match the NPT.dat (see below) information
2. The NPT.dat file contains the conditions at which each simulation should run. The first column refers to the temperature, the second one to the pressure and the third one to the number of molecules. The number of molecules should match with the name of a initial box in your `<rot_dir>` folder.
2. Create a folder to save all the results (This is the `<target_dir>`)
3. Run the `sizeEffect.sh` script as:
```
bash sizeEffect.sh -s <root_dir> -o <target_dir> -npt <NPT_filename> -nstep <Number of steps in the production stage> -nsim <Number of tries per point>
```
### Required files in the root_dir
- equil.mdp: The equilibration config file
- final.mdp: The production config file
- topol_mod.top: The file with all the topologic information
- init_box_size.gro: The initial box for the MD simulation

### Output files 
This script generates a tree with all the results. You can process them with the msdAverage.py script
## Stuff to do


## Contributors

Luis Salas lfsalasg@unal.edu.co
