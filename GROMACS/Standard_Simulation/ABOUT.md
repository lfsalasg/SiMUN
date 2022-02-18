#ABOUT

The current script runs a standard simulation following these steps:
1. A minimization step using the `minim.mdp` file.
2. A equilibration step using the `equil.mdp` file.
3. Scaling the box obtained during the equilibration stage to the average density. The script calculates the density starting at the 800ps. You may need to change this!
4. A production step using the `final.mdp` file. Be aware that this script was thought to produce on the NVT ensemble.

## How to use?

1. Prepare all your files (topologic, structure, and the mdp files) on a folder (this folder is the `<root_file>`). You may need to edit the given mdp files.
2. Create a folder to save all the results (This is the `<target_file>`). Results will have the same name as the mdp files.
3. Run the `startSimulation.sh` script as:
```
bash startSimulation.sh <root_file> <target_file>
```

## Stuff to do

- Add `help` option
- Let choose NVT or NPT

## Contributors

Luis Salas lfsalasg@unal.edu.co
