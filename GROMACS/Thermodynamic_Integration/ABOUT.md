# ABOUT

With this script you can run a thermodynamic integration and segment it to run on different nodes. **BE AWARE** This script is not fully finished so some options may not work properly. It would be nice to optimice this script, see the _stuff to do_ below.

The current script works as follows:

1. For each point listed on the mdp files, the script runs a NPT equilibriation and an NPT production.
2. The first point will take the initial configuration given by argument in the script. This starting configuration must be minimized, because the script does not make any minimizations step.
3. The next point will take the configuration generated by the production of the previous step and will use it as starting configuration.

## How to use it?

1. Organize all your files on a folder, this will be the `<root_folder>`. See the _listing lambda points_ for more information about handling thermodynamic integrations.
2. Create the `<target_folder>`
3. Run the script. The script has the following arguments (optionally you can do `bash thermoInt.sh -h` to get some help)
-	s: The source folder with all the topologic and mdp files
-	c: The initial configuration file (a .gro file)
-	o: The target folder
-	startVdW: The lambda state where the system starts (using the VdW mdp files)
-	endVdW: The lambda state where the system ends (using the VdW mdp files)
-	startCoulomb: The lambda state where the system starts (using the charge mdp files)
-	endCoulomb: The lambda state where the system ends (using the charge mdp files)

So for example, you want to run on a node the first five states of the VdW integration, the script may look as follows:
```
bash thermoInt.sh -s <root_file> -o <target_file> -c <init_box.gro> -startVdW 0 -endVdW 4 -startCoulomb 0 -endCoulomb 0
```

4. The script creates a directory for each point and a `progress.log` file where the execution progress of each point is printed.

## Listing lambda points

The thermodynamic integration run each simulation on the specified lambda state (for more information about what thermodynamic integration works see [this paper]() ). You need to specify the thermodynamic states in the mdp files. See the corresponding section on the mdp files for more info. Some key parameters are:
- The vdw\_lambdas and coul\_lambdas: Here we list all the states. Remember that each point shall be between 0 and 1.
- init\_lambda\_state: Specifies the state from the states' table (the first element of the table is indexed as the 0th element)
- couple-moltype: It is the residue over which the activation occurs.
- couple-lambda0: Initial state
- couple-lambda1: Final state

## Stuff to do

**MAJOR CHANGES ARE PLANNED FOR THIS SCRIPT!**

### Some known bugs
- Specifying the ranges to operate is somewhat tricky. If you want to integrate modifying the VdW state, run all with startCoulomb and endCoulomb set to 0. If you want to integrate modifying the Coulomb state, run all states with startVdW and endVdW set to 0. On the future is planned to change the way the integration range (and states) work.


