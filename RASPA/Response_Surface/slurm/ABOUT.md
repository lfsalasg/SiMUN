# ABOUT

The current script runs a response surface used to vary the temperature and pressure for a
GCMC simulation. Further adaptations of the script can modify other variables of the 
configuration files.

## Requirements

- Python >=3.8.
- RASPA installed either as a module or in the node
- The slurm service running with enough permissions to queue jobs

## How to use?

1. Prepare all your files on a folder (this is the `<root_dir>`) and include the `launch.py` and `run.sh` files in this directory.
2. Edit the `fugacity.dat` file to include the points of the response surface. For now, each point consists of the pressure (in Pa), the temperature (in K) and the fugacity coefficient.
3. Run the script `launch.py` from the `<root_dir>`. It can be run as:

```terminal
python launch.py <pname> <walltime>
```

where `<pname>` is the prefix of the `<output_dirs>` and `<walltime>` is the walltime for the job (in minutes). Further information of this script can be obtained by running `python launch.py --help`.

### Output files

The script will generate a series of directories with the name `<output_dirs>_#` where `#` corresponds to the index of the point. Each generated folder contains the result of the simulation.

## Contributors

Luis Salas lfsalasg@unal.edu.co