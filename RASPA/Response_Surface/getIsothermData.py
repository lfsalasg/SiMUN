import argparse
import glob
import os
import re
import shutil
from datetime import datetime

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# 1. FUNCTIONS

#Internal functions

def _seek(filename,expr,starts_with=False,ends_with=False):
    #start = '\A' if starts_with else '.*'
    #end = '\Z' if ends_with else '.*'
    #exp = start + expr + end
    exp = '.*' + expr + '.*'
    with open(filename) as file:
        txt = file.read()
        matches = re.findall(exp,txt)

    result = list()
    for row in matches:
        if (not starts_with and not ends_with) or (starts_with and row.startswith(expr)) or (ends_with and row.endswith(expr)):
            result.append(re.split(r'\s+',row))

    return result

def _create_run_folder(dir_name,source_files,param):
    """
    Description
    ------------
    Create the folder with the corresponding files to run a single RASPA simulation

    Parameters
    ------------
    name: (str) Name of the folder
    source_files: ([str]) Path where the template files are available
    param:(dict) Dictionary with the keys and values to replace in the simulation.input

    Returns
    -----------
    (bool)
    """

    try:
        os.mkdir(dir_name)
    except FileExistsError:
        print("WARNING: Directory already exists")

    for file in source_files:
        shutil.copyfile(file,dir_name)
        
    return 1

def _replace_keys(filename,param):
    """
    Description
    ------------
    Replace keys according to the information given in param

    Parameters
    ------------
    param:(dict) Dictionary with the keys and values to replace in the simulation.input

    Returns
    -----------
    (bool)
    """    

# CLI linked functions
def cli_status(args):
    root_dir = args.root_dir + '/' if args.root_dir[-1] != '/' else args.root_dir
    pname = root_dir[:-1] if args.pname == 'same' else args.pname
    pname += '_'
    points = glob.glob(root_dir+pname+'*')
    print(f"{len(points)} Points found")
    for p in points:
        files = glob.glob(p+'/Output/System_0/*.data')
        r = _seek(files[0],'Simulation finished,')
        if len(r):
            print(f"Point {p}: Done")
        else:
            print(f"Point {p}: Running")

def cli_control(args):
    """
    Description
    ------------
    An utility command to check information about multiple simulations at once
    
    Parameters
    -----------
    - args: A Namespace containing the following properties:
        - root_dir: (str) Root directory for each simulation
        - pname: (str) Project name or "same". If "same", the variable will take the same name as the root folder
        - point: ([int]) A list of the points to check. By default it will check every point
        - evolution: (bool) Print the ensamble average against the number of cycles
        - plot: (bool) Plot the results (if available)
        - save: (bool) Save the dataset
    """
    root_dir = args.root_dir + '/' if args.root_dir[-1] != '/' else args.root_dir
    pname = root_dir[:-1] if args.pname == 'same' else args.pname
    pname += '_'
    points = list()
    if not len(args.point):
        points = glob.glob(root_dir+pname+'*')
    else:
        for i in args.point:
            points.append(root_dir+pname+str(i))
    
    if args.evolution:
        data = dict()
        for p in points:
            files = glob.glob(p+'/Output/System_0/*.data')
            cycle = _seek(files[0],'Current cycle:',True)
            evolution = _seek(files[0],'absolute adsorption') 
            evolution = evolution[len(evolution) - len(cycle):]
            history = list()
            for c, e in zip(cycle,evolution):
                history.append([c[2],e[5][:-1]])
            df = pd.DataFrame(history,columns=['Cycle','Loading'],dtype=np.float32)
            data[p] = df
        for p in data:
            print(data[p])
            if args.save:
                df.to_csv(args.save, index=False)
            plt.plot(data[p]['Cycle'],data[p]['Loading'])
        plt.savefig("test.png")
        plt.show()
        return 1
    data = list()
    for p in points:
        files = glob.glob(p+'/Output/System_0/*.data')
        done = True if len(_seek(files[0],'Simulation finished')) else False
        if done:
            data.append([p,'Done'])
            continue
        runtime = _seek(files[0],'Current cycle')
        data.append([p,str(runtime[-1][2]) + '/' + str(runtime[-1][5])])
        #print(f"Point {p}: Cycle "+runtime[-1][2]+" of "+runtime[-1][5])

    df = pd.DataFrame(data,
        columns=['point','cycle'],
    )

    df.sort_values(by='point',inplace=True)
    print(df)
    return 1

def cli_adsorbed(args):
    """
    Description
    -----------
    Generates the data for the isotherm using a set of simulated points

    Parameters
    ----------
    - args: A Namespace containing the following properties
        - root_dir: (str) Root directory for each simulation
        - pname: (str) Project name or "same". If "same", the variable will take the same name as the root folder
        - units: (str) Units to show for the adsorbed molecules
        Other options are:
        - save: (bool) Wether or not save the dataframe with results to a file
        - name: (str)  Name of the saved file
        - plot: (bool) Plot the results using the builtin pyplot library
    Returns
    ---------
    df: (Dataframe) A Dataframe containing the adsorption isotherm
    """
    root_dir = args.root_dir + '/' if args.root_dir[-1] != '/' else args.root_dir
    pname = root_dir[:-1] if args.pname == 'same' else args.pname
    pname += '_'
    points = glob.glob(root_dir+pname+'*')
    units = {
        'mgg':{
            'abs':'Average loading absolute \[milligram\/gram',
            'exc':'Average loading excess \[milligram\/gram',
            'index':6
        },
        'molkg':{
            'abs':'Average loading absolute \[mol\/kg',
            'exc':'Average loading excess \[mol\/kg',
            'index':6
        },
        'mcuc':{
            'abs':'Average loading absolute \[molecules\/unit',
            'exc':'Average loading excess \[molecules\/unit',            
            'index':6
        },
        'ccg':{
            'abs':'Average loading absolute \[cm\^3 \(STP\)\/gr',
            'exc':'Average loading excess \[cm\^3 \(STP\)\/gr',
            'index':7
        }
    }
    index = units[args.units]['index'] 
    data = list()

    for p in points:
        files = glob.glob(p+'/Output/System_0/*.data')
        t = _seek(files[0],'External temperature')
        p = _seek(files[0],'External Pressure')
        abs = _seek(files[0],units[args.units]['abs'])
        exc = _seek(files[0],units[args.units]['exc'])
        if len(abs) and len(exc):
            data.append([t[0][2],p[0][2],abs[0][index],abs[0][index+2],exc[0][index],exc[0][index+2]])
        else:
            print(f'WARNING: Point {p} does not have enough information')
    df = pd.DataFrame(data,
        columns=['Temperature [K]','Pressure [Pa]','Loading Abs. '+args.units,'Error','Loading Exc. '+args.units,'Error'],
        dtype=np.float32)
    df.sort_values(by=['Pressure [Pa]'],inplace=True)
    #Output options
    print(df)
    if args.save:
        if not args.name:
            filename = pname + "_" + datetime.now().strftime('%Y%m%d_%H%M%S')
        else:
            filename = args.name
        
        df.to_csv(filename+'.dat',sep='\t')

    if args.plot:
        x_values = df['Pressure [Pa]'].tolist()
        abs_values = df['Loading Abs. '+args.units].tolist()
        exc_values = df['Loading Exc. '+args.units].tolist()
        plt.figure(figsize=(9,3))
        plt.subplot(121)
        plt.plot(x_values,abs_values,'k--o')
        plt.ylabel('Loading Abs. '+args.units)
        plt.xlabel('Pressure [Pa]')
        plt.subplot(122)
        plt.plot(x_values,exc_values,'k--o')
        plt.ylabel('Loading Exc. '+args.units)
        plt.xlabel('Pressure [Pa]')        
        plt.show()
    return df

def cli_run(args):
    pass

    
#points = [1,2,3,4]
#for i in points:
#    files = glob.glob(f'sample/{i}/*.data')
#    r = _seek(files[0],'Average loading absolute \[milligram')
#    print(r[0][6])


# 1. ENTRY POINT
# Top level parser
parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='Sub-command help')

# Modules handling

# Subparser status
parser_status = subparsers.add_parser('status')
parser_status.add_argument('root_dir',type=str,help='Root directory where the isotherms are running')
#parser_status.add_argument('-b',type=int,default=1,help='Check points starting from this')
#parser_status.add_argument('-e',type=int,default=-1,help='Check points up to this')
parser_status.add_argument('--pname',type=str,default='same',help='Project name. It will be the same as the folder unless specified')
parser_status.set_defaults(func=cli_status)

# Subparser adsorbed
parser_adsorbed = subparsers.add_parser('adsorbed')
parser_adsorbed.add_argument('root_dir',type=str,help='Root directory where the isotherms are running')
parser_adsorbed.add_argument('--name',type=str,help='Save the results into this location. If not specified, the name is guessed and file is saved in current directory')
parser_adsorbed.add_argument('--plot',action='store_true',help='Plot the result using the pyplot window')
parser_adsorbed.add_argument('--pname',type=str,default='same',help='Project name. It will be the same as the folder unless specified')
parser_adsorbed.add_argument('--save',action='store_true',help='Save calculations to a .dat file')
parser_adsorbed.add_argument('--units',choices=['mgg','molkg','mcuc','ccg'],default='mcuc',help='Selec the units for the results')
parser_adsorbed.set_defaults(func=cli_adsorbed)

# Subparser control
parser_control = subparsers.add_parser('control')
parser_control.add_argument('root_dir',type=str,help='Root directory where the isotherms are running')
parser_control.add_argument('point',type=int,nargs='*',help='Points to be checked')
parser_control.add_argument('--evolution',action='store_true',help='Show how the average absolute loading changed over time.')
parser_control.add_argument('--plot',action='store_true',help='Plot results')
parser_control.add_argument('--pname',type=str,default='same',help='Project name. It will be the same as the folder unless specified')
parser_control.add_argument('--save', type=str, help='Save the cycle on a text file with the given name')
parser_control.set_defaults(func=cli_control)

# Subparser run
parser_run = subparsers.add_parser('run')
parser_run.add_argument('root_dir',type=str,help='Root directory where the isotherms will be stored')
parser_run.add_argument('-b',type=int,default=1,help='Start from this point of the file')
parser_run.add_argument('-e',type=int,default=-1,help='Run up to this point of the file')
parser_run.add_argument('-t',type=int,default=240,help='Simulation wall time (in minutes)')
parser_run.add_argument('--cluster',type=str,choices=['local','slurm'], default='local',help='If you are running in a clustesr, use this option followed by the corresponding cluster (gridengine, torque, slurm) in order to use the proper job manager. Use local to run on local machine')
parser_run.add_argument('--fugacity',type=str,default='fugacity.dat',help='File with the points to simulate')
parser_run.add_argument('--pname',type=str,default='same',help='Project name. It will be the same as the folder unless specified')
parser_run.add_argument('--template',type=str,default='same',help='Folder with the template files. If not specified, the script search the files inside the root directory')

# Subparser continue


args = parser.parse_args()
args.func(args)
exit()
