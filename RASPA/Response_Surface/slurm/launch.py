import argparse
import datetime
import os
import shutil
import subprocess
import sys

VERSION = '0.2.0'

def replace_in_place(file:str, keyval:dict):
    with open(file, 'r') as f:
        content = f.read()
    
    for k, v in keyval.items():
        content = content.replace(k, str(v))
    
    with open(file, 'w') as f:
        f.write(content)

def timedelta_to_str(delta:datetime.timedelta):
    days = delta.days
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    _time = datetime.time(hours, minutes, seconds)

    return str(days) + '-' + _time.strftime('%H:%M:%S')

def parse_conditions(source:str)->list:
    data = []
    with open(source, 'r') as f:
        for i, line in enumerate(f):
            raw_data = line.split()
            if len(raw_data) < 3:
                raise ValueError('At least three values are expected for each line of the conditions file')

            data.append({
                'pressure': raw_data[0],
                'temperature': raw_data[1],
                'phi': raw_data[2],
                'job_number': i + 1
            })
    
    return data


def create_dirs(target_name:str, conditions:list, walltime:datetime.timedelta):
    content = os.listdir('./')
    for condition in conditions:
        target_dir = target_name + '_' + str(condition['job_number'])
        os.mkdir(target_dir)
        
        for file in content:
            shutil.copy(
                file,
                os.path.join(target_dir, file)
            )
        replace_in_place(os.path.join(target_dir, 'simulation.input'), condition)
        replace_in_place(
            os.path.join(target_dir, 'run.sh'), 
            {
                'maxTime': timedelta_to_str(walltime),
                'filename': target_dir
            }
        )

def run_simulations(target_name:str, conditions:list):
    for condition in conditions:
        target_dir = target_name + '_' + str(condition['job_number'])
        print(f'Running point {condition["job_number"]} with P={condition["pressure"]} and T={condition["temperature"]}')
        subprocess.run(['sbatch', 'run.sh'], cwd=target_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python script to launch multiple adsorption points on RASPA using slurm as worload manager')
    parser.add_argument('pname', type=str, help='Name of the simulation')
    parser.add_argument('walltime', type=int, help='Max time to run each simulation')
    parser.add_argument('-b', type=int, default=1, help='start at this line of the conditions file')
    parser.add_argument('-e', type=int, default=0, help='end at this line of the conditions file')
    parser.add_argument('-d', type=str, default='fugacity.dat', help='name of the conditions file. Default is fugacity.dat')
    parser.add_argument('-v', action='version', version=VERSION)

    args = parser.parse_args()

    walltime = datetime.timedelta(minutes=args.walltime)
    conditions = parse_conditions(args.d)
    start = args.b - 1
    end = len(conditions) if args.e == 0 else args.e
    conditions = conditions[start:end]
    create_dirs(args.pname, conditions, walltime)
    run_simulations(args.pname, conditions)