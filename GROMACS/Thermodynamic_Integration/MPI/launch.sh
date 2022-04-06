#!/bin/bash

dir=$1
start=$2
end=$3
type=$4

for i in {$start..$end}
do
	runSimulation.sh $dir/$type/L_$i $type
done
