control(){
    source=$1
    if ! test -f "$source/equil.mdp"; then
        echo "equil.mdp file not found in $source"
        exit 1
    fi
    if ! test -f "$source/final.mdp"; then
        echo "final.mdp file not found in $source"
        exit 1
    fi
    if ! test -f "$source/topol.top"; then
        echo "topol.top file not found in $source"
        exit 1
    fi
}

setUp(){
    source=$1
    target=$2
    genVel=$3
    productionSteps=$4
    numSimulations=$5

    if [ $numSimulations -ge 0 ]; then
        i=0
        #echo "Starting Diffusion calculation" > ${target}'../progress.log'
        #echo "Beginnig at " $(date) >> ${target}'progress.log'
        while [ $i -le $numSimulations ]
        do
            #Create the folders and copy the content in the template folder
            mkdir ${target}'/'${i}
            cp -a ${source}'.' ${target}'/'${i}'/'
            #Generate the random initial temperature
            lowBond=$(($genVel-30))
            upBound=$(($genVel+30))
            initTemp=$(shuf -i $lowBond-$upBound -n 1)
            #Edit the mdp template
            
            
            sed -i "s/&initTemp/$initTemp/" ${target}'/'${i}'/equil.mdp'
            sed -i "s/&genVel/$genVel/" ${target}'/'${i}'/equil.mdp'
            sed -i "s/&genVel/$genVel/" ${target}'/'${i}'/final.mdp'
            sed -i "s/&productionSteps/$productionSteps/" ${target}'/'${i}'/final.mdp'  
            echo "Run ${i} set with initial temperature ${initTemp}" >> ${target}'../progress.log'
            i=$(($i+1))
        done 
        echo $i' folders were created.' >> ${target}'../progress.log'
    fi
}

runSimulations(){
    target=$1
    numSimulations=$2
    i=0
    echo $i
    while [ $i -le $numSimulations ]
    do
        #Run equilibration
        echo $(date) "Starting equilibration for run $i" >> ${target}'../progress.log'
        gmx grompp -f ${target}/${i}/equil.mdp -c ${target}/${i}/init_box.gro -p ${target}/${i}/topol_mod.top -o ${target}/${i}/equil.tpr
        gmx mdrun -deffnm ${target}/${i}/equil -s ${target}/${i}/equil.tpr
        #Scale box
        printf '19 20 21 23 0' | gmx energy -f ${target}/${i}/equil.edr -o ${target}/${i}/box.xvg
        gmx analyze -f ${target}/${i}/box.xvg -w -b 800 > ${target}/${i}/box.log
        n=0
        while IFS= read -r line
        do
            if [ $n -eq 9 ]
            then
                resultado=$line
            fi
            n=$((n+1))
        done <${target}/${i}/"box.log"

        read -a fields <<< $resultado
        gmx editconf -f ${target}/${i}/equil.gro -density ${fields[1]} -o ${target}/${i}/scaled.gro

        #Run production
        echo $(date) "Starting production for run $i" >> ${target}'../progress.log'
        gmx grompp -f ${target}/${i}/final.mdp -c ${target}/${i}/scaled.gro -p ${target}/${i}/topol_mod.top -o ${target}/${i}/final.tpr
        gmx mdrun -deffnm ${target}/${i}/final -s ${target}/${i}/final.tpr
        #Run MSD
        printf "2" | gmx msd -f ${target}/${i}/final.xtc -s ${target}/${i}/final.tpr -o ${target}/${i}/diff.xvg
        #Remove big files
        rm ${target}/${i}/*.trr 
        rm ${target}/${i}/*.xtc
        i=$(($i+1))
    done
    

}
checkVariable(){
    var=$OUTPUT
    if [ -z ${var+x} ]; then result=0; else result=1; fi
}

#PARSE THE INPUT ARGUMENTS
if [ $# -le 0 ];
then
    echo "No input arguments where provided. Use -h or --help for more information"
    exit 0
fi
while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo " "
      echo "[arguments]"
      echo " "
      echo "arguments:"
      echo "-h, --help                  show brief help"
      echo "-s,                         Specify the directory containing the required files"
      echo "-o,                         Specify the target folder"
      echo "-npt,                       Specify the file with the NPT tuples"
      echo "-nstep                      Specify the number of steps"
      echo "-nsim                       Number of simulations"                  
      exit 0
      ;;
    -s)
      shift
      if test $# -gt 0; then
        SOURCE=$1
      else
        echo "No source specified"
        exit 1
      fi
      shift
      ;;

    -o)
      shift
      if test $# -gt 0; then
        OUTPUT=$1
      else
        echo "No output dir specified"
        exit 1
      fi
      shift
      ;;
    -npt)
        shift
        if test $# -gt 0; then
        NPT=$1
        else
            echo "No NPT data file specified"
            exit 1
        fi
        shift
        ;;
    -nstep)
      shift
      if test $# -gt 0; then
        NSTEPS=$1
      else
        echo "Please specify the number of steps for each run"
        exit 1
      fi
      shift
      ;;
    -nsim)
      shift
      if test $# -gt 0; then
        NSIMS=$(($1-1))
      else
        echo "Please specify the number of simulations to run"
        exit 1
      fi
      shift
      ;;
    *)
      break
      ;;
  esac
done

control $SOURCE

# MAIN PROGRAM
echo "Starting Diffusion calculation" > $OUTPUT'/progress.log'
echo "Beginnig at " $(date) >> $OUTPUT'/progress.log'
while IFS= read -r line
do
    if [[ $line != \#* ]]; then
        temp=$(echo $line | awk '{print $1}') 
        pressure=$(echo $line | awk '{print $2}') 
        nmol=$(echo $line | awk '{print $3}') 
        #Create and copy the folders
        echo "Running dataset with T=$temp P=$pressure AND Nmol=$nmol" >> $OUTPUT/progress.log
        dir=$OUTPUT/Dataset_T_$temp\_P_$pressure\_nmol_$nmol
        mkdir $dir
        mkdir $dir/building
        
        cp $SOURCE/* $dir/building/
        rm $dir/building/*.gro
        cp $SOURCE/init_box_$nmol.gro $dir/building/init_box.gro

        mv $dir/building/topol.top $dir/building/topol_mod.top
        sed -i 's/&nmol/'$nmol'/' $dir/building/topol_mod.top
        sed -i 's/&press/'$pressure'/' $dir/building/equil.mdp
        
        # Run the  functions
        setUp $dir/building/ $dir/ $temp $NSTEPS $NSIMS
        runSimulations $dir/ $NSIMS
    fi
done < $NPT
