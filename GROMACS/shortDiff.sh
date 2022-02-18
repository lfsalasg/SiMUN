setUp(){
    source=$1
    target=$2
    genVel=$3
    productionSteps=$4
    numSimulations=$5

    if [ $numSimulations -ge 0 ]; then
        i=0
        echo "Starting Diffusion calculation" > ${target}'progress.log'
        echo "Beginnig at " $(date) >> ${target}'progress.log'
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
            echo "Run ${i} set with initial temperature ${initTemp}" >> ${target}'progress.log'
            i=$(($i+1))
        done 
        echo $i' folders were created.' >> ${target}'progress.log'
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
        echo $(date) "Starting equilibration for run $i" >> ${target}'progress.log'
        gmx grompp -f ${target}/${i}/equil.mdp -c ${target}/${i}/init_box.gro -p ${target}/${i}/topol_mod.top -o ${target}/${i}/equil.tpr
        gmx mdrun -deffnm ${target}/${i}/equil -s ${target}/${i}/equil.tpr
        #Run production
        echo $(date) "Starting production for run $i" >> ${target}'progress.log'
        gmx grompp -f ${target}/${i}/final.mdp -c ${target}/${i}/equil.gro -p ${target}/${i}/topol_mod.top -o ${target}/${i}/final.tpr
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
      echo "-t,                         Specify the target temperature"
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
    -t)
      shift
      if test $# -gt 0; then
        TEMPERATURE=$1
      else
        echo "No target temperature specified"
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

setUp $SOURCE $OUTPUT $TEMPERATURE $NSTEPS $NSIMS
runSimulations $OUTPUT $NSIMS

