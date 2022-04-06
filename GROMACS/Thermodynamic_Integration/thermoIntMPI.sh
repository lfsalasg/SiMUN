#!/bin/bash

writeResult(){
  step=$1
  fsource=$2
  writeTo=$3

  n=0
  while IFS= read -r line
  do
        if [ $n -eq 6 ]
        then
            resultado=$line
        fi
        n=$((n+1))
    done <$fsource
    read -a fields <<< $resultado
  echo "$step   ${fields[1]}" >> $writeTo
}

setUp(){
    source=$1
    aim=$2
    countVdW=$3
    countCoulomb=$4
    startVdW=$5
    startCoulomb=$6
    #Create the log file
    echo "Starting Thermodynamic integration"
    echo "Starting Thermodynamic integration" > ${aim}'progress.log'
    echo "Beginnig at " $(date)
    echo "Beginnig at " $(date) >> ${aim}'progress.log'
    if [ $countVdW -gt 0 ]; then
        mkdir ${aim}'VdW'
        i=$startVdW
        while [ $i -le $countVdW ]
        do
            #Create the folders and copy the content in the template folder
            mkdir ${aim}'VdW/L_'${i}
            cp -a ${source}'.' ${aim}'VdW/L_'${i}'/'
            #Delete unnecessary files 
            if test -f ${aim}'VdW/L_'${i}'/equil_charge.mdp'; then
                rm ${aim}'VdW/L_'${i}'/equil_charge.mdp';
            fi

            if test -f ${aim}'VdW/L_'${i}'/final_charge.mdp'; then
                rm ${aim}'VdW/L_'${i}'/final_charge.mdp';
            fi
            #Edit the mdp template
            sed -i "s/&state/$i/" ${aim}'VdW/L_'${i}'/equil_vdw.mdp'
            sed -i "s/&state/$i/" ${aim}'VdW/L_'${i}'/final_vdw.mdp'  
            i=$(($i+1))
        done 
        echo "$(($i-$startVdW)) folders for the Van der Waals integration were created in ${aim}"
        echo "$(($i-$startVdW)) folders for the Van der Waals integration were created in ${aim}" >> ${aim}'progress.log'

    fi
    
    if [ $countCoulomb -gt 0 ]; then
        mkdir ${aim}'Coulomb'
        i=$startCoulomb
        while [ $i -le $countCoulomb ]
        do
            mkdir ${aim}'Coulomb/L_'${i}
            cp -a ${source}'.' ${aim}'Coulomb/L_'${i}'/'
            
            if test -f ${aim}'Coulomb/L_'${i}'/equil_vdw.mdp'; then
                rm ${aim}'Coulomb/L_'${i}'/equil_vdw.mdp';
            fi

            if test -f ${aim}'Coulomb/L_'${i}'/final_vdw.mdp'; then
                rm ${aim}'Coulomb/L_'${i}'/final_vdw.mdp';
            fi
            sed -i "s/&state/$i/" ${aim}'Coulomb/L_'${i}'/equil_charge.mdp'
            sed -i "s/&state/$i/" ${aim}'Coulomb/L_'${i}'/final_charge.mdp'
	    i=$(($i+1))
        done
        echo "$(($i-$startVdW)) folders for the Coulomb integration were created in ${aim}" 
        echo "$(($i-$startVdW)) folders for the Coulomb integration were created in ${aim}" >> ${aim}'progress.log'
    fi
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
      echo "-o,                         Specify a directory to store the output in"
      echo "-c                          Specify the initial structure (.gro file)"
      echo "-startVdW                    Initial state for the VdW lambda"
      echo "-endVdW                      Specify the last lambda state" 
      echo "-startCoulomb                Initial state for the Coulomb lambda"
      echo "-endCoulomb                  Specify the last lambda state"                  
      exit 0
      ;;
    -s)
      shift
      if test $# -gt 0; then
        TEMPLATE=$1
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
    -c)
      shift
      if test $# -gt 0; then
        INITCONF=$1
      else
        echo "No initial structure specified"
        exit 1
      fi
      shift
      ;;
    -startVdW)
      shift
      if test $# -gt 0; then
        STARTVDW=$1
      else
        echo "Please specify the initial state for lambda VdW"
        exit 1
      fi
      shift
      ;;
    -endVdW)
      shift
      if test $# -gt 0; then
        ENDVDW=$1
      else
        echo "Please specify the last state for lambda VdW"
        exit 1
      fi
      shift
      ;;
    -startCoulomb)
      shift
      if test $# -gt 0; then
        STARTCOL=$1
      else
        echo "Please specify the initial state for lambda Coulomb"
        exit 1
      fi
      shift
      ;;
    -endCoulomb)
      shift
      if test $# -gt 0; then
        ENDCOL=$1
      else
        echo "Please specify the last state for lambda Coulomb"
        exit 1
      fi
      shift
      ;;
    *)
      break
      ;;
  esac
done


setUp $TEMPLATE $OUTPUT $ENDVDW $ENDCOL $STARTVDW $STARTCOL
#runIntegration $INITCONF $ENDVDW $ENDCOL $OUTPUT $STARTVDW $STARTCOL
