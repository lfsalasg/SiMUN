while IFS= read -r line
do
    if [[ $line != \#* ]]; then
        temp=$(echo $line | awk '{print $1}') 
        pressure=$(echo $line | awk '{print $2}') 
        nmol=$(echo $line | awk '{print $3}') 
        #Create and copy the folders
        echo $temp $pressure $nmol
        dir=$SOURCE/Dataset_T_$temp\_P_$pressure\_nmol_$nmol
        mkdir $dir
        cp $SOURCE/* $dir/
        mv $SOURCE/topol.top $SOURCE/topol_mod.top
        sed -i 's/&nmol/'$nmol'/' $dir/topol_mod.top
        # Run the  functions
        setUp $dir $dir $temp $NSTEPS $NSIMS
        #runSimulations $OUTPUT $NSIMS

    fi
done < $SOURCE"/NPT.dat"