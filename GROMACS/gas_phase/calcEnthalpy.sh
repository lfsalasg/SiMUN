#source /usr/local/gromacs/gmx-2018.7/bin/GMXRC

rm -r tmp/
mkdir tmp/
mkdir results
for $t in "$@"
do
    cp one_molecule.mdp THC.itp one_molecule.gro topol.top tmp/
    cd tmp
    sed -i "s/&temperature/$t/" 'one_molecule.mdp'
    gmx grompp -f one_molecule.mdp -c one_molecule.gro -p topol.top -o one_molecule
    gmx mdrun -deffnm one_molecule -s one_molecule.tpr
    printf '9 10 11' | gmx energy -f one_molecule.edr -o T_$t
    cp T$t.xvg ../results/
    cd ../
    rm tmp/*
done

for t in "$@"
do
    gmx analyze -f 
    awk '/SS1/{print $0}' re
done
echo "Done!"