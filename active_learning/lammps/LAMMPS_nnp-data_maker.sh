#!/bin/bash

index=0


new_committee_LAMMPS='ninth_committee'          #name of committee
committee_leader='4'                            #is the NNP that guides the MD, its weights and scaling are in the principal folder "new_committee_LAMMPS", the other's weights are in nnp-data-i
weights=(  22 25 14  17  14  )                    #change the weights from learning curve!!


for i in {1..5}; do
        mkdir "nnp-data-$i"
        cp run_$i/train/weights.001.0000${weights[$index]}.out   nnp-data-$i/weights.001.data
        cp run_$i/train/weights.008.0000${weights[$index]}.out   nnp-data-$i/weights.008.data
        cp run_$i/train/weights.026.0000${weights[$index]}.out   nnp-data-$i/weights.026.data
        ((index++))
done


mkdir "$new_committee_LAMMPS"

cp "run_$committee_leader/scaling.data"  $new_committee_LAMMPS
cp "run_$committee_leader/input.nn"  $new_committee_LAMMPS

mv nnp-data* $new_committee_LAMMPS

cd $new_committee_LAMMPS

mv nnp-data-$committee_leader/* ./
rm -r "nnp-data-$committee_leader"
mv nnp-data-5 "nnp-data-$committee_leader"

pwd > info.txt

mkdir MD
cd ../
