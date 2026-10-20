#!/bin/bash

index=0


new_committee_LAMMPS='tenth_committee'          #name of committee
committee_size='3'                              #number of NNP in the committee
committee_leader='1'                            #its NNP-weights and scaling are in the principal folder "new_committee_LAMMPS" and the other's weights are in nnp-data-i (in validation mode is that NNP that guides the MD)

weights=(  16 21 21  )                    #change the weights from learning curve!!

for (( i=1; i<=$committee_size; i++ )); do
        mkdir "nnp-data-$i"
        echo $i
        cp run_$i/train/weights.001.0000${weights[$index]}.out   nnp-data-$i/weights.001.data
        cp run_$i/train/weights.008.0000${weights[$index]}.out   nnp-data-$i/weights.008.data
        cp run_$i/train/weights.026.0000${weights[$index]}.out   nnp-data-$i/weights.026.data
        ((index++))
done


mkdir "$new_committee_LAMMPS"

cp "run_$committee_leader/scaling.data"  $new_committee_LAMMPS
cp input.data input.nn $new_committee_LAMMPS



mv nnp-data* $new_committee_LAMMPS

cd $new_committee_LAMMPS

mv nnp-data-$committee_leader/* ./
rm -r "nnp-data-$committee_leader"
mv "nnp-data-$committee_size" "nnp-data-$committee_leader"

mkdir MD
cd ../

