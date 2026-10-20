#!/bin/bash
echo 'log.lammps'
#for i in $(find . -mindepth 1 -type d | grep -v act |  sort )  ; do   echo  $i; tail -n 1 $i/log.lammps ; done
for i in $(find . -mindepth 1 -type d | grep -v act |  sort )  ; do   echo  $i; tail -n 5010 $i/log.lammps | grep NN -B1  | tail -n 1 ; done
echo 'disagreements'
for i in $(find . -mindepth 1 -type d | grep -v act | sort )  ; do  if [ -f "$i/disagreement.txt" ]; then echo $i; cat $i/disagreement.txt| awk '{print $1 } ' | uniq | wc -l; fi; done

echo 'pruned disagreement'
for i in $(find . -mindepth 1 -type d     | grep -v act |sort )  ; do  if [ -f "$i/pruned_disagreement.txt" ]; then echo $i; cat $i/pruned_disagreement.txt| awk '{print $1 } ' | uniq | wc -l; fi; done
echo 'active_learning'
for i in $(find . -mindepth 1 -type d     |  sort)  ; do    if [ -d "$i/active_learn_POSCARs" ];  then echo $i; ls $i/active_lea* | wc -l; fi; done

