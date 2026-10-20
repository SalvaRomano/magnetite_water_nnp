#!/bin/bash
for i in */; do cd $i;   if grep -q "Total wall" log.lammps; then  echo 'FINISHED log.lammps'; cat log.lammps | grep -v WAR  | grep -v NNP | awk '/Step Temp/,0' |  tail -n +2 | awk '/Loop time/{p=1} !p' > table_log; else echo 'NOT FINISHED log.lammps'; cat log.lammps | grep -v WAR  | grep -v NNP | awk '/Step Temp/,0' | tail -n +2 | head -n -1 > table_log; fi;       cd ../; done

