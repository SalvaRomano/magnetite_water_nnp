#/bin/bash
for i in */; do
#if [ -f " $i/equi.dump" ]; then
last_occurrence_line=$(grep -n 'TIMESTEP' $i/equi.dump | tail -n 1 | cut -d ':' -f 1)
tail -n +"$((last_occurrence_line-1))"  $i/equi.dump   | tail -n +2 > $i/last.equi;
#fi;
done

