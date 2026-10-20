#!/bin/bash
 for i in */; do    i="${i%/}";  grep -Fwf <(sort pressure_good_structures_$i) $i/disagreement.txt > $i/reduced_disagreement.txt; done
