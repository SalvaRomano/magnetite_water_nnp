#!/bin/bash


#creates 5 (you can change) run_i folders corresponding to 5 different trainings (5 committee members)

#inside ./ there has to be: input.nn input.data.   "1307" has to be the random seed in input.nn
#check if the random seeds are matching in input.nn and sed command down!!!

for k in {1..4}; do
mkdir run_$k
cp input.data validation.data input.nn run_$k
cd run_$k
random_number=$RANDOM
sed -i "s/1307/$RANDOM/" input.nn
cd ../
done

