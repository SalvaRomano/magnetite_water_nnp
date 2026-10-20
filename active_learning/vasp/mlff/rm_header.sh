#!bin/bash



lines_rm=$((-50 - $NTOT -$NTOT ))                              #Number of lines of the last configuration

Last_Config=$(head -n 5 ML_ABN | tail -n 1 | awk '{print $1}')  #NUMBER OF CONFIG
Config_line=$(head -n 5 ML_ABN | tail -n 1)                     #HEADER LINE WHERE THE CONFIG INFO IS

Basis_line=$(head -n 33 ML_ABN | tail -n 1)                      #HEADER LINE WHERE THE BASIS FUNCTION INFO IS
Length_H=$(head -n 33 ML_ABN | tail -n 1 | awk '{print $1}')     #LENGTH OF THE BASIS FUNCTION OF H
Length_O=$(head -n 33 ML_ABN | tail -n 1 | awk '{print $2}')
Length_Fe=$(head -n 33 ML_ABN | tail -n 1 | awk '{print $3}')

Start_H=37                                  #EXACT STARTING AND ENDING LINE NUMBER OF H BASIS FUNC BLOCK
End_H=$(($Start_H+ $Length_H -1))
Start_O=$(( $End_H + 4  ))
End_O=$(( $Start_O + $Length_O -1   ))
Start_Fe=$(( $End_O +4  ))
End_Fe=$(( $Start_Fe + $Length_Fe -1    ))

#NUMBER OF LINES WHICH REFER IN THE BASIS FUNC BLOCKS TO LAST CONFIG THAT WILL BE REMOVED 
Removed_H=$(head -n $End_H ML_ABN | tail -n $Length_H    | awk "BEGIN {count=0} {if(\$1 == $Last_Config) { count+=1 }} END {print count}"  )
Removed_O=$(head -n $End_O ML_ABN | tail -n $Length_O    | awk "BEGIN {count=0} {if(\$1 == $Last_Config) { count+=1 }} END {print count}" )
Removed_Fe=$(head -n $End_Fe ML_ABN | tail -n $Length_Fe | awk "BEGIN {count=0} {if(\$1 == $Last_Config) { count+=1 }} END {print count}" )

#SURVIVED LINES
Survived_H=$(($Length_H-$Removed_H))
Survived_O=$(($Length_O-$Removed_O))
Survived_Fe=$(($Length_Fe-$Removed_Fe))

if [[ $Last_Config -ge 1 ]]; then
    #REMOVING BASIS FUNCTIONS LINES REFERRING TO THE LAST CONFIG
    sed -i "/ $Last_Config \+/d"  ML_ABN
    #UPDATE BASIS FUNCTION NUMBERS
    sed -i "s/$Basis_line\$/       $Survived_H   $Survived_O   $Survived_Fe/" ML_ABN
    #UPDATE TOTAL NUMBER OF CONFIG
    sed -i "s/$Config_line\$/         $(($Last_Config-1))/" ML_ABN
    #FINALLY REMOVING LAST CONFIGURATION
    head -n $lines_rm ML_ABN > temp_ML_ABN
    mv temp_ML_ABN ML_ABN
else echo "//////////////ERROR NO CONFIGURATION IN ML_ABN/////// "
fi
