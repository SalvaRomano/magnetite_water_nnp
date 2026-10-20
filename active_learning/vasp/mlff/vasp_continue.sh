#!/bin/bash                                                                     
#SBATCH -J SF                                           
#SBATCH -N 3                                                                    
#SBATCH --partition=zen3_0512                                                  
#SBATCH --qos=zen3_0512_devel
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-core=1                                                     

module purge

module load intel-oneapi-compilers/2022.0.2-gcc-11.2.0-yzi4tsu intel-oneapi-mpi\
/2021.4.0-gcc-11.2.0-2e7zm7z intel-oneapi-mkl/2022.0.2-gcc-11.2.0-2ugv2zs
module load hdf5/1.12.1-gcc-11.2.0-523hswj

export I_MPI_PIN_RESPECT_CPUSET=0
export OMP_NUM_THREADS=1



# Next to where this script is executed there has to be placed "setup", a folder which contains the initial files for VASP:
# INCAR  KPOINTS  ML_ABN (THIS COULD BE CHANHGED)  POSCAR  POTCAR  vdw_kernel.bindat 


T_RANGE_L=50
T_RANGE_H=100    # TEMP RANGE AS (T_LOW -T_RANGE_L, T_HIGH + T_RANGE_H) WHERE T_LOW AND T_HIGH ARE THE INITIAL TEMP (TEBEG) AND FINAL (TEND) TEMP OR VICE VERSA

M_RANGE=1       # MAG RANGE

N_RUNS=100
delta=10



# Checking which is the last directory and create the new working directory run_j

j=0 #j is the counter of folders

for ((i=1; i <= $N_RUNS; i++)); do
DIR=run_$i
if [ -d "$DIR" ]; then
    echo "$DIR  exists already."
else
    echo "$DIR  does NOT exist. Creating $DIR, copying files, and going inside it"
    if [ $i == 1 ]; then 
        mkdir run_1
        cp setup/* run_1
        cd run_1
        
        #Extract TEMPERATURE and N ATOMS info from INCAR and POSCAR 
        TEBEG=$(grep TEBEG INCAR | awk '{print $3}'  )
        TEND=$(grep TEBEG INCAR | awk '{print $7}'   )
        NH=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $1}')
        NO=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $2}' )
        NFe=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $3}' )          
        
        NTOT=$(($NH + $NO + $NFe))
        
        echo "TBEG= $TEBEG"
        echo "TEND= $TEND"
        echo "H=$NH, O= $NO, Fe=$NFe"

        j=$(($i+1))
        break

    else
        mkdir run_$i
        
        #going into previous folder
        cd run_$(($i-1)) 
        
        #Extract TEMP and N atom from previous POSCAR and INCAR
        TEBEG=$(grep TEBEG INCAR | awk '{print $3}'  )
        TEND=$(grep TEBEG INCAR | awk '{print $7}'   )
        
        
        NH=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $1}')
        NO=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $2}' )
        NFe=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $3}' )          
        
        NTOT=$(($NH + $NO + $NFe))
        
        echo "TBEG= $TEBEG"
        echo "TEND= $TEND"
        echo "H=$NH, O= $NO, Fe=$NFe"

        #CONTCAR TO POSCAR (Remove Velocity and Correction Prediction DATA from CONCTCAR)
        head -n $((8 + $NH + $NO + $NFe )) CONTCAR > temp 
        mv temp CONTCAR
        
        cp CONTCAR ../run_$i/POSCAR   #PRODUCTION
        #cp POSCAR ../run_$i/POSCAR     #DEVEL TESTING

        cp ML_ABN ../run_$i/ML_AB
        cp INCAR  POTCAR KPOINTS vdw_kernel.bindat    ../run_$i
        cd ../run_$i
        j=$(($i+1))
        break
    fi
fi
done


#set limit of temperature (including the T_RANGE)
if [[ $TEBEG < $TEND ]]; then
    TL=$(($TEBEG -$T_RANGE_L))
    TH=$(($TEND + $T_RANGE_H))
else
    TL=$(($TEND - $T_RANGE_L))
    TH=$(($TEBEG + $T_RANGE_H))
fi


echo "TLOW=$TL"
echo "THIGH=$TH"

#Setting MAGMOM range accoding to the number of irons: for the 44Fe and for 38Fe (48 and 56 respectively)
if [[ "$NFe" -eq 38 ]]; then
    MH=$((48 + $M_RANGE))   
    ML=$((48 - $M_RANGE))
    echo " ////SMALL//// "
elif [[ "$NFe" == 44 ]]; then
    MH=$((56 + $M_RANGE))
    ML=$((56 - $M_RANGE))
    echo " ////BIG//// " 
else                        
    echo "UNEXPECTED NUMBER OF IRON"
fi







#Repeating the procedure N_RUNS times
for ((i=$j; i <= $N_RUNS; i++))
do
	mpirun -np 192 /gpfs/data/fs71796/salvrom/vasp-conda_wo_conda_intel/bin/vasp_std	&
    
    #Check block begins
    lambda=true 
    while $lambda
	do  
    #Check if the first magnetization was printed
    if  ( grep mag OSZICAR > /dev/null); then
        #getting mag and temp. ! IMPORTANT awk '{sub(/.$/,"")}1' remove a weird . at the end of the #T in OSZICAR
    	mag=$( grep mag OSZICAR  | tail -1 | awk '{print $17}' )
        temp=$( grep mag OSZICAR  | tail -1 | awk '{print $3}' |awk '{sub(/.$/,"")}1' )
        echo "TEMPERATURE = $temp"
        echo "MAGNETIZATION =  $mag" 
        #MAGMOM CHECK
    	if  [[ "$mag" > $MH    ||  "$mag" < $ML ]]; then
            killall -9 vasp_std 
            lambda=false 
            echo "////MAG OUT OF RANGE////"
        #TEMPERATURE CHECK
        elif [[ "$temp" -le $TL    ||  "$temp" -ge $TH  ]]  ; then
            killall -9 vasp_std
			lambda=false
            echo "////TEMP OUT OF RANGE////" 
    	else
            echo "////ALLES IN ORDNUNG////" 
        fi 
   	    sleep $delta
    fi
	done
    #Check block ends
    
    mkdir ../run_$i
    
    #Remove Velocity and Correction/Prediction DATA from CONCTCAR (it leaves only Cartesian coordinates)
    head -n $((8 + $NTOT)) CONTCAR > temp_CONT 
    mv temp_CONT CONTCAR

    
    #Removing last ML config from ML_ABN and copy as new ML_AB for the next run
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
    
    
    #Copy
    cp INCAR POTCAR KPOINTS vdw_kernel.bindat    ../run_$i
    
    #DEVEL OPTION OR PRODUCTION COMMENT ONE OF NEXT LINES
    #cp CONTCAR ../run_$i/POSCAR   #PRODUCTION
    cp POSCAR ../run_$i/POSCAR     #DEVEL
    
    cp ML_ABN ../run_$i/ML_AB

    #get 5th final temperature if exists or 4th or 3th or 2nd or 1st
    #T=$(grep mag OSZICAR | tail -n5 | head -1   | awk '{print $3}' | awk '{sub(/.$/,"")}1')
    
    #echo "CURRENT TEMPERATURE IS $T"
    #CHECK IF IT IS OUTSIDE THE TEMP RANGE (FLEXIBLE) AND IF SO RESTART FROM TEBEG OF THE PREVIOUS RUN
    #if  [[ "$T" -le $(( $TL -50 ))    ||  "$T" -ge $(($TH + $T_RANGE)) ]]; then
    #    T=$TEBEG    
    #    echo "T OUTSIDE RANGE, RESTARTING FROM TEBEG"
    #else
    #    echo "T INSIDE RANGE, TEBEG IS CURRENT T"   
    #fi 
    #
    
    cd ../run_$i

    #changing temperature and update
    #sed -i "s/ $TEBEG / $T /" INCAR
    #TEBEG=$T
done













