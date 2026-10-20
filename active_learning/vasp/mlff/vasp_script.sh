#!/bin/bash                                                                     
#SBATCH -J SF                                           
#SBATCH -N 3                                                                    
#SBATCH --partition=zen3_0512                                                  
#SBATCH --qos=zen3_0512
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-core=1                                                     

module purge

module load intel-oneapi-compilers/2022.0.2-gcc-11.2.0-yzi4tsu intel-oneapi-mpi\
/2021.4.0-gcc-11.2.0-2e7zm7z intel-oneapi-mkl/2022.0.2-gcc-11.2.0-2ugv2zs
module load hdf5/1.12.1-gcc-11.2.0-523hswj

export I_MPI_PIN_RESPECT_CPUSET=0
export OMP_NUM_THREADS=1

#SET TIME WAITING

delta=10

#setup has to be a folder next to wher this script runs, containing the intial files for VASP:
#INCAR  KPOINTS  ML_ABN  POSCAR  POTCAR  vdw_kernel.bindat 

mkdir run_1
cp setup/* run_1
cd run_1

#Extract usefull info from INCAR and POSCAR 
TEBEG=$(grep TEBEG INCAR | awk '{print $3}'  )
TEND=$(grep TEBEG INCAR | awk '{print $7}'   )

NH=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $1}')
NO=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $2}' )
NFe=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $3}' )          


echo "TBEG= $TEBEG"
echo "TEND= $TEND"

#Setting MAGMOM range accoding to the number of irons of the big 44 Fe or small 38 Fe systems
if [[ "$NFe" -eq 38 ]]; then
    M=49   
    L=47
    echo " ////SMALL//// "
elif [[ "$NFe" == 44 ]]; then
    M=57
    L=55
    echo " ////BIG//// " 
else                        
    echo "UNEXPECTED NUMBER OF IRON"
fi


#This for loop is just for repeating the procedure N times
for i in {2..100}
do
	mpirun -np 192 /gpfs/data/fs71796/salvrom/vasp-conda_wo_conda_intel/bin/vasp_std	&
    
    #Check block begins
    lambda=true 
    while $lambda
	do  
        #getting mag and temp every delta. !!!!!!! IMPORTANT awk '{sub(/.$/,"")}1' remove a weird . at the end of the T in OSZICAR
    	mag=$( grep mag OSZICAR  | tail -1 | awk '{print $17}' )
        temp=$( grep mag OSZICAR  | tail -1 | awk '{print $3}' |awk '{sub(/.$/,"")}1' )
        echo "temp= $temp"
        echo "mag=  $mag" 
        #MAGMOM 
    	if  ( grep mag OSZICAR > /dev/null && [[ "$mag" > $M    ||  "$mag" < $L ]]); then
            killall -9 vasp_std 
            lambda=false 
               echo "////MAG OUT OF RANGE////"
        #TEMPERATURE CHECK
        elif  ( grep mag OSZICAR > /dev/null  &&  [[ "$temp" -le $(( $TEBEG -10 ))    ||  "$temp" -ge $(($TEND + 50)) ]] )  ; then
            killall -9 vasp_std
			lambda=false
            echo "////TEMP OUT OF RANGE////" 
    	else
            echo "////ALLES IN ORDUNG////" 
        fi 
   	    sleep $delta
	done
    #Check block ends
    
    mkdir ../run_$i
    
    #Remove Velocity and Correction Prediction DATA from CONCTCAR (it leaves only Cartesian Coordinates). 
    head -n $((8 + $NH + $NO + $NFe )) CONTCAR > temp 
    temp > CONTCAR
    rm temp
    #Copy
    cp INCAR  POTCAR KPOINTS vdw_kernel.bindat    ../run_$i
     
    #DEVEL OPTION OR PRODUCTION COMMENT ONE OF NEXT LINES
    cp CONTCAR ../run_$i/POSCAR   #PRODUCTION
    #cp POSCAR ../run_$i/POSCAR     #DEVELL

    cp ML_ABN ../run_$i/ML_AB

    #get 5th final temperature if exists or 4th or 3th or 2nd or 1st
    T=$(grep mag OSZICAR | tail -n5 | head -1   | awk '{print $3}' | awk '{sub(/.$/,"")}1')
    
    cd ../run_$i

    #changing temperature and update
    sed -i "s/$TEBEG/$T/" INCAR
    TEBEG=$T
done













