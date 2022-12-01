#!/bin/bash

echo "running batch scripts to creat job files from template.job"
cd /lfs/home/ychen/scripts/R/Rscripts/LAI_STUDY_EA/
# LOAD R MODULE 
#module load r/3.1.1

#echo "removing previous jobs!"
#rm -f R_* 
# /home/orchidee03/jryder/cru_regridded2 
# submit the job on que 
# chmod u+x $CONFIG_FILE
# qsub -jeo -q short $CONFIG_FILE
CONT=0

# argunment 1 to 6 are parameters for the function lai_change_analysis.R
#fun.dyn.lai <- function ( Ref.wind=3.0, Ref.rain=100, Win.size=40, offdays=30, pdays=60)
 
# arg1 : EXPERINMENTS 
# arg2 : start year & end year
# arg3 : start day & end day 
# arg4 :

#declare -A arg_arr
#case 00
#arg_arr=()
#arg_arr+=('12 30 50 0 50 TRACK_DATA_2D')
#arg_arr+=('12 30 50 0 50 TRACK_DATA_2D')




#for iwind in 12 14 16  

#do 
#arg1="5"   # wind speed
#arg1="${iwind}"

#for irainf in 20 
#do 

#arg2="80" #rainf 
#arg2="${irainf}"  
#arg3="50" #window size 
#arg4="0" #minimum offset 

# post days
#for i in 40 50 70 90 120 150 180    
#for i in 80 100 120 150 180      
#for i in  110 130 150 170 190 200   
#for i in  90 100 120 140 160 180 210   
#for ipost in   50    
#for i in 20 30 40 50 60  
#----------------
#do 
#
#arg5="50"   #"${ipost}"

#arg6="TRACK_DATA_2D"
#

for irun in {5..5}

do
#combine based
#case 1 
if [ "${irun}" == "1" ]; then 
   arg1="8"; arg2="60"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_2D"
fi
#case 2
if [ "${irun}" == "2" ]; then 
   arg1="10"; arg2="80"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_3D"
fi
#case 3
if [ "${irun}" == "3" ]; then 
   arg1="12"; arg2="100"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_4D"
fi
#wind based
#case 4 
if [ "${irun}" == "4" ]; then 
   arg1="8"; arg2="0"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_2D"
fi
#case 5
if [ "${irun}" == "5" ]; then 
   arg1="10"; arg2="0"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_3D"
fi
#case 6
if [ "${irun}" == "6" ]; then 
   arg1="12"; arg2="0"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_4D"
fi

#rainfall based
#case 7
if [ "${irun}" == "7" ]; then 
   arg1="0"; arg2="60"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_2D"
fi
#case 8
if [ "${irun}" == "8" ]; then 
   arg1="0"; arg2="80"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_3D"
fi
#case 9
if [ "${irun}" == "9" ]; then 
   arg1="0"; arg2="100"; arg3="60"; arg4="0"; arg5="60"; arg6="TRACK_DATA_4D"
fi


#counter for the iteration or loop 
CONT=$(($CONT+1))
echo "The No. of for loop: $CONT ."


#copy the template from template.job (script for submitting R job with some argunments)
TEMPLATE_FILE="template.submit"
CONFIG_FILE="R_${arg1}_${arg2}_${arg3}_${arg4}_${arg5}_${arg6}.job"
cp $TEMPLATE_FILE $CONFIG_FILE

#dynamic text for the argunments, input and output files and directory
log_file="w_${arg1}_p_${arg2}_window_${arg3}_offset_${arg4}_pos_${arg5}_track_${arg6}_tc.occ.spei_02_qc.log"
#
TARGET_KEY="work_dir"
REPLACEMENT_VALUE="\/lfs\/home\/ychen\/scripts\/R\/Rscripts\/LAI_STUDY_EA\/"
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

TARGET_KEY="R_filename"
REPLACEMENT_VALUE="lai_change_analysis.R"
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

# RE-SET THE PARAMETER VALUES FOR EACH JOB
TARGET_KEY="arg1"
REPLACEMENT_VALUE=${arg1}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

TARGET_KEY="arg2"
REPLACEMENT_VALUE=${arg2}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

TARGET_KEY="arg3"
REPLACEMENT_VALUE=${arg3}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

TARGET_KEY="arg4"
REPLACEMENT_VALUE=${arg4}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

TARGET_KEY="arg5"
REPLACEMENT_VALUE=${arg5}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE

TARGET_KEY="arg6"
REPLACEMENT_VALUE=${arg6}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE


TARGET_KEY="log_file"
REPLACEMENT_VALUE=${log_file}
sed -i -c "s/$TARGET_KEY/$REPLACEMENT_VALUE/" $CONFIG_FILE


echo " arg1:$arg1, arg2:$arg2, arg3:$arg3, arg4:$arg4 , arg5:$arg5, arg6:$arg6 "  
# submit the job on queuing system  
chmod u+x $CONFIG_FILE

#run on curie by tesing Queue
#ccc_msub -q standard -T 1800 -A gen6328 -Q TEST $CONFIG_FILE 

#normal run on CURIE
#ccc_msub -q standard -T 86400 -A gen6328 $CONFIG_FILE 
#run on airain
#ccc_msub -q ivybridge -T 86400 -A dsm $CONFIG_FILE 
# excute the job 
./$CONFIG_FILE &
echo $CONFIG_FILE

#done
#done


done 
