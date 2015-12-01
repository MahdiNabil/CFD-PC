#!/bin/bash
# Continuous integration + validation code for interThermalPhaseChangeFoam

#Initialize a summary log file for all the operations
rm -rf CFD-PC
rm Summary.log
rm errors.log
touch Summary.log
SummaryFile=`readlink -f Summary.log`

touch errors.log

# First check, clone the case
git clone https://github.com/MahdiNabil/CFD-PC.git 2> errors.log
#Check if successful
FATAL_ERROR=`grep fatal errors.log`
if [ -n "$FATAL_ERROR" ]
then
  MSG="Clone repository check: FAIL"
  echo $MSG
  echo $MSG >> $SummaryFile
  exit
else
  MSG="Clone repository check: PASS"
  echo $MSG
  echo $MSG >> $SummaryFile
fi


# Next, do the build check
cd CFD-PC/interThermalPhaseFoam
touch errors.log
export WM_NCOMPPROCS=12
./Allwmake.sh 2> errors.log 1> compile.log

#Check if successful
FATAL_ERROR=`grep Error errors.log`
if [ -n "$FATAL_ERROR" ]
then
  MSG="Build master branch check: FAIL"
  echo $MSG
  echo $MSG >> $SummaryFile
  exit
else
  MSG="Build master branch check: PASS"
  echo $MSG
  echo $MSG >> $SummaryFile
fi

cat $SummaryFile | mail -s "test" nabil.mehdi1@gmail.com
