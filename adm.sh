#!/bin/bash
#PBS -q parallel12
#PBS -l select=1:ncpus=12:mem=45GB
#PBS -l walltime=720:00:00
#PBS -j oe

cd /home/svu/dbstq/ADM
for i in *_adm.ped
 do 
 for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  do 
    /home/svu/dbstq/ADM/admixture -j12 --cv $i $K | tee log${K}.out 
  done
 grep -h CV log*.out > ${i}_cv.txt
done
  