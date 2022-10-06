#!/bin/bash
#$ -cwd
#PBS -l select=1:ncpus=12:mem=45GB
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -q parallel_test

cd /home/svu/dbstq/SCENARIONAME/bestrun
./fsc2702 -c 12 -i SCENARIONAME_boot.par -n100 -j -m -s0 -x -I -q --foldedSFS

