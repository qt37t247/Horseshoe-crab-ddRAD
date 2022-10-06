#!/bin/bash
#$ -cwd
#PBS -l select=1:ncpus=12:mpiprocs=12:mem=45GB
#PBS -l walltime=720:00:00
#PBS -j oe
#PBS -q parallel12
#chmod +x ./fsc2702

cd /home/svu/dbstq/SCENARIONAME/bestrun/SCENARIONAME_boot/
for i in {1..100}
do
cd SCENARIONAME_boot_$i
cp ../../fsc2702 fsc2702
cp ../../SCENARIONAME.tpl SCENARIONAME_boot.tpl
cp ../../SCENARIONAME.est SCENARIONAME_boot.est
./fsc2702 -c 12 -t SCENARIONAME_boot.tpl -n 500000 -m -e SCENARIONAME_boot.est -M -L 50 -q --foldedSFS
cd ..
done
