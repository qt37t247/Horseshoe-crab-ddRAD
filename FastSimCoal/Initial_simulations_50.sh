#!/bin/bash
cd $PBS_O_WORKDIR;
fsc=fsc2702
jobcount=0
msgs=conOutputs_SCENARIONAME

#-------- Number of different runs per data set ------
numRuns=50
runBase=1 
#-----------------------------

mkdir $msgs 2>/dev/null

#-------- Default run values ------
numSims=500000                    #-n command line option
numCycles=40                      #-L command line option
minValidSFSEntry=1                #-C command line option

#-------- Ascertainment ------
withAscertainment=0
ascPop=0                          #-a command line option
ascSize=2                         #-A command line option
#-----------------------------
useMonoSites="" 
#useMonoSites="-0"                #-0 command line option
#----------multiSF------------
multiSFS=""
#multiSFS="--multiSFS"            #--multiSFS command line option
#-----------------------------

#-------- Generic Name ------
genericName=SCENARIONAME
tplGenericName=SCENARIONAME
estGenericName=SCENARIONAME
#-----------------------------

for dirs in $genericName
do
	#Check that dirs is a directory
	if [ -d "$dirs" ];	then
		cd $dirs
			echo "Main directory : $dirs"
			estFile=$estGenericName.est
			tplFile=$tplGenericName.tpl     
			for (( runsDone=$runBase; runsDone<=$numRuns; runsDone++ ))
			do
				runDir="run$runsDone"
				mkdir $runDir 2>/dev/null
				echo "--------------------------------------------------------------------"
				echo ""
				echo "Currrent file: $subDirs $runDir"
				echo ""
				cd $runDir
				#Copying necessary files
				cp ../$fsc .
				cp ../$tplFile .
				cp ../$estFile .
				cp ../*.obs .
				#Renaming files for consistency
				mv $tplFile ${genericName}.tpl 2>/dev/null
				mv $estFile ${genericName}.est 2>/dev/null

				let jobcount=jobcount+1
				jobName=${genericName}${jobcount}.sh

				#Creating bash file on the fly
				(
				echo "#!/bin/bash"
				echo ""
				echo "#$ -cwd"
				echo ""
				echo "#PBS -l select=1:ncpus=12:mem=45GB"
				echo "#PBS -l walltime=500:00:00"
				echo "#PBS -j oe"
				echo ""
				echo "#$ -N j1P_${jobcount}"
				echo "#$ -o ../../../$msgs/1P_$runsDone.out"
				echo "#$ -e ../../../$msgs/1P_$runsDone.err"
				echo "#$ -m a"
				echo "#PBS -q parallel12"
				echo ""
				echo "#chmod +x ./$fsc"
				echo ""
				echo "echo \"Analysis of file $files\""
				echo "#Computing likelihood of the parameters using the ECM-Brent algorithm"
				echo "echo \"\""
				echo "cd $PBS_O_WORKDIR/$genericName/$runDir"
				if [ $withAscertainment -eq 1 ] ; then									  
					echo "./$fsc -c 12 -t ${genericName}.tpl -n $numSims -m -e ${genericName}.est -M $stopCrit -L $numCycles -a${ascPop} -A${ascSize} -q ${useMonoSites} ${multiSFS}"
				else
					echo "./$fsc -c 12 -t ${genericName}.tpl -n $numSims -m -e ${genericName}.est -M $stopCrit -L $numCycles -q ${useMonoSites} ${multiSFS} --foldedSFS"
				fi
				echo ""
				echo "echo \"\""
				echo "echo \"Job $jobcount terminated\""
				) > $jobName
				chmod +x $jobName

				echo "Bash file $jobName created"
				qsub ./${jobName}
				#./${jobName}
				cd .. #$runDir
			done
		cd .. #dirs
	fi
done



