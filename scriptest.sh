# clean the previous output files
function removematchedfiles {
    result=`ls $1 | grep "No such file or directory"`
    echo result
    if [ -z $result ]
    then
	ls $1 | cut -d" " -f1 | while read LINE
	do	
	    onefile=$LINE
	    if [ -e $onefile ]
	    then    
		echo rm $onefile
		rm $onefile
	    fi
	done
    fi   
}
removematchedfiles "*.out.scwrl4out"
removematchedfiles "*.dna"
removematchedfiles "*.dna.cmd*"
if [ -d results ]
then
    rm -r results
fi
mkdir results
# we need to wait until the job is finished. 
# When the job is finished, when we repeatedly check
# condor_q yzheng (Here I will set it to every 20 minutes)

# we type condor_q, and we read the first line of last line of the result of this command, if it is 0, then okay, we go to the next step

# Lots of *.dna, *.dha.cmd files generated, we handle all files end with .cmd
ls *.dna.cmd | cut -d" " -f1 | while read LINE
do
    currentfile=$LINE
    condor_submit $currentfile
done

# the result are kept in *.dna.cmd.out
mkdir results
ls *.dna.cmd.out | cut -d" " -f1 | while read LINE
do
    thefile=$LINE
    mv $thefile results
done
if [ -d cmderrlog ]
then 
    echo # do nothing
else
    mkdir cmderrlog
fi
ls *.dna.cmd* | cut -d" " -f1 | while read LINE
do
    currfile=$LINE
    mv $currfile cmderrlog
done


# Note that there may be some errrs in the form of *.out, *.out.scwrl4out. Those *.out keep those dnas whose robustness was not successfully computed
if [ -d failedjobs ]
then
    echo # do nothing
else
    mkdir failedjobs
fi
ls -l *.dna.cmd* | cut -d" " -f12 | while read LINE
do
    currfile=$LINE
    mv $currfile cmderrlog
done

if [ -e *.out ]
then
    mv *.out failedjobs
fi
ls -l *.out | cut -d" " -f12 | while read LINE
do
    currfile=$LINE
    mv $currfile failedjobs
done
ls -l *.out.scwrl4out | cut -d" " -f12 | while read LINE
do
    currfile=$LINE
    mv $currfile failedjobs
done


