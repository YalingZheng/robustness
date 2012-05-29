rm *.out
rm *.out.scwrl4out
rm *.dna
rm *.dna.cmd
rm *.cmd.out
rm *.cmd.err
rm *.cmd.log

if [ -d results ]
then
    rm -r results
fi
if [ -d cmderrlog ]
then 
    rm -r cmderrlog
fi
if [ -d failedjobs ]
then
    rm -r failedjobs
fi
# compile the c source files
g++ red_black_tree.c stack.c misc.c neutralrobustness.cpp userintf.cpp stoc1.cpp mersenne.cpp -o neutralrobustness
# submit a job 
condor_submit runpop1000gen500.cmd 
# we need to wait until the job is finished. 
# When the job is finished, when we repeatedly check
# condor_q yzheng (Here I will set it to every 1 minute)

# we type condor_q, and we read the first line of last line of the result of this command, if it is 0, then okay, we go to the next step
# wait for the job of list all the dnas during the generations to finish
stopflag=0
while [ $stopflag = 0 ]
do
    sleep 60
    jobstotal=`condor_q yzheng | grep jobs | grep idle | grep running | grep held | cut -d" " -f1`
    echo $jobstotal jobs still ...
    if [ $jobstotal = 0 ]
    then
	stopflag=1
    fi
done

# Lots of *.dna, *.dha.cmd files generated, we handle all files end with .cmd
ls *.dna.cmd | cut -d" " -f1 | while read LINE
do
    currentfile=$LINE
    condor_submit $currentfile
done

# wait for all the jobs of computing the robustnesses of dnas to finish
stopflag=0
while [ $stopflag = 0 ]
do
    sleep 60
    jobstotal=`condor_q yzheng | grep jobs | grep idle | grep running | grep held | cut -d" " -f1`
    echo $jobstotal jobs still ...
    if [ $jobstotal = 0 ]
    then
	stopflag=1
    fi
done


# the result are kept in *.dna.cmd.out
if [ -d results ]
then
    echo directory results already exists ...
else
    mkdir results
fi

mv *.dna.cmd.out results
mv *.dna results

if [ -d cmderrlog ]
then 
    echo directory cmderrlog already exists ... # do nothing
else
    mkdir cmderrlog
fi

mv *.dna.cmd* cmderrlog
mv *.cmd.err cmderrlog
mv *.cmd.log cmderrlog
mv *.cmd.out cmderrlog


# Note that there may be some errrs in the form of *.out, *.out.scwrl4out. Those *.out keep those dnas whose robustness was not successfully computed
if [ -d failedjobs ]
then
    echo directory failedjobs already exists ... # do nothing
else
    mkdir failedjobs
fi


mv *.out failedjobs

mv *.out.scwrl4out failedjobs

