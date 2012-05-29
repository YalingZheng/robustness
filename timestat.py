import fileinput

generation=0
robustness=False
count=0
stablesize=0
exectime=0.0
lastexectime=0.0
for line in fileinput.input():
    line=line.strip()
    if line == "":
        continue

    if line.find("Generation #") != -1:
        generation=generation+1
        if lastexectime != 0.0:
            exectime=float((line.split(", ")[1]).split(" = ")[1][:-1])
            elapsed=exectime-lastexectime
            print "~%3.2f energy calculations per second." % (((1000)+(float(stablesize)*630))/elapsed)
            lastexectime=exectime
        else:
            lastexectime=float((line.split(", ")[1]).split(" = ")[1][:-1])
        print line
        if robustness:
            robustness=False
            count=0
            stablesize=0
    elif line.find("stablesize") != -1:
        print line
        stablesize = int((line.split(",")[0]).split(" = ")[1])
    elif line.find("Robustness") != -1:
        robustness=True
    elif robustness:
        count=count+1


if stablesize != 0:
    print "%3.2f%%" % (100*((float(count)/2)/float(stablesize)))
else:
    print "0.0%"

