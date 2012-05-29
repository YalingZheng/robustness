all: 
	g++ red_black_tree.c stack.c misc.c neutralrobustness.cpp userintf.cpp stoc1.cpp mersenne.cpp -o neutralrobustness 
clean:
	rm -rf mutatedprotein nurefined.pdb refined.pdb hot.grp logfile logfile1 aa2 *.out *.out.pdb *.dSYM *.scwrl4out *.o *.cmd *.cmd.* cmd
distclean:
	rm -rf neutralrobustness mutatedprotein nurefined.pdb refined.pdb hot.grp logfile logfile1 aa2 *.out *.out.pdb *.scwrl4out *.dSYM *.o *.cmd *.cmd.* cmd
