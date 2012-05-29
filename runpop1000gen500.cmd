####################
##
## Robustness command file
##
####################

universe        = vanilla
+ProjectName = "TG-STA120004S"
executable      = neutralrobustness
output          = neutralrobustnessp1000g500.cmd.out
error           = neutralrobustnessp1000g500.cmd.err
log             = neutralrobustnessp1000g500.cmd.log
arguments       = -p 1000 -g 500 
requirements = (Memory>=1)&&(Arch=="X86_64")
transfer_input_files = dna,Scwrl4,bbDepRotLib.bin,Scwrl4.ini,refined.pdb, temp.pdb
should_transfer_files = IF_NEEDED
when_to_transfer_output = ON_EXIT
queue
