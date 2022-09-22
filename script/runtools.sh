#!/bin/bash
declare -a Domain=("PF05175" "PF08241" "PF08242" "PF13489" "PF13649" "PF13847")
for group in "${Domain[@]}"; do
        cmd1=`hmmalign --amino ../autodata/align/${group}.hmm ../autodata/sequences/${group}_rm.fasta > ../autodata/align/${group}_hmmalign_out`
        cmd2=`hmmbuild ../autodata/align/${group}_hmmalign_out.hmm ../autodata/align/${group}_hmmalign_out`
        echo "start"
        wait
        echo $cmd1
        wait
        echo $cmd2
done

#declare -a StringArray=("S" "O" "C" "N" )
#for group in "${StringArray[@]}"; do
#	cmd1=`mmseqs createindex `$group_seq_db tmp`
#	echo $cmd1
#done
