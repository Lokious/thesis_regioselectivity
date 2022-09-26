#!/bin/bash
echo Please input the domains: ("domain1" "domian2"...)
read domains
echo It\'s nice to meet you $varname

declare -a Domain=domains
:'
for group in "${Domain[@]}"; do
	wait
	cmd=`hmmsearch -T 15 --domT 15 --domtblout ../autodata/align/${group}_domout.tsv --tblout ../autodata/align/${group}_tbout.tsv --cpu 4 ../autodata/align/${group}.hmm ../autodata/rawdata/uniprot_ec2.1.1.fasta`
	echo $cmd
done
echo `python3 hmm.py`
for group in "${Domain[@]}"; do
	cmd1=`hmmalign --amino --trim ../autodata/align/${group}.hmm ../autodata/sequences/${group}.fasta > ../autodata/align/${group}_hmmalign_out_trim`
	cmd2=`hmmbuild ../autodata/align/${group}_hmmalign_out_trim.hmm ../autodata/align/${group}_hmmalign_out_trim`
	cmd3=`hmmsearch ../autodata/align/${group}_hmmalign_out_trim.hmm pdbaa > pdb_${group}_trim.hmmer`
	echo "start"
	wait
	echo $cmd1
	wait
	echo $cmd2
	wait
	echo $cmd3
done
 '
