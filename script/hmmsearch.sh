#!/bin/bash
echo hmmsearch for different version of Pfam

declare -a Version=("Pfam25.0" "Pfam29.0" "Pfam35.0")

for group in "${Version[@]}"; do
	wait
	cmd0=`mkdir ../autodata/align/${group}/`
	echo $cmd0
	wait

	FILE=../autodata/align/${group}/Pfam-A.hmm.gz
	if [ -f "$FILE" ]; then
		echo "$FILE exists."
	else 
	        echo "$FILE does not exist."
		cmd1=`wget -P ../autodata/align/${group}/ http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/${group}/Pfam-A.hmm.gz`
		echo $cmd1
	fi
	wait
	cmd2=`gzip -d ../autodata/align/${group}/Pfam-A.hmm.gz`
	echo $cmd2

done
for group in "${Version[@]}"; do 
	cmd3=`hmmpress ../autodata/align/${group}/Pfam-A.hmm`
	#-T 15 --domT 15 
	cmd4=`hmmsearch --domtblout ../autodata/align/${group}uniprot_2_1_1_domout.tsv --tblout ../autodata/align/${group}uniprot_2_1_1_tbout.tsv --cpu 4 ../autodata/align/${group}/Pfam-A.hmm ../autodata/rawdata/uniprot_ec2.1.1.fasta`
	wait
	echo $cmd3
	wait
	echo $cmd4
done