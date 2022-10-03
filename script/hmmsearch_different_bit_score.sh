#!/bin/bash
echo hmmsearch for different version of Pfam

declare -a BitScore=("3" "5" "7" "9" "11" "13" "15" "17" "19" "21")
FILE=../autodata/align/different_version_pfam/Pfam35.0/Pfam-A.hmm.gz
if [ -f "$FILE" ]; then
    echo "$FILE exists."
else
    echo "$FILE does not exist."
cmd1=`wget -P ../autodata/align/different_version_pfam/Pfam35.0/ http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/${group}/Pfam-A.hmm.gz`
echo $cmd1
fi
wait

cmd2=`gzip -d ../autodata/align/different_version_pfam/Pfam35.0/Pfam-A.hmm.gz`
echo $cmd2
cmd3=`hmmpress ../autodata/align/different_version_pfam/Pfam35.0/Pfam-A.hmm`
wait
echo $cmd3

for group in "${BitScore[@]}"; do
	wait
	cmd0=`mkdir ../autodata/align/different_version_pfam/Pfam35.0/Bit_Score_${group}/`
	echo $cmd0
	wait
done

for group in "${BitScore[@]}"; do 

	#-T 15 --domT 15 
	cmd4=`hmmsearch -T ${group} --domT ${group} --domtblout ../autodata/align/different_version_pfam/Pfam35.0/Bit_Score_${group}/uniprot_2_1_1_domout.tsv --tblout ../autodata/align/different_version_pfam/Pfam35.0/Bit_Score_${group}/uniprot_2_1_1_tbout.tsv --cpu 4 ../autodata/align/different_version_pfam/Pfam35.0/Pfam-A.hmm ../autodata/rawdata/uniprot_ec2.1.1.fasta`
	echo "bitscore ${group}" 
	wait
	echo $cmd4
done
