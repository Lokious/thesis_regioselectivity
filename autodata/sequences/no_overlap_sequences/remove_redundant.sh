#!/bin/bash
fasta=`ls ./*.fasta`
for entry in $fasta; do
	echo "$entry"
	echo `mmseqs easy-cluster ${entry} ${entry}cluster_result tmp`
done
