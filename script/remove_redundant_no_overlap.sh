#!/bin/bash
fasta=`ls ../autodata/sequences/no_overlap_sequences/*.fasta`
for entry in $fasta; do
	echo `python3 get_no_overlap_sequences.py`
	wait
	echo "$entry"
	echo `mmseqs easy-cluster ${entry} ${entry}cluster_result tmp`
done
