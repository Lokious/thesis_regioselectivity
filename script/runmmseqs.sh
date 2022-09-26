#declare -a StringArray=("S" "O" "C" "N" )
#for group in "${StringArray[@]}"; do
#	cmd1=`mmseqs createindex `$group_seq_db tmp`
#	echo $cmd1
#done

mmseqs easy-cluster ../autodata/sequences/uniprot_ec2.1.1.fasta uniprot_ec2.1.1.fasta uniprot——2.1.1cluster_result tmp
