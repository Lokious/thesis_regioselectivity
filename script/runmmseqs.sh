#declare -a StringArray=("S" "O" "C" "N" )
#for group in "${StringArray[@]}"; do
#	cmd1=`mmseqs createindex `$group_seq_db tmp`
#	echo $cmd1
#done

mmseqs easy-cluster ../autodata/sequences/${group}.fasta ${group}cluster_result tmp
