#This script counts the number of reads per sample for MiSeq, AVITI, and NovaSeq sequencing technologies

myseq=("MiSeq" "AVITI" "NovaSeq")

for seq in ${myseq[@]}
do 
	for i in /path/to/sequencing/data/$seq/*.gz
	do
		name=$(basename $i)
		reads=$(zcat ${i}|wc -l)
		reads=$(($reads/4))
		echo $name $reads
	done > /path/to/results/${seq}_reads.txt
done