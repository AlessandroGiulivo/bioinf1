echo "Run salmon..."

mkdir -p salmon_results

for filename in *.fastq.gz
do
   base=$(basename $filename .fastq.gz)
   
   echo "Align sample ${base}..."
   
   salmon quant -i yeast_index \
		--libType A \
        	-r ${base}.fastq.gz \
        	-o salmon_results/${base}
done
