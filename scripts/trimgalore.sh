mkdir -p trimmed_fastq

echo "Run TrimGalore..."

for fastq in *.fastq.gz
do trim_galore --quality 25 --stringency 5 --length 35 --fastqc $fastq
done

echo "trimming is done..."

mv *fq.gz trimmed_fastq
