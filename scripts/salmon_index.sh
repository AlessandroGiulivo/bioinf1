# download transcriptome sequences
wget http://ftp.ensembl.org/pub/release-104/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# create salmon index
salmon index -t Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i yeast_index
