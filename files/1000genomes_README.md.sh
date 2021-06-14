wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -O NA12878.chr20.bam
samtools index NA12878.chr20.bam

samtools sort -n -o NA12878.chr20.namesorted.bam NA12878.chr20.bam

bedtools bamtofastq -i NA12878.chr20.namesorted.bam -fq NA12878_R1.fq -fq2 NA12878_R2.fq

mkdir -p sub_reads

head NA12878_R1.fq -n 120000 | gzip > sub_reads/NA12878_R1.fq.gz
head NA12878_R2.fq -n 120000 | gzip > sub_reads/NA12878_R2.fq.gz

wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/NA12874/alignment/NA12874.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam -O NA12874.chr20.bam
samtools index NA12874.chr20.bam

samtools sort -n -o NA12874.chr20.namesorted.bam NA12874.chr20.bam

bedtools bamtofastq -i NA12874.chr20.namesorted.bam -fq NA12874_R1.fq -fq2 NA12874_R2.fq

head NA12874_R1.fq -n 120000 | gzip > sub_reads/NA12874_R1.fq.gz
head NA12874_R2.fq -n 120000 | gzip > sub_reads/NA12874_R2.fq.gz


wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/NA12873/alignment/NA12873.chrom20.ILLUMINA.bwa.CEU.low_coverage.20130415.bam -O NA12873.chr20.bam
samtools index NA12873.chr20.bam

samtools sort -n -o NA12873.chr20.namesorted.bam NA12873.chr20.bam

bedtools bamtofastq -i NA12873.chr20.namesorted.bam -fq NA12873_R1.fq -fq2 NA12873_R2.fq

head NA12873_R1.fq -n 120000 | gzip > sub_reads/NA12873_R1.fq.gz
head NA12873_R2.fq -n 120000 | gzip > sub_reads/NA12873_R2.fq.gz
