set -e
cd ~/dc_workshop/results

genome=~/dc_workshop/data/ref_genome/chr20.fa

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in ~/dc_workshop/data/trimmed_fastq/*_1.trim.fq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.fq.gz)
    echo "base name is $base"

    fq1=~/dc_workshop/data/trimmed_fastq_small/${base}_1.trim.fq.gz
    fq2=~/dc_workshop/data/trimmed_fastq_small/${base}_2.trim.fq.gz
    sam=~/dc_workshop/results/sam/${base}.aligned.sam
    bam=~/dc_workshop/results/bam/${base}.aligned.bam
    sorted_bam=~/dc_workshop/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/dc_workshop/results/bcf/${base}_raw.bcf
    variants=~/dc_workshop/results/bcf/${base}_variants.vcf

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    freebayes -f $genome $bam > $variants

    done
