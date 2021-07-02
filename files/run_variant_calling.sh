set -e
cd ~/dc_workshop/results

genome=~/dc_workshop/data/ref_genome/chr20.fa
## Load the required modules
module load BWA
module load freebayes
module load SAMtools

bwa index $genome

mkdir -p sam bam vcf

for fq1 in ~/dc_workshop/data/trimmed_fastq/*_R1.trim.fq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _R1.trim.fq.gz)
    echo "base name is $base"

    fq1=~/dc_workshop/data/trimmed_fastq/${base}_R1.trim.fq.gz
    fq2=~/dc_workshop/data/trimmed_fastq/${base}_R2.trim.fq.gz
    sam=~/dc_workshop/results/sam/${base}.aligned.sam
    bam=~/dc_workshop/results/bam/${base}.aligned.bam
    sorted_bam=~/dc_workshop/results/bam/${base}.aligned.sorted.bam
    variants=~/dc_workshop/results/vcf/${base}_chr20.vcf

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    freebayes -f $genome $sorted_bam > $variants

    done
