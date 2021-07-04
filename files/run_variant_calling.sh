set -e
cd ~/dc_workshop/results

genome=~/dc_workshop/data/ref_genome/chr20.fa
## Load the required modules
module load BWA
module load freebayes
module load SAMtools
module load VCFtools
module load annovar

bwa index $genome

mkdir -p sam bam vcf vcf_annotated

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
    variants_filtered=~/dc_workshop/results/vcf/${base}_chr20_filtered.vcf 
    
    annovar_input=~/dc_workshop/results/vcf_annotated/${base}_avinput 
    annovar_db=~/dc_workshop/results/vcf_annotated/humandb
    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    
    echo "Running freebayes..."
    
    freebayes -f $genome $sorted_bam > $variants

    vcftools --vcf $variants --minQ 20 --recode --recode-INFO-all --out $variants_filtered
    convert2annovar.pl -format vcf4 $variants_filtered > $annovar_input
    
    ## Need to change directory as annovar will create output in the working directory
    cd ~/dc_workshop/results/vcf_annotated/
    
    table_annovar.pl $annovar_input $annovar_db -buildver hg38 -out ${base}_final -remove -protocol refGene,1000g2015aug_all,cosmic70,dbnsfp30a -operation g,f,f,f -nastring NA -csvout

    done
