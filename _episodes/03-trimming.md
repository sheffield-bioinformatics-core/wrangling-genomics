---
title: "Trimming and Filtering"
teaching: 30
exercises: 25
questions:
- "How can I get rid of sequence data that doesn't meet my quality standards?"
objectives:
- "Clean FASTQ reads using Trimmomatic."
- "Select and set multiple options for command-line bioinformatic tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is an essential step in a genomics workflow."
---
<img src="../img/logo-sm.png" align=right>

# Cleaning Reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We visualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This doesn't mean,
though, that our samples should be thrown out! It's very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to
filter poor quality reads and trim poor quality bases from our samples.

## Trimmomatic Options

Trimmomatic has a variety of options to trim your reads. If we run the following command, we can see some of our options.

~~~
$ trimmomatic
~~~
{: .bash}

Which will give you the following output:
~~~
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version
~~~
{: .output}

This output shows us that we must first specify whether we have paired end (`PE`) or single end (`SE`) reads.
Next, we specify what flag we would like to run. For example, you can specify `threads` to indicate the number of
processors on your computer that you want Trimmomatic to use. In most cases using multiple threads (processors) can help to run the trimming faster. These flags are not necessary, but they can give you more control over the command. The flags are followed by positional arguments, meaning the order in which you specify them is important. 
In paired end mode, Trimmomatic expects the two input files, and then the names of the output files. These files are described below. While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file.

| option    | meaning |
| ------- | ---------- |
|  \<inputFile1>  | Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.|
| \<inputFile2> | Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.|
|  \<outputFile1P> | Output file that contains surviving pairs from the `_1` file. |
|  \<outputFile1U> | Output file that contains orphaned reads from the `_1` file. |
|  \<outputFile2P> | Output file that contains surviving pairs from the `_2` file.|
|  \<outputFile2U> | Output file that contains orphaned reads from the `_2` file.|

The last thing trimmomatic expects to see is the trimming parameters:

| step   | meaning |
| ------- | ---------- |
| `ILLUMINACLIP` | Perform adapter removal. |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |

We will use only a few of these options and trimming steps in our
analysis. It is important to understand the steps you are using to
clean your data. For more information about the Trimmomatic arguments
and options, see [the Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

However, a complete command for Trimmomatic will look something like the command below. This command is an example and will not work, as we do not have the files it refers to:

~~~
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~
{: .bash}

In this example, we've told Trimmomatic:

| code   | meaning |
| ------- | ---------- |
| `PE` | that it will be taking a paired end file as input |
| `-threads 4` | to use four computing threads to run (this will spead up our run) |
| `SRR_1056_1.fastq` | the first input file name |
| `SRR_1056_2.fastq` | the second input file name |
| `SRR_1056_1.trimmed.fastq` | the output file for surviving pairs from the `_1` file |
| `SRR_1056_1un.trimmed.fastq` | the output file for orphaned reads from the `_1` file |
| `SRR_1056_2.trimmed.fastq` | the output file for surviving pairs from the `_2` file |
| `SRR_1056_2un.trimmed.fastq` | the output file for orphaned reads from the `_2` file |
| `ILLUMINACLIP:SRR_adapters.fa`| to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa` |
|`SLIDINGWINDOW:4:20` | to use a sliding window of size 4 that will remove bases if their phred score is below 20 |



> ## Multi-line commands 
> Some of the commands we ran in this lesson are long! When typing a long 
> command into your terminal, you can use the `\` character
> to separate code chunks onto separate lines. This can make your code more readable.
{: .callout}



## Running Trimmomatic

Now we will run Trimmomatic on our data. To begin, navigate to your `untrimmed_fastq` data directory:

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
~~~
{: .bash}

We are going to run Trimmomatic on one of our paired-end samples. 
While using FastQC we saw that Nextera adapters were present in our samples. 
The adapter sequences came with the installation of trimmomatic, so we will first copy these sequences into our current directory.

~~~
$ cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .
~~~
{: .bash}

We will also use a sliding window of size 4 that will remove bases if their
phred score is below 20 (like in our example above). We will also
discard any reads that do not have at least 25 bases remaining after
this trimming step. This command will take a few minutes to run.

~~~
$ trimmomatic PE -phred33 NA12873_R1.fq NA12873_R2.fq.gz \
                NA12873_R1.trim.fq.gz NA12873_R1un.trim.fq.gz \
                NA12873_R2.trim.fq.gz NA12873_R2un.trim.fq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
~~~
{: .bash}


~~~
TrimmomaticPE: Started with arguments:
 -phred33 NA12873_R1.fq NA12873_R2.fq NA12873_R1.trim.fq.gz NA12873_R1un.trim.fq.gz NA12873_R2.trim.fq.gz NA12873_R2un.trim.fq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
Multiple cores found: Using 4 threads
Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 30000 Both Surviving: 26344 (87.81%) Forward Only Surviving: 259 (0.86%) Reverse Only Surviving: 2920 (9.73%) Dropped: 477 (1.59%)
TrimmomaticPE: Completed successfully
~~~
{: .output}

> ## Exercise
>
> Use the output from your Trimmomatic command to answer the
> following questions.
>
> 1) What percent of reads did we discard from our sample?
> 2) What percent of reads did we keep both pairs?
>
>> ## Solution
>> 1) 1.59%
>> 2) 87.81%
> {: .solution}
{: .challenge}

We can confirm that we have our output files:

~~~
$ ls NA12873_R*
~~~
{: .bash}

~~~
NA12873_R1.fq           NA12873_R1_fastqc.zip    NA12873_R2.trim.fq.gz   NA12873_R2un.trim.fq.gz
NA12873_R1.trim.fq.gz   NA12873_R1un.trim.fq.gz  NA12873_R2_fastqc.html
NA12873_R1_fastqc.html  NA12873_R2.fq.gz         NA12873_R2_fastqc.zip
~~~
{: .output}

The output files are also FASTQ files. It should be smaller than our
input file, because we've removed reads. We can confirm this:

~~~
$ ls NA12873_R* -l -h
~~~
{: .bash}

~~~
-rwxrwxrwx 1 root root 6.5M Jun 11 13:57 NA12873_R1.fq
-rw-r--r-- 1 root root 2.4M Jun 14 14:28 NA12873_R1.trim.fq.gz
-rw-r--r-- 1 root root 585K Jun 14 11:05 NA12873_R1_fastqc.html
-rw-r--r-- 1 root root 329K Jun 14 11:05 NA12873_R1_fastqc.zip
-rw-r--r-- 1 root root  21K Jun 14 14:28 NA12873_R1un.trim.fq.gz
-rwxrwxrwx 1 root root 6.5M Jun 11 13:57 NA12873_R2.fq.gz
-rw-r--r-- 1 root root 2.4M Jun 14 14:28 NA12873_R2.trim.fq.gz
-rw-r--r-- 1 root root 660K Jun 14 11:05 NA12873_R2_fastqc.html
-rw-r--r-- 1 root root 354K Jun 14 11:05 NA12873_R2_fastqc.zip
-rw-r--r-- 1 root root 221K Jun 14 14:28 NA12873_R2un.trim.fq.gz
~~~
{: .output}


We've just successfully run Trimmomatic on one of our FASTQ files!
However, there is some bad news. Trimmomatic can only operate on
one sample at a time and we have more than one sample. The good news
is that we can use a `for` loop to iterate through our sample files
quickly! 

We unzipped one of our files before to work with it, let's compress it again before we run our for loop.

~~~
gzip NA12873_R1.fq 
~~~
{: .bash}

~~~
$ for infile in *_R1.fq.gz
> do
>   base=$(basename ${infile} _R1.fq.gz)
>   trimmomatic PE ${infile} ${base}_2.fq.gz \
>                ${base}_1.trim.fq.gz ${base}_1un.trim.fq.gz \
>                ${base}_2.trim.fq.gz ${base}_2un.trim.fq.gz \
>                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
> done
~~~
{: .bash}


Go ahead and run the for loop. It should take a few minutes for
Trimmomatic to run for each of our six input files. Once it's done
running, take a look at your directory contents. You'll notice that even though we ran Trimmomatic on file `NA12973` before running the for loop, there is only one set of files for it. Because we matched the ending `_1.fq.gz`, we re-ran Trimmomatic on this file, overwriting our first results. That's ok, but it's good to be aware that it happened.

~~~
$ ls
~~~
{: .bash}

~~~
NA12873_1.trim.fq.gz     NA12873_R2_fastqc.zip           NA12874_R1.fq.gz        NA12878_R1_fastqc.zip
NA12873_1untrim.fq.gz    NA12873_R2un.trim.fq.gz         NA12874_R1_fastqc.html  NA12878_R2.fq.gz
NA12873_2.trim.fq.gz     NA12873_trimmed_R1.fq           NA12874_R1_fastqc.zip   NA12878_R2_fastqc.html
NA12873_2untrim.fq.gz    NA12873_trimmed_R1_fastqc.html  NA12874_R2.fq.gz        NA12878_R2_fastqc.zip
NA12873_R1.fq.gz         NA12873_trimmed_R1_fastqc.zip   NA12874_R2_fastqc.html  NexteraPE-PE.fa
NA12873_R1.trim.fq.gz    NA12873_trimmed_R2.fq           NA12874_R2_fastqc.zip   SRR622461_1.filt.fastq.gz
NA12873_R1_fastqc.html   NA12873_untrimmed_R1.fq         NA12878_1.trim.fq.gz    multiqc_data
NA12873_R1_fastqc.zip    NA12873_untrimmed_R2.fq         NA12878_1untrim.fq.gz   multiqc_report.html
NA12873_R1un.trim.fq.gz  NA12874_1.trim.fq.gz            NA12878_2.trim.fq.gz    qc_report
NA12873_R2.fq.gz         NA12874_1untrim.fq.gz           NA12878_2untrim.fq.gz
NA12873_R2.trim.fq.gz    NA12874_2.trim.fq.gz            NA12878_R1.fq.gz
NA12873_R2_fastqc.html   NA12874_2untrim.fq.gz           NA12878_R1_fastqc.html

~~~
{: .output}

> ## Exercise
> We trimmed our fastq files with Nextera adapters, 
> but there are other adapters that are commonly used.
> What other adapter files came with Trimmomatic?
>
>
>> ## Solution
>> ~~~
>> $ ls ~/miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/
>> ~~~
>> {: .bash}
>>
>> ~~~
>> NexteraPE-PE.fa  TruSeq2-SE.fa    TruSeq3-PE.fa
>> TruSeq2-PE.fa    TruSeq3-PE-2.fa  TruSeq3-SE.fa
>> ~~~
>> {: .output}
>>
> {: .solution}
{: .challenge}

We've now completed the trimming and filtering steps of our quality
control process! Before we move on, let's move our trimmed FASTQ files
to a new subdirectory within our `data/` directory.

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
$ mkdir ../trimmed_fastq
$ mv *.trim* ../trimmed_fastq
$ cd ../trimmed_fastq
$ ls
~~~
{: .bash}

~~~
SRR2584863_1.trim.fastq.gz    SRR2584866_1.trim.fastq.gz    SRR2589044_1.trim.fastq.gz
SRR2584863_2.trim.fastq.gz    SRR2584866_2.trim.fastq.gz    SRR2589044_2.trim.fastq.gz
~~~
{: .output}

> ## Bonus Exercise (Advanced)
>
> Now that our samples have gone through quality control, they should perform
> better on the quality tests run by FastQC. Go ahead and re-run
> FastQC on your trimmed FASTQ files and visualize the HTML files
> to see whether your per base sequence quality is higher after
> trimming.
>
>> ## Solution
>>
>> In your AWS terminal window do:
>>
>> ~~~
>> $ fastqc ~/dc_workshop/data/trimmed_fastq/*.fastq*
>> ~~~
>> {: .bash}
>>
>> In a new tab in your terminal do:
>>
>> ~~~
>> $ mkdir ~/Desktop/fastqc_html/trimmed
>> $ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/trimmed_fastq/*.html ~/Desktop/fastqc_html/trimmed
>> ~~~
>> {: .bash}
>> 
>> Then take a look at the html files in your browser.
>> 
>> Remember to replace everything between the `@` and `:` in your scp
>> command with your AWS instance number.
>>
>> After trimming and filtering, our overall quality is much higher, 
>> we have a distribution of sequence lengths, and more samples pass 
>> adapter content. However, quality trimming is not perfect, and some
>> programs are better at removing some sequences than others. Because our
>> sequences still contain 3' adapters, it could be important to explore
>> other trimming tools like [cutadapt](http://cutadapt.readthedocs.io/en/stable/) to remove these, depending on your
>> downstream application. Trimmomatic did pretty well though, and its performance
>> is good enough for our workflow.
> {: .solution}
{: .challenge}
