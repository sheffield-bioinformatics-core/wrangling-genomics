---
title: "Running Genomics workflows on HPC"
teaching: 30
exercises: 15
questions:
- "How can I run my analyses in parallel using HPC?"
objectives:
- "Submitting a analysis script to a queue."
- "Identifying the opportunities for parallelisation."
- "Using a job array"
keypoints:
- "Job arrays can make our lives easier"
- "Some tools are able to use multiple threads"

---

<img src="../img/logo-sm.png" align=right>


# Submitting a script as a batch job

We successfully written a script to automate our analysis. However, we are not taking full advantage of the features that a HPC system offers us. As saw previously, a batch script may be submitted to a scheduler for execution on one of the computing nodes.

We have modified the script slightly to include lines of code required for scheduling. th

~~~
curl -O https://raw.githubusercontent.com/sheffield-bioinformatics-core/wrangling-genomics/gh-pages/files/run_variant_calling_batch.sh

~~~
{: .bash}

The lines we have added are as follows:-

~~~
#!/usr/bin/env bash
#SBATCH -J variant_calling_batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=6G

~~~
{: .bash}

As we have seen previously, this can submitted to the queue with the following command. We could now work on other tasks whilst the analysis runs in the background.

~~~
sbatch run_variant_calling_batch.sh
~~~
{: .bash}


# Introducing the task array

We have automated our analysis and it should hopefully work on larger sample sizes than we using in this workshop. We have chosen a small, unrealistic, number of reads for each individual. In a real-life experiment each individuals may have many millions of reads and many Gbs of data to process. The data from each individual could take a day, or more, to run by itself, meaning the pipeline we have created could take many days, or weeks, to run to completion on a more realistic dataset. 

Our workflow processes individuals in a sequential manner; i.e. NA12874 is aligned only after all the steps for NA12873 have completed. This is clearly inefficient as the processing of NA12874 is **completely independent** of NA12873 so there is no reason why we should wait to complete NA12873 before starting the analysis of NA12874. The term for this is [Embarrassing parallel](https://en.wikipedia.org/wiki/Embarrassingly_parallel), or pleasingly parallel. Such a situation can be solved using a technique called a job array or task array 



