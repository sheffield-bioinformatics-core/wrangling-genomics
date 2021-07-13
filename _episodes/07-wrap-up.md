---
title: "Workshop wrap-up"
teaching: 30
exercises: 5
questions:
- "What should I learn about next?"
objectives:
- "Describing workflow managers"
- "Make a copy of everything you've done in the workshop"
keypoints:
- "Many reproducible pipelines and workflows are already available"
- "No need to re-invent the wheel"

---

<img src="../img/logo-sm.png" align=right>

# Next steps: Workflow Managers and Reproducible Pipelines

At this point in the course we hope that you have gained an understanding of what it means to write and run a Bioinformatics script on HPC infrastructure; and how to troubleshoot if things goes wrong. 

However, in practice we would **not** recommend that you now start to write and develop your own pipelines from scratch. In all researchers wrote their own pipelines from scratch, this would lead to a huge amount of redundancy and discrepancies in the manner in which data are processed. 

Many Bioinformaticians instead prefer to use tried-and-tested *best practice* pipelines for data processing; especially for steps such as alignment where the tools are established, well-understood and don't usually require much intervention from the user. 


A few *workflow manager* tools also exist, the purpose of which is to manage to monitor the submission of jobs without requiring the user to monitor the workflow and manage scripts that are dependent on the outputs of others. In other words, alignment of samples will automatically be triggered once the workflow manager detects that trimming has been completed; even if the commands to perform trimming and alignment are found in different scripts.

A few tools that the Sheffield Bioinformatics Core recommends are:-

- [bcbio](https://bcbio-nextgen.readthedocs.io/en/latest/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [nextflow](https://www.nextflow.io/)
  + [nf-core (pipelines built using nextflow)](https://nf-co.re/)
  
# Making a copy of the scripts and data from this workshop

Due to the way in which we are running this workshop, the **cluster instance will not be available** after today. If you want to look back and review anything you have written, it would be a good idea to make a local copy


**Make sure you run this from your own machine (Desktop / laptop)**

~~~
$ mkdir -p ~/Desktop/n8_workshop/
$ scp -r <your_username>@ephemeron.n8cir.org.uk:~/dc_workshop ~/Desktop/n8_workshop/
~~~
{: .bash}


> ## Exercise
> Use the command listed above to copy everything you have created in the `~/dc_workshop` folder to your home directory.
{: .challenge}