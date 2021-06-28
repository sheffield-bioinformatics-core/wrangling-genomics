---
title: "Background and Metadata"
teaching: 10
exercises: 5
questions:
- "What data are we using?"
- "Why is this experiment important?"
objectives:
- "Why study Human Genomes?"
- "Understand the data set."
keypoints:
- "It's important to record and understand your experiment's metadata."
output:
  html_notebook:
    toc: yes
    toc_float: yes
    css: ../stylesheets/styles.css
  html_document:
    df_print: paged
    toc: yes
editor_options:
  chunk_output_type: inline
---
<img src="../img/logo-sm.png" align=right>


# Background

We are going to use a sequencing dataset from healthy humans. 

 - **What is the 1000 genomes project**
    - The 1000 genomes project was established in 2008 to study variation in the *human genome* and provide a solid foundation on which to build an understanding of genetic variation in the human population.
    

 - **Why is the 1000 genomes project important**
    - The data generated for the 1000 genomes project can be incorporated into many healthcare studies. For example, when identifying mutations in a diseased individual we can use mutations identified among healthy individuals to narrow-down our search of potential disease-causing mutations.
    
# The Data

 - We have selected three individuals from 1000 genomes and will be working with a subset of the data for these individuals. This is to make the tools and workflows run in a reasonable amount of time. When analysing your own data the same steps can be applied, although they will take much longer to complete
 
 
## View the Metadata

The metadata file associated with this lesson can be [downloaded directly here]("https://raw.githubusercontent.com/sheffield-bioinformatics-core/wrangling-genomics/gh-pages/files/igsr-1000 genomes phase 3 release.tsv.tsv") or [viewed in Github]("https://github.com/sheffield-bioinformatics-core/wrangling-genomics/blob/gh-pages/files/igsr-1000 genomes phase 3 release.tsv.tsv"). If you would like to know details of how the file was created, you can look at [some notes and sources here](https://github.com/datacarpentry/wrangling-genomics/blob/gh-pages/files/1000genomes_README.md).



This metadata describes information on the samples sequences as part of the dataset and the columns represent:

| Column           | Description                                |
|------------------|--------------------------------------------|
| Sample name           | Sample name					|
| Sex       | Sex	|
| Biosample ID            | 		|
| Population code        | Short code for the population			|
| Population name       | Longer, descriptive name for the population |
| Superpopulation code          | Grouping of populations from a similar geographic area (e.g. continent)|
| Population elastic ID       |  |
| Data collections            | Which datasets the sample belongs to		|



> ## Challenge
> 
> Based on the metadata, can you answer the following questions?
> 
> 1. How many different individuals exist in the data?
> 2. How many rows and how many columns are in this data?
> 3. How many different super populations are there?
> 4. How many different populations exist with European origin?
>
>> ## Solution
>> 
>> 1. 25 different generations
>> 2. 62 rows, 12 columns
>> 3. Nine different sub-populations
>> 4. Five populations within Europe
> {: .solution}
{: .challenge}



