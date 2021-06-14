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
    - 
    

 - **Why is the 1000 genomes project important**
    - 
    
# The Data

 - The data we are going to use is part of the 1000 genomes project. We have selected three individuals and will be working with a subset of the data for these individuals. This is to make the tools and workflows run in a reasonable amount of time. When analysing your own data the same steps can be applied, although they will take much longer to complete
 
 
## View the Metadata

We will be working with three individuals from the 1000 genomes project 


The metadata file associated with this lesson can be [downloaded directly here]("https://raw.githubusercontent.com/sheffield-bioinformatics-core/wrangling-genomics/gh-pages/files/igsr-1000 genomes phase 3 release.tsv.tsv") or [viewed in Github]("https://github.com/sheffield-bioinformatics-core/wrangling-genomics/blob/gh-pages/files/igsr-1000 genomes phase 3 release.tsv.tsv"). If you would like to know details of how the file was created, you can look at [some notes and sources here](https://github.com/datacarpentry/wrangling-genomics/blob/gh-pages/files/1000genomes_README.md).



This metadata describes information on the samples sequences as part of the dataset and the columns represent:

| Column           | Description                                |
|------------------|--------------------------------------------|
| Sample name           | Sample name					|
| Sex       | Sex	|
| Biosample ID            | Something		|
| Population code        | Short code for the population			|
| Population name       | Longer, descriptive name for the population |
| Superpopulation code          | e.g. the continent |
| Population elastic ID       | Something |
| Data collections            | Which datasets the sample belongs to		|


<div class="exercise">

> ## Challenge
> 
> Based on the metadata, can you answer the following questions?
> 
> 1. How many different populations exist in the data?
> 2. How many rows and how many columns are in this data?
> 3. How many citrate+ mutants have been recorded in **Ara-3**?
> 4. How many hypermutable mutants have been recorded in **Ara-3**?
>
> > ## Solution
>> 
> > 1. 25 different generations
> > 2. 62 rows, 12 columns
> > 3. 10 citrate+ mutants
> > 4. 6 hypermutable mutants
> {: .solution}
{: .challenge}

</div>

<!-- can add some additional info relevant to interplay of hypermutability and Cit+ adaptations, but keep it simple for now -->

