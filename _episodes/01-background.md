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


- **What is Genome sequencing**
  - Genome sequencing (sometimes called *next-generation sequencing* (NGS) or high throughput sequencing) is the process by which small stretches of an individuals' DNA are "read" to see which bases (A, T, C or G) they are comprised of. These reads are they compared to a reference genome to see where they originated from and what mutations are present. Mutations, differences in DNA sequence between individuals, are not always harmful and can be responsible in normal variations such as eye colour. However, some mutations have the potential to lead to the progression of disease. The most popular method of sequencing is that employed by Illumina, which is demonstrated in this short video.
  
<iframe width="560" height="315" src="https://www.youtube.com/embed/fCd6B5HRaZ8" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

 - **What is the 1000 genomes project**
    - The 1000 genomes project was established in 2008 to study variation in the *human genome* and provide a solid foundation on which to build an understanding of genetic variation in the human population.
    

 - **Why is the 1000 genomes project important**
    - The data generated for the 1000 genomes project can be incorporated into many healthcare studies. For example, when identifying mutations in a diseased individual we can use mutations identified among healthy individuals to narrow-down our search of potential disease-causing mutations.
    
# The Data

 - We have selected three individuals from 1000 genomes and will be working with a subset of the data for these individuals. This is to make the tools and workflows run in a reasonable amount of time. When analysing your own data the same steps can be applied, although they will take much longer to complete
 
 
## View the Metadata

The metadata file associated with this lesson can be [downloaded directly here](https://raw.githubusercontent.com/sheffield-bioinformatics-core/wrangling-genomics/gh-pages/files/1000_genomes_meta.tsv) (right-click and Save Link as) or [viewed in Github](https://github.com/sheffield-bioinformatics-core/wrangling-genomics/blob/gh-pages/files/1000_genomes_meta.tsv). If you would like to know details of how the file was created, you can look at [some notes and sources here](https://github.com/sheffield-bioinformatics-core/wrangling-genomics/blob/gh-pages/files//1000genomes_README.md).



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
> Based on the metadata, can you answer the following questions using a spreadsheet such as Excel?
> 
> 1. How many rows and how many columns are in this data?
> 2. How many different super populations are there?
> 3. How many different populations exist with European origin?
>
>> ## Solution
>> 
>> 1. 3116 rows and 9 columns
>> 2. Nine different sub-populations
>> 3. Five populations within Europe
> {: .solution}
{: .challenge}


> ## Creating and editing metadata
> 
> The metadata for a project is usually entered *by-hand* using software such as Microsoft Excel. When creating such metadata it would be good to bear in mind some common errors that can be inadvertently introduced that complicate computational analysis. These materials from Data Carpentry can be consulted if you are not sure about this:-[https://datacarpentry.org/spreadsheet-ecology-lesson/02-common-mistakes/index.html](https://datacarpentry.org/spreadsheet-ecology-lesson/02-common-mistakes/index.html)
{: .callout}


