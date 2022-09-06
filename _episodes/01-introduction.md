---
title: "Introduction to single cell RNA-seq"
teaching: 30
exercises: 10
questions:
- "What is single-cell RNA-seq?"
- "What is the difference between bulk RNA-seq and single-cell RNA-seq?"
- "How do I choose between bulk RNA-seq and single-cell RNA-seq"
objectives:
- "Describe what single-cell RNA-seq is."
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

## Brief overview of single cell transcriptomics technology 

Single cell RNA-sequencing (scRNA-Seq) is a method of quantifying transcript expression levels in individual cells. scRNA-Seq technology can take on many different forms and this area of research is rapidly evolving. In 2022, the most widely used systems for performing scRNA-Seq involve separating cells and introducing them into a microfluidic system which performs the chemistry on each cell individually (droplet-based scRNA-Seq).

In this workshop we will primarily focus on the 10X Genomics technology. 10X Genomics is a market leader in the single cell space and was among the first technologies that made it feasible to profile thousands of cells simultaneously. Single cell technology is changing rapidly and it is not clear whether any other companies will be able to successfully challenge 10X's dominance in this space. 

The steps in droplet scRNA-Seq are:

1.  Cell isolation: 
    * If the cells are part of a tissue, the cells are disaggregated using collagenase or other reagents. The specifics of this protocol can vary greatly due to differences between tissues that are biological in nature. If the cells are in culture or suspension, they may be used as is.
2. Assess cell viability.
    * If scRNA-Seq is being performed on fresh tissue, the cells are usually checked for viability. We want "happy" cells loaded into the machine. We might hope for >90% viable and set a minimum threshold of >70%, although these numbers can vary depending on the experiment.
3. Cell suspension:
    * Using a microfluidic system, each cell is suspended in a nanoliter-size droplet along with a barcoded primer bead. The cells are kept separate from each other in an oil/water emulsion.
4. Cell lysis, generating cDNA:
    * The cells are lysed in each droplet. Each cell was already encapsulated with a barcoded primer bead which has a primer specific to that cell. Often a poly-d(T) primer is used to prime the poly(A) tail of mRNA. 
7. Library generation:
    * [DAS check this] Oil-water emulsion is "de-emulsified". Amplify cDNA and add correct primers for Illuina sequencing. Sequence on any old Illumina machine. Sequencing should be paired-end, one read contains cell and molecule barcodes, other read contains the bit of transcript that was captured.

## Comparing and contrasting scRNA-Seq with bulk RNA-Seq 

Bulk RNA-Seq and single cell RNA-Seq are related in that they both assess transcription levels by sequencing short reads, but these two technologies have a variety of differences. Neither technology is always better. The approach that one might use should depend upon the information one hopes to gather.

Consider the following points when assessing the differences between the technologies and choosing which to utilize for your own experiment:
 * Tissues are heterogeneous mixtures of diverse cell types. Bulk RNA-Seq data consists of average measures of transcripts expressed across many different cell types, while scRNA-Seq data is cell-type resolved.
 * Bulk RNA-Seq data may obscure the processes of changes in gene expression or tissue cell composition.
 * Bulk RNA-Seq allows for much higher sequencing coverage for each gene.
 * Bulk RNA-Seq allows for better isoform detection due to the higher sequencing depth and relatively uniform coverage across transcripts (vs. a typical 3' bias in scRNA-Seq).
 * Genes without poly-A tail (e.g. some noncoding RNAs) might not be detected in scRNA-Seq, but can be reliably assessed using bulk RNA-Seq.

> ## Challeng
> For each of these scenarios, choose between using bulk RNA-Seq and scRNA-Seq to address your problem.
>
> Differentiation of embryonic stem cells to another cell type
> > ## Solution
> > You would likely find single cell RNA-Seq most powerful in this situation since the cells are differentiating along a continuous transcriptional gradient.
> {: .solution}
> 
> Aging and Cdkn2a
> > ## Solution
> > bulk
> {: .solution}
> 
> PBMCs -- 
> > ## Solution
> > single cell
> {: .solution}
> 
> New non-model species
> > ## Solution
> > bulk
> {: .solution}
> 
> Studying miRNAs
> > ## Solution
> > bulk
> {: .solution}
> 
> Very heterogeneous tissue
> > ## Solution
> > You would likely find single cell RNA-Seq most powerful in this situation since ...
> {: .solution}
> 
> eQTL mapping
> > ## Solution
> > probably bulk, but both could be informative!
> {: .solution}
{: .challenge}


## What is scRNA-Seq useful for? 

Single cell RNA-Seq is a new technology and its uses are limited only by your imagination! A few examples of problems that have been addressed using scRNA-Seq include:
 * developmental studies & studies of cellular trajectories.
 * detailed tissue atlases.
 * characterization of tumor clonality.
 * definition of cell-type specific transcriptional responses (e.g. T-cell response to infection)
 * profiling of changes in cell state (i.e. homeostasis vs. response state)
 * a variety of different types of CRISPR screens

## Emphasize focus of this course: 10X Genomics, mouse

While there are several commercially available scRNA-seq technologies, this course will focus on data generated by the 10X Genomics platform. This technology is widely used at JAX, and it is flexible, reliable, and relatively cost-efficient. We will focus on a mouse data set, although most of the techniques that you will learn in this workshop apply equally well to other species.  

{% include links.md %}

