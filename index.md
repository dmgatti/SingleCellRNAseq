---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---
Single cell RNA-sequencing (scRNA-Seq) is a method of quantifying transcript expression levels in individual cells. scRNA-Seq technology can take on many different forms and this area of research is rapidly evolving. In 2022, the most widely used systems for performing scRNA-Seq involved separating cells and introducing them into a microfluidic system which performs the chemistry on each cell individually (droplet-based scRNA-Seq).

In this workshop we will primarily focus on the 10X Genomics technology. We will review the steps in single cell isolation and library preparation. We will then work with gene count data produced by CellRanger and will proceed through filtering, normalization, selection of variable genes, clustering, and finding marker genes for each cluster. Each student will analyze the data on their laptop.

<!-- this is an html comment -->

{% comment %} This is a comment in Liquid {% endcomment %}

> ## Prerequisites
>
> * Basic knowledge of R.
> * Basic understanding of cellular and molecular biology.
> * Understanding of bulk RNA-sequencing.
> 
{: .prereq}

{% include links.md %}
