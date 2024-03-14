# Transcriptional control of the *Cryptosporidium* lifecycle

## ABSTRACT 

The parasite *Cryptosporidium* is a leading agent of diarrheal disease in young children and a cause and consequence of chronic malnutrition. There are no vaccines and only limited treatment options. The parasite infects enterocytes, where it engages in asexual and sexual replication, both of which are essential to continued infection and transmission. However, their molecular mechanisms remain largely unknown. Here we use single-cell RNA sequencing to reveal the gene expression program of the entire *C. parvum* life cycle in culture and infected animals. Diverging from the prevailing model, we find support for only three intracellular stages: asexual type-I meronts, male gamonts, and female gametes. We uncover a highly organized program for the assembly of components at each stage. Dissecting the underlying regulatory network, we identify the transcription factor Myb-M as the earliest determinant of male fate, in an organism that lacks genetic sex determination. Conditional expression of this factor overrides the developmental program and induces widespread maleness, while conditional deletion ablates male development. Both have a profound impact on the infection. A large set of stage-specific genes now provides the opportunity to understand, engineer, and disrupt parasite sex and life cycle progression to advance the development of vaccines and treatments.

## DIRECTORY OUTLINE

/ANALYSIS - contains the Rmarkdown (.rmd) file that includes all code for the paper, now published in mBio. Also contains the R file with code. Outputs used for sequencing analysis are in this folder, as well as the study design file. The table for differentially expressed genes (color-coded for AP2-F and crystalloid body) is also in this folder.

/ANALYSIS/GSEA - contains the files used for gene set enrichment analysis of oocyst wall and crystalloid body genes.

/QC - contains the outputs from running fastqc on all raw fastq files.
