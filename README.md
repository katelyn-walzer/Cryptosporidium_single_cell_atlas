# Transcriptional control of the *Cryptosporidium* lifecycle

## ABSTRACT 

The parasite *Cryptosporidium* is a leading agent of diarrheal disease in young children and a cause and consequence of chronic malnutrition. There are no vaccines and only limited treatment options. The parasite infects enterocytes, where it engages in asexual and sexual replication, both of which are essential to continued infection and transmission. However, their molecular mechanisms remain largely unknown. Here we use single-cell RNA sequencing to reveal the gene expression program of the entire *C. parvum* life cycle in culture and infected animals. Diverging from the prevailing model, we find support for only three intracellular stages: asexual type-I meronts, male gamonts, and female gametes. We uncover a highly organized program for the assembly of components at each stage. Dissecting the underlying regulatory network, we identify the transcription factor Myb-M as the earliest determinant of male fate, in an organism that lacks genetic sex determination. Conditional expression of this factor overrides the developmental program and induces widespread maleness, while conditional deletion ablates male development. Both have a profound impact on the infection. A large set of stage-specific genes now provides the opportunity to understand, engineer, and disrupt parasite sex and life cycle progression to advance the development of vaccines and treatments.

## DIRECTORY OUTLINE

**/Bulk_RNA_Seq_Analysis -** contains the code and results for bulk RNA sequencing of sorted asexual, male, and female populations. Asexual and female datasets were obtained from GEO:GSE129267 while the male in vitro dataset was produced in this study. Supplementary_Data_1.R contains the code for QC, read mapping, filtering, normalization, and differential expression analysis of these stage-specific datasets. Outputs used for sequencing analysis are in this folder, as well as the study design file. Supplementary_Table_3 is a multi-tab Excel file containing differentially expressed genes between all samples while Supplementary_Table_4 is a multi-tab Excel file listing male and female marker genes identified in this study. The folder **QC** contains the outputs from running fastqc on all raw fastq files.

**/Bulk_RNA_Seq_Myb-M_Gain-of-Function_Analysis -** 

**/scRNA_Seq_Analysis -**

**/scRNA_Seq_Analysis_with_lncRNAs -** 
