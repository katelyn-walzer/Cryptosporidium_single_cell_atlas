# Transcriptional control of the *Cryptosporidium* lifecycle

## ABSTRACT 

The parasite *Cryptosporidium* is a leading agent of diarrheal disease in young children and a cause and consequence of chronic malnutrition. There are no vaccines and only limited treatment options. The parasite infects enterocytes, where it engages in asexual and sexual replication, both of which are essential to continued infection and transmission. However, their molecular mechanisms remain largely unknown. Here we use single-cell RNA sequencing to reveal the gene expression program of the entire *C. parvum* life cycle in culture and infected animals. Diverging from the prevailing model, we find support for only three intracellular stages: asexual type-I meronts, male gamonts, and female gametes. We uncover a highly organized program for the assembly of components at each stage. Dissecting the underlying regulatory network, we identify the transcription factor Myb-M as the earliest determinant of male fate, in an organism that lacks genetic sex determination. Conditional expression of this factor overrides the developmental program and induces widespread maleness, while conditional deletion ablates male development. Both have a profound impact on the infection. A large set of stage-specific genes now provides the opportunity to understand, engineer, and disrupt parasite sex and life cycle progression to advance the development of vaccines and treatments.

## DIRECTORY OUTLINE

**/Bulk_RNA_Seq_Analysis -** contains the code and results for bulk RNA sequencing of sorted asexual, male, and female populations. Asexual and female datasets were obtained from GEO:GSE129267 while the male in vitro dataset was produced in this study. Supplementary_Data_1.R contains the code for QC, read mapping, filtering, normalization, and differential expression analysis of these stage-specific datasets. Outputs used for sequencing analysis are in this folder, as well as the study design file. Supplementary_Table_3 is a multi-tab Excel file containing differentially expressed genes between all samples while Supplementary_Table_4 is a multi-tab Excel file listing male and female marker genes identified in this study. The folder **QC** contains the outputs from running fastqc on all raw fastq files.

**/Bulk_RNA_Seq_Myb-M_Gain-of-Function_Analysis -** contains the code and results for bulk RNA sequencing of samples collected at 30 hours post invasion after 6 hours of vehicle or Shield-1 treatment. Shield-1 drives Myb-M protein expression in a gain-of-function mutant. Supplementary_Data_2.R contains the code for QC, read mapping, filtering, normalization, and differential expression analysis across each strain and condition. Outputs used for sequencing analysis are in this folder, as well as the study design file. Supplementary_Table_12 is a multi-tab Excel file containing differentially expressed genes between all samples (asexual- and male-specific genes are highlighted in green and blue, respectively). The folder **QC** contains the outputs from running fastqc on all raw fastq files.

**/scRNA_Seq_Analysis -** contains the code and results for single-cell RNA sequencing (scRNA-seq) of parasites captured from HCT-8 culture or an IFN-gamma knockout mouse. In vitro samples were collected at 24, 36, 42, 46, and 54 hours post invasion while the in vivo sample was collected from an acutely infected mouse on day 6. The filtered_feature_bc_matrix (Cell Ranger 3.1.0 output) consists of features (tsv), barcodes (tsv), and matrix (mtx) files which are all in this folder. The features (tsv) file is the same for all samples. Supplementary_Table_1 contains a summary of parasites captured for each sample and their sequencing information. Analyses were carried out on asexual samples (24 and 36 hours) or on all samples.

&nbsp;&nbsp;&nbsp;&nbsp;<ins>Asexual Single-Cell Atlas</ins> <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Supplementary_Data_3.R contains the code used for alignment, filtering, normalization, clustering, and differential gene expression &nbsp;&nbsp;&nbsp;&nbsp;analysis. <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Supplementary_Data_5.R contains the code for RNA velocity and pseudotime analyses of asexual single-cell transcriptomes. <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Asexual_metadata includes the metadata from this analysis. <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Supplementary_Table_2 is a multi-tab Excel file containing differentially expressed genes across asexual clusters (first tab) and lists of 
&nbsp;&nbsp;&nbsp;&nbsp;genes that were plotted across pseudotime (all other tabs). <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Supplementary_Table_8 is a multi-tab Excel file containing the results of differential gene expression analysis between in vitro 24 and 36 
&nbsp;&nbsp;&nbsp;&nbsp;hour samples by cluster. <br />

&nbsp;&nbsp;&nbsp;&nbsp;<ins>*Cryptosporidium* Single-Cell Atlas (All Samples)</ins> <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Supplementary_Data_4.R contains the code used for alignment, filtering, normalization, clustering, and differential gene expression &nbsp;&nbsp;&nbsp;&nbsp;analysis. <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Supplementary_Data_5.R contains the code for pseudotime analyses of asexual, male, and female single-cell transcriptomes. <br />
&nbsp;&nbsp;&nbsp;&nbsp;* Metadata

**/scRNA_Seq_Analysis_with_lncRNAs -** 
