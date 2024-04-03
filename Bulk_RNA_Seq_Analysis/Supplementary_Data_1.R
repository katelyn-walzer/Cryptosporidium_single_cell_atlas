#First do QC and read mapping on the command line using python
#Quality control of raw reads was carried out using fastqc
#Raw reads were mapped to the version 46 Cryptosporidium parvum Iowa II reference transcriptome available on VEuPathDB using Kallisto, version 0.46.0. 

#Kallisto make new reference using CryptoDB version 46
#Use annotated transcripts file
kallisto index -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.fasta

#Run fastqc on male sorted samples
#Note: fastq.gz filenames for GEO accession number GSE232438 have been updated to reflect male replicates 1-4
#13-MALE-KWJT_S13_mergedLanes_R1_001.fastq.gz equals Male_in_vitro_1_mergedLanes_R1_001.fastq.gz
#14-MALE-KWJT_S14_mergedLanes_R1_001.fastq.gz equals Male_in_vitro_2_mergedLanes_R1_001.fastq.gz
#15-MALE-KWJT_S15_mergedLanes_R1_001.fastq.gz equals Male_in_vitro_3_mergedLanes_R1_001.fastq.gz
#16-MALE-KWJT_S16_mergedLanes_R1_001.fastq.gz equals Male_in_vitro_4_mergedLanes_R1_001.fastq.gz
fastqc 13-MALE-KWJT_S13_mergedLanes_R1_001.fastq.gz -t 24 -o /venice/striepenlab/Katelyn_Walzer_data/Jayesh_manuscript/Male_Sort_CryptoDB46/fastqc
fastqc 14-MALE-KWJT_S14_mergedLanes_R1_001.fastq.gz -t 24 -o /venice/striepenlab/Katelyn_Walzer_data/Jayesh_manuscript/Male_Sort_CryptoDB46/fastqc
fastqc 15-MALE-KWJT_S15_mergedLanes_R1_001.fastq.gz -t 24 -o /venice/striepenlab/Katelyn_Walzer_data/Jayesh_manuscript/Male_Sort_CryptoDB46/fastqc
fastqc 16-MALE-KWJT_S16_mergedLanes_R1_001.fastq.gz -t 24 -o /venice/striepenlab/Katelyn_Walzer_data/Jayesh_manuscript/Male_Sort_CryptoDB46/fastqc

#Male data and data from Tandel et al., Nature Microbiology 2019 aligned to CryptoDB46
#Average length of fragments from male replicates is 400 bp from TapeStation
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o male_CryptoDB_46_rep_1 -t 24 -b 60 --single -l 400 -s 100 13-MALE-KWJT_S13_mergedLanes_R1_001.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o male_CryptoDB_46_rep_2 -t 24 -b 60 --single -l 400 -s 100 14-MALE-KWJT_S14_mergedLanes_R1_001.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o male_CryptoDB_46_rep_3 -t 24 -b 60 --single -l 400 -s 100 15-MALE-KWJT_S15_mergedLanes_R1_001.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o male_CryptoDB_46_rep_4 -t 24 -b 60 --single -l 400 -s 100 16-MALE-KWJT_S16_mergedLanes_R1_001.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vitro_CryptoDB_46_rep_1 -t 24 -b 60 --single -l 500 -s 100 Female_in_vitro_1_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vitro_CryptoDB_46_rep_2 -t 24 -b 60 --single -l 500 -s 100 Female_in_vitro_2_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vitro_CryptoDB_46_rep_3 -t 24 -b 60 --single -l 500 -s 100 Female_in_vitro_3_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vitro_CryptoDB_46_rep_4 -t 24 -b 60 --single -l 500 -s 100 Female_in_vitro_4_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vivo_CryptoDB_46_rep_1 -t 24 -b 60 --single -l 500 -s 100 Female_in_vivo_1_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vivo_CryptoDB_46_rep_2 -t 24 -b 60 --single -l 500 -s 100 Female_in_vivo_2_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vivo_CryptoDB_46_rep_3 -t 24 -b 60 --single -l 500 -s 100 Female_in_vivo_3_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o female_in_vivo_CryptoDB_46_rep_4 -t 24 -b 60 --single -l 500 -s 100 Female_in_vivo_4_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o asexual_CryptoDB_46_rep_1 -t 24 -b 60 --single -l 500 -s 100 Asexual_1_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o asexual_CryptoDB_46_rep_2 -t 24 -b 60 --single -l 500 -s 100 Asexual_2_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o asexual_CryptoDB_46_rep_3 -t 24 -b 60 --single -l 500 -s 100 Asexual_3_mergedLanes.fastq.gz
kallisto quant -i CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.index -o asexual_CryptoDB_46_rep_4 -t 24 -b 60 --single -l 500 -s 100 Asexual_4_mergedLanes.fastq.gz
#End python on terminal


# load packages ----
library(tidyverse) 
library(tximport) 
library(hrbrthemes) 
library(RColorBrewer) 
library(reshape2)
library(genefilter) 
library(edgeR) 
library(matrixStats)
library(DT) 
library(gt) 
library(plotly) 
library(skimr)
library(limma)

setwd("/Volumes/KATIE/2023_1_3_Sorted_samples_bulk_Crypto46")
getwd()

# read in study design ----
targets <- read_tsv("Study_Design_Sorted_Samples_Crypto46.txt")

# create file paths to the abundance files generated by Kallisto using the 'file.path' function
path <- file.path(targets$sample, "abundance.h5")

# check to make sure this path is correct by seeing if the files exist
all(file.exists(path)) 

# use dplyr to modify study design to include these file paths as a new column.
targets <- mutate(targets, path)

# import Kallisto transcript counts into R using Tximport ----
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     txOut = TRUE,
                     countsFromAbundance = "lengthScaledTPM")

#take a look at the type of object
class(Txi_gene)
names(Txi_gene)

myCPM <- as_tibble(Txi_gene$abundance, rownames = "geneSymbol") # counts after adjusting for transcript length
myCounts <- as_tibble(Txi_gene$counts, rownames = "geneSymbol") # counts per million (CPM) 

# Identify variables of interest in study design file ----
stage <- as.factor(targets$stage)
sex <- as.factor(targets$sex)
origin <- as.factor(targets$origin)
host <- as.factor(targets$host)

#Capture sample labels for later use
sampleLabels <-targets$sample

# Examine data up to this point ----
myCPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts

# graph both matrices
colSums(myCPM)
colSums(myCounts)

# Take a look at the heteroskedasticity of the data ----
# first, calculate row means and standard deviations for each transcript or gene 
# and add these to your data matrix
myCPM.stats <- transform(myCPM, 
                         SD=rowSds(myCPM), 
                         AVG=rowMeans(myCPM),
                         MED=rowMedians(myCPM)
)

head(myCPM.stats)

#produce a scatter plot of the transformed data
ggplot(myCPM.stats, aes(x=SD, y=MED)) +
  geom_point(shape=16, size=2)


# Make a DGElist from the counts, and plot ----
myDGEList <- DGEList(Txi_gene$counts)

save(myDGEList, file = "myDGEList")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList) 
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)

# Take a look at the distribution of the Log2 CPM
nsamples <- ncol(log2.cpm)
# now select colors from a single palette
myColors <- brewer.pal(nsamples, "Paired")

# 'coerce' the data matrix to a dataframe to use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm)

colnames(log2.cpm.df) <- sampleLabels

# use the reshape2 package to 'melt' dataframe (from wide to tall)
log2.cpm.df.melt <- melt(log2.cpm.df)

ggplot(log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
  caption=paste0("produced on ", Sys.time())) +
  coord_flip() 

# Filter the data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==16)

# set some cut-off to get rid of genes/transcripts with low counts
keepers <- rowSums(cpm>10)>=4 #last number is replicates, min replicates, at least expressed ten times
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered) 
colnames(log2.cpm.filtered.df) <- sampleLabels
log2.cpm.filtered.df.melt <- melt(log2.cpm.filtered.df)

ggplot(log2.cpm.filtered.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
  caption=paste0("produced on ", Sys.time())) +
  coord_flip() 

# Normalize the data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)
colnames(log2.cpm.filtered.norm.df) <- sampleLabels
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)

ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.y = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
  caption=paste0("produced on ", Sys.time())) +
  coord_flip() 

# Need to convert the datamatrix to a dataframe, while preserving the rownames as a new column in the dataframe
mydata.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneSymbol")
colnames(mydata.df) <- c("geneSymbol", sampleLabels)
skim(mydata.df)
mydata.melt <- as_tibble(melt(mydata.df))

write_tsv(mydata.df, "normData_CryptoDB_46.txt") # Note: this is the data to use as input .gct file for GSEA analysis

#Label so we can see which samples affect PCs
colnames(log2.cpm.filtered.norm) <- sampleLabels

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
x <- pca.res$rotation 
pc.var<-pca.res$sdev^2 
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize the PCA result ------------------
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=targets$stage)) +
  geom_point(size=4) +
  theme_bw() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot of Cryptosporidium stages")


# Set up design matrix ----
stage <- relevel(stage, "male")
design <- model.matrix(~0 + stage)
colnames(design) <- levels(stage)

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DGEList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DGEList.filtered.norm, design)

# Contrast matrix ----
# how are males different from other stages?
contrast.matrix <- makeContrasts(male_vs_female_in_vitro = male - female_in_vitro,
                                 male_vs_female_in_vivo = male - female_in_vivo,
                                 male_vs_asexual = male - asexual,
                                 female_in_vitro_vs_asexual = female_in_vitro - asexual,
                                 female_in_vivo_vs_asexual = female_in_vivo - asexual,
                                 female_in_vivo_vs_female_in_vitro = female_in_vivo - female_in_vitro,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for linear model fit
ebFit <- eBayes(fits)

#Compare male versus female in vitro
M_vs_F_invitro_top_hits <- topTable(ebFit, adjust ="BH", coef=1, number=10000, sort.by="logFC")

# convert to a tibble
M_vs_F_invitro_top_hits <- as_tibble(M_vs_F_invitro_top_hits, rownames = "geneSymbol")


#Compare male versus female in vivo
M_vs_F_invivo_top_hits <- topTable(ebFit, adjust ="BH", coef=2, number=10000, sort.by="logFC")

# convert to a tibble
M_vs_F_invivo_top_hits <- as_tibble(M_vs_F_invivo_top_hits, rownames = "geneSymbol")


#Compare male versus asexual
M_vs_asexual_top_hits <- topTable(ebFit, adjust ="BH", coef=3, number=10000, sort.by="logFC")

# convert to a tibble
M_vs_asexual_top_hits <- as_tibble(M_vs_asexual_top_hits, rownames = "geneSymbol")


#Make volcano plots for males versus asexual or female in vivo

#Do male versus female in vivo first
#Subset male genes
AP2M_HAP2_GGC_subset_female_invivo <- subset(M_vs_F_invivo_top_hits,
                                             geneSymbol=="cgd8_2220-RA" |
                                               geneSymbol=="cgd6_2670-RA" |
                                               geneSymbol=="cgd2_2610-RA" |
                                               geneSymbol=="cgd5_3570-RA" |
                                               geneSymbol=="cgd7_5500-RA" |
                                               geneSymbol=="cgd8_1740-RA")

#Plot, Figure 2d
ggplot(M_vs_F_invivo_top_hits, aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneSymbol))) +
  geom_point(size=4) +
  geom_point(mapping = NULL, AP2M_HAP2_subset_female_invivo, size = 4, colour = "#0095FF", inherit.aes = TRUE) +
  ylim(-0.5,12) +
  xlim(-15,15) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", colour="grey", size=0.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="grey", size=1) +
  theme_bw() +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size = 14),
        axis.title = element_text(size=16),
        plot.title=element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size=16)) +
  labs(title="Male versus female in vivo") 


#Now do male versus asexual
#Subset male genes
AP2M_HAP2_GGC_subset_asexual <- subset(M_vs_asexual_top_hits,
                                       geneSymbol=="cgd8_2220-RA" |
                                         geneSymbol=="cgd6_2670-RA" |
                                         geneSymbol=="cgd2_2610-RA" |
                                         geneSymbol=="cgd5_3570-RA" |
                                         geneSymbol=="cgd7_5500-RA" |
                                         geneSymbol=="cgd8_1740-RA")

#Plot, Figure 2c
ggplot(M_vs_asexual_top_hits, aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneSymbol))) +
  geom_point(size=4) +
  geom_point(mapping = NULL, AP2M_HAP2_GGC_subset_asexual, size = 4, colour = "#0095FF", inherit.aes = TRUE) +
  ylim(-0.5,12) +
  xlim(-15,15) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", colour="grey", size=0.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="grey", size=1) +
  theme_bw() +
  theme(axis.text=element_text(size = 14),
        axis.title = element_text(size=16),
        plot.title=element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size=16)) +
  labs(title="Male versus asexual") 


# decideTests to pull out the DEGs ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1)

# retrieve expression data for DEGs ----
head(v.DGEList.filtered.norm$E)
colnames(v.DGEList.filtered.norm$E) <- sampleLabels

#This gives all differentially expressed genes (only needs to be significant in at least one comparison)
diffGenes <- v.DGEList.filtered.norm$E[results[,1] !=0 | results[,2] !=0 | results[,3] !=0 | results[,4] !=0 | results[,5] !=0 | results[,6] !=0,]
head(diffGenes)
dim(diffGenes)

#convert DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneSymbol")

#Adding averages to Diff Genes
Differential_genes_avg_all.df <- mutate(diffGenes.df,
                                        asexual_Crypto.AVG = (asexual_CryptoDB_46_rep_1 + asexual_CryptoDB_46_rep_2 + asexual_CryptoDB_46_rep_3 + asexual_CryptoDB_46_rep_4)/4,
                                        female_in_vitro_Crypto.AVG = (female_in_vitro_CryptoDB_46_rep_1 + female_in_vitro_CryptoDB_46_rep_2 + female_in_vitro_CryptoDB_46_rep_3 + female_in_vitro_CryptoDB_46_rep_4)/4,
                                        female_in_vivo_Crypto.AVG = (female_in_vivo_CryptoDB_46_rep_1 + female_in_vivo_CryptoDB_46_rep_2 + female_in_vivo_CryptoDB_46_rep_3 + female_in_vivo_CryptoDB_46_rep_4)/4,
                                        male_Crypto.AVG = (male_CryptoDB_46_rep_1 + male_CryptoDB_46_rep_2 + male_CryptoDB_46_rep_3 + male_CryptoDB_46_rep_4)/4,
                                        
                                        
                                        #now make columns comparing each of the averages above that you're interested in
                                        LogFC.maleCrypto_vs_asexualCrypto = (male_Crypto.AVG - asexual_Crypto.AVG),
                                        LogFC.female_invitro_Crypto_vs_asexualCrypto = (female_in_vitro_Crypto.AVG - asexual_Crypto.AVG),
                                        LogFC.female_invivo_Crypto_vs_asexualCrypto = (female_in_vivo_Crypto.AVG - asexual_Crypto.AVG),
                                        LogFC.maleCrypto_vs_female_invitro_Crypto = (male_Crypto.AVG - female_in_vitro_Crypto.AVG),
                                        LogFC.maleCrypto_vs_female_invivo_Crypto = (male_Crypto.AVG - female_in_vivo_Crypto.AVG),
                                        LogFC.female_invivo_Crypto_vs_female_invitro_Crypto = (female_in_vivo_Crypto.AVG - female_in_vitro_Crypto.AVG)) %>%
  mutate_if(is.numeric, round, 3)                                        

#This was re-named to supplementary table 3, this is the first tab
write_tsv(Differential_genes_avg_all.df, "Differential_genes_CryptoDB_ver46_LFC1.txt")


#Do individually for each spreadsheet of final supplementary file
#Male versus female in vitro, tab 2
diffGenes_male_female_invitro <- v.DGEList.filtered.norm$E[results[,1] !=0,]
head(diffGenes_male_female_invitro)
dim(diffGenes_male_female_invitro)

#convert DEGs to a dataframe using as_tibble
diffGenes_male_female_invitro.df <- as_tibble(diffGenes_male_female_invitro, rownames = "geneSymbol")

Differential_genes_avg_male_female_invitro.df <- dplyr::select(diffGenes_male_female_invitro.df, geneSymbol, 
                                                               male_CryptoDB_46_rep_1:male_CryptoDB_46_rep_4, 
                                                               female_in_vitro_CryptoDB_46_rep_1:female_in_vitro_CryptoDB_46_rep_4) %>%
  
  dplyr::mutate(
    male_Crypto.AVG = (male_CryptoDB_46_rep_1 + male_CryptoDB_46_rep_2 + male_CryptoDB_46_rep_3 + male_CryptoDB_46_rep_4)/4,
    female_in_vitro_Crypto.AVG = (female_in_vitro_CryptoDB_46_rep_1 + female_in_vitro_CryptoDB_46_rep_2 + female_in_vitro_CryptoDB_46_rep_3 + female_in_vitro_CryptoDB_46_rep_4)/4,
    
    LogFC.maleCrypto_vs_female_invitro_Crypto = (male_Crypto.AVG - female_in_vitro_Crypto.AVG)) %>%
  
  mutate_if(is.numeric, round, 3)  


write_tsv(Differential_genes_avg_male_female_invitro.df, "2023_1_3_Differential_genes_CryptoDB_ver46_LFC1_male_female_invitro.txt")


#Male versus female in vivo, tab 3
diffGenes_male_female_invivo <- v.DGEList.filtered.norm$E[results[,2] !=0,]
head(diffGenes_male_female_invivo)
dim(diffGenes_male_female_invivo)

#convert DEGs to a dataframe using as_tibble
diffGenes_male_female_invivo.df <- as_tibble(diffGenes_male_female_invivo, rownames = "geneSymbol")

Differential_genes_avg_male_female_invivo.df <- dplyr::select(diffGenes_male_female_invivo.df, geneSymbol, 
                                                              male_CryptoDB_46_rep_1:male_CryptoDB_46_rep_4, 
                                                              female_in_vivo_CryptoDB_46_rep_1:female_in_vivo_CryptoDB_46_rep_4) %>%
  
  dplyr::mutate(
    male_Crypto.AVG = (male_CryptoDB_46_rep_1 + male_CryptoDB_46_rep_2 + male_CryptoDB_46_rep_3 + male_CryptoDB_46_rep_4)/4,
    female_in_vivo_Crypto.AVG = (female_in_vivo_CryptoDB_46_rep_1 + female_in_vivo_CryptoDB_46_rep_2 + female_in_vivo_CryptoDB_46_rep_3 + female_in_vivo_CryptoDB_46_rep_4)/4,
    
    LogFC.maleCrypto_vs_female_invivo_Crypto = (male_Crypto.AVG - female_in_vivo_Crypto.AVG)) %>%
  
  mutate_if(is.numeric, round, 3)  


write_tsv(Differential_genes_avg_male_female_invivo.df, "2023_1_3_Differential_genes_CryptoDB_ver46_LFC1_male_female_invivo.txt")



#Male versus asexual, tab 4
diffGenes_male_asexual <- v.DGEList.filtered.norm$E[results[,3] !=0,]
head(diffGenes_male_asexual)
dim(diffGenes_male_asexual)

#convert DEGs to a dataframe using as_tibble
diffGenes_male_asexual.df <- as_tibble(diffGenes_male_asexual, rownames = "geneSymbol")

Differential_genes_avg_male_asexual.df <- dplyr::select(diffGenes_male_asexual.df, geneSymbol, 
                                                        male_CryptoDB_46_rep_1:male_CryptoDB_46_rep_4, 
                                                        asexual_CryptoDB_46_rep_1:asexual_CryptoDB_46_rep_4) %>%
  
  dplyr::mutate(
    male_Crypto.AVG = (male_CryptoDB_46_rep_1 + male_CryptoDB_46_rep_2 + male_CryptoDB_46_rep_3 + male_CryptoDB_46_rep_4)/4,
    asexual_Crypto.AVG = (asexual_CryptoDB_46_rep_1 + asexual_CryptoDB_46_rep_2 + asexual_CryptoDB_46_rep_3 + asexual_CryptoDB_46_rep_4)/4,
    
    LogFC.maleCrypto_vs_asexualCrypto = (male_Crypto.AVG - asexual_Crypto.AVG)) %>%
  
  mutate_if(is.numeric, round, 3)  


write_tsv(Differential_genes_avg_male_asexual.df, "2023_1_3_Differential_genes_CryptoDB_ver46_LFC1_male_asexual.txt")



#Female in vitro versus asexual, tab 5
diffGenes_female_in_vitro_asexual <- v.DGEList.filtered.norm$E[results[,4] !=0,]
head(diffGenes_female_in_vitro_asexual)
dim(diffGenes_female_in_vitro_asexual)

#convert DEGs to a dataframe using as_tibble
diffGenes_female_in_vitro_asexual.df <- as_tibble(diffGenes_female_in_vitro_asexual, rownames = "geneSymbol")

Differential_genes_avg_female_in_vitro_asexual.df <- dplyr::select(diffGenes_female_in_vitro_asexual.df, geneSymbol, 
                                                                   female_in_vitro_CryptoDB_46_rep_1:female_in_vitro_CryptoDB_46_rep_4, 
                                                                   asexual_CryptoDB_46_rep_1:asexual_CryptoDB_46_rep_4) %>%
  
  dplyr::mutate(
    female_in_vitro_Crypto.AVG = (female_in_vitro_CryptoDB_46_rep_1 + female_in_vitro_CryptoDB_46_rep_2 + female_in_vitro_CryptoDB_46_rep_3 + female_in_vitro_CryptoDB_46_rep_4)/4,
    asexual_Crypto.AVG = (asexual_CryptoDB_46_rep_1 + asexual_CryptoDB_46_rep_2 + asexual_CryptoDB_46_rep_3 + asexual_CryptoDB_46_rep_4)/4,
    
    LogFC.female_invitro_Crypto_vs_asexualCrypto = (female_in_vitro_Crypto.AVG - asexual_Crypto.AVG)) %>%
  
  mutate_if(is.numeric, round, 3)  


write_tsv(Differential_genes_avg_female_in_vitro_asexual.df, "2023_1_3_Differential_genes_CryptoDB_ver46_LFC1_female_invitro_asexual.txt")



#Female in vivo versus asexual, tab 6
diffGenes_female_in_vivo_asexual <- v.DGEList.filtered.norm$E[results[,5] !=0,]
head(diffGenes_female_in_vivo_asexual)
dim(diffGenes_female_in_vivo_asexual)

#convert DEGs to a dataframe using as_tibble
diffGenes_female_in_vivo_asexual.df <- as_tibble(diffGenes_female_in_vivo_asexual, rownames = "geneSymbol")

Differential_genes_avg_female_in_vivo_asexual.df <- dplyr::select(diffGenes_female_in_vivo_asexual.df, geneSymbol, 
                                                                  female_in_vivo_CryptoDB_46_rep_1:female_in_vivo_CryptoDB_46_rep_4, 
                                                                  asexual_CryptoDB_46_rep_1:asexual_CryptoDB_46_rep_4) %>%
  
  dplyr::mutate(
    female_in_vivo_Crypto.AVG = (female_in_vivo_CryptoDB_46_rep_1 + female_in_vivo_CryptoDB_46_rep_2 + female_in_vivo_CryptoDB_46_rep_3 + female_in_vivo_CryptoDB_46_rep_4)/4,
    asexual_Crypto.AVG = (asexual_CryptoDB_46_rep_1 + asexual_CryptoDB_46_rep_2 + asexual_CryptoDB_46_rep_3 + asexual_CryptoDB_46_rep_4)/4,
    
    LogFC.female_invivo_Crypto_vs_asexualCrypto = (female_in_vivo_Crypto.AVG - asexual_Crypto.AVG)) %>%
  
  mutate_if(is.numeric, round, 3)  


write_tsv(Differential_genes_avg_female_in_vivo_asexual.df, "2023_1_3_Differential_genes_CryptoDB_ver46_LFC1_female_invivo_asexual.txt")



#Female in vivo versus female in vitro, tab 7
diffGenes_female_in_vivo_in_vitro <- v.DGEList.filtered.norm$E[results[,6] !=0,]
head(diffGenes_female_in_vivo_in_vitro)
dim(diffGenes_female_in_vivo_in_vitro)

#convert your DEGs to a dataframe using as_tibble
diffGenes_female_in_vivo_in_vitro.df <- as_tibble(diffGenes_female_in_vivo_in_vitro, rownames = "geneSymbol")

Differential_genes_avg_female_in_vivo_in_vitro.df <- dplyr::select(diffGenes_female_in_vivo_in_vitro.df, geneSymbol, 
                                                                   female_in_vivo_CryptoDB_46_rep_1:female_in_vivo_CryptoDB_46_rep_4, 
                                                                   female_in_vitro_CryptoDB_46_rep_1:female_in_vitro_CryptoDB_46_rep_4) %>%
  
  dplyr::mutate(
    female_in_vivo_Crypto.AVG = (female_in_vivo_CryptoDB_46_rep_1 + female_in_vivo_CryptoDB_46_rep_2 + female_in_vivo_CryptoDB_46_rep_3 + female_in_vivo_CryptoDB_46_rep_4)/4,
    female_in_vitro_Crypto.AVG = (female_in_vitro_CryptoDB_46_rep_1 + female_in_vitro_CryptoDB_46_rep_2 + female_in_vitro_CryptoDB_46_rep_3 + female_in_vitro_CryptoDB_46_rep_4)/4,
    
    LogFC.female_invivo_Crypto_vs_female_invitro_Crypto = (female_in_vivo_Crypto.AVG - female_in_vitro_Crypto.AVG)) %>%
  
  mutate_if(is.numeric, round, 3)  


write_tsv(Differential_genes_avg_female_in_vivo_in_vitro.df, "2023_1_3_Differential_genes_CryptoDB_ver46_LFC1_female_invivo_invitro.txt")


#Get highly expressed genes as bulk markers for male and female
#This is supplementary table 4
#Re-run decideTest to get genes with logFC of 2 or higher
results_high_expressed <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

# take a look at what the results of decideTests looks like
head(results_high_expressed)
summary(results_high_expressed)

#Male versus asexual and female in vivo
#Supplementary table 4 tab 1
diffGenes_male_expressers <- v.DGEList.filtered.norm$E[results_high_expressed[,2] >0 & results_high_expressed[,3] >0,]
head(diffGenes_male_expressers)
dim(diffGenes_male_expressers)

#convert your DEGs to a dataframe using as_tibble
diffGenes_male_expressers.df <- as_tibble(diffGenes_male_expressers, rownames = "geneSymbol")

#Adding averages to Diff Genes
#This is where you would make new files for gene expression log-fold changes
Male_highly_expressed.df <- mutate(diffGenes_male_expressers.df,
                                   asexual_Crypto.AVG = (asexual_CryptoDB_46_rep_1 + asexual_CryptoDB_46_rep_2 + asexual_CryptoDB_46_rep_3 + asexual_CryptoDB_46_rep_4)/4,
                                   female_in_vivo_Crypto.AVG = (female_in_vivo_CryptoDB_46_rep_1 + female_in_vivo_CryptoDB_46_rep_2 + female_in_vivo_CryptoDB_46_rep_3 + female_in_vivo_CryptoDB_46_rep_4)/4,
                                   male_Crypto.AVG = (male_CryptoDB_46_rep_1 + male_CryptoDB_46_rep_2 + male_CryptoDB_46_rep_3 + male_CryptoDB_46_rep_4)/4,
                                   
                                   
                                   #now make columns comparing each of the averages above that you're interested in
                                   LogFC.maleCrypto_vs_asexualCrypto = (male_Crypto.AVG - asexual_Crypto.AVG),
                                   LogFC.maleCrypto_vs_female_invivo_Crypto = (male_Crypto.AVG - female_in_vivo_Crypto.AVG)) %>%
  mutate_if(is.numeric, round, 3)                                        


write_tsv(Male_highly_expressed.df, "2023_1_30_male_expressers_LFC1.txt")


#Female in vivo versus asexual and male
#Supplementary table 4 tab 2
diffGenes_female_expressers <- v.DGEList.filtered.norm$E[results_high_expressed[,2] <0 & results_high_expressed[,5] >0,]
head(diffGenes_female_expressers)
dim(diffGenes_female_expressers)

#convert your DEGs to a dataframe using as_tibble
diffGenes_female_expressers.df <- as_tibble(diffGenes_female_expressers, rownames = "geneSymbol")

#Adding averages to Diff Genes
#This is where you would make new files for gene expression log-fold changes
Female_highly_expressed.df <- mutate(diffGenes_female_expressers.df,
                                     asexual_Crypto.AVG = (asexual_CryptoDB_46_rep_1 + asexual_CryptoDB_46_rep_2 + asexual_CryptoDB_46_rep_3 + asexual_CryptoDB_46_rep_4)/4,
                                     female_in_vivo_Crypto.AVG = (female_in_vivo_CryptoDB_46_rep_1 + female_in_vivo_CryptoDB_46_rep_2 + female_in_vivo_CryptoDB_46_rep_3 + female_in_vivo_CryptoDB_46_rep_4)/4,
                                     male_Crypto.AVG = (male_CryptoDB_46_rep_1 + male_CryptoDB_46_rep_2 + male_CryptoDB_46_rep_3 + male_CryptoDB_46_rep_4)/4,
                                     
                                     
                                     #now make columns comparing each of the averages above that you're interested in
                                     LogFC.female_invivo_Crypto_vs_asexualCrypto = (female_in_vivo_Crypto.AVG - asexual_Crypto.AVG),
                                     LogFC.female_invivo_Crypto_vs_maleCrypto = (female_in_vivo_Crypto.AVG - male_Crypto.AVG)) %>%
  mutate_if(is.numeric, round, 3)                                        

write_tsv(Female_highly_expressed.df, "2023_1_30_female_expressers_LFC1.txt")


