#Alignment, filtering, normalization, clustering, visualization, and differential gene expression analysis of asexual in vitro samples

#Run Cell Ranger 3.1.0 from the command line to build a reference genome for C. parvum Iowa II version 46 and process sequencing reads

#Making genome reference for Cryptosporidium
#Code I ran to get CryptoDB gff file to gtf file
gffread -E CryptoDB-46_CparvumIowaII.gff -T -o CryptoDB-46_CparvumIowaII.gtf

#Basic command to make reference genome: cellranger mkref --genome=output_genome --fasta=input.fa --genes=input.gtf
cellranger mkref --genome=Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref --fasta=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_ii_EuPathDB46/CryptoDB-46_CparvumIowaII_Genome.fasta --genes=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_ii_EuPathDB46/CryptoDB-46_CparvumIowaII.gtf

#Demultiplex data and align to transcriptome
cellranger mkfastq --id=10x_fastq_files \
--run=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/2020_2_27_10x_BCL_files \
--csv=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/cellranger-KAW-bcl-simple-2020-2-27.csv \
--localcores=18

#The data is in the 10x_fastq_files folder
#Need to go to outs, then fastq_path, then flow cell, then all samples are listed
#Each sample has the fastq files listed

#Run cellranger count on samples for Crypto only using CryptoDB reference
#Performing single cell gene expression quantification
cellranger count --id=crypto24_CryptoDB_only_1 \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto24hrs \
--expect-cells=4000 \
--localcores=18

cellranger count --id=crypto36_CryptoDB_only_1 \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto36hrs \
--expect-cells=6000 \
--localcores=18
#Can see Cell Ranger output in Supplementary Table 1
#End code in terminal


#Load R packages
library(tidyverse)
library(ggthemes)
library(patchwork)
library(cowplot)
library(Seurat)
library(reticulate)
library(plotly)
library(ggrepel)
library(made4)
library(ggpubr)
library(scater)

#Load the Crypto 24 hours crypto only dataset
crypto_24_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_February_2020/24hrs_CryptoDB_only_analysis/filtered_feature_bc_matrix")
crypto_24_hrs_crypto_only <- CreateSeuratObject(counts = crypto_24_hrs_crypto_only.data, project = "crypto_24hrs")

#Load the Crypto 36 hours crypto only dataset
crypto_36_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_February_2020/36hrs_CryptoDB_only_analysis/filtered_feature_bc_matrix")
crypto_36_hrs_crypto_only <- CreateSeuratObject(counts = crypto_36_hrs_crypto_only.data, project = "crypto_36hrs")

#View features versus counts
FeatureScatter(crypto_24_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_36_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)

crypto_24_hrs_crypto_only$sample <- "24hrs"
crypto_36_hrs_crypto_only$sample <- "36hrs"

#Calculate percentage rRNA gene in each seurat object based on CryptoDB rRNA genes 9/9/20 from Omar
crypto_24_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_24_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
crypto_36_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_36_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))

#Create list of seurats
list.of.seurats <- c(crypto_24_hrs_crypto_only, crypto_36_hrs_crypto_only)

#Merge objects to visualize violin plot of features
merge.seurat.objects <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)], add.cell.ids = c("24hrs", "36hrs")
)

VlnPlot(merge.seurat.objects, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0.1)

#Set filtering criteria for seurats
min.features <- 100
max.features.1200 <- 1200
max.counts.4000 <- 4000
max.rRNA <- 60

crypto_24_hrs_crypto_only <- subset(crypto_24_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)
crypto_36_hrs_crypto_only <- subset(crypto_36_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)

#Create list of seurats after filtering
list.of.seurats <- c(crypto_24_hrs_crypto_only, crypto_36_hrs_crypto_only)

#Merge objects to visualize violin plot of features
merge.seurat.objects <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)], add.cell.ids = c("24hrs", "36hrs")
)

VlnPlot(merge.seurat.objects, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0.1)

all.genes.crypto <- rownames(merge.seurat.objects)

for (i in 1:length(list.of.seurats)) {
  list.of.seurats[[i]] <- NormalizeData(list.of.seurats[[i]], verbose = FALSE)
  list.of.seurats[[i]] <- FindVariableFeatures(list.of.seurats[[i]], selection.method = "vst",
                                               nfeatures = 2000, verbose = FALSE)
}

VariableFeaturePlot(list.of.seurats[[i]])

#Identify anchors and integrate the data 
seurat.anchors <- FindIntegrationAnchors(object.list = list.of.seurats, anchor.features = all.genes.crypto, dims = 1:30)

seurat.integrated <- IntegrateData(anchorset = seurat.anchors, features.to.integrate = all.genes.crypto, dims = 1:30)

#Switch to integrated assay, not RNA
DefaultAssay(seurat.integrated) <- "integrated"

#Run standard workflow
seurat.integrated <- ScaleData(seurat.integrated, features = all.genes.crypto, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE, npcs = 60)
DimPlot(seurat.integrated, reduction = "pca", group.by = "orig.ident")

seurat.integrated <- JackStraw(seurat.integrated, dims = 60)
seurat.integrated <- ScoreJackStraw(seurat.integrated, dims = 1:60)
JackStrawPlot(seurat.integrated, dims = 1:60)
#PC10 has ribosomal genes as drivers but still highly significant p-value
#PC21 has ribosomal genes as drivers then more of a dropoff, use PC20
#PC20: 3.47e−32
#PC21: 1.24e−14

ElbowPlot(seurat.integrated, ndims = 50)

#Find markers for resolution 0.4 without re-numbering and combining
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:20)
seurat.integrated_0.4 <- FindClusters(seurat.integrated, resolution = 0.4)
head(Idents(seurat.integrated_0.4), 5)

seurat.integrated_0.4 <- RunUMAP(seurat.integrated_0.4, dims = 1:20)
DimPlot(seurat.integrated_0.4, reduction = "umap")
DimPlot(seurat.integrated_0.4, reduction = "umap", group.by = "orig.ident")

#Find all markers integration 30, 20 PCs resolution 0.4
DefaultAssay(seurat.integrated_0.4) <- "RNA"
seurat.integrated.res0.4.markers <- FindAllMarkers(seurat.integrated_0.4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.res0.4.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write_tsv(seurat.integrated.res0.4.markers, "2022_6_22_invitro_markers_CryptoDBver46_20PCs_res04")

#Rename and combine clusters
#Asexual clusters, 9, combine two ribosomal clusters since many top markers are the same
new.cluster.numbers <- c("3", "2", "7", "9", "2", "6", "4", "8", "1", "5")

names(new.cluster.numbers) <- levels(seurat.integrated_0.4)
seurat.integrated_0.4 <- RenameIdents(seurat.integrated_0.4, new.cluster.numbers)

#Define the order of cluster identities
my_levels_combined = c(1, 2, 3, 4, 5, 6, 7, 8, 9)

#Relevel seurat object
seurat.integrated_0.4@active.ident <- factor(x = seurat.integrated_0.4@active.ident, levels = my_levels_combined)

#Label new cluster identities and add column to metadata
cell.labels.combined <-seurat.integrated_0.4@active.ident
seurat.integrated_0.4 <- AddMetaData(seurat.integrated_0.4, metadata = cell.labels.combined, col.name = "cluster_labels")

#After more analyses and examination, a new start of the asexual cycle was defined
#Add new start for cluster 1 and cell cycle
new.cluster.numbers.start <- c("9", "1", "2", "3", "4", "5", "6", "7", "8")

names(new.cluster.numbers.start) <- levels(seurat.integrated_0.4)
seurat.integrated_0.4 <- RenameIdents(seurat.integrated_0.4, new.cluster.numbers.start)

#Define the order of cluster identities
my_levels_combined = c(1, 2, 3, 4, 5, 6, 7, 8, 9)

#Relevel seurat object
seurat.integrated_0.4@active.ident <- factor(x = seurat.integrated_0.4@active.ident, levels = my_levels_combined)

#Label new cluster identities and add column to metadata
cell.labels.combined <-seurat.integrated_0.4@active.ident
seurat.integrated_0.4 <- AddMetaData(seurat.integrated_0.4, metadata = cell.labels.combined, col.name = "cluster_labels")

#Plot asexual UMAP with clusters, Figure 1d
#No legend
DimPlot(seurat.integrated_0.4, group.by = "cluster_labels", cols = c("lightgreen", "seagreen1", "seagreen3",
                                                                     "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                                                                     "mediumseagreen"), pt.size = 2) +
                                                            FontSize(x.title = 20, y.title = 20) + NoLegend() + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
                                                            coord_fixed(ratio = 1)

#With legend
DimPlot(seurat.integrated_0.4, group.by = "cluster_labels", cols = c("lightgreen", "seagreen1", "seagreen3",
                                                                     "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                                                                     "mediumseagreen"), pt.size = 2) +
                                                            FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
                                                            coord_fixed(ratio = 1)


#Re-do find markers with new combined clusters in numerical order, for Supplementary Table 2
DefaultAssay(seurat.integrated_0.4) <- "RNA"
seurat.integrated.final.markers <- FindAllMarkers(seurat.integrated_0.4, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.final.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is Supplementary Table 2
write_tsv(seurat.integrated.final.markers, "2022_12_20_invitro_markers_24hrs_36hrs_9_clusters")


#Run marker analysis on 24 versus 36 hours, this does it by entire sample
Idents(seurat.integrated_0.4) <- "orig.ident"

seurat.integrated.markers_24_36 <- FindAllMarkers(seurat.integrated_0.4, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.markers_24_36 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write_tsv(seurat.integrated.markers_24_36, "2022_12_20_invitro_markers_24hrs_vs_36hrs")

#Be sure to also scale data for RNA assay and not just integrated data
seurat.integrated_0.4 <- ScaleData(seurat.integrated_0.4, features = all.genes.crypto, verbose = FALSE)

#Subset the atlas by cluster, then compare between 24 and 36 hours
#This analysis is in Supplementary Table 8
#Cluster 1
cluster1 <- subset(seurat.integrated_0.4, idents = "1")

DimPlot(cluster1, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster1) <- "orig.ident"

DimPlot(cluster1, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster1 <- FindAllMarkers(cluster1, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster1 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster1, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c1")


#Cluster 2
cluster2 <- subset(seurat.integrated_0.4, idents = "2")

DimPlot(cluster2, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster2) <- "orig.ident"

DimPlot(cluster1, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster2 <- FindAllMarkers(cluster2, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster2 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster2, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c2")


#Cluster 3
cluster3 <- subset(seurat.integrated_0.4, idents = "3")

DimPlot(cluster3, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster3) <- "orig.ident"

DimPlot(cluster3, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster3 <- FindAllMarkers(cluster3, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster3 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster3, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c3")


#Cluster 4
cluster4 <- subset(seurat.integrated_0.4, idents = "4")

DimPlot(cluster4, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster4) <- "orig.ident"

DimPlot(cluster4, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster4 <- FindAllMarkers(cluster4, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster4 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster4, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c4")


#Cluster 5
cluster5 <- subset(seurat.integrated_0.4, idents = "5")

DimPlot(cluster5, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster5) <- "orig.ident"

DimPlot(cluster5, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster5 <- FindAllMarkers(cluster5, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster5 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster5, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c5")


#Cluster 6
cluster6 <- subset(seurat.integrated_0.4, idents = "6")

DimPlot(cluster6, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster6) <- "orig.ident"

DimPlot(cluster6, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster6 <- FindAllMarkers(cluster6, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster6 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster6, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c6")


#Cluster 7
cluster7 <- subset(seurat.integrated_0.4, idents = "7")

DimPlot(cluster7, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster7) <- "orig.ident"

DimPlot(cluster7, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster7 <- FindAllMarkers(cluster7, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster7 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster7, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c7")


#Cluster 8
cluster8 <- subset(seurat.integrated_0.4, idents = "8")

DimPlot(cluster8, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster8) <- "orig.ident"

DimPlot(cluster8, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster8 <- FindAllMarkers(cluster8, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster8 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster8, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c8")


#Cluster 9
cluster9 <- subset(seurat.integrated_0.4, idents = "9")

DimPlot(cluster9, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

Idents(cluster9) <- "orig.ident"

DimPlot(cluster9, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.markers.24_36_cluster9 <- FindAllMarkers(cluster9, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers.24_36_cluster9 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is part of Supplementary Table 8
write_tsv(seurat.markers.24_36_cluster9, "2023_12_15_invitro_markers_24hrs_vs_36hrs_c9")



#Save R object, this file was used in RNA velocity and pseudotime analyses (see Supplementary Data 5 for R script)
save(seurat.integrated_0.4, file = "seurat.integrated_0.4.Robj")

#RDS file called seurat.object.crypto.asexual.rds containing all metadata, including pseudotime, is available through GEO under accession number GSE232438
