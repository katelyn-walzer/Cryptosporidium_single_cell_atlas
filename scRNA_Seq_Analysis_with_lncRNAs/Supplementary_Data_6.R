#Alignment, filtering, normalization, clustering, visualization, and differential gene expression analysis of entire Cryptosporidium parvum life cycle including lncRNAs

#10X analysis with genome from Jessie Kissinger's lab
#Contains lncRNA annotations
#Added cryspovirus

#Run Cell Ranger 3.1.0 from the command line to build a reference genome for C. parvum IOWA ATCC and process sequencing reads

#Making genome reference for Cryptosporidium and cryspovirus
cellranger mkref --genome=Cparvum_IOWA_ATCC_with_cryspovirus_April_2021 --fasta=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_genome_with_cryspovirus_March_2021.fa --genes=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/CPATCC_Oct_2020_CSPV_March_2021_fixedIDs.gtf

#Data was already demultiplexed

#Run cellranger count on samples for Crypto only using reference from Jessie and also cryspovirus
cellranger count --id=crypto24hrs_IOWA_ATCC \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_with_cryspovirus_April_2021/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto24hrs \
--expect-cells=4000 \
--chemistry=SC3Pv3 \
--localcores=18

cellranger count --id=crypto36hrs_IOWA_ATCC \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_with_cryspovirus_April_2021/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto36hrs \
--expect-cells=6000 \
--localcores=18

cellranger count --id=crypto42hrs_IOWA_ATCC \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_with_cryspovirus_April_2021/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto42hrs \
--expect-cells=2000 \
--localcores=18

cellranger count --id=crypto46hrs__IOWA_ATCC \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_with_cryspovirus_April_2021/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_July_2020/10x_files/outs/fastq_path/H23F5BGXG/ \
--sample=crypto46hrs \
--expect-cells=3000 \
--localcores=18

cellranger count --id=crypto54hrs_IOWA_ATCC \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_with_cryspovirus_April_2021/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto54hrs \
--expect-cells=800 \
--chemistry=SC3Pv3 \
--localcores=18

cellranger count --id=in_vivo_parasite_IOWA_ATCC \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/10x_with_lncRNAs/Cparvum_IOWA_ATCC_with_cryspovirus_April_2021/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_July_2020/10x_files/outs/fastq_path/H23F5BGXG/ \
--sample=infected-in-vivo \
--expect-cells=8000 \
--localcores=18
#End code in terminal


#Integrate all datasets aligned to Rodrigo and Jessie's IOWA ATCC with lncRNAs and cryspovirus
#Packages:
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

#Realized after starting analysis that gene names (which are actually descriptions) are imported
#Change these all manually in features.tsv and gzip the file OR gene.column specify 1


#Load the in vivo infected data
in_vivo_parasite_crypto.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/Analysis_with_lncRNAs_IOWA_ATCC_reference_March_2021/in_vivo_parasite_IOWA_ATCC/filtered_feature_bc_matrix", gene.column = 1)
in_vivo_parasite_crypto <- CreateSeuratObject(counts = in_vivo_parasite_crypto.data, project = "in_vivo_crypto")

#Load the Crypto 24 hours crypto only dataset from second replicate (February 2020 run)
crypto_24_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/Analysis_with_lncRNAs_IOWA_ATCC_reference_March_2021/24hrs_IOWA_ATCC/filtered_feature_bc_matrix", gene.column = 1)
crypto_24_hrs_crypto_only <- CreateSeuratObject(counts = crypto_24_hrs_crypto_only.data, project = "crypto_24hrs")

#Load the Crypto 36 hours crypto only dataset
crypto_36_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/Analysis_with_lncRNAs_IOWA_ATCC_reference_March_2021/36hrs_IOWA_ATCC/filtered_feature_bc_matrix", gene.column = 1)
crypto_36_hrs_crypto_only <- CreateSeuratObject(counts = crypto_36_hrs_crypto_only.data, project = "crypto_36hrs")

#Load the Crypto 42 hours crypto only dataset
crypto_42_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/Analysis_with_lncRNAs_IOWA_ATCC_reference_March_2021/42hrs_IOWA_ATCC/filtered_feature_bc_matrix", gene.column = 1)
crypto_42_hrs_crypto_only <- CreateSeuratObject(counts = crypto_42_hrs_crypto_only.data, project = "crypto_42hrs")

#Load the Crypto 46 hours crypto only dataset from July 2020 run
crypto_46_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/Analysis_with_lncRNAs_IOWA_ATCC_reference_March_2021/46hrs_IOWA_ATCC/filtered_feature_bc_matrix", gene.column = 1)
crypto_46_hrs_crypto_only <- CreateSeuratObject(counts = crypto_46_hrs_crypto_only.data, project = "crypto_46hrs")

#Load the Crypto 54 hours crypto only dataset
crypto_54_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/Analysis_with_lncRNAs_IOWA_ATCC_reference_March_2021/54hrs_IOWA_ATCC/filtered_feature_bc_matrix", gene.column = 1)
crypto_54_hrs_crypto_only <- CreateSeuratObject(counts = crypto_54_hrs_crypto_only.data, project = "crypto_54hrs")

#View features versus counts
FeatureScatter(in_vivo_parasite_crypto, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_24_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_36_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_42_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_46_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_54_hrs_crypto_only, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)

in_vivo_parasite_crypto$sample <- "in_vivo"
crypto_24_hrs_crypto_only$sample <- "24hrs"
crypto_36_hrs_crypto_only$sample <- "36hrs"
crypto_42_hrs_crypto_only$sample <- "42hrs"
crypto_46_hrs_crypto_only$sample <- "46hrs"
crypto_54_hrs_crypto_only$sample <- "54hrs"

#Calculate percentage rRNA gene in each seurat object based on CryptoDB rRNA search 4/5/2021
#Once I ran this, there are a lot less rRNA reads than other alignment
crypto_24_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_24_hrs_crypto_only , features = c("CPATCC-0010890", "CPATCC-0010900", "CPATCC-0010910", "CPATCC-0029440", "CPATCC-0029470", "CPATCC-0031560", "CPATCC-0031570", "CPATCC-0031580", "CPATCC-0031600", "CPATCC-0039410", "CPATCC-0039420"))
crypto_36_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_36_hrs_crypto_only , features = c("CPATCC-0010890", "CPATCC-0010900", "CPATCC-0010910", "CPATCC-0029440", "CPATCC-0029470", "CPATCC-0031560", "CPATCC-0031570", "CPATCC-0031580", "CPATCC-0031600", "CPATCC-0039410", "CPATCC-0039420"))
crypto_42_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_42_hrs_crypto_only , features = c("CPATCC-0010890", "CPATCC-0010900", "CPATCC-0010910", "CPATCC-0029440", "CPATCC-0029470", "CPATCC-0031560", "CPATCC-0031570", "CPATCC-0031580", "CPATCC-0031600", "CPATCC-0039410", "CPATCC-0039420"))
crypto_46_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_46_hrs_crypto_only , features = c("CPATCC-0010890", "CPATCC-0010900", "CPATCC-0010910", "CPATCC-0029440", "CPATCC-0029470", "CPATCC-0031560", "CPATCC-0031570", "CPATCC-0031580", "CPATCC-0031600", "CPATCC-0039410", "CPATCC-0039420"))
crypto_54_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_54_hrs_crypto_only , features = c("CPATCC-0010890", "CPATCC-0010900", "CPATCC-0010910", "CPATCC-0029440", "CPATCC-0029470", "CPATCC-0031560", "CPATCC-0031570", "CPATCC-0031580", "CPATCC-0031600", "CPATCC-0039410", "CPATCC-0039420"))
in_vivo_parasite_crypto[['percent.rRNA']] <- PercentageFeatureSet(object = in_vivo_parasite_crypto, features = c("CPATCC-0010890", "CPATCC-0010900", "CPATCC-0010910", "CPATCC-0029440", "CPATCC-0029470", "CPATCC-0031560", "CPATCC-0031570", "CPATCC-0031580", "CPATCC-0031600", "CPATCC-0039410", "CPATCC-0039420"))

#Create list of seurats
list.of.seurats <- c(in_vivo_parasite_crypto, crypto_24_hrs_crypto_only, crypto_36_hrs_crypto_only,
                     crypto_42_hrs_crypto_only, crypto_46_hrs_crypto_only, crypto_54_hrs_crypto_only)


#Merge objects to visualize violin plot of features
merge.seurat.objects <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)], add.cell.ids = c("in_vivo", "24hrs", "36hrs", "42hrs", "46hrs", "54hrs")
)

VlnPlot(merge.seurat.objects, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0.1)

#Set filtering criteria for seurats
#Do everything the same as for CryptoDB version 46 except rRNA since the numbers are low here (maybe those cells were already eliminated?)
#Even minimum features of 100 probably isn't necessary
min.features <- 100
max.features.400 <- 400
max.features.500 <- 500
max.features.1200 <- 1200
max.features.1800 <- 1800
max.counts.1000 <- 1000
max.counts.4000 <- 4000
max.counts.7500 <- 7500

in_vivo_parasite_crypto <- subset(in_vivo_parasite_crypto, nFeature_RNA > min.features & nFeature_RNA < max.features.1800 & nCount_RNA < max.counts.7500)

crypto_24_hrs_crypto_only <- subset(crypto_24_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000)
crypto_36_hrs_crypto_only <- subset(crypto_36_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000)
crypto_42_hrs_crypto_only <- subset(crypto_42_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000)
crypto_46_hrs_crypto_only <- subset(crypto_46_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.500 & nCount_RNA < max.counts.1000)
crypto_54_hrs_crypto_only <- subset(crypto_54_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.400 & nCount_RNA < max.counts.1000)

#Create list of seurats after filtering
list.of.seurats <- c(in_vivo_parasite_crypto, crypto_24_hrs_crypto_only, crypto_36_hrs_crypto_only,
                     crypto_42_hrs_crypto_only, crypto_46_hrs_crypto_only, crypto_54_hrs_crypto_only)

#Merge objects to visualize violin plot of features
merge.seurat.objects <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)], add.cell.ids = c("in_vivo", "24hrs", "36hrs", "42hrs", "46hrs", "54hrs")
)

VlnPlot(merge.seurat.objects, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0.1)

rownames(merge.seurat.objects)

#Gets all gene IDs
all.genes <- rownames(merge.seurat.objects)

#Gives only Crypto gene IDs and gets rid of cryspovirus
all.genes.crypto <- all.genes[1:4426]

#Normalize data and find variable features
for (i in 1:length(list.of.seurats)) {
  list.of.seurats[[i]] <- NormalizeData(list.of.seurats[[i]], verbose = FALSE)
  list.of.seurats[[i]] <- FindVariableFeatures(list.of.seurats[[i]], selection.method = "vst",
                                               nfeatures = 2000, verbose = FALSE)
}

VariableFeaturePlot(list.of.seurats[[i]])

#Identify anchors and integrate the data to correct for batch effects and platform differences
#Skip reference list for now
#Do based on all crypto genes here but can also use variable features and/or cryspovirus

seurat.anchors.crypto <- FindIntegrationAnchors(object.list = list.of.seurats, anchor.features = all.genes.crypto, dims = 1:50)

seurat.integrated.crypto <- IntegrateData(anchorset = seurat.anchors.crypto, features.to.integrate = all.genes.crypto, dims = 1:50)

#Switch to integrated assay, not RNA
DefaultAssay(seurat.integrated.crypto) <- "integrated"

#Run standard workflow
seurat.integrated.crypto <- ScaleData(seurat.integrated.crypto, features = all.genes.crypto, verbose = FALSE)


seurat.integrated.crypto <- RunPCA(seurat.integrated.crypto, verbose = FALSE, npcs = 100)
DimPlot(seurat.integrated.crypto, reduction = "pca", group.by = "orig.ident")

seurat.integrated.crypto <- JackStraw(seurat.integrated.crypto, dims = 100)
seurat.integrated.crypto <- ScoreJackStraw(seurat.integrated.crypto, dims = 1:100)
JackStrawPlot(seurat.integrated.crypto, dims = 1:100)
#PC 44 2.79e−13

ElbowPlot(seurat.integrated.crypto, ndims = 50)

seurat.integrated.crypto <- FindNeighbors(seurat.integrated.crypto, dims = 1:44)
seurat.integrated.crypto <- FindClusters(seurat.integrated.crypto, resolution = 1.0)
#res 1.0 is good, dense granules never separate and over 1 is overclustered
head(Idents(seurat.integrated.crypto), 5)

seurat.integrated.crypto <- RunUMAP(seurat.integrated.crypto, dims = 1:44)
DimPlot(seurat.integrated.crypto, reduction = "umap")
DimPlot(seurat.integrated.crypto, reduction = "umap", group.by = "orig.ident")

#Find all markers integrated 44 PCs resolution 1.0
DefaultAssay(seurat.integrated.crypto) <- "RNA"
seurat.integrated.crypto.1.0.markers <- FindAllMarkers(seurat.integrated.crypto, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.crypto.1.0.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write_tsv(seurat.integrated.crypto.1.0.markers, "2023_11_30_invivo_invitro_markers_IOWA_ATCC_44PCs_res1.0")

#Combine and re-order clusters
#Use same order as CryptoDB 46 atlas
seurat.integrated.crypto.combined <- seurat.integrated.crypto
new.cluster.numbers <- c("2", "2", "1", "15", "3", "5", "18", "14", "1", "13", "12",
                         "8", "16", "9", "1", "17", "7", "4", "10", "11", "1", "6", "8", "5")

names(new.cluster.numbers) <- levels(seurat.integrated.crypto.combined)
seurat.integrated.crypto.combined <- RenameIdents(seurat.integrated.crypto.combined, new.cluster.numbers)

DimPlot(seurat.integrated.crypto.combined, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8) + NoLegend()
DimPlot(seurat.integrated.crypto.combined, reduction = "umap")

#Define the order of cluster identities
my_levels_combined = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)

#Relevel seurat object
seurat.integrated.crypto.combined@active.ident <- factor(x = seurat.integrated.crypto.combined@active.ident, levels = my_levels_combined)

#Label new cluster identities and add column to metadata
cell.labels.combined <- seurat.integrated.crypto.combined@active.ident
seurat.integrated.crypto.combined <- AddMetaData(seurat.integrated.crypto.combined, metadata = cell.labels.combined, col.name = "cluster_labels")

DimPlot(seurat.integrated.crypto.combined, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8) + NoLegend()
DimPlot(seurat.integrated.crypto.combined, reduction = "umap")

#Find all markers combined and re-ordered clusters
seurat.integrated.crypto.combined.markers <- FindAllMarkers(seurat.integrated.crypto.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.crypto.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is Supplementary Table 6
write_tsv(seurat.integrated.crypto.combined.markers, "2023_12_1_invivo_invitro_markers_IOWA_ATCC_18_clusters")

#Now plot
#Crypto atlas all, Extended Data Figure 5a
#No legend
DimPlot(seurat.integrated.crypto.combined, group.by = "cluster_labels", 
        cols = c("lightgreen", "seagreen1", "seagreen3",
                 "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                 "mediumseagreen", "#00B9E3", "#00ADFA", "#619CFF", "orchid3", "orchid2", "orchid1", "hotpink2",
                 "hotpink1", "lightpink2"), pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
DimPlot(seurat.integrated.crypto.combined, group.by = "cluster_labels", 
        cols = c("lightgreen", "seagreen1", "seagreen3",
                 "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                 "mediumseagreen", "#00B9E3", "#00ADFA", "#619CFF", "orchid3", "orchid2", "orchid1", "hotpink2",
                 "hotpink1", "lightpink2"), pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#Scale data for RNA assay
seurat.integrated.crypto.combined <- ScaleData(seurat.integrated.crypto.combined, do.scale = TRUE, do.center = TRUE, features = all.genes.crypto, verbose = FALSE)

#Make feature plots for top lncRNA gene
#Extended Data Figure 5b
#CPATCC_0004745
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0004745", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0004745", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0031553
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0031553", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0031553", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0007353
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0007353", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0007353", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0030843
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0030843", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0030843", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0016893
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0016893", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0016893", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0002913
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0002913", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0002913", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0003653
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0003653", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0003653", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0003973
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0003973", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0003973", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0020213
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0020213", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0020213", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0012253
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0012253", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0012253", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0021643
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0021643", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0021643", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0034083
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0034083", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0034083", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0008643
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0008643", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0008643", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0001633
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0001633", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0001633", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0031483
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0031483", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0031483", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0019943
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0019943", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0019943", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0022135
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0022135", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0022135", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#CPATCC_0032625
#No legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0032625", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated.crypto.combined, slot = "data", features = "CPATCC-0032625", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 14)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

seurat.object.crypto.atlas.all.lncRNAs <- seurat.integrated.crypto.combined


#Add additional labels and categories to samples for comparision and visualization

#In vitro versus in vivo
#This is annotated "Host Model" in final RDS file
#Change cluster identity to orig.ident to then mark in vitro or in vivo
Idents(seurat.object.crypto.atlas.all.lncRNAs) <- "orig.ident"

in_vivo_in_vitro_labels <- c("in_vivo", "in_vitro", "in_vitro", "in_vitro", "in_vitro", "in_vitro")

names(in_vivo_in_vitro_labels) <- levels(seurat.object.crypto.atlas.all.lncRNAs)
seurat.object.crypto.atlas.all.lncRNAs <- RenameIdents(seurat.object.crypto.atlas.all.lncRNAs, in_vivo_in_vitro_labels)
in_vivo_in_vitro_labels <-seurat.object.crypto.atlas.all.lncRNAs@active.ident
seurat.object.crypto.atlas.all.lncRNAs <- AddMetaData(seurat.object.crypto.atlas.all.lncRNAs, metadata = in_vivo_in_vitro_labels, col.name = "in_vivo_in_vitro_labels")

DimPlot(seurat.object.crypto.atlas.all.lncRNAs, split.by = "in_vivo_in_vitro_labels")

seurat.object.crypto.atlas.all.lncRNAs$in_vivo_in_vitro_labels <- as.character(seurat.object.crypto.atlas.all.lncRNAs$in_vivo_in_vitro_labels)


#Asexual, Male, and Female labels
#This is annotated "Stage" in final RDS file
#Change identity back to cluster to mark each cluster as asexual, male, or female
Idents(seurat.object.crypto.atlas.all.lncRNAs) <- "cluster_labels"

Asex_M_F_labels <- c("Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Male", "Male", "Male",
                     "Female", "Female", "Female", "Female", "Female", "Female")

names(Asex_M_F_labels) <- levels(seurat.object.crypto.atlas.all.lncRNAs)
seurat.object.crypto.atlas.all.lncRNAs <- RenameIdents(seurat.object.crypto.atlas.all.lncRNAs, Asex_M_F_labels)
Asex_M_F_labels <-seurat.object.crypto.atlas.all.lncRNAs@active.ident
seurat.object.crypto.atlas.all.lncRNAs <- AddMetaData(seurat.object.crypto.atlas.all.lncRNAs, metadata = Asex_M_F_labels, col.name = "Asex_M_F_labels")

DimPlot(seurat.object.crypto.atlas.all.lncRNAs)

seurat.object.crypto.atlas.all.lncRNAs$Asex_M_F_labels <- as.character(seurat.object.crypto.atlas.all.lncRNAs$Asex_M_F_labels)

#Change identity back to clusters
Idents(seurat.object.crypto.atlas.all.lncRNAs) <- "cluster_labels"


#Save .rds file for GitHub and GEO submission
saveRDS(seurat.object.crypto.atlas.all.lncRNAs, file = "seurat.object.crypto.atlas.all.lncRNAs.rds")


#RDS file called seurat.object.crypto.atlas.all.lncRNAs.rds containing all metadata is available through GEO under accession number GSE232438

