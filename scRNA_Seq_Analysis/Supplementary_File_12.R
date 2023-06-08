#Alignment, filtering, normalization, clustering, visualization, and differential gene expression analysis of entire Cryptosporidium parvum life cycle

#Run Cell Ranger 3.1.0 from the command line to build a reference genome for C. parvum Iowa II version 46 and process sequencing reads

#Making genome reference for Cryptosporidium
#Code I ran to get CryptoDB gff file to gtf file
gffread -E CryptoDB-46_CparvumIowaII.gff -T -o CryptoDB-46_CparvumIowaII.gtf

#Basic command to make reference genome: cellranger mkref --genome=output_genome --fasta=input.fa --genes=input.gtf
cellranger mkref --genome=Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref --fasta=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_ii_EuPathDB46/CryptoDB-46_CparvumIowaII_Genome.fasta --genes=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_ii_EuPathDB46/CryptoDB-46_CparvumIowaII.gtf

#Demultiplex data and align to transcriptome
#For in vitro samples 24 hrs, 36 hrs, 42 hrs, and 54 hrs
cellranger mkfastq --id=10x_fastq_files \
--run=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/2020_2_27_10x_BCL_files \
--csv=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/cellranger-KAW-bcl-simple-2020-2-27.csv \
--localcores=18

#For in vitro sample 46 hrs and in vivo sample
cellranger mkfastq --id=10x_files \
--run=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_July_2020/2020_7_24_10x_BCL_files \
--csv=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_July_2020/cellranger-KAW-bcl-simple-2020-7-24.csv \
--localcores=18

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

cellranger count --id=crypto42_CryptoDB_only_1 \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto42hrs \
--expect-cells=2000 \
--localcores=18

cellranger count --id=crypto54_CryptoDB_only_1 \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/10x_fastq_files/outs/fastq_path/HV2NVBGXC/ \
--sample=crypto54hrs \
--expect-cells=800 \
--localcores=18

cellranger count --id=crypto46hrs_CryptoDB_1 \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_July_2020/10x_files/outs/fastq_path/H23F5BGXG/ \
--sample=crypto46hrs \
--expect-cells=3000 \
--localcores=18

cellranger count --id=in_vivo_parasite_CryptoDB_1 \
--transcriptome=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/ \
--fastqs=/data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_July_2020/10x_files/outs/fastq_path/H23F5BGXG/ \
--sample=infected-in-vivo \
--expect-cells=8000 \
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

#Load the in vivo infected data
in_vivo_parasite_crypto.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_July_2020/In_vivo_parasite_CryptoDB_July_2020/filtered_feature_bc_matrix")
in_vivo_parasite_crypto <- CreateSeuratObject(counts = in_vivo_parasite_crypto.data, project = "in_vivo_crypto")

#Load the Crypto 24 hours crypto only dataset 
crypto_24_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_February_2020/24hrs_CryptoDB_only_analysis/filtered_feature_bc_matrix")
crypto_24_hrs_crypto_only <- CreateSeuratObject(counts = crypto_24_hrs_crypto_only.data, project = "crypto_24hrs")

#Load the Crypto 36 hours crypto only dataset
crypto_36_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_February_2020/36hrs_CryptoDB_only_analysis/filtered_feature_bc_matrix")
crypto_36_hrs_crypto_only <- CreateSeuratObject(counts = crypto_36_hrs_crypto_only.data, project = "crypto_36hrs")

#Load the Crypto 42 hours crypto only dataset
crypto_42_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_February_2020/42hrs_CryptoDB_only_analysis/filtered_feature_bc_matrix")
crypto_42_hrs_crypto_only <- CreateSeuratObject(counts = crypto_42_hrs_crypto_only.data, project = "crypto_42hrs")

#Load the Crypto 46 hours crypto only dataset 
crypto_46_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_July_2020/In_vitro_parasite_CryptoDB_August_2020/46hrs_data/filtered_feature_bc_matrix")
crypto_46_hrs_crypto_only <- CreateSeuratObject(counts = crypto_46_hrs_crypto_only.data, project = "crypto_46hrs")

#Load the Crypto 54 hours crypto only dataset
crypto_54_hrs_crypto_only.data <- Read10X(data.dir = "/Users/katelynwalzer/Desktop/Single_cell_transcriptomics/10x_February_2020/54hrs_CryptoDB_only_analysis/filtered_feature_bc_matrix")
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

#Calculate percentage rRNA gene in each seurat object based on CryptoDB rRNA genes 9/9/20 from Omar
crypto_24_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_24_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
crypto_36_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_36_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
crypto_42_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_42_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
crypto_46_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_46_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
crypto_54_hrs_crypto_only[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_54_hrs_crypto_only , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
in_vivo_parasite_crypto[['percent.rRNA']] <- PercentageFeatureSet(object = in_vivo_parasite_crypto, features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))

#Create list of seurats
list.of.seurats <- c(in_vivo_parasite_crypto, crypto_24_hrs_crypto_only, crypto_36_hrs_crypto_only,
                     crypto_42_hrs_crypto_only, crypto_46_hrs_crypto_only, crypto_54_hrs_crypto_only)


#Merge objects to visualize violin plot of features
merge.seurat.objects <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)], add.cell.ids = c("in_vivo", "24hrs", "36hrs", "42hrs", "46hrs", "54hrs")
)

VlnPlot(merge.seurat.objects, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0)


#Set filtering criteria for seurats
min.features <- 100
max.features.400 <- 400
max.features.500 <- 500
max.features.1200 <- 1200
max.features.1800 <- 1800
max.counts.1000 <- 1000
max.counts.4000 <- 4000
max.counts.7500 <- 7500
max.rRNA <- 60

in_vivo_parasite_crypto <- subset(in_vivo_parasite_crypto, nFeature_RNA > min.features & nFeature_RNA < max.features.1800 & nCount_RNA < max.counts.7500 & percent.rRNA < max.rRNA)

crypto_24_hrs_crypto_only <- subset(crypto_24_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)
crypto_36_hrs_crypto_only <- subset(crypto_36_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)
crypto_42_hrs_crypto_only <- subset(crypto_42_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)
crypto_46_hrs_crypto_only <- subset(crypto_46_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.500 & nCount_RNA < max.counts.1000 & percent.rRNA < max.rRNA)
crypto_54_hrs_crypto_only <- subset(crypto_54_hrs_crypto_only, nFeature_RNA > min.features & nFeature_RNA < max.features.400 & nCount_RNA < max.counts.1000 & percent.rRNA < max.rRNA)

#Create list of seurats after filtering
list.of.seurats <- c(in_vivo_parasite_crypto, crypto_24_hrs_crypto_only, crypto_36_hrs_crypto_only,
                     crypto_42_hrs_crypto_only, crypto_46_hrs_crypto_only, crypto_54_hrs_crypto_only)


#Merge objects to visualize violin plot of features
merge.seurat.objects <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)], add.cell.ids = c("in_vivo", "24hrs", "36hrs", "42hrs", "46hrs", "54hrs")
)

#Plot after filtering
#Extended Data Figure 3
VlnPlot(merge.seurat.objects, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0)

all.genes.crypto <- rownames(merge.seurat.objects)

for (i in 1:length(list.of.seurats)) {
  list.of.seurats[[i]] <- NormalizeData(list.of.seurats[[i]], verbose = FALSE)
  list.of.seurats[[i]] <- FindVariableFeatures(list.of.seurats[[i]], selection.method = "vst",
                                               nfeatures = 2000, verbose = FALSE)
}

VariableFeaturePlot(list.of.seurats[[i]])

#Identify anchors and integrate the data to correct for batch effects and platform differences
seurat.anchors <- FindIntegrationAnchors(object.list = list.of.seurats, anchor.features = all.genes.crypto, dims = 1:30)

seurat.integrated <- IntegrateData(anchorset = seurat.anchors, features.to.integrate = all.genes.crypto, dims = 1:30)

#Switch to integrated assay, not RNA
DefaultAssay(seurat.integrated) <- "integrated"

#Run standard workflow
seurat.integrated <- ScaleData(seurat.integrated, features = all.genes.crypto, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE, npcs = 60)
DimPlot(seurat.integrated, reduction = "pca", group.by = "orig.ident")

DimHeatmap(seurat.integrated, dims = 33:35, cells = 500, balanced = TRUE)
DimHeatmap(seurat.integrated, dims = 35:37, cells = 500, balanced = TRUE)
DimHeatmap(seurat.integrated, dims = 38:40, cells = 500, balanced = TRUE)

seurat.integrated <- JackStraw(seurat.integrated, dims = 60)
seurat.integrated <- ScoreJackStraw(seurat.integrated, dims = 1:60)
JackStrawPlot(seurat.integrated, dims = 1:60)
#I used all crypto genes, integration dims 30, FindNeighbors and UMAP 33 PCs
#PC 33: 4e-61
#PC 34: 1.1e-20
#Used resolution 1.0
#The top gene in PC34 is a ribosomal gene, so I used up to 33 PCs

ElbowPlot(seurat.integrated, ndims = 50)

#Find markers for resolution 1.0 without re-numbering and combining
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:33)
seurat.integrated_1.0 <- FindClusters(seurat.integrated, resolution = 1.0)
head(Idents(seurat.integrated_1.0), 5)

seurat.integrated_1.0 <- RunUMAP(seurat.integrated_1.0, dims = 1:33)
DimPlot(seurat.integrated_1.0, reduction = "umap")
DimPlot(seurat.integrated_1.0, reduction = "umap", group.by = "orig.ident")

DimPlot(seurat.integrated_1.0, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.8) + NoLegend()

#Find all markers integrated 33 PCs resolution 1.0
DefaultAssay(seurat.integrated_1.0) <- "RNA"
seurat.integrated.res1.0.markers <- FindAllMarkers(seurat.integrated_1.0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.res1.0.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write_tsv(seurat.integrated.res1.0.markers, "2020_9_13_invivo_invitro_markers_CryptoDBver46_33PCs_res1.0")

#Rename and combine clusters
#Asexual first, then male, then female
new.cluster.numbers <- c("3", "3", "2", "2", "4", "15", "18", "5", "13", "9", "14", "6",
                         "7", "1", "16", "18", "17", "8", "11", "12", "10", "9", "6")

names(new.cluster.numbers) <- levels(seurat.integrated_1.0)
seurat.integrated_1.0 <- RenameIdents(seurat.integrated_1.0, new.cluster.numbers)

#Define the order of cluster identities
my_levels_combined = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)

#Relevel seurat object
seurat.integrated_1.0 @active.ident <- factor(x = seurat.integrated_1.0@active.ident, levels = my_levels_combined)

#Label new cluster identities and add column to metadata
cell.labels.combined <-seurat.integrated_1.0@active.ident
seurat.integrated_1.0 <- AddMetaData(seurat.integrated_1.0, metadata = cell.labels.combined, col.name = "cluster_labels")

#Rename clusters for start of infection, determined after further investigation
#Asexual first, then male, then female
new.cluster.numbers.infection.start <- c("9", "1", "2", "3", "4", "5", "6", "7", "8", "10", "11", "12",
                                         "13", "14", "15", "16", "17", "18")

names(new.cluster.numbers.infection.start) <- levels(seurat.integrated_1.0)
seurat.integrated_1.0_new_start <- RenameIdents(seurat.integrated_1.0, new.cluster.numbers.infection.start)

#Define the order of cluster identities
my_levels_combined = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)

#Relevel seurat object
seurat.integrated_1.0_new_start @active.ident <- factor(x = seurat.integrated_1.0_new_start@active.ident, levels = my_levels_combined)

#Label new cluster identities and add column to metadata
cell.labels.combined.new <-seurat.integrated_1.0_new_start@active.ident
seurat.integrated_1.0_new_start <- AddMetaData(seurat.integrated_1.0_new_start, metadata = cell.labels.combined.new, col.name = "cluster_labels")

#Find all markers integrated 33 PCs resolution 1.0 new start for 18 clusters, for Supplementary Table 5
DefaultAssay(seurat.integrated_1.0_new_start) <- "RNA"
seurat.integrated.res1.0.markers.new <- FindAllMarkers(seurat.integrated_1.0_new_start, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.integrated.res1.0.markers.new %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#This is Supplementary Table 5
write_tsv(seurat.integrated.res1.0.markers.new, "2022_12_18_invivo_invitro_markers_CryptoDBver46_33PCs_RNA")

#Scale data for RNA assay
seurat.integrated_1.0_new_start_scaled <- ScaleData(seurat.integrated_1.0_new_start, do.scale = TRUE, do.center = TRUE, features = all.genes.crypto, verbose = FALSE)

#Now plot
#Crypto atlas all, Figure 2f
#No legend
DimPlot(seurat.integrated_1.0_new_start_scaled, group.by = "cluster_labels", 
        cols = c("lightgreen", "seagreen1", "seagreen3",
                 "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                 "mediumseagreen", "#00B9E3", "#00ADFA", "#619CFF", "orchid3", "orchid2", "orchid1", "hotpink2",
                 "hotpink1", "lightpink2"), pt.size = 2) +
        FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
        coord_fixed(ratio = 1)

#With legend
DimPlot(seurat.integrated_1.0_new_start_scaled, group.by = "cluster_labels", 
        cols = c("lightgreen", "seagreen1", "seagreen3",
                 "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                 "mediumseagreen", "#00B9E3", "#00ADFA", "#619CFF", "orchid3", "orchid2", "orchid1", "hotpink2",
                 "hotpink1", "lightpink2"), pt.size = 2) +
       FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
       coord_fixed(ratio = 1)

#Now plot each sample separately 
Idents(seurat.integrated_1.0_new_start_scaled) <- "orig.ident"

#Subset each original identity separately to plot UMAP by sample
subset_invivo <- subset(seurat.integrated_1.0_new_start_scaled, idents = "in_vivo_crypto")
subset_24hrs <- subset(seurat.integrated_1.0_new_start_scaled, idents = "crypto_24hrs")
subset_36hrs <- subset(seurat.integrated_1.0_new_start_scaled, idents = "crypto_36hrs")
subset_42hrs <- subset(seurat.integrated_1.0_new_start_scaled, idents = "crypto_42hrs")
subset_46hrs <- subset(seurat.integrated_1.0_new_start_scaled, idents = "crypto_46hrs")
subset_54hrs <- subset(seurat.integrated_1.0_new_start_scaled, idents = "crypto_54hrs")

#Plot each sample separately for Figure 2
#Figure 2a
DimPlot(subset_invivo, reduction = "umap", group.by = "orig.ident", pt.size = 2, cols = "black") + 
  FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + NoLegend() + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#Figure 2b
DimPlot(subset_24hrs, reduction = "umap", group.by = "orig.ident", pt.size = 2, cols = "black") + 
  FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + NoLegend() + ggtitle(NULL) +
  coord_fixed(ratio = 1)

DimPlot(subset_36hrs, reduction = "umap", group.by = "orig.ident", pt.size = 2, cols = "black") + 
  FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + NoLegend() + ggtitle(NULL) +
  coord_fixed(ratio = 1)

DimPlot(subset_42hrs, reduction = "umap", group.by = "orig.ident", pt.size = 2, cols = "black") + 
  FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + NoLegend() + ggtitle(NULL) +
  coord_fixed(ratio = 1)

DimPlot(subset_46hrs, reduction = "umap", group.by = "orig.ident", pt.size = 2, cols = "black") + 
  FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + NoLegend() + ggtitle(NULL) +
  coord_fixed(ratio = 1)

DimPlot(subset_54hrs, reduction = "umap", group.by = "orig.ident", pt.size = 2, cols = "black") + 
  FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + NoLegend() + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#Change identities back
Idents(seurat.integrated_1.0_new_start_scaled) <- "cluster_labels"

#Now paint male and female signatures (Supplementary Table 4) onto UMAP
#Figure 2e
#Plot male signature logFC 2
male.signature.logFC2 <- list(c('cgd1-1140',
                                'cgd1-1150',
                                'cgd1-1260',
                                'cgd1-1360',
                                'cgd1-1370',
                                'cgd1-1830',
                                'cgd1-1860',
                                'cgd1-1880',
                                'cgd1-1890',
                                'cgd1-2200',
                                'cgd1-2240',
                                'cgd1-2470',
                                'cgd1-2550',
                                'cgd1-2980',
                                'cgd1-3150',
                                'cgd1-3200',
                                'cgd1-3280',
                                'cgd1-3290',
                                'cgd1-3830',
                                'cgd1-400',
                                'cgd1-420',
                                'cgd1-710',
                                'cgd1-720',
                                'cgd1-730',
                                'cgd1-960',
                                'cgd2-1230',
                                'cgd2-1300',
                                'cgd2-1410',
                                'cgd2-1610',
                                'cgd2-1620',
                                'cgd2-1920',
                                'cgd2-2030',
                                'cgd2-2390',
                                'cgd2-2490',
                                'cgd2-2533',
                                'cgd2-2800',
                                'cgd2-2860',
                                'cgd2-3000',
                                'cgd2-3210',
                                'cgd2-3560',
                                'cgd2-3870',
                                'cgd2-3880',
                                'cgd2-480',
                                'cgd2-660',
                                'cgd2-670',
                                'cgd2-810',
                                'cgd3-1100',
                                'cgd3-1240',
                                'cgd3-1370',
                                'cgd3-1830',
                                'cgd3-2140',
                                'cgd3-2150',
                                'cgd3-2260',
                                'cgd3-2540',
                                'cgd3-2650',
                                'cgd3-3050',
                                'cgd3-3060',
                                'cgd3-3130',
                                'cgd3-3450',
                                'cgd3-3550',
                                'cgd3-3700',
                                'cgd3-3710',
                                'cgd3-3960',
                                'cgd3-4200',
                                'cgd3-4340',
                                'cgd3-4350',
                                'cgd3-453',
                                'cgd3-520',
                                'cgd3-830',
                                'cgd3-90',
                                'cgd4-1103',
                                'cgd4-1133',
                                'cgd4-123',
                                'cgd4-1230',
                                'cgd4-1283',
                                'cgd4-130',
                                'cgd4-1300',
                                'cgd4-1410',
                                'cgd4-1450',
                                'cgd4-1460',
                                'cgd4-1543',
                                'cgd4-1620',
                                'cgd4-1650',
                                'cgd4-1830',
                                'cgd4-2030',
                                'cgd4-230',
                                'cgd4-2430',
                                'cgd4-2553',
                                'cgd4-2690',
                                'cgd4-2820',
                                'cgd4-2900',
                                'cgd4-3240',
                                'cgd4-3560',
                                'cgd4-3760',
                                'cgd4-3900',
                                'cgd4-400',
                                'cgd4-4370',
                                'cgd4-4380',
                                'cgd5-1423',
                                'cgd5-1510',
                                'cgd5-1530',
                                'cgd5-1650',
                                'cgd5-1810',
                                'cgd5-1830',
                                'cgd5-210',
                                'cgd5-2163',
                                'cgd5-2510',
                                'cgd5-2650',
                                'cgd5-2680',
                                'cgd5-2690',
                                'cgd5-2800',
                                'cgd5-2870',
                                'cgd5-2890',
                                'cgd5-300',
                                'cgd5-3030',
                                'cgd5-3170',
                                'cgd5-3380',
                                'cgd5-3570',
                                'cgd5-3720',
                                'cgd5-3740',
                                'cgd5-4080',
                                'cgd5-4140',
                                'cgd5-4160',
                                'cgd5-4190',
                                'cgd5-4380',
                                'cgd5-4570',
                                'cgd5-580',
                                'cgd5-590',
                                'cgd5-80',
                                'cgd5-810',
                                'cgd5-910',
                                'cgd5-940',
                                'cgd6-1040',
                                'cgd6-1083',
                                'cgd6-1100',
                                'cgd6-1350',
                                'cgd6-1360',
                                'cgd6-1363',
                                'cgd6-160',
                                'cgd6-1623',
                                'cgd6-2120',
                                'cgd6-2670',
                                'cgd6-2680',
                                'cgd6-2710',
                                'cgd6-2790',
                                'cgd6-2913',
                                'cgd6-3060',
                                'cgd6-3600',
                                'cgd6-440',
                                'cgd6-4620',
                                'cgd6-5130',
                                'cgd6-5340',
                                'cgd6-570',
                                'cgd6-900',
                                'cgd6-970',
                                'cgd6-980',
                                'cgd6-990',
                                'cgd7-1150',
                                'cgd7-1290',
                                'cgd7-1350',
                                'cgd7-2280',
                                'cgd7-2320',
                                'cgd7-2340',
                                'cgd7-2350',
                                'cgd7-2460',
                                'cgd7-2470',
                                'cgd7-2483',
                                'cgd7-2490',
                                'cgd7-2510',
                                'cgd7-2820',
                                'cgd7-2860',
                                'cgd7-2870',
                                'cgd7-2960',
                                'cgd7-2970',
                                'cgd7-3090',
                                'cgd7-3400',
                                'cgd7-3690',
                                'cgd7-380',
                                'cgd7-4310',
                                'cgd7-4670',
                                'cgd7-4700',
                                'cgd7-5',
                                'cgd7-5090',
                                'cgd7-5253',
                                'cgd7-5350',
                                'cgd7-540',
                                'cgd7-5500',
                                'cgd7-5510',
                                'cgd7-660',
                                'cgd7-750',
                                'Cgd8-1200',
                                'cgd8-10',
                                'cgd8-1160',
                                'cgd8-1540',
                                'cgd8-1740',
                                'cgd8-1750',
                                'cgd8-2180',
                                'cgd8-2190',
                                'cgd8-2220',
                                'cgd8-2310',
                                'cgd8-240',
                                'cgd8-2680',
                                'cgd8-2710',
                                'cgd8-293',
                                'cgd8-2950',
                                'cgd8-2970',
                                'cgd8-3190',
                                'cgd8-3200',
                                'cgd8-3290',
                                'cgd8-3860',
                                'cgd8-4490',
                                'cgd8-4743',
                                'cgd8-5010',
                                'cgd8-5060',
                                'cgd8-5323',
                                'cgd8-5400',
                                'cgd8-560',
                                'cgd8-820',
                                'cgd8-830',
                                'cgd8-950',
                                'cgd8-960'))

seurat.integrated_1.0_new_start_scaled <- AddModuleScore(
  object = seurat.integrated_1.0_new_start_scaled,
  features = male.signature.logFC2,
  name = 'male_signature_logFC2_'
)

#No legend, Figure 2e
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'male_signature_logFC2_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "blue")) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'male_signature_logFC2_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "blue")) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#Figure 2e
#Plot female signature logFC 2
female.signature.logFC2 <- list(c('cgd1-1100',
                                  'cgd1-1220',
                                  'cgd1-1910',
                                  'cgd1-1920',
                                  'cgd1-2360',
                                  'cgd1-2370',
                                  'cgd1-3060',
                                  'cgd1-3120',
                                  'cgd1-320',
                                  'cgd1-3630',
                                  'cgd1-470',
                                  'cgd1-810',
                                  'cgd1-970',
                                  'cgd1-980',
                                  'cgd1-990',
                                  'cgd2-1090',
                                  'cgd2-1590',
                                  'cgd2-210',
                                  'cgd2-2510',
                                  'cgd2-3250',
                                  'cgd2-3490',
                                  'cgd2-3660',
                                  'cgd2-420',
                                  'cgd2-490',
                                  'cgd2-510',
                                  'cgd2-540',
                                  'cgd2-780',
                                  'cgd2-790',
                                  'cgd3-1400',
                                  'cgd3-150',
                                  'cgd3-1770',
                                  'cgd3-1860',
                                  'cgd3-1880',
                                  'cgd3-330',
                                  'cgd3-3430',
                                  'cgd3-440',
                                  'cgd3-580',
                                  'cgd4-1310',
                                  'cgd4-2110',
                                  'cgd4-2390',
                                  'cgd4-3280',
                                  'cgd4-3350',
                                  'cgd4-3410',
                                  'cgd4-3440',
                                  'cgd4-3830',
                                  'cgd4-4280',
                                  'cgd4-4290',
                                  'cgd4-4400',
                                  'cgd4-4450',
                                  'cgd4-670',
                                  'cgd4-910',
                                  'cgd5-1220',
                                  'cgd5-1240',
                                  'cgd5-2020',
                                  'cgd5-2050',
                                  'cgd5-2130',
                                  'cgd5-2180',
                                  'cgd5-2570',
                                  'cgd5-350',
                                  'cgd5-4440',
                                  'cgd5-700',
                                  'cgd5-743',
                                  'cgd5-760',
                                  'cgd5-790',
                                  'cgd6-1160',
                                  'cgd6-1270',
                                  'cgd6-1450',
                                  'cgd6-2090',
                                  'cgd6-2450',
                                  'cgd6-2610',
                                  'cgd6-3280',
                                  'cgd6-3460',
                                  'cgd6-3610',
                                  'cgd6-3960',
                                  'cgd6-4440',
                                  'cgd6-4510',
                                  'cgd6-4640',
                                  'cgd6-4650',
                                  'cgd6-4800',
                                  'cgd6-4840',
                                  'cgd6-4880',
                                  'cgd6-4950',
                                  'cgd6-4970',
                                  'cgd6-5000',
                                  'cgd6-5040',
                                  'cgd6-5050',
                                  'cgd6-5240',
                                  'cgd6-5500',
                                  'cgd8-20230',
                                  'cgd7-1370',
                                  'cgd7-1730',
                                  'cgd7-1800',
                                  'cgd7-2310',
                                  'cgd7-270',
                                  'cgd7-2740',
                                  'cgd7-2780',
                                  'cgd7-300',
                                  'cgd7-3540',
                                  'cgd7-370',
                                  'cgd7-4030',
                                  'cgd7-4440',
                                  'cgd7-4450',
                                  'cgd7-4530',
                                  'cgd7-4780',
                                  'cgd7-4790',
                                  'cgd7-4810',
                                  'cgd7-4920',
                                  'cgd7-5040',
                                  'cgd7-5150',
                                  'cgd7-5390',
                                  'cgd7-5400',
                                  'cgd7-5440',
                                  'cgd7-5450',
                                  'cgd7-80',
                                  'cgd8-1700',
                                  'cgd8-1720',
                                  'cgd8-1980',
                                  'cgd8-2080',
                                  'cgd8-2360',
                                  'cgd8-2800',
                                  'cgd8-2810',
                                  'cgd8-3270',
                                  'cgd8-3350',
                                  'cgd8-3910',
                                  'cgd8-3950',
                                  'cgd8-40',
                                  'cgd8-4860',
                                  'cgd8-5380',
                                  'cgd8-780',
                                  'cgd8-920'))

seurat.integrated_1.0_new_start_scaled <- AddModuleScore(
  object = seurat.integrated_1.0_new_start_scaled,
  features = female.signature.logFC2,
  name = 'female_signature_logFC2_'
)

#No legend, Figure 2e
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'female_signature_logFC2_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "magenta")) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'female_signature_logFC2_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "magenta")) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#Extended Data Figure 7, plot gene signatures on the life cycle
#Plot glideosome genes on life cycle
glideosome.genes <- list(c("cgd7-3790",
                           "cgd6-2220",
                           "cgd2-640",
                           "cgd6-2210",
                           "cgd6-1500",
                           "cgd7-550",
                           "cgd6-4460",
                           "cgd7-4420"))

seurat.integrated_1.0_new_start_scaled <- AddModuleScore(
  object = seurat.integrated_1.0_new_start_scaled,
  features = glideosome.genes,
  name = 'glideosome_'
)

#No legend, Extended Data Figure 7
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'glideosome_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "green3")) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'glideosome_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "green3")) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Plot rhoptry genes on life cycle
rhoptry.genes <- list(c("cgd4-1790", "cgd8-540", "cgd8-2530", "cgd4-2420", "cgd1-1870", "cgd2-370", 
                        "cgd3-2010",  "cgd1-1380", "cgd1-950", "cgd2-2900", "cgd3-1330", "cgd3-1710", 
                        "cgd3-1730", "cgd3-1770", "cgd3-1780", "cgd3-2750", "cgd3-910", "cgd5-20", 
                        "cgd5-2760", "cgd6-1000", "cgd6-3630", "cgd6-3930", "cgd6-3940", "cgd8-520"))

seurat.integrated_1.0_new_start_scaled <- AddModuleScore(
  object = seurat.integrated_1.0_new_start_scaled,
  features = rhoptry.genes,
  name = 'rhoptry_'
)

#No legend, Extended Data Figure 7
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'rhoptry_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "green3")) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'rhoptry_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "green3")) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Plot microneme genes on life cycle
microneme.genes <- list(c("cgd4-2350",
                          "cgd8-2330",
                          "cgd3-3100",
                          "cgd2-4270",
                          "cgd2-3080",
                          "cgd3-3370",
                          "cgd4-1970",
                          "cgd4-2340",
                          "cgd4-30",
                          "cgd4-32",
                          "cgd4-3620",
                          "cgd4-4503",
                          "cgd5-1520",
                          "cgd5-2820",
                          "cgd5-4470",
                          "cgd6-1080",
                          "cgd6-1660",
                          "cgd6-2330",
                          "cgd6-780",
                          "cgd6-800",
                          "cgd7-1960",
                          "cgd7-2850",
                          "cgd7-4020",
                          "cgd7-4330",
                          "cgd7-5520",
                          "cgd7-5530",
                          "cgd8-3830",
                          "cgd8-530"))

seurat.integrated_1.0_new_start_scaled <- AddModuleScore(
  object = seurat.integrated_1.0_new_start_scaled,
  features = microneme.genes,
  name = 'microneme_'
)

#No legend, Extended Data Figure 7
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'microneme_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "green3")) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, features = 'microneme_1', min.cutoff = 0, order = TRUE, pt.size = 2, cols = c("lightgray", "green3")) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)  



#Make heatmap of all cluster markers, Extended Data Figure 5
DoHeatmap(seurat.integrated_1.0_new_start_scaled, features = seurat.integrated.res1.0.markers.new$gene, 
          group.bar = TRUE, group.colors = c("lightgreen", "seagreen1", "seagreen3",
                                             "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                                             "mediumseagreen", "#00B9E3", "#00ADFA", "#619CFF", "orchid3", 
                                             "orchid2", "orchid1", "hotpink2", "hotpink1", "lightpink2"))

#How many unique markers, even across multiple clusters
length(unique(seurat.integrated.res1.0.markers.new$gene))
#2880 unique genes

#Can also check with total sum and number of duplicates
length(seurat.integrated.res1.0.markers.new$gene)
#8721
sum(duplicated(seurat.integrated.res1.0.markers.new$gene))
#5841



#Make heatmap of all AP2 transcription factors, Extended Data Figure 9
DoHeatmap(seurat.integrated_1.0_new_start_scaled, slot = "data", assay = "RNA", disp.min = 0, disp.max = 3, features = c("cgd8-3230", "cgd6-5320", "cgd3-1980", "cgd5-4250",
                                                                                                                         "cgd5-2570", "cgd4-2950", "cgd3-2970", "cgd6-2670",
                                                                                                                         "cgd2-3490", "cgd4-1110", "cgd8-810", "cgd6-1140",
                                                                                                                         "cgd8-3130", "cgd4-600", "cgd1-3520", "cgd4-3820"),
                                                                                                  group.bar = TRUE, group.colors = c("lightgreen", "seagreen1", "seagreen3",
                                                                                                           "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                                                                                                           "mediumseagreen", "#00B9E3", "#00ADFA", "#619CFF", "orchid3", "orchid2", "orchid1", "hotpink2",
                                                                                                           "hotpink1", "lightpink2")) +
                                                                                                 scale_fill_gradient2(mid = "black", high = rev('yellow'), guide = "colourbar", aesthetics = "fill")



#Make feature plots for marker genes
#MybM cgd6_2250, Figure 4a and Extended Data Figure 10a
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd6-2250", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd6-2250", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Plot the AP2 transcription factors, Extended Data Figure 10a
#AP2 cgd3_2970
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd3-2970", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd3-2970", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#AP2 cgd6_2670 (AP2_M mBio)
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd6-2670", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd6-2670", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#AP2 cgd2_3490
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd2-3490", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd2-3490", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#AP2 cgd4_1110 (AP2_F mBio)
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd4-1110", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd4-1110", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#AP2 cgd8_810
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-810", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-810", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Plot the GGCs, Extended Data Figure 8a
#cgd7_5500 - the one I tagged
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd7-5500", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd7-5500", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd2_2610
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd2-2610", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd2-2610", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd5_3570
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd5-3570", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd5-3570", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd8_1740
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-1740", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-1740", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)



#Plot the female-specific meiosis genes, Extended Data Figure 6
#DMC1 cgd7_1690
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd7-1690", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd7-1690", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd1_60
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd1-60", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd1-60", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd2_510
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd2-510", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd2-510", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd5_1750
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd5-1750", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd5-1750", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd6_4420
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd6-4420", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd6-4420", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd5_2560
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd5-2560", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd5-2560", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Spo11 cgd8_1350
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-1350", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-1350", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd3_4050
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd3-4050", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd3-4050", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Additional meiosis genes identified through inspection of female-specific genes in this study
#cgd1_840
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd1-840", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd1-840", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd7_1000
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd7-1000", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd7-1000", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)


#cgd8_1020
#No legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-1020", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
FeaturePlot(seurat.integrated_1.0_new_start_scaled, slot = "data", features = "cgd8-1020", min.cutoff = 0, order = TRUE, pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + FontSize(x.title = 20, y.title = 20) + xlim(-11, 11)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16)) + ggtitle(NULL) +
  coord_fixed(ratio = 1)



#Add additional labels and categories to samples for comparision and visualization

#In vitro versus in vivo
#This is annotated "Host Model" in final RDS file
#Change cluster identity to orig.ident to then mark in vitro or in vivo
Idents(seurat.integrated_1.0_new_start_scaled) <- "orig.ident"

in_vivo_in_vitro_labels <- c("in_vivo", "in_vitro", "in_vitro", "in_vitro", "in_vitro", "in_vitro")

names(in_vivo_in_vitro_labels) <- levels(seurat.integrated_1.0_new_start_scaled)
seurat.integrated_1.0_new_start_scaled <- RenameIdents(seurat.integrated_1.0_new_start_scaled, in_vivo_in_vitro_labels)
in_vivo_in_vitro_labels <-seurat.integrated_1.0_new_start_scaled@active.ident
seurat.integrated_1.0_new_start_scaled <- AddMetaData(seurat.integrated_1.0_new_start_scaled, metadata = in_vivo_in_vitro_labels, col.name = "in_vivo_in_vitro_labels")

DimPlot(seurat.integrated_1.0_new_start_scaled)

seurat.integrated_1.0_new_start_scaled$in_vivo_in_vitro_labels <- as.character(seurat.integrated_1.0_new_start_scaled$in_vivo_in_vitro_labels)


#Asexual, Male, and Female labels
#This is annotated "Stage" in final RDS file
#Change identity back to cluster to mark each cluster as asexual, male, or female
Idents(seurat.integrated_1.0_new_start_scaled) <- "cluster_labels"

Asex_M_F_labels <- c("Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Asexual", "Male", "Male", "Male",
                     "Female", "Female", "Female", "Female", "Female", "Female")

names(Asex_M_F_labels) <- levels(seurat.integrated_1.0_new_start_scaled)
seurat.integrated_1.0_new_start_scaled <- RenameIdents(seurat.integrated_1.0_new_start_scaled, Asex_M_F_labels)
Asex_M_F_labels <-seurat.integrated_1.0_new_start_scaled@active.ident
seurat.integrated_1.0_new_start_scaled <- AddMetaData(seurat.integrated_1.0_new_start_scaled, metadata = Asex_M_F_labels, col.name = "Asex_M_F_labels")

DimPlot(seurat.integrated_1.0_new_start_scaled)

seurat.integrated_1.0_new_start_scaled$Asex_M_F_labels <- as.character(seurat.integrated_1.0_new_start_scaled$Asex_M_F_labels)

#Change identity back to clusters
Idents(seurat.integrated_1.0_new_start_scaled) <- "cluster_labels"


#Figure out the median gene count in vitro in vivo
#This is included in Supplementary Table 1
#This gives the counts (UMIs)
summary(colSums(seurat.integrated_1.0_new_start_scaled))
summary(colSums(in_vivo_parasite_crypto))
summary(colSums(crypto_24_hrs_crypto_only))
summary(colSums(crypto_36_hrs_crypto_only))
summary(colSums(crypto_42_hrs_crypto_only))
summary(colSums(crypto_46_hrs_crypto_only))
summary(colSums(crypto_54_hrs_crypto_only))

#This gives the genes detected
summary(seurat.integrated_1.0_new_start_scaled@meta.data$nFeature_RNA)
summary(in_vivo_parasite_crypto@meta.data$nFeature_RNA)
summary(crypto_24_hrs_crypto_only@meta.data$nFeature_RNA)
summary(crypto_36_hrs_crypto_only@meta.data$nFeature_RNA)
summary(crypto_42_hrs_crypto_only@meta.data$nFeature_RNA)
summary(crypto_46_hrs_crypto_only@meta.data$nFeature_RNA)
summary(crypto_54_hrs_crypto_only@meta.data$nFeature_RNA)

#Do in vitro only for genes detected
seurat.integrated_1.0_new_start_scaled_invitro <- subset(seurat.integrated_1.0_new_start_scaled, idents = "in_vitro")

summary(seurat.integrated_1.0_new_start_scaled_invitro@meta.data$nFeature_RNA)



#Read in CSV files of male and female cluster markers with p-value < e-50 to determine number of unique enriched genes
#This refers to Supplementary Table 8, male genes
male_enriched_genes <- read.csv("Male_enriched_cluster_markers_e_neg50.csv")

#Gives number of unique genes from all male clusters
nrow(unique(male_enriched_genes))
#389 male-enriched genes

#How many genes are duplicated?
sum(duplicated(male_enriched_genes$Genes))
#69 are duplicated in list of 458 genes

#Alternative way that also lists the 69 duplicated genes
male_enriched_genes$Genes[duplicated(male_enriched_genes$Genes)]


#This refers to Supplementary Table 7, female genes
female_enriched_genes <- read.csv("Female_enriched_cluster_markers_e_neg50.csv")

#Gives number of unique genes from all female clusters
nrow(unique(female_enriched_genes))
#773 female-enriched genes

#How many genes are duplicated?
sum(duplicated(female_enriched_genes$Genes))
#353 are duplicated in list of 1126 genes

#Alternative way that also lists the 353 duplicated genes
female_enriched_genes$Genes[duplicated(female_enriched_genes$Genes)]



#Cluster 18 has some late females along the asexual circle, which affects branching of female pseudotime trajectory in Monocle
#Subset data to only cluster 18 and then filter these cells out
seurat.integrated_1.0_new_start_scaled_c18 <- subset(seurat.integrated_1.0_new_start_scaled, idents = "18")

DimPlot(seurat.integrated_1.0_new_start_scaled_c18)

UMAP_plot <- DimPlot(seurat.integrated_1.0_new_start_scaled_c18)
select.cells <- CellSelector(plot = UMAP_plot)

#Make new seurat object for monocle
seurat.for.monocle <- seurat.integrated_1.0_new_start_scaled

#Identify cells to filter for trajectory analysis from along the circle in cluster 18
#This sets this group of cells as a new active identity
Idents(seurat.for.monocle, cells = select.cells) <- 'FilterCircleCells'

#Now subset only on the active identity of the cluster labels to get rid of "select.cells"
seurat.for.monocle <- subset(seurat.for.monocle, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                                            "10", "11", "12", "13", "14", "15", "16", 
                                                            "17", "18"))

#Double check that cells were filtered, also samples went from 9098 to 8739
table(seurat.for.monocle@meta.data$cluster_labels, seurat.for.monocle@meta.data$orig.ident)
DimPlot(seurat.for.monocle)

#Now save this as an R object to use in Monocle
save(seurat.for.monocle, file = "seurat.for.monocle.Robj")

#RDS file called seurat.object.crypto.atlas.all.rds containing all metadata, including male and female pseudotime, is available through GEO under accession number GSE232438


