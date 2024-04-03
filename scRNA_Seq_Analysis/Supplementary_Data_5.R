#RNA Velocity and Pseudotime Analyses

#RNA velocity of Asexual Cycle
#Run velocyto on the command line to generate loom files with splicing information

#24hrs in vitro 
velocyto run10x /data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/crypto24_CryptoDB_only_1 /data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/genes/genes.gtf

#36hrs in vitro
velocyto run10x /data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Feb_2020/crypto36_CryptoDB_only_1 /data/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_Dec_2019/Cryptosporidium_parvum_iowa_ii_EuPathDB46_Feb2020_ref/genes/genes.gtf
#End python code in terminal


#Packages for RNA velocity analysis
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)

#Read in the loom data for velocyto
#This reads in the spliced/unspliced reads
crypto_24hrs_vel <- ReadVelocity(file = "/venice/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_RNA_velocity/crypto24_CryptoDB_only_1.loom")
crypto_36hrs_vel <- ReadVelocity(file = "/venice/striepenlab/Katelyn_Walzer_data/Single_Cell_Transcriptomics/Striepen_10x_RNA_velocity/crypto36_CryptoDB_only_1.loom")

#This creates a total RNA value from the three bins of RNA
crypto_24hrs_vel$RNA <- (crypto_24hrs_vel$spliced + crypto_24hrs_vel$unspliced + crypto_24hrs_vel$ambiguous)
crypto_36hrs_vel$RNA <- (crypto_36hrs_vel$spliced + crypto_36hrs_vel$unspliced + crypto_36hrs_vel$ambiguous)

#Create Seurat object
crypto_24_hrs_crypto_only_vel <- as.Seurat(x = crypto_24hrs_vel, project = "crypto_24hrs")
crypto_36_hrs_crypto_only_vel <- as.Seurat(x = crypto_36hrs_vel, project = "crypto_36hrs")

#View features versus counts
FeatureScatter(crypto_24_hrs_crypto_only_vel, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)
FeatureScatter(crypto_36_hrs_crypto_only_vel, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)

crypto_24_hrs_crypto_only_vel$sample <- "24hrs"
crypto_36_hrs_crypto_only_vel$sample <- "36hrs"

#Change default assay to RNA
DefaultAssay(crypto_24_hrs_crypto_only_vel) <- "RNA"
DefaultAssay(crypto_36_hrs_crypto_only_vel) <- "RNA"

#Calculate percentage rRNA gene in each seurat object based on CryptoDB rRNA genes 9/9/20 from Omar
crypto_24_hrs_crypto_only_vel[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_24_hrs_crypto_only_vel , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))
crypto_36_hrs_crypto_only_vel[['percent.rRNA']] <- PercentageFeatureSet(object = crypto_36_hrs_crypto_only_vel , features = c("cgd2-1372", "cgd2-1373", "cgd3-665", "cgd3-666", "cgd3-667"))

#Create list of seurats
list.of.seurats.asex.velocity <- c(crypto_24_hrs_crypto_only_vel, crypto_36_hrs_crypto_only_vel)

#Merge objects to visualize violin plot of features
merge.seurat.objects.asex.velocity <- merge(
  x = list.of.seurats.asex.velocity[[1]],
  y = list.of.seurats.asex.velocity[2:length(list.of.seurats.asex.velocity)], add.cell.ids = c("24hrs", "36hrs")
)


VlnPlot(merge.seurat.objects.asex.velocity, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0.1)


#Set filtering criteria for seurats
min.features <- 100
max.features.1200 <- 1200
max.counts.4000 <- 4000
max.rRNA <- 60

crypto_24_hrs_crypto_only_vel <- subset(crypto_24_hrs_crypto_only_vel, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)
crypto_36_hrs_crypto_only_vel <- subset(crypto_36_hrs_crypto_only_vel, nFeature_RNA > min.features & nFeature_RNA < max.features.1200 & nCount_RNA < max.counts.4000 & percent.rRNA < max.rRNA)

#Create list of seurats after filtering
list.of.seurats.asex.velocity <- c(crypto_24_hrs_crypto_only_vel, crypto_36_hrs_crypto_only_vel)


#Merge objects to visualize violin plot of features
merge.seurat.objects.asex.velocity <- merge(
  x = list.of.seurats.asex.velocity[[1]],
  y = list.of.seurats.asex.velocity[2:length(list.of.seurats.asex.velocity)], add.cell.ids = c("24hrs", "36hrs")
)

VlnPlot(merge.seurat.objects.asex.velocity, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0.1)

#Gets all Cryptosporidium gene IDs
all.genes.crypto <- rownames(merge.seurat.objects.asex.velocity)

#Normalize data and find variable features
for (i in 1:length(list.of.seurats.asex.velocity)) {
  list.of.seurats.asex.velocity[[i]] <- NormalizeData(list.of.seurats.asex.velocity[[i]], verbose = FALSE)
  list.of.seurats.asex.velocity[[i]] <- FindVariableFeatures(list.of.seurats.asex.velocity[[i]], selection.method = "vst",
                                                             nfeatures = 2000, verbose = FALSE)
}

VariableFeaturePlot(list.of.seurats.asex.velocity[[i]])

#Identify anchors and integrate the data 
seurat.anchors.asex.vel <- FindIntegrationAnchors(object.list = list.of.seurats.asex.velocity, anchor.features = all.genes.crypto, dims = 1:30)

seurat.integrated.asex.vel <- IntegrateData(anchorset = seurat.anchors.asex.vel, features.to.integrate = all.genes.crypto, dims = 1:30)

#Switch to integrated assay, not RNA
DefaultAssay(seurat.integrated.asex.vel) <- "integrated"

#Run standard workflow
seurat.integrated.asex.vel <- ScaleData(seurat.integrated.asex.vel, features = all.genes.crypto, verbose = FALSE)

seurat.integrated.asex.vel <- RunPCA(seurat.integrated.asex.vel, verbose = FALSE, npcs = 60)
DimPlot(seurat.integrated.asex.vel, reduction = "pca", group.by = "orig.ident")


seurat.integrated.asex.vel <- FindNeighbors(seurat.integrated.asex.vel, dims = 1:20)
seurat.integrated.asex.vel_0.4 <- FindClusters(seurat.integrated.asex.vel, resolution = 0.4)
head(Idents(seurat.integrated.asex.vel_0.4), 5)

seurat.integrated.asex.vel_0.4 <- RunUMAP(seurat.integrated.asex.vel_0.4, dims = 1:20)
DimPlot(seurat.integrated.asex.vel_0.4, reduction = "umap")
DimPlot(seurat.integrated.asex.vel_0.4, reduction = "umap", group.by = "orig.ident")

DimPlot(seurat.integrated.asex.vel_0.4, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")

cell.labels <- seurat.integrated.asex.vel_0.4@active.ident
seurat.integrated.asex.vel_0.4 <- AddMetaData(seurat.integrated.asex.vel_0.4, metadata = cell.labels, col.name = "cluster_labels")


#Now make new Seurat object to run velocity
seurat.integrated.asex.vel_0.4_velocity <- seurat.integrated.asex.vel_0.4

DefaultAssay(seurat.integrated.asex.vel_0.4_velocity) <- "spliced"
seurat.integrated.asex.vel_0.4_velocity <- NormalizeData(object = seurat.integrated.asex.vel_0.4_velocity)
DefaultAssay(seurat.integrated.asex.vel_0.4_velocity) <- "unspliced"
seurat.integrated.asex.vel_0.4_velocity <- NormalizeData(object = seurat.integrated.asex.vel_0.4_velocity)
DefaultAssay(seurat.integrated.asex.vel_0.4_velocity) <- "ambiguous"
seurat.integrated.asex.vel_0.4_velocity <- NormalizeData(object = seurat.integrated.asex.vel_0.4_velocity)
DefaultAssay(seurat.integrated.asex.vel_0.4_velocity) <-"RNA"

#Load the original Seurat object
load("seurat.integrated_0.4.Robj")
DefaultAssay(seurat.integrated_0.4) <-"RNA"
DimPlot(seurat.integrated_0.4, group.by = "sample")


p1 <- DimPlot(seurat.integrated.asex.vel_0.4_velocity, group.by = "sample", pt.size = 2) + ggtitle('New object (with splice info)')
p2 <- DimPlot(seurat.integrated_0.4, group.by = "sample", pt.size = 2) + ggtitle('Original object')
plot_grid(p2, p1)

old.list <- SplitObject(seurat.integrated_0.4, split.by = "sample")
new.list <- SplitObject(seurat.integrated.asex.vel_0.4_velocity, split.by = "sample")

for(current.number in 1:length(old.list)){
  # current.number <- 1
  current.old.seurat <- old.list[[current.number]]
  current.new.seurat <- new.list[[current.number]]
  
  
  #Pull the cell names from seurat.integrated_1.0
  old.cells <- Cells(x = current.old.seurat)
  head(old.cells)
  
  #remove the last digit
  old.cells <- (substr(old.cells, start = 1, stop = 16))
  head(old.cells)
  
  #Pull the cell names from new.integrated.vel
  new.cells <- Cells(x = current.new.seurat)
  head(new.cells)
  
  # cut out the first 28 and the last 1
  new.cells <- (substr(new.cells, start = start.number[current.number], stop = (start.number[current.number] + 15)))
  head(new.cells)
  
  length(old.cells)
  length(new.cells)
  length(intersect(old.cells, new.cells))
  
  head(Cells(x = current.old.seurat))
  old.list[[current.number]] <- RenameCells(current.old.seurat,
                                            new.names = old.cells)
  # head(Cells(x = old.list[[current.number]]))
  
  head(Cells(x = current.new.seurat))
  new.list[[current.number]] <- RenameCells(current.new.seurat,
                                            new.names = new.cells)
  # head(Cells(x = new.list[[current.number]]))
  
  dim(old.list[[current.number]])
  old.list[[current.number]] <- subset(x = old.list[[current.number]], cells = intersect(old.cells, new.cells))
  dim(old.list[[current.number]])
  
  DimPlot(old.list[[current.number]])
  
  
  dim(new.list[[current.number]])
  new.list[[current.number]] <- subset(x = new.list[[current.number]], cells = intersect(old.cells, new.cells))
  dim(new.list[[current.number]])
  
  
  head(Cells(new.list[[current.number]]))
  head(Cells(old.list[[current.number]]))
  
  old.list[[current.number]]@reductions$umap -> new.list[[current.number]]@reductions$umap
  
  old.list[[current.number]]@meta.data[["cluster_labels"]] -> new.list[[current.number]]@meta.data[["cluster_labels"]]
  
  DimPlot(old.list[[current.number]])
  DimPlot(new.list[[current.number]])
  
}

new.object.vel <- merge(x = new.list[[1]], y = new.list[2], merge.dr = c("umap", "pca"), merge.data = T)
p1 <- DimPlot(new.object.vel, group.by = "cluster_labels") + ggtitle('New object (with splice info)')

old.object.seurat <- merge(x = old.list[[1]], y = old.list[2],  merge.dr = c("umap", "pca"), merge.data = T)
p2 <- DimPlot(old.object.seurat, group.by = "cluster_labels") + ggtitle('Original object')
plot_grid(p2, p1)

#Now run velocity
new.object.vel <- RunVelocity(object = new.object.vel, deltaT = 1, kCells = 25, fit.quantile = 0.02)

#Re-label to correct clusters
Idents(object = new.object.vel) <- "cluster_labels"

#Was changed to character vector along the way, want as factor to relevel
new.object.vel$cluster_labels <- as.factor(new.object.vel$cluster_labels)

my_levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9)

#Relevel seurat object
new.object.vel@active.ident <- factor(x = new.object.vel@active.ident, levels = my_levels)


#Make lightgray
lightgray_color <- c("lightgray", "lightgray", "lightgray",
                     "lightgray", "lightgray", "lightgray", "lightgray", "lightgray",
                     "lightgray")


ident.colors <- lightgray_color
names(x = ident.colors) <- levels(x = new.object.vel)
cell.colors <- ident.colors[Idents(object = new.object.vel)]
names(x = cell.colors) <- colnames(x = new.object.vel)

#Plot velocity, Figure 1e
show.velocity.on.embedding.cor(emb = Embeddings(object = new.object.vel, reduction = "umap"), vel = Tool(object = new.object.vel, 
                                                                                                         slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors), 
                               cex = 1, arrow.scale = 1.2, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0, frame = FALSE)




#Packages for pseudotime analysis
library(monocle3)
library(Seurat)
library(ggplot2)
library(Signac)
library(SeuratWrappers)
library(Matrix)
library(patchwork)
library(dplyr)
library(data.table)


#Determine pseudotime for asexual

load("seurat.integrated_0.4.Robj")

seurat.integrated_0.4.updated = UpdateSeuratObject(object = seurat.integrated_0.4)

seurat.integrated_0.4.updated[["UMAP"]] <- seurat.integrated_0.4.updated[["umap"]]

#Convert to cell dataset
cds <- as.cell_data_set(seurat.integrated_0.4.updated)

cds <- cluster_cells(cds, reduction_method = c("UMAP"), k = 20, resolution = 5e-4)

plot_cells(cds, show_trajectory_graph = FALSE)

cds <- learn_graph(cds, 
                   close_loop = T,
                   learn_graph_control=list(ncenter=200, minimal_branch_len=12),
                   use_partition = T)

#Order cells and pick start of trajectory
cds <- order_cells(cds)

#Plot trajectory, Figure 1f
#No legend
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 2, trajectory_graph_color = "red", alpha = 1, trajectory_graph_segment_size = 1.2, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16)) +
  viridis::scale_color_viridis(begin = 0.05, end = 1, option = "D") +
  theme(axis.line.y = element_line(size=0.5, color="black")) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  coord_fixed(ratio = 1) 

#With legend
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 2, trajectory_graph_color = "red", alpha = 1, trajectory_graph_segment_size = 1.2, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE) + FontSize(x.title = 20, y.title = 20) + xlim(-10, 10)+ylim(-10, 10) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16), legend.title = element_blank()) +
  viridis::scale_color_viridis(begin = 0.05, end = 1, option = "D") +
  theme(axis.line.y = element_line(size=0.5, color="black")) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  coord_fixed(ratio = 1) 


pseudotime_start <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
seurat.integrated_0.4.updated <- AddMetaData(seurat.integrated_0.4.updated, metadata = pseudotime_start, col.name = "Pseudotime")

data_to_write_out <- as.data.frame(as.matrix(seurat.integrated_0.4.updated@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "2023_1_19_Asexual_24hrs_36hrs_with_pseudotime.csv")


#Pull out normalized expression data for LOPIT groups and gene lists

DefaultAssay(seurat.integrated_0.4.updated) <- "RNA"

all.genes.crypto <- rownames(seurat.integrated_0.4.updated)

all.parasites <- colnames(seurat.integrated_0.4.updated)


#Rhoptry LOPIT first, Extended Data Figure 3a
list.of.rhoptries <- c("cgd4-1790", "cgd8-540", "cgd8-2530", "cgd4-2420", "cgd1-1870", "cgd2-370", 
                       "cgd3-2010",  "cgd1-1380", "cgd1-950", "cgd2-2900", "cgd3-1330", "cgd3-1710", 
                       "cgd3-1730", "cgd3-1770", "cgd3-1780", "cgd3-2750", "cgd3-910", "cgd5-20", 
                       "cgd5-2760", "cgd6-1000", "cgd6-3630", "cgd6-3930", "cgd6-3940", "cgd8-520")

#Add barcodes and pseudotime as row
rhoptry_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 24))

rownames(rhoptry_pseudotime_results) <- c("cgd4-1790", "cgd8-540", "cgd8-2530", "cgd4-2420", "cgd1-1870", "cgd2-370", 
                                          "cgd3-2010",  "cgd1-1380", "cgd1-950", "cgd2-2900", "cgd3-1330", "cgd3-1710", 
                                          "cgd3-1730", "cgd3-1770", "cgd3-1780", "cgd3-2750", "cgd3-910", "cgd5-20", 
                                          "cgd5-2760", "cgd6-1000", "cgd6-3630", "cgd6-3930", "cgd6-3940", "cgd8-520")

#Run loop to get expression data, I only need to run pseudotime once though
for (current.gene in 1:length(list.of.rhoptries)) {
  gene.name <- list.of.rhoptries[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  all_pseudotime <- as.data.frame(as.matrix(seurat.integrated_0.4.updated@meta.data[["Pseudotime"]]), row.names = colnames(seurat.integrated_0.4.updated))
  colnames(all_pseudotime) <- "Pseudotime"
  print(gene.name)
  as.vector(all_expression_data)
  rhoptry_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(rhoptry_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
rhoptry_pseudotime_results_transposed <- as.data.frame(t(rhoptry_pseudotime_results))

#Scale the data
rhoptry_pseudotime_results_transposed <- as.data.frame(scale(rhoptry_pseudotime_results_transposed))

#Get the mean of the gene group
rhoptry_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(rhoptry_pseudotime_results_transposed))

#Rename first column for rhoptry
colnames(rhoptry_pseudotime_results_transposed_mean) <- "Rhoptry"

#Now add this data to new data frame where I will append the rest of the means as a new column, do pseudotime first
mean_gene_groups_pseudotime <- all_pseudotime
colnames(mean_gene_groups_pseudotime) <- "Pseudotime"

mean_gene_groups_pseudotime$Rhoptry <- rhoptry_pseudotime_results_transposed_mean$Rhoptry

#Add pseudotime to column after finding mean, for plotting
rhoptry_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime


#Rename columns because ggplot does not like dashes
colnames(rhoptry_pseudotime_results_transposed) <- c("cgd4_1790", "cgd8_540", "cgd8_2530", "cgd4_2420", "cgd1_1870", "cgd2_370", 
                                                     "cgd3_2010",  "cgd1_1380", "cgd1_950", "cgd2_2900", "cgd3_1330", "cgd3_1710", 
                                                     "cgd3_1730", "cgd3_1770", "cgd3_1780", "cgd3_2750", "cgd3_910", "cgd5_20", 
                                                     "cgd5_2760", "cgd6_1000", "cgd6_3630", "cgd6_3930", "cgd6_3940", "cgd8_520", "Pseudotime")

#Now plot, Extended Data Figure 3a
ggplot(rhoptry_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd4_1790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_540), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_370), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2010), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_950), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2900), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1710), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1730), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1770), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2750), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_910), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_20), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2760), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1000), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3630), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3930), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3940), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_520), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Rhoptry), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Microneme LOPIT, for Extended Data Figure 3a
list.of.micronemes <- c("cgd4-2350",
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
                        "cgd8-530")

microneme_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 28))

rownames(microneme_pseudotime_results) <- c("cgd4-2350",
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
                                            "cgd8-530")

#Run loop to get expression data
for (current.gene in 1:length(list.of.micronemes)) {
  gene.name <- list.of.micronemes[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  microneme_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(microneme_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
microneme_pseudotime_results_transposed <- as.data.frame(t(microneme_pseudotime_results))

#Scale the data
microneme_pseudotime_results_transposed <- as.data.frame(scale(microneme_pseudotime_results_transposed))

#Get the mean of the gene group
microneme_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(microneme_pseudotime_results_transposed))

#Rename first column for microneme
colnames(microneme_pseudotime_results_transposed_mean) <- "Microneme"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Microneme <- microneme_pseudotime_results_transposed_mean$Microneme

#Add pseudotime to column after finding mean, for plotting
microneme_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(microneme_pseudotime_results_transposed) <- c("cgd4_2350",
                                                       "cgd8_2330",
                                                       "cgd3_3100",
                                                       "cgd2_4270",
                                                       "cgd2_3080",
                                                       "cgd3_3370",
                                                       "cgd4_1970",
                                                       "cgd4_2340",
                                                       "cgd4_30",
                                                       "cgd4_32",
                                                       "cgd4_3620",
                                                       "cgd4_4503",
                                                       "cgd5_1520",
                                                       "cgd5_2820",
                                                       "cgd5_4470",
                                                       "cgd6_1080",
                                                       "cgd6_1660",
                                                       "cgd6_2330",
                                                       "cgd6_780",
                                                       "cgd6_800",
                                                       "cgd7_1960",
                                                       "cgd7_2850",
                                                       "cgd7_4020",
                                                       "cgd7_4330",
                                                       "cgd7_5520",
                                                       "cgd7_5530",
                                                       "cgd8_3830",
                                                       "cgd8_530", 
                                                       "Pseudotime")

#Now plot, Extended Data Figure 3a
ggplot(microneme_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd4_2350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3100), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_4270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3370), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_30), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_32), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_4503), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1520), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_4470), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1660), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_800), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1960), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4020), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5520), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Microneme), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Dense granule LOPIT, Figure 4j, Extended Data 3a
list.of.dense.granules <- c("cgd1-1230",
                            "cgd1-1250",
                            "cgd1-1433",
                            "cgd1-3780",
                            "cgd1-3790",
                            "cgd1-590",
                            "cgd1-620",
                            "cgd2-2530",
                            "cgd2-2580",
                            "cgd2-3030",
                            "cgd2-3210",
                            "cgd2-3360",
                            "cgd2-430",
                            "cgd3-1690",
                            "cgd3-1700",
                            "cgd3-1750",
                            "cgd3-530",
                            "cgd3-560",
                            "cgd3-600",
                            "cgd3-650",
                            "cgd3-693",
                            "cgd4-230",
                            "cgd4-2510",
                            "cgd4-3460",
                            "cgd4-3530",
                            "cgd4-3580",
                            "cgd4-3970",
                            "cgd4-850",
                            "cgd5-1440",
                            "cgd5-1450",
                            "cgd5-1480",
                            "cgd5-1490",
                            "cgd5-1530",
                            "cgd5-3210",
                            "cgd6-1030",
                            "cgd6-1110",
                            "cgd6-1170",
                            "cgd6-1180",
                            "cgd6-3050",
                            "cgd6-3080",
                            "cgd6-5280",
                            "cgd6-5310",
                            "cgd6-5330",
                            "cgd6-5360",
                            "cgd6-5420",
                            "cgd6-5440",
                            "cgd6-60",
                            "cgd7-1270",
                            "cgd7-1300",
                            "cgd7-1340",
                            "cgd7-2070",
                            "Cgd7-2213",
                            "cgd7-2340",
                            "cgd7-2620",
                            "cgd7-3800",
                            "cgd7-3810",
                            "cgd7-3820",
                            "cgd7-3860",
                            "cgd7-3870",
                            "cgd7-4280",
                            "cgd7-4340",
                            "cgd7-4420",
                            "cgd7-4490",
                            "cgd7-4500",
                            "cgd7-4530",
                            "cgd7-4680",
                            "cgd8-1780",
                            "cgd8-1800",
                            "cgd8-2160",
                            "cgd8-3030",
                            "cgd8-3600",
                            "cgd8-3610",
                            "cgd8-5300",
                            "cgd8-5310",
                            "cgd8-5350",
                            "cgd8-5360",
                            "cgd8-5380",
                            "cgd8-680",
                            "cgd8-690")


dense_granule_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 79))

rownames(dense_granule_pseudotime_results) <- c("cgd1-1230",
                                                "cgd1-1250",
                                                "cgd1-1433",
                                                "cgd1-3780",
                                                "cgd1-3790",
                                                "cgd1-590",
                                                "cgd1-620",
                                                "cgd2-2530",
                                                "cgd2-2580",
                                                "cgd2-3030",
                                                "cgd2-3210",
                                                "cgd2-3360",
                                                "cgd2-430",
                                                "cgd3-1690",
                                                "cgd3-1700",
                                                "cgd3-1750",
                                                "cgd3-530",
                                                "cgd3-560",
                                                "cgd3-600",
                                                "cgd3-650",
                                                "cgd3-693",
                                                "cgd4-230",
                                                "cgd4-2510",
                                                "cgd4-3460",
                                                "cgd4-3530",
                                                "cgd4-3580",
                                                "cgd4-3970",
                                                "cgd4-850",
                                                "cgd5-1440",
                                                "cgd5-1450",
                                                "cgd5-1480",
                                                "cgd5-1490",
                                                "cgd5-1530",
                                                "cgd5-3210",
                                                "cgd6-1030",
                                                "cgd6-1110",
                                                "cgd6-1170",
                                                "cgd6-1180",
                                                "cgd6-3050",
                                                "cgd6-3080",
                                                "cgd6-5280",
                                                "cgd6-5310",
                                                "cgd6-5330",
                                                "cgd6-5360",
                                                "cgd6-5420",
                                                "cgd6-5440",
                                                "cgd6-60",
                                                "cgd7-1270",
                                                "cgd7-1300",
                                                "cgd7-1340",
                                                "cgd7-2070",
                                                "Cgd7-2213",
                                                "cgd7-2340",
                                                "cgd7-2620",
                                                "cgd7-3800",
                                                "cgd7-3810",
                                                "cgd7-3820",
                                                "cgd7-3860",
                                                "cgd7-3870",
                                                "cgd7-4280",
                                                "cgd7-4340",
                                                "cgd7-4420",
                                                "cgd7-4490",
                                                "cgd7-4500",
                                                "cgd7-4530",
                                                "cgd7-4680",
                                                "cgd8-1780",
                                                "cgd8-1800",
                                                "cgd8-2160",
                                                "cgd8-3030",
                                                "cgd8-3600",
                                                "cgd8-3610",
                                                "cgd8-5300",
                                                "cgd8-5310",
                                                "cgd8-5350",
                                                "cgd8-5360",
                                                "cgd8-5380",
                                                "cgd8-680",
                                                "cgd8-690")

#Run loop to get expression data
for (current.gene in 1:length(list.of.dense.granules)) {
  gene.name <- list.of.dense.granules[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  dense_granule_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(dense_granule_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
dense_granule_pseudotime_results_transposed <- as.data.frame(t(dense_granule_pseudotime_results))

#Scale the data
dense_granule_pseudotime_results_transposed <- as.data.frame(scale(dense_granule_pseudotime_results_transposed))

#Get the mean of the gene group
dense_granule_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(dense_granule_pseudotime_results_transposed))

#Rename first column for dense granule
colnames(dense_granule_pseudotime_results_transposed_mean) <- "Dense_Granule"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Dense_Granule <- dense_granule_pseudotime_results_transposed_mean$Dense_Granule

#Add pseudotime to column after finding mean, for plotting
dense_granule_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(dense_granule_pseudotime_results_transposed) <- c("cgd1_1230",
                                                           "cgd1_1250",
                                                           "cgd1_1433",
                                                           "cgd1_3780",
                                                           "cgd1_3790",
                                                           "cgd1_590",
                                                           "cgd1_620",
                                                           "cgd2_2530",
                                                           "cgd2_2580",
                                                           "cgd2_3030",
                                                           "cgd2_3210",
                                                           "cgd2_3360",
                                                           "cgd2_430",
                                                           "cgd3_1690",
                                                           "cgd3_1700",
                                                           "cgd3_1750",
                                                           "cgd3_530",
                                                           "cgd3_560",
                                                           "cgd3_600",
                                                           "cgd3_650",
                                                           "cgd3_693",
                                                           "cgd4_230",
                                                           "cgd4_2510",
                                                           "cgd4_3460",
                                                           "cgd4_3530",
                                                           "cgd4_3580",
                                                           "cgd4_3970",
                                                           "cgd4_850",
                                                           "cgd5_1440",
                                                           "cgd5_1450",
                                                           "cgd5_1480",
                                                           "cgd5_1490",
                                                           "cgd5_1530",
                                                           "cgd5_3210",
                                                           "cgd6_1030",
                                                           "cgd6_1110",
                                                           "cgd6_1170",
                                                           "cgd6_1180",
                                                           "cgd6_3050",
                                                           "cgd6_3080",
                                                           "cgd6_5280",
                                                           "cgd6_5310",
                                                           "cgd6_5330",
                                                           "cgd6_5360",
                                                           "cgd6_5420",
                                                           "cgd6_5440",
                                                           "cgd6_60",
                                                           "cgd7_1270",
                                                           "cgd7_1300",
                                                           "cgd7_1340",
                                                           "cgd7_2070",
                                                           "Cgd7_2213",
                                                           "cgd7_2340",
                                                           "cgd7_2620",
                                                           "cgd7_3800",
                                                           "cgd7_3810",
                                                           "cgd7_3820",
                                                           "cgd7_3860",
                                                           "cgd7_3870",
                                                           "cgd7_4280",
                                                           "cgd7_4340",
                                                           "cgd7_4420",
                                                           "cgd7_4490",
                                                           "cgd7_4500",
                                                           "cgd7_4530",
                                                           "cgd7_4680",
                                                           "cgd8_1780",
                                                           "cgd8_1800",
                                                           "cgd8_2160",
                                                           "cgd8_3030",
                                                           "cgd8_3600",
                                                           "cgd8_3610",
                                                           "cgd8_5300",
                                                           "cgd8_5310",
                                                           "cgd8_5350",
                                                           "cgd8_5360",
                                                           "cgd8_5380",
                                                           "cgd8_680",
                                                           "cgd8_690",
                                                           "Pseudotime")

#Now plot, Figure 1j, Extended Data 3a
ggplot(dense_granule_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1433), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_590), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1700), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1750), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_560), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_600), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_650), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_693), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2510), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3460), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1480), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1180), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3050), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_60), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=Cgd7_2213), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3800), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3810), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3860), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4500), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1800), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3600), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3610), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Dense_Granule), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Small Granule LOPIT, Extended Data Figure 3a
list.of.small.granules <- c("cgd7-200",
                            "cgd7-3390",
                            "cgd1-3640",
                            "cgd1-3650",
                            "cgd1-3810",
                            "cgd1-640",
                            "cgd2-3110",
                            "cgd2-340",
                            "cgd2-3730",
                            "cgd2-4020",
                            "cgd3-4250",
                            "cgd4-2460",
                            "cgd4-3050",
                            "cgd5-2720",
                            "cgd5-2730",
                            "cgd5-2740",
                            "cgd5-320",
                            "cgd6-3070",
                            "cgd6-5380",
                            "cgd7-3830",
                            "cgd8-1770",
                            "cgd8-2490",
                            "cgd8-3590")

small_granule_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 23))

rownames(small_granule_pseudotime_results) <- c("cgd7-200",
                                                "cgd7-3390",
                                                "cgd1-3640",
                                                "cgd1-3650",
                                                "cgd1-3810",
                                                "cgd1-640",
                                                "cgd2-3110",
                                                "cgd2-340",
                                                "cgd2-3730",
                                                "cgd2-4020",
                                                "cgd3-4250",
                                                "cgd4-2460",
                                                "cgd4-3050",
                                                "cgd5-2720",
                                                "cgd5-2730",
                                                "cgd5-2740",
                                                "cgd5-320",
                                                "cgd6-3070",
                                                "cgd6-5380",
                                                "cgd7-3830",
                                                "cgd8-1770",
                                                "cgd8-2490",
                                                "cgd8-3590")

#Run loop to get expression data
for (current.gene in 1:length(list.of.small.granules)) {
  gene.name <- list.of.small.granules[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  small_granule_pseudotime_results[current.gene,] <- all_expression_data
  
}


#Get barcode cell labels as column names
colnames(small_granule_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
small_granule_pseudotime_results_transposed <- as.data.frame(t(small_granule_pseudotime_results))

#Scale the data
small_granule_pseudotime_results_transposed <- as.data.frame(scale(small_granule_pseudotime_results_transposed))

#Get the mean of the gene group
small_granule_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(small_granule_pseudotime_results_transposed))

#Rename column for small granule
colnames(small_granule_pseudotime_results_transposed_mean) <- "Small_Granule"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Small_Granule <- small_granule_pseudotime_results_transposed_mean$Small_Granule

#Add pseudotime to column after finding mean, for plotting
small_granule_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(small_granule_pseudotime_results_transposed) <- c("cgd7_200",
                                                           "cgd7_3390",
                                                           "cgd1_3640",
                                                           "cgd1_3650",
                                                           "cgd1_3810",
                                                           "cgd1_640",
                                                           "cgd2_3110",
                                                           "cgd2_340",
                                                           "cgd2_3730",
                                                           "cgd2_4020",
                                                           "cgd3_4250",
                                                           "cgd4_2460",
                                                           "cgd4_3050",
                                                           "cgd5_2720",
                                                           "cgd5_2730",
                                                           "cgd5_2740",
                                                           "cgd5_320",
                                                           "cgd6_3070",
                                                           "cgd6_5380",
                                                           "cgd7_3830",
                                                           "cgd8_1770",
                                                           "cgd8_2490",
                                                           "cgd8_3590",
                                                           "Pseudotime")


#Now plot, Extended Data Figure 3a
ggplot(small_granule_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd7_200), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3390), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3640), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3650), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3810), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_640), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3730), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_4020), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_4250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2460), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3050), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2730), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2740), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_320), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1770), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3590), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Small_Granule), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#20S Proteasome LOPIT, Extended Data Figure 3a
list.of.20S.proteasome <- c("cgd1-2490",
                            "cgd1-420",
                            "cgd2-1440",
                            "cgd2-2050",
                            "cgd2-530",
                            "cgd2-860",
                            "cgd3-2170",
                            "cgd3-2200",
                            "cgd3-2530",
                            "cgd4-250",
                            "cgd5-1820",
                            "cgd5-3220",
                            "cgd5-4210",
                            "cgd7-3660")

proteasome20S_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 14))

rownames(proteasome20S_pseudotime_results) <- c("cgd1-2490",
                                                "cgd1-420",
                                                "cgd2-1440",
                                                "cgd2-2050",
                                                "cgd2-530",
                                                "cgd2-860",
                                                "cgd3-2170",
                                                "cgd3-2200",
                                                "cgd3-2530",
                                                "cgd4-250",
                                                "cgd5-1820",
                                                "cgd5-3220",
                                                "cgd5-4210",
                                                "cgd7-3660")

#Run loop to get expression data
for (current.gene in 1:length(list.of.20S.proteasome)) {
  gene.name <- list.of.20S.proteasome[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  proteasome20S_pseudotime_results[current.gene,] <- all_expression_data
  
}


#Get barcode cell labels as column names
colnames(proteasome20S_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
proteasome20S_pseudotime_results_transposed <- as.data.frame(t(proteasome20S_pseudotime_results))

#Scale the data
proteasome20S_pseudotime_results_transposed <- as.data.frame(scale(proteasome20S_pseudotime_results_transposed))

#Get the mean of the gene group
proteasome20S_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(proteasome20S_pseudotime_results_transposed))

#Rename first column for 20S Proteasome
colnames(proteasome20S_pseudotime_results_transposed_mean) <- "Proteasome_20S"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Proteasome_20S <- proteasome20S_pseudotime_results_transposed_mean$Proteasome_20S

#Add pseudotime to column after finding mean, for plotting
proteasome20S_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(proteasome20S_pseudotime_results_transposed) <- c("cgd1_2490",
                                                           "cgd1_420",
                                                           "cgd2_1440",
                                                           "cgd2_2050",
                                                           "cgd2_530",
                                                           "cgd2_860",
                                                           "cgd3_2170",
                                                           "cgd3_2200",
                                                           "cgd3_2530",
                                                           "cgd4_250",
                                                           "cgd5_1820",
                                                           "cgd5_3220",
                                                           "cgd5_4210",
                                                           "cgd7_3660",
                                                           "Pseudotime")


#Now plot, Extended Data Figure 3a
ggplot(proteasome20S_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_2490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2050), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_860), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2200), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3220), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_4210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3660), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Proteasome_20S), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#40S Ribosome LOPIT, Extended Data Figure 3a
list.of.40S.ribosome <- c("cgd1-850",
                          "cgd2-1070",
                          "cgd2-170",
                          "cgd2-3000",
                          "cgd2-4260",
                          "cgd3-2440",
                          "cgd4-3160",
                          "cgd4-4020",
                          "cgd5-2210",
                          "cgd5-3040",
                          "cgd5-3720",
                          "cgd6-1390",
                          "cgd6-3180",
                          "cgd6-3710",
                          "cgd6-4320",
                          "cgd6-4630",
                          "cgd7-2250",
                          "cgd7-4760",
                          "cgd8-1840",
                          "cgd8-4360")

ribosome40S_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 20))

rownames(ribosome40S_pseudotime_results) <- c("cgd1-850",
                                              "cgd2-1070",
                                              "cgd2-170",
                                              "cgd2-3000",
                                              "cgd2-4260",
                                              "cgd3-2440",
                                              "cgd4-3160",
                                              "cgd4-4020",
                                              "cgd5-2210",
                                              "cgd5-3040",
                                              "cgd5-3720",
                                              "cgd6-1390",
                                              "cgd6-3180",
                                              "cgd6-3710",
                                              "cgd6-4320",
                                              "cgd6-4630",
                                              "cgd7-2250",
                                              "cgd7-4760",
                                              "cgd8-1840",
                                              "cgd8-4360")

#Run loop to get expression data
for (current.gene in 1:length(list.of.40S.ribosome)) {
  gene.name <- list.of.40S.ribosome[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  ribosome40S_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(ribosome40S_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
ribosome40S_pseudotime_results_transposed <- as.data.frame(t(ribosome40S_pseudotime_results))

#Scale the data
ribosome40S_pseudotime_results_transposed <- as.data.frame(scale(ribosome40S_pseudotime_results_transposed))

#Get the mean of the gene group
ribosome40S_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(ribosome40S_pseudotime_results_transposed))

#Rename first column for 40S Ribosome
colnames(ribosome40S_pseudotime_results_transposed_mean) <- "Ribosome_40S"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Ribosome_40S <- ribosome40S_pseudotime_results_transposed_mean$Ribosome_40S

#Add pseudotime to column after finding mean, for plotting
ribosome40S_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(ribosome40S_pseudotime_results_transposed) <- c("cgd1_850",
                                                         "cgd2_1070",
                                                         "cgd2_170",
                                                         "cgd2_3000",
                                                         "cgd2_4260",
                                                         "cgd3_2440",
                                                         "cgd4_3160",
                                                         "cgd4_4020",
                                                         "cgd5_2210",
                                                         "cgd5_3040",
                                                         "cgd5_3720",
                                                         "cgd6_1390",
                                                         "cgd6_3180",
                                                         "cgd6_3710",
                                                         "cgd6_4320",
                                                         "cgd6_4630",
                                                         "cgd7_2250",
                                                         "cgd7_4760",
                                                         "cgd8_1840",
                                                         "cgd8_4360",
                                                         "Pseudotime")


#Now plot, Extended Data Figure 3a
ggplot(ribosome40S_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3000), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_4260), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_4020), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3040), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1390), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3180), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3710), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4320), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4630), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4760), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1840), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Ribosome_40S), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#60S Ribosome LOPIT, Figure 1i, Extended Data Figure 3a
list.of.60S.ribosome <- c("cgd1-1660",
                          "cgd1-3000",
                          "cgd2-120",
                          "cgd2-2200",
                          "cgd2-280",
                          "Cgd2-2990",
                          "cgd3-1250",
                          "cgd3-1300",
                          "cgd3-2250",
                          "cgd3-300",
                          "cgd3-3790",
                          "cgd3-3890",
                          "cgd3-3930",
                          "cgd3-830",
                          "cgd4-1230",
                          "cgd4-2260",
                          "cgd4-2400",
                          "cgd4-470",
                          "cgd4-840",
                          "cgd5-1580",
                          "cgd5-3823",
                          "cgd5-970",
                          "cgd6-2460",
                          "cgd6-3190",
                          "cgd6-3340",
                          "cgd6-4190",
                          "cgd6-4620",
                          "cgd6-570",
                          "cgd7-1873",
                          "cgd7-2110",
                          "cgd7-2420",
                          "cgd7-2540",
                          "cgd7-320",
                          "cgd7-4050",
                          "cgd7-4460",
                          "cgd8-2340",
                          "cgd8-2870",
                          "cgd8-3450",
                          "cgd8-3480",
                          "cgd8-430",
                          "cgd8-4350",
                          "cgd8-440")


ribosome60S_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 42))

rownames(ribosome60S_pseudotime_results) <- c("cgd1-1660",
                                              "cgd1-3000",
                                              "cgd2-120",
                                              "cgd2-2200",
                                              "cgd2-280",
                                              "Cgd2-2990",
                                              "cgd3-1250",
                                              "cgd3-1300",
                                              "cgd3-2250",
                                              "cgd3-300",
                                              "cgd3-3790",
                                              "cgd3-3890",
                                              "cgd3-3930",
                                              "cgd3-830",
                                              "cgd4-1230",
                                              "cgd4-2260",
                                              "cgd4-2400",
                                              "cgd4-470",
                                              "cgd4-840",
                                              "cgd5-1580",
                                              "cgd5-3823",
                                              "cgd5-970",
                                              "cgd6-2460",
                                              "cgd6-3190",
                                              "cgd6-3340",
                                              "cgd6-4190",
                                              "cgd6-4620",
                                              "cgd6-570",
                                              "cgd7-1873",
                                              "cgd7-2110",
                                              "cgd7-2420",
                                              "cgd7-2540",
                                              "cgd7-320",
                                              "cgd7-4050",
                                              "cgd7-4460",
                                              "cgd8-2340",
                                              "cgd8-2870",
                                              "cgd8-3450",
                                              "cgd8-3480",
                                              "cgd8-430",
                                              "cgd8-4350",
                                              "cgd8-440")

#Run loop to get expression data
for (current.gene in 1:length(list.of.60S.ribosome)) {
  gene.name <- list.of.60S.ribosome[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  ribosome60S_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(ribosome60S_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
ribosome60S_pseudotime_results_transposed <- as.data.frame(t(ribosome60S_pseudotime_results))

#Scale the data
ribosome60S_pseudotime_results_transposed <- as.data.frame(scale(ribosome60S_pseudotime_results_transposed))

#Get the mean of the gene group
ribosome60S_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(ribosome60S_pseudotime_results_transposed))

#Rename first column for 60S Ribosome
colnames(ribosome60S_pseudotime_results_transposed_mean) <- "Ribosome_60S"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Ribosome_60S <- ribosome60S_pseudotime_results_transposed_mean$Ribosome_60S

#Add pseudotime to column after finding mean, for plotting
ribosome60S_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(ribosome60S_pseudotime_results_transposed) <- c("cgd1_1660",
                                                         "cgd1_3000",
                                                         "cgd2_120",
                                                         "cgd2_2200",
                                                         "cgd2_280",
                                                         "Cgd2_2990",
                                                         "cgd3_1250",
                                                         "cgd3_1300",
                                                         "cgd3_2250",
                                                         "cgd3_300",
                                                         "cgd3_3790",
                                                         "cgd3_3890",
                                                         "cgd3_3930",
                                                         "cgd3_830",
                                                         "cgd4_1230",
                                                         "cgd4_2260",
                                                         "cgd4_2400",
                                                         "cgd4_470",
                                                         "cgd4_840",
                                                         "cgd5_1580",
                                                         "cgd5_3823",
                                                         "cgd5_970",
                                                         "cgd6_2460",
                                                         "cgd6_3190",
                                                         "cgd6_3340",
                                                         "cgd6_4190",
                                                         "cgd6_4620",
                                                         "cgd6_570",
                                                         "cgd7_1873",
                                                         "cgd7_2110",
                                                         "cgd7_2420",
                                                         "cgd7_2540",
                                                         "cgd7_320",
                                                         "cgd7_4050",
                                                         "cgd7_4460",
                                                         "cgd8_2340",
                                                         "cgd8_2870",
                                                         "cgd8_3450",
                                                         "cgd8_3480",
                                                         "cgd8_430",
                                                         "cgd8_4350",
                                                         "cgd8_440",
                                                         "Pseudotime")

#Now plot, Figure 1i, Extended Data Figure 3a
ggplot(ribosome60S_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1660), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3000), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_120), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2200), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=Cgd2_2990), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3890), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3930), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2260), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2400), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_470), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_840), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3823), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2460), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3190), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4190), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_570), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1873), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2540), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_320), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4050), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4460), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3480), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Ribosome_60S), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Apical LOPIT, Extended Data Figure 3a
list.of.apical <- c("cgd1-1210",
                    "cgd1-3170",
                    "cgd1-510",
                    "cgd3-2860",
                    "cgd4-2610",
                    "cgd6-1690",
                    "cgd7-4350",
                    "cgd8-4550")

apical_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 8))

rownames(apical_pseudotime_results) <- c("cgd1-1210",
                                         "cgd1-3170",
                                         "cgd1-510",
                                         "cgd3-2860",
                                         "cgd4-2610",
                                         "cgd6-1690",
                                         "cgd7-4350",
                                         "cgd8-4550")

#Run loop to get expression data
for (current.gene in 1:length(list.of.apical)) {
  gene.name <- list.of.apical[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  apical_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(apical_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
apical_pseudotime_results_transposed <- as.data.frame(t(apical_pseudotime_results))

#Scale the data
apical_pseudotime_results_transposed <- as.data.frame(scale(apical_pseudotime_results_transposed))

#Get the mean of the gene group
apical_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(apical_pseudotime_results_transposed))

#Rename first column for apical
colnames(apical_pseudotime_results_transposed_mean) <- "Apical"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$Apical <- apical_pseudotime_results_transposed_mean$Apical

#Add pseudotime to column after finding mean, for plotting
apical_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(apical_pseudotime_results_transposed) <- c("cgd1_1210",
                                                    "cgd1_3170",
                                                    "cgd1_510",
                                                    "cgd3_2860",
                                                    "cgd4_2610",
                                                    "cgd6_1690",
                                                    "cgd7_4350",
                                                    "cgd8_4550",
                                                    "Pseudotime")


#Now plot, Extended Data Figure 3a
ggplot(apical_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_510), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2860), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2610), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=Apical), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#ER LOPIT, Extended Data Figure 3a
list.of.ER <- c("cgd1-1110",
                "cgd1-1770",
                "cgd1-1830",
                "cgd1-1890",
                "cgd1-230",
                "cgd1-2400",
                "cgd1-2870",
                "cgd1-2970",
                "cgd1-3200",
                "cgd1-3850",
                "cgd1-440",
                "cgd1-570",
                "cgd1-780",
                "cgd2-1090",
                "cgd2-1200",
                "cgd2-1520",
                "cgd2-1650",
                "cgd2-1690",
                "cgd2-2160",
                "cgd2-2660",
                "cgd2-3100",
                "cgd2-3290",
                "cgd2-3750",
                "cgd2-550",
                "cgd2-940",
                "cgd3-1580",
                "cgd3-1900",
                "cgd3-2150",
                "cgd3-2503",
                "cgd3-2630",
                "cgd3-2690",
                "cgd3-3160",
                "cgd3-510",
                "cgd3-550",
                "cgd4-1080",
                "cgd4-1120",
                "cgd4-1150",
                "cgd4-2270",
                "cgd4-2580",
                "cgd4-2790",
                "cgd4-300",
                "cgd4-3610",
                "cgd4-3760",
                "cgd4-4140",
                "cgd4-4200",
                "cgd4-4400",
                "cgd4-80",
                "cgd5-1080",
                "cgd5-1270",
                "cgd5-2230",
                "cgd5-2300",
                "cgd5-2580",
                "cgd5-2590",
                "cgd5-3430",
                "cgd5-3500",
                "cgd5-3690",
                "cgd5-3730",
                "cgd5-3820",
                "cgd6-1070",
                "cgd6-110",
                "cgd6-1270",
                "cgd6-1340",
                "cgd6-1490",
                "cgd6-1960",
                "cgd6-2040",
                "cgd6-2280",
                "cgd6-2503",
                "cgd6-260",
                "cgd6-2620",
                "cgd6-2720",
                "cgd6-2980",
                "cgd6-3170",
                "cgd6-3530",
                "cgd6-4530",
                "cgd6-470",
                "cgd6-4780",
                "cgd6-5070",
                "cgd6-5140",
                "cgd6-660",
                "cgd6-70",
                "cgd6-720",
                "cgd6-840",
                "cgd6-850",
                "cgd7-10",
                "cgd7-1060",
                "cgd7-1310",
                "cgd7-1430",
                "cgd7-1810",
                "cgd7-2013",
                "cgd7-2260",
                "cgd7-3010",
                "cgd7-3070",
                "cgd7-3250",
                "cgd7-3280",
                "cgd7-3290",
                "cgd7-330",
                "cgd7-3880",
                "cgd7-4120",
                "cgd7-4160",
                "cgd7-4190",
                "cgd7-5080",
                "cgd7-5120",
                "cgd7-700",
                "cgd7-950",
                "cgd8-20",
                "cgd8-240",
                "cgd8-2410",
                "cgd8-2550",
                "cgd8-2820",
                "cgd8-293",
                "cgd8-3170",
                "cgd8-3190",
                "cgd8-3520",
                "cgd8-3540",
                "cgd8-3850",
                "cgd8-4850",
                "cgd8-4880",
                "cgd8-5140",
                "cgd8-70")

ER_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 119))

rownames(ER_pseudotime_results) <- c("cgd1-1110",
                                     "cgd1-1770",
                                     "cgd1-1830",
                                     "cgd1-1890",
                                     "cgd1-230",
                                     "cgd1-2400",
                                     "cgd1-2870",
                                     "cgd1-2970",
                                     "cgd1-3200",
                                     "cgd1-3850",
                                     "cgd1-440",
                                     "cgd1-570",
                                     "cgd1-780",
                                     "cgd2-1090",
                                     "cgd2-1200",
                                     "cgd2-1520",
                                     "cgd2-1650",
                                     "cgd2-1690",
                                     "cgd2-2160",
                                     "cgd2-2660",
                                     "cgd2-3100",
                                     "cgd2-3290",
                                     "cgd2-3750",
                                     "cgd2-550",
                                     "cgd2-940",
                                     "cgd3-1580",
                                     "cgd3-1900",
                                     "cgd3-2150",
                                     "cgd3-2503",
                                     "cgd3-2630",
                                     "cgd3-2690",
                                     "cgd3-3160",
                                     "cgd3-510",
                                     "cgd3-550",
                                     "cgd4-1080",
                                     "cgd4-1120",
                                     "cgd4-1150",
                                     "cgd4-2270",
                                     "cgd4-2580",
                                     "cgd4-2790",
                                     "cgd4-300",
                                     "cgd4-3610",
                                     "cgd4-3760",
                                     "cgd4-4140",
                                     "cgd4-4200",
                                     "cgd4-4400",
                                     "cgd4-80",
                                     "cgd5-1080",
                                     "cgd5-1270",
                                     "cgd5-2230",
                                     "cgd5-2300",
                                     "cgd5-2580",
                                     "cgd5-2590",
                                     "cgd5-3430",
                                     "cgd5-3500",
                                     "cgd5-3690",
                                     "cgd5-3730",
                                     "cgd5-3820",
                                     "cgd6-1070",
                                     "cgd6-110",
                                     "cgd6-1270",
                                     "cgd6-1340",
                                     "cgd6-1490",
                                     "cgd6-1960",
                                     "cgd6-2040",
                                     "cgd6-2280",
                                     "cgd6-2503",
                                     "cgd6-260",
                                     "cgd6-2620",
                                     "cgd6-2720",
                                     "cgd6-2980",
                                     "cgd6-3170",
                                     "cgd6-3530",
                                     "cgd6-4530",
                                     "cgd6-470",
                                     "cgd6-4780",
                                     "cgd6-5070",
                                     "cgd6-5140",
                                     "cgd6-660",
                                     "cgd6-70",
                                     "cgd6-720",
                                     "cgd6-840",
                                     "cgd6-850",
                                     "cgd7-10",
                                     "cgd7-1060",
                                     "cgd7-1310",
                                     "cgd7-1430",
                                     "cgd7-1810",
                                     "cgd7-2013",
                                     "cgd7-2260",
                                     "cgd7-3010",
                                     "cgd7-3070",
                                     "cgd7-3250",
                                     "cgd7-3280",
                                     "cgd7-3290",
                                     "cgd7-330",
                                     "cgd7-3880",
                                     "cgd7-4120",
                                     "cgd7-4160",
                                     "cgd7-4190",
                                     "cgd7-5080",
                                     "cgd7-5120",
                                     "cgd7-700",
                                     "cgd7-950",
                                     "cgd8-20",
                                     "cgd8-240",
                                     "cgd8-2410",
                                     "cgd8-2550",
                                     "cgd8-2820",
                                     "cgd8-293",
                                     "cgd8-3170",
                                     "cgd8-3190",
                                     "cgd8-3520",
                                     "cgd8-3540",
                                     "cgd8-3850",
                                     "cgd8-4850",
                                     "cgd8-4880",
                                     "cgd8-5140",
                                     "cgd8-70")


#Run loop to get expression data
for (current.gene in 1:length(list.of.ER)) {
  gene.name <- list.of.ER[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  ER_pseudotime_results[current.gene,] <- all_expression_data
  
}


#Get barcode cell labels as column names
colnames(ER_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
ER_pseudotime_results_transposed <- as.data.frame(t(ER_pseudotime_results))

#Scale the data
ER_pseudotime_results_transposed <- as.data.frame(scale(ER_pseudotime_results_transposed))

#Get the mean of the gene group
ER_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(ER_pseudotime_results_transposed))

#Rename first column for ER
colnames(ER_pseudotime_results_transposed_mean) <- "ER"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$ER <- ER_pseudotime_results_transposed_mean$ER

#Add pseudotime to column after finding mean, for plotting
ER_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(ER_pseudotime_results_transposed) <- c("cgd1_1110",
                                                "cgd1_1770",
                                                "cgd1_1830",
                                                "cgd1_1890",
                                                "cgd1_230",
                                                "cgd1_2400",
                                                "cgd1_2870",
                                                "cgd1_2970",
                                                "cgd1_3200",
                                                "cgd1_3850",
                                                "cgd1_440",
                                                "cgd1_570",
                                                "cgd1_780",
                                                "cgd2_1090",
                                                "cgd2_1200",
                                                "cgd2_1520",
                                                "cgd2_1650",
                                                "cgd2_1690",
                                                "cgd2_2160",
                                                "cgd2_2660",
                                                "cgd2_3100",
                                                "cgd2_3290",
                                                "cgd2_3750",
                                                "cgd2_550",
                                                "cgd2_940",
                                                "cgd3_1580",
                                                "cgd3_1900",
                                                "cgd3_2150",
                                                "cgd3_2503",
                                                "cgd3_2630",
                                                "cgd3_2690",
                                                "cgd3_3160",
                                                "cgd3_510",
                                                "cgd3_550",
                                                "cgd4_1080",
                                                "cgd4_1120",
                                                "cgd4_1150",
                                                "cgd4_2270",
                                                "cgd4_2580",
                                                "cgd4_2790",
                                                "cgd4_300",
                                                "cgd4_3610",
                                                "cgd4_3760",
                                                "cgd4_4140",
                                                "cgd4_4200",
                                                "cgd4_4400",
                                                "cgd4_80",
                                                "cgd5_1080",
                                                "cgd5_1270",
                                                "cgd5_2230",
                                                "cgd5_2300",
                                                "cgd5_2580",
                                                "cgd5_2590",
                                                "cgd5_3430",
                                                "cgd5_3500",
                                                "cgd5_3690",
                                                "cgd5_3730",
                                                "cgd5_3820",
                                                "cgd6_1070",
                                                "cgd6_110",
                                                "cgd6_1270",
                                                "cgd6_1340",
                                                "cgd6_1490",
                                                "cgd6_1960",
                                                "cgd6_2040",
                                                "cgd6_2280",
                                                "cgd6_2503",
                                                "cgd6_260",
                                                "cgd6_2620",
                                                "cgd6_2720",
                                                "cgd6_2980",
                                                "cgd6_3170",
                                                "cgd6_3530",
                                                "cgd6_4530",
                                                "cgd6_470",
                                                "cgd6_4780",
                                                "cgd6_5070",
                                                "cgd6_5140",
                                                "cgd6_660",
                                                "cgd6_70",
                                                "cgd6_720",
                                                "cgd6_840",
                                                "cgd6_850",
                                                "cgd7_10",
                                                "cgd7_1060",
                                                "cgd7_1310",
                                                "cgd7_1430",
                                                "cgd7_1810",
                                                "cgd7_2013",
                                                "cgd7_2260",
                                                "cgd7_3010",
                                                "cgd7_3070",
                                                "cgd7_3250",
                                                "cgd7_3280",
                                                "cgd7_3290",
                                                "cgd7_330",
                                                "cgd7_3880",
                                                "cgd7_4120",
                                                "cgd7_4160",
                                                "cgd7_4190",
                                                "cgd7_5080",
                                                "cgd7_5120",
                                                "cgd7_700",
                                                "cgd7_950",
                                                "cgd8_20",
                                                "cgd8_240",
                                                "cgd8_2410",
                                                "cgd8_2550",
                                                "cgd8_2820",
                                                "cgd8_293",
                                                "cgd8_3170",
                                                "cgd8_3190",
                                                "cgd8_3520",
                                                "cgd8_3540",
                                                "cgd8_3850",
                                                "cgd8_4850",
                                                "cgd8_4880",
                                                "cgd8_5140",
                                                "cgd8_70",
                                                "Pseudotime")

#Now plot, Extended Data Figure 3a
ggplot(ER_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1770), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1890), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_2400), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_2870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_2970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3200), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_570), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1090), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1200), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1520), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1650), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2660), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3100), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3290), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3750), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_940), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1900), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2150), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2503), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2630), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_2690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_510), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1120), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1150), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3610), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3760), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_4140), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_4200), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_4400), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_80), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2590), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3500), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3730), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1960), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2040), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2503), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_260), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2980), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_470), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5140), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_660), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_70), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_840), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_10), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1060), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1810), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2013), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2260), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3010), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3290), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3880), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4120), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4190), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5120), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_700), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_950), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_20), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_240), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2410), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_293), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3190), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3520), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3540), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4880), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5140), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_70), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=ER), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#IMC/PM LOPIT, Extended Data Figure 3a
list.of.IMC.PM <- c("cgd1-1970",
                    "cgd1-2060",
                    "cgd1-3210",
                    "cgd1-450",
                    "cgd1-820",
                    "cgd2-1140",
                    "cgd2-1290",
                    "cgd2-1450",
                    "cgd2-1910",
                    "cgd2-2000",
                    "cgd2-2310",
                    "cgd2-3350",
                    "cgd2-4150",
                    "cgd2-4180",
                    "cgd2-540",
                    "cgd2-640",
                    "cgd3-3210",
                    "cgd3-3740",
                    "cgd3-3830",
                    "cgd3-4040",
                    "cgd3-590",
                    "cgd3-710",
                    "cgd3-880",
                    "cgd4-2080",
                    "cgd4-2280",
                    "cgd4-2560",
                    "cgd4-3290",
                    "cgd4-3300",
                    "cgd4-3430",
                    "cgd4-4340",
                    "cgd4-700",
                    "cgd5-1020",
                    "cgd5-1060",
                    "cgd5-1590",
                    "cgd5-1663",
                    "cgd5-2530",
                    "cgd5-2900",
                    "cgd5-3130",
                    "cgd5-3420",
                    "cgd5-4080",
                    "cgd5-640",
                    "cgd5-690",
                    "cgd6-1083",
                    "cgd6-1500",
                    "cgd6-2210",
                    "cgd6-2220",
                    "cgd6-2913",
                    "cgd6-3130",
                    "cgd6-3220",
                    "cgd6-3470",
                    "cgd6-3510",
                    "cgd6-3830",
                    "cgd6-4160",
                    "cgd6-4390",
                    "cgd6-4540",
                    "cgd6-4940",
                    "cgd6-5350",
                    "cgd6-5450",
                    "cgd6-5460",
                    "cgd6-680",
                    "cgd6-760",
                    "cgd7-120",
                    "cgd7-1713",
                    "cgd7-1900",
                    "cgd7-3640",
                    "cgd7-3790",
                    "cgd7-4300",
                    "cgd7-4380",
                    "cgd7-4540",
                    "cgd7-4720",
                    "cgd7-4830",
                    "cgd7-4880",
                    "cgd7-550",
                    "cgd7-580",
                    "cgd8-1180",
                    "cgd8-2270",
                    "cgd8-2790",
                    "cgd8-3430",
                    "cgd8-3730",
                    "cgd8-4570",
                    "cgd8-5030")

IMC_PM_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 81))

rownames(IMC_PM_pseudotime_results) <- c("cgd1-1970",
                                         "cgd1-2060",
                                         "cgd1-3210",
                                         "cgd1-450",
                                         "cgd1-820",
                                         "cgd2-1140",
                                         "cgd2-1290",
                                         "cgd2-1450",
                                         "cgd2-1910",
                                         "cgd2-2000",
                                         "cgd2-2310",
                                         "cgd2-3350",
                                         "cgd2-4150",
                                         "cgd2-4180",
                                         "cgd2-540",
                                         "cgd2-640",
                                         "cgd3-3210",
                                         "cgd3-3740",
                                         "cgd3-3830",
                                         "cgd3-4040",
                                         "cgd3-590",
                                         "cgd3-710",
                                         "cgd3-880",
                                         "cgd4-2080",
                                         "cgd4-2280",
                                         "cgd4-2560",
                                         "cgd4-3290",
                                         "cgd4-3300",
                                         "cgd4-3430",
                                         "cgd4-4340",
                                         "cgd4-700",
                                         "cgd5-1020",
                                         "cgd5-1060",
                                         "cgd5-1590",
                                         "cgd5-1663",
                                         "cgd5-2530",
                                         "cgd5-2900",
                                         "cgd5-3130",
                                         "cgd5-3420",
                                         "cgd5-4080",
                                         "cgd5-640",
                                         "cgd5-690",
                                         "cgd6-1083",
                                         "cgd6-1500",
                                         "cgd6-2210",
                                         "cgd6-2220",
                                         "cgd6-2913",
                                         "cgd6-3130",
                                         "cgd6-3220",
                                         "cgd6-3470",
                                         "cgd6-3510",
                                         "cgd6-3830",
                                         "cgd6-4160",
                                         "cgd6-4390",
                                         "cgd6-4540",
                                         "cgd6-4940",
                                         "cgd6-5350",
                                         "cgd6-5450",
                                         "cgd6-5460",
                                         "cgd6-680",
                                         "cgd6-760",
                                         "cgd7-120",
                                         "cgd7-1713",
                                         "cgd7-1900",
                                         "cgd7-3640",
                                         "cgd7-3790",
                                         "cgd7-4300",
                                         "cgd7-4380",
                                         "cgd7-4540",
                                         "cgd7-4720",
                                         "cgd7-4830",
                                         "cgd7-4880",
                                         "cgd7-550",
                                         "cgd7-580",
                                         "cgd8-1180",
                                         "cgd8-2270",
                                         "cgd8-2790",
                                         "cgd8-3430",
                                         "cgd8-3730",
                                         "cgd8-4570",
                                         "cgd8-5030")


#Run loop to get expression data
for (current.gene in 1:length(list.of.IMC.PM)) {
  gene.name <- list.of.IMC.PM[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  IMC_PM_pseudotime_results[current.gene,] <- all_expression_data
  
}


#Get barcode cell labels as column names
colnames(IMC_PM_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
IMC_PM_pseudotime_results_transposed <- as.data.frame(t(IMC_PM_pseudotime_results))

#Scale the data
IMC_PM_pseudotime_results_transposed <- as.data.frame(scale(IMC_PM_pseudotime_results_transposed))

#Get the mean of the gene group
IMC_PM_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(IMC_PM_pseudotime_results_transposed))

#Rename first column for IMC/PM
colnames(IMC_PM_pseudotime_results_transposed_mean) <- "IMC_PM"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$IMC_PM <- IMC_PM_pseudotime_results_transposed_mean$IMC_PM

#Add pseudotime to column after finding mean, for plotting
IMC_PM_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(IMC_PM_pseudotime_results_transposed) <- c("cgd1_1970",
                                                    "cgd1_2060",
                                                    "cgd1_3210",
                                                    "cgd1_450",
                                                    "cgd1_820",
                                                    "cgd2_1140",
                                                    "cgd2_1290",
                                                    "cgd2_1450",
                                                    "cgd2_1910",
                                                    "cgd2_2000",
                                                    "cgd2_2310",
                                                    "cgd2_3350",
                                                    "cgd2_4150",
                                                    "cgd2_4180",
                                                    "cgd2_540",
                                                    "cgd2_640",
                                                    "cgd3_3210",
                                                    "cgd3_3740",
                                                    "cgd3_3830",
                                                    "cgd3_4040",
                                                    "cgd3_590",
                                                    "cgd3_710",
                                                    "cgd3_880",
                                                    "cgd4_2080",
                                                    "cgd4_2280",
                                                    "cgd4_2560",
                                                    "cgd4_3290",
                                                    "cgd4_3300",
                                                    "cgd4_3430",
                                                    "cgd4_4340",
                                                    "cgd4_700",
                                                    "cgd5_1020",
                                                    "cgd5_1060",
                                                    "cgd5_1590",
                                                    "cgd5_1663",
                                                    "cgd5_2530",
                                                    "cgd5_2900",
                                                    "cgd5_3130",
                                                    "cgd5_3420",
                                                    "cgd5_4080",
                                                    "cgd5_640",
                                                    "cgd5_690",
                                                    "cgd6_1083",
                                                    "cgd6_1500",
                                                    "cgd6_2210",
                                                    "cgd6_2220",
                                                    "cgd6_2913",
                                                    "cgd6_3130",
                                                    "cgd6_3220",
                                                    "cgd6_3470",
                                                    "cgd6_3510",
                                                    "cgd6_3830",
                                                    "cgd6_4160",
                                                    "cgd6_4390",
                                                    "cgd6_4540",
                                                    "cgd6_4940",
                                                    "cgd6_5350",
                                                    "cgd6_5450",
                                                    "cgd6_5460",
                                                    "cgd6_680",
                                                    "cgd6_760",
                                                    "cgd7_120",
                                                    "cgd7_1713",
                                                    "cgd7_1900",
                                                    "cgd7_3640",
                                                    "cgd7_3790",
                                                    "cgd7_4300",
                                                    "cgd7_4380",
                                                    "cgd7_4540",
                                                    "cgd7_4720",
                                                    "cgd7_4830",
                                                    "cgd7_4880",
                                                    "cgd7_550",
                                                    "cgd7_580",
                                                    "cgd8_1180",
                                                    "cgd8_2270",
                                                    "cgd8_2790",
                                                    "cgd8_3430",
                                                    "cgd8_3730",
                                                    "cgd8_4570",
                                                    "cgd8_5030",
                                                    "Pseudotime")

#Now plot, Extended Data Figure 3a
ggplot(IMC_PM_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_2060), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1140), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1290), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1910), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2000), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_4150), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_4180), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_540), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_640), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3740), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_4040), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_590), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_710), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_880), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2560), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3290), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_4340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_700), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1020), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1060), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1590), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1663), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_2900), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3130), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_4080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_640), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1083), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1500), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2220), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2913), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3130), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3220), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3470), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3510), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4390), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4540), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4940), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5460), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_760), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_120), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1713), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1900), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3640), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4540), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4880), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1180), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3730), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4570), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=IMC_PM), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#DNA Replication GO term, Figure 3h
list.of.DNA.replication <- c("cgd2-1100",
                             "cgd2-1250",
                             "cgd2-1550",
                             "cgd2-1600",
                             "cgd2-2500",
                             "cgd2-3180",
                             "cgd3-1450",
                             "cgd3-3170",
                             "cgd3-3470",
                             "cgd3-3820",
                             "cgd3-4290",
                             "cgd4-1283",
                             "cgd4-1490",
                             "cgd4-1930",
                             "cgd4-430",
                             "cgd4-970",
                             "cgd5-410",
                             "cgd5-4173",
                             "cgd6-1710",
                             "cgd6-1940",
                             "cgd6-1950",
                             "cgd6-240",
                             "cgd6-4410",
                             "cgd7-1930",
                             "cgd7-2140",
                             "cgd7-2920",
                             "cgd7-5280",
                             "cgd7-710",
                             "cgd8-1240",
                             "cgd8-1410",
                             "cgd8-1620",
                             "cgd8-1630",
                             "cgd8-2150",
                             "cgd8-2380",
                             "cgd8-2940",
                             "cgd8-3820",
                             "cgd8-4650",
                             "cgd8-610",
                             "cgd8-870")

#This gives the right thing, add barcodes and pseudotime as row
DNA_replication_pseudotime_results <- data.frame(matrix(ncol = 2989, nrow = 39))

rownames(DNA_replication_pseudotime_results) <- c("cgd2-1100",
                                                  "cgd2-1250",
                                                  "cgd2-1550",
                                                  "cgd2-1600",
                                                  "cgd2-2500",
                                                  "cgd2-3180",
                                                  "cgd3-1450",
                                                  "cgd3-3170",
                                                  "cgd3-3470",
                                                  "cgd3-3820",
                                                  "cgd3-4290",
                                                  "cgd4-1283",
                                                  "cgd4-1490",
                                                  "cgd4-1930",
                                                  "cgd4-430",
                                                  "cgd4-970",
                                                  "cgd5-410",
                                                  "cgd5-4173",
                                                  "cgd6-1710",
                                                  "cgd6-1940",
                                                  "cgd6-1950",
                                                  "cgd6-240",
                                                  "cgd6-4410",
                                                  "cgd7-1930",
                                                  "cgd7-2140",
                                                  "cgd7-2920",
                                                  "cgd7-5280",
                                                  "cgd7-710",
                                                  "cgd8-1240",
                                                  "cgd8-1410",
                                                  "cgd8-1620",
                                                  "cgd8-1630",
                                                  "cgd8-2150",
                                                  "cgd8-2380",
                                                  "cgd8-2940",
                                                  "cgd8-3820",
                                                  "cgd8-4650",
                                                  "cgd8-610",
                                                  "cgd8-870")

#Run loop to get expression data
for (current.gene in 1:length(list.of.DNA.replication)) {
  gene.name <- list.of.DNA.replication[current.gene]
  FindMyGene <- which(all.genes.crypto==gene.name)
  all_expression_data <- GetAssayData(seurat.integrated_0.4.updated, slot = 'data')[FindMyGene, all.parasites]
  print(gene.name)
  as.vector(all_expression_data)
  DNA_replication_pseudotime_results[current.gene,] <- all_expression_data
  
}

#Get barcode cell labels as column names
colnames(DNA_replication_pseudotime_results) <- seurat.integrated_0.4.updated@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
DNA_replication_pseudotime_results_transposed <- as.data.frame(t(DNA_replication_pseudotime_results))

#Scale the data
DNA_replication_pseudotime_results_transposed <- as.data.frame(scale(DNA_replication_pseudotime_results_transposed))

#Get the mean of the gene group
DNA_replication_pseudotime_results_transposed_mean <- as.data.frame(rowMeans(DNA_replication_pseudotime_results_transposed))

#Rename first column for DNA replication
colnames(DNA_replication_pseudotime_results_transposed_mean) <- "DNA_Replication"

#Now add this data to existing data frame as a new column
mean_gene_groups_pseudotime$DNA_Replication <- DNA_replication_pseudotime_results_transposed_mean$DNA_Replication

#Add pseudotime to column after finding mean, for plotting
DNA_replication_pseudotime_results_transposed$Pseudotime <- all_pseudotime$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(DNA_replication_pseudotime_results_transposed) <- c("cgd2_1100",
                                                             "cgd2_1250",
                                                             "cgd2_1550",
                                                             "cgd2_1600",
                                                             "cgd2_2500",
                                                             "cgd2_3180",
                                                             "cgd3_1450",
                                                             "cgd3_3170",
                                                             "cgd3_3470",
                                                             "cgd3_3820",
                                                             "cgd3_4290",
                                                             "cgd4_1283",
                                                             "cgd4_1490",
                                                             "cgd4_1930",
                                                             "cgd4_430",
                                                             "cgd4_970",
                                                             "cgd5_410",
                                                             "cgd5_4173",
                                                             "cgd6_1710",
                                                             "cgd6_1940",
                                                             "cgd6_1950",
                                                             "cgd6_240",
                                                             "cgd6_4410",
                                                             "cgd7_1930",
                                                             "cgd7_2140",
                                                             "cgd7_2920",
                                                             "cgd7_5280",
                                                             "cgd7_710",
                                                             "cgd8_1240",
                                                             "cgd8_1410",
                                                             "cgd8_1620",
                                                             "cgd8_1630",
                                                             "cgd8_2150",
                                                             "cgd8_2380",
                                                             "cgd8_2940",
                                                             "cgd8_3820",
                                                             "cgd8_4650",
                                                             "cgd8_610",
                                                             "cgd8_870",
                                                             "Pseudotime")


#Now plot, Figure 3h
ggplot(DNA_replication_pseudotime_results_transposed, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd2_1100), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1600), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2500), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3180), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3470), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_4290), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1283), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_1930), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_410), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_4173), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1710), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1940), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1950), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_240), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4410), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1930), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2140), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2920), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_710), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1240), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1410), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1630), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2150), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2940), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4650), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_610), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(data = mean_gene_groups_pseudotime, aes(y=DNA_Replication), level = 0.95, colour="black", fill = "red1", show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#Plot all means together, Extended Data Figure 3b
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Rhoptry), fill=TRUE, colour="royalblue1", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Microneme), fill=TRUE, colour="purple1", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Dense_Granule), fill=TRUE, colour="seagreen3", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Small_Granule), fill=TRUE, colour="violet", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Proteasome_20S), fill=TRUE, colour="sienna3", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Ribosome_40S), fill=TRUE, colour="coral", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Ribosome_60S), fill=TRUE, colour="firebrick1", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=Apical), fill=TRUE, colour="turquoise3", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=ER), fill=TRUE, colour="tan2", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=IMC_PM), fill=TRUE, colour="chartreuse3", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) +
  coord_cartesian(xlim=c(0, 60), ylim=c(-1, 2))

#Plot each LOPIT group mean on a min/max axis, Figure 1k-m
#Rhoptry
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Rhoptry), fill=TRUE, colour="royalblue1", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Rhoptry") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Microneme
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Microneme), fill=TRUE, colour="purple1", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Microneme") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#Dense Granule
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Dense_Granule), fill=TRUE, colour="seagreen3", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Dense_Granule") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#Small_Granule
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Small_Granule), fill=TRUE, colour="violet", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Small_Granule") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#Proteasome_20S
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Proteasome_20S), fill=TRUE, colour="sienna3", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Proteasome_20S") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#40S Ribosome
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Ribosome_40S), fill=TRUE, colour="coral", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("40S Ribosome") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16))

#60S Ribosome
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Ribosome_60S), fill=TRUE, colour="firebrick1", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("60S Ribosome") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16))

#Apical
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=Apical), fill=TRUE, colour="turquoise3", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Apical") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#ER
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=ER), fill=TRUE, colour="tan2", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("ER") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#IMC_PM
ggplot(mean_gene_groups_pseudotime, aes(Pseudotime)) + 
  geom_smooth(aes(y=IMC_PM), fill=TRUE, colour="chartreuse3", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("IMC/PM") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 

#Make heatmap of scaled expression of genes from LOPIT groups (those plotted in Extended Data Figure 3a) across asexual clusters
#Order of genes corresponds to order of LOPIT groups in Extended Data Figure 3a
#This is Extended Data Figure 3c
DoHeatmap(seurat.integrated_0.4.updated,  
          features = c("cgd1-850",
                       "cgd2-1070",
                       "cgd2-170",
                       "cgd2-3000",
                       "cgd2-4260",
                       "cgd3-2440",
                       "cgd4-3160",
                       "cgd4-4020",
                       "cgd5-2210",
                       "cgd5-3040",
                       "cgd5-3720",
                       "cgd6-1390",
                       "cgd6-3180",
                       "cgd6-3710",
                       "cgd6-4320",
                       "cgd6-4630",
                       "cgd7-2250",
                       "cgd7-4760",
                       "cgd8-1840",
                       "cgd8-4360",
                       "cgd1-1660",
                       "cgd1-3000",
                       "cgd2-120",
                       "cgd2-2200",
                       "cgd2-280",
                       "Cgd2-2990",
                       "cgd3-1250",
                       "cgd3-1300",
                       "cgd3-2250",
                       "cgd3-300",
                       "cgd3-3790",
                       "cgd3-3890",
                       "cgd3-3930",
                       "cgd3-830",
                       "cgd4-1230",
                       "cgd4-2260",
                       "cgd4-2400",
                       "cgd4-470",
                       "cgd4-840",
                       "cgd5-1580",
                       "cgd5-3823",
                       "cgd5-970",
                       "cgd6-2460",
                       "cgd6-3190",
                       "cgd6-3340",
                       "cgd6-4190",
                       "cgd6-4620",
                       "cgd6-570",
                       "cgd7-1873",
                       "cgd7-2110",
                       "cgd7-2420",
                       "cgd7-2540",
                       "cgd7-320",
                       "cgd7-4050",
                       "cgd7-4460",
                       "cgd8-2340",
                       "cgd8-2870",
                       "cgd8-3450",
                       "cgd8-3480",
                       "cgd8-430",
                       "cgd8-4350",
                       "cgd8-440",
                       "cgd1-2490",
                       "cgd1-420",
                       "cgd2-1440",
                       "cgd2-2050",
                       "cgd2-530",
                       "cgd2-860",
                       "cgd3-2170",
                       "cgd3-2200",
                       "cgd3-2530",
                       "cgd4-250",
                       "cgd5-1820",
                       "cgd5-3220",
                       "cgd5-4210",
                       "cgd7-3660",
                       "cgd1-1110",
                       "cgd1-1770",
                       "cgd1-1830",
                       "cgd1-1890",
                       "cgd1-230",
                       "cgd1-2400",
                       "cgd1-2870",
                       "cgd1-2970",
                       "cgd1-3200",
                       "cgd1-3850",
                       "cgd1-440",
                       "cgd1-570",
                       "cgd1-780",
                       "cgd2-1090",
                       "cgd2-1200",
                       "cgd2-1520",
                       "cgd2-1650",
                       "cgd2-1690",
                       "cgd2-2160",
                       "cgd2-2660",
                       "cgd2-3100",
                       "cgd2-3290",
                       "cgd2-3750",
                       "cgd2-550",
                       "cgd2-940",
                       "cgd3-1580",
                       "cgd3-1900",
                       "cgd3-2150",
                       "cgd3-2503",
                       "cgd3-2630",
                       "cgd3-2690",
                       "cgd3-3160",
                       "cgd3-510",
                       "cgd3-550",
                       "cgd4-1080",
                       "cgd4-1120",
                       "cgd4-1150",
                       "cgd4-2270",
                       "cgd4-2580",
                       "cgd4-2790",
                       "cgd4-300",
                       "cgd4-3610",
                       "cgd4-3760",
                       "cgd4-4140",
                       "cgd4-4200",
                       "cgd4-4400",
                       "cgd4-80",
                       "cgd5-1080",
                       "cgd5-1270",
                       "cgd5-2230",
                       "cgd5-2300",
                       "cgd5-2580",
                       "cgd5-2590",
                       "cgd5-3430",
                       "cgd5-3500",
                       "cgd5-3690",
                       "cgd5-3730",
                       "cgd5-3820",
                       "cgd6-1070",
                       "cgd6-110",
                       "cgd6-1270",
                       "cgd6-1340",
                       "cgd6-1490",
                       "cgd6-1960",
                       "cgd6-2040",
                       "cgd6-2280",
                       "cgd6-2503",
                       "cgd6-260",
                       "cgd6-2620",
                       "cgd6-2720",
                       "cgd6-2980",
                       "cgd6-3170",
                       "cgd6-3530",
                       "cgd6-4530",
                       "cgd6-470",
                       "cgd6-4780",
                       "cgd6-5070",
                       "cgd6-5140",
                       "cgd6-660",
                       "cgd6-70",
                       "cgd6-720",
                       "cgd6-840",
                       "cgd6-850",
                       "cgd7-10",
                       "cgd7-1060",
                       "cgd7-1310",
                       "cgd7-1430",
                       "cgd7-1810",
                       "cgd7-2013",
                       "cgd7-2260",
                       "cgd7-3010",
                       "cgd7-3070",
                       "cgd7-3250",
                       "cgd7-3280",
                       "cgd7-3290",
                       "cgd7-330",
                       "cgd7-3880",
                       "cgd7-4120",
                       "cgd7-4160",
                       "cgd7-4190",
                       "cgd7-5080",
                       "cgd7-5120",
                       "cgd7-700",
                       "cgd7-950",
                       "cgd8-20",
                       "cgd8-240",
                       "cgd8-2410",
                       "cgd8-2550",
                       "cgd8-2820",
                       "cgd8-293",
                       "cgd8-3170",
                       "cgd8-3190",
                       "cgd8-3520",
                       "cgd8-3540",
                       "cgd8-3850",
                       "cgd8-4850",
                       "cgd8-4880",
                       "cgd8-5140",
                       "cgd8-70",
                       "cgd1-1970",
                       "cgd1-2060",
                       "cgd1-3210",
                       "cgd1-450",
                       "cgd1-820",
                       "cgd2-1140",
                       "cgd2-1290",
                       "cgd2-1450",
                       "cgd2-1910",
                       "cgd2-2000",
                       "cgd2-2310",
                       "cgd2-3350",
                       "cgd2-4150",
                       "cgd2-4180",
                       "cgd2-540",
                       "cgd2-640",
                       "cgd3-3210",
                       "cgd3-3740",
                       "cgd3-3830",
                       "cgd3-4040",
                       "cgd3-590",
                       "cgd3-710",
                       "cgd3-880",
                       "cgd4-2080",
                       "cgd4-2280",
                       "cgd4-2560",
                       "cgd4-3290",
                       "cgd4-3300",
                       "cgd4-3430",
                       "cgd4-4340",
                       "cgd4-700",
                       "cgd5-1020",
                       "cgd5-1060",
                       "cgd5-1590",
                       "cgd5-1663",
                       "cgd5-2530",
                       "cgd5-2900",
                       "cgd5-3130",
                       "cgd5-3420",
                       "cgd5-4080",
                       "cgd5-640",
                       "cgd5-690",
                       "cgd6-1083",
                       "cgd6-1500",
                       "cgd6-2210",
                       "cgd6-2220",
                       "cgd6-2913",
                       "cgd6-3130",
                       "cgd6-3220",
                       "cgd6-3470",
                       "cgd6-3510",
                       "cgd6-3830",
                       "cgd6-4160",
                       "cgd6-4390",
                       "cgd6-4540",
                       "cgd6-4940",
                       "cgd6-5350",
                       "cgd6-5450",
                       "cgd6-5460",
                       "cgd6-680",
                       "cgd6-760",
                       "cgd7-120",
                       "cgd7-1713",
                       "cgd7-1900",
                       "cgd7-3640",
                       "cgd7-3790",
                       "cgd7-4300",
                       "cgd7-4380",
                       "cgd7-4540",
                       "cgd7-4720",
                       "cgd7-4830",
                       "cgd7-4880",
                       "cgd7-550",
                       "cgd7-580",
                       "cgd8-1180",
                       "cgd8-2270",
                       "cgd8-2790",
                       "cgd8-3430",
                       "cgd8-3730",
                       "cgd8-4570",
                       "cgd8-5030",
                       "cgd1-1230",
                       "cgd1-1250",
                       "cgd1-1433",
                       "cgd1-3780",
                       "cgd1-3790",
                       "cgd1-590",
                       "cgd1-620",
                       "cgd2-2530",
                       "cgd2-2580",
                       "cgd2-3030",
                       "cgd2-3210",
                       "cgd2-3360",
                       "cgd2-430",
                       "cgd3-1690",
                       "cgd3-1700",
                       "cgd3-1750",
                       "cgd3-530",
                       "cgd3-560",
                       "cgd3-600",
                       "cgd3-650",
                       "cgd3-693",
                       "cgd4-230",
                       "cgd4-2510",
                       "cgd4-3460",
                       "cgd4-3530",
                       "cgd4-3580",
                       "cgd4-3970",
                       "cgd4-850",
                       "cgd5-1440",
                       "cgd5-1450",
                       "cgd5-1480",
                       "cgd5-1490",
                       "cgd5-1530",
                       "cgd5-3210",
                       "cgd6-1030",
                       "cgd6-1110",
                       "cgd6-1170",
                       "cgd6-1180",
                       "cgd6-3050",
                       "cgd6-3080",
                       "cgd6-5280",
                       "cgd6-5310",
                       "cgd6-5330",
                       "cgd6-5360",
                       "cgd6-5420",
                       "cgd6-5440",
                       "cgd6-60",
                       "cgd7-1270",
                       "cgd7-1300",
                       "cgd7-1340",
                       "cgd7-2070",
                       "Cgd7-2213",
                       "cgd7-2340",
                       "cgd7-2620",
                       "cgd7-3800",
                       "cgd7-3810",
                       "cgd7-3820",
                       "cgd7-3860",
                       "cgd7-3870",
                       "cgd7-4280",
                       "cgd7-4340",
                       "cgd7-4420",
                       "cgd7-4490",
                       "cgd7-4500",
                       "cgd7-4530",
                       "cgd7-4680",
                       "cgd8-1780",
                       "cgd8-1800",
                       "cgd8-2160",
                       "cgd8-3030",
                       "cgd8-3600",
                       "cgd8-3610",
                       "cgd8-5300",
                       "cgd8-5310",
                       "cgd8-5350",
                       "cgd8-5360",
                       "cgd8-5380",
                       "cgd8-680",
                       "cgd8-690",
                       "cgd1-1210",
                       "cgd1-3170",
                       "cgd1-510",
                       "cgd3-2860",
                       "cgd4-2610",
                       "cgd6-1690",
                       "cgd7-4350",
                       "cgd8-4550",
                       "cgd4-1790", "cgd8-540", "cgd8-2530", "cgd4-2420", "cgd1-1870", "cgd2-370", 
                       "cgd3-2010",  "cgd1-1380", "cgd1-950", "cgd2-2900", "cgd3-1330", "cgd3-1710", 
                       "cgd3-1730", "cgd3-1770", "cgd3-1780", "cgd3-2750", "cgd3-910", "cgd5-20", 
                       "cgd5-2760", "cgd6-1000", "cgd6-3630", "cgd6-3930", "cgd6-3940", "cgd8-520",
                       "cgd4-2350",
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
                       "cgd8-530",
                       "cgd7-200",
                       "cgd7-3390",
                       "cgd1-3640",
                       "cgd1-3650",
                       "cgd1-3810",
                       "cgd1-640",
                       "cgd2-3110",
                       "cgd2-340",
                       "cgd2-3730",
                       "cgd2-4020",
                       "cgd3-4250",
                       "cgd4-2460",
                       "cgd4-3050",
                       "cgd5-2720",
                       "cgd5-2730",
                       "cgd5-2740",
                       "cgd5-320",
                       "cgd6-3070",
                       "cgd6-5380",
                       "cgd7-3830",
                       "cgd8-1770",
                       "cgd8-2490",
                       "cgd8-3590"), 
          group.bar = TRUE, group.colors = c("lightgreen", "seagreen1", "seagreen3",
                                             "chartreuse3", "yellowgreen", "chartreuse1", "green1", "green3",
                                             "mediumseagreen")) 




#Now determine pseudotime for males and females

load("seurat.for.monocle.Robj")
DimPlot(seurat.for.monocle, group.by = "sample")

seurat.for.monocle.updated = UpdateSeuratObject(object = seurat.for.monocle)

seurat.for.monocle[["UMAP"]] <- seurat.for.monocle[["umap"]]

#Need to subset data for sexual male and female
#Create male subset, new seurat object
male_for_monocle <- subset(seurat.for.monocle, idents = c("10", "11", "12"))

#Check to make sure it is male trajectory
#No legend, Figure 3i
DimPlot(male_for_monocle, group.by = "cluster_labels", cols = c("#00B9E3", "#00ADFA", "#619CFF"), pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + 
  NoLegend() + xlim(-3, 4)+ylim(-4, 7) + theme(axis.text=element_text(size = 16)) +
  ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
DimPlot(male_for_monocle, group.by = "cluster_labels", cols = c("#00B9E3", "#00ADFA", "#619CFF"), pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + 
  xlim(-3, 4)+ylim(-4, 7) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) +
  ggtitle(NULL) +
  coord_fixed(ratio = 1)


#Create female subset,new seurat object
female_for_monocle <- subset(seurat.for.monocle, idents = c("13", "14", "15", "16", "17", "18"))

#Check to make sure it is female trajectory
#No legend, Figure 3b
DimPlot(female_for_monocle, group.by = "cluster_labels", cols = c("orchid3", "orchid2", "orchid1", "hotpink2", "hotpink1", "lightpink2"), pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + 
  NoLegend() + xlim(0, 11.5)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) +
  ggtitle(NULL) +
  coord_fixed(ratio = 1)

#With legend
DimPlot(female_for_monocle, group.by = "cluster_labels", cols = c("orchid3", "orchid2", "orchid1", "hotpink2", "hotpink1", "lightpink2"), pt.size = 2) +
  FontSize(x.title = 20, y.title = 20) + 
  xlim(0, 11.5)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(0.6, 'cm'), legend.text=element_text(size=16)) +
  ggtitle(NULL) +
  coord_fixed(ratio = 1)

#Do male analysis
#Convert to cell dataset
cds_male <- as.cell_data_set(male_for_monocle)

#Pick a trajectory
cds_male <- cluster_cells(cds_male, reduction_method = c("UMAP"), k = 20, resolution = 5e-3)

plot_cells(cds_male, show_trajectory_graph = FALSE)

cds_male <- learn_graph(cds_male, 
                        close_loop = T,
                        learn_graph_control=list(ncenter=200, minimal_branch_len=12),
                        use_partition = T)

cds_male <- order_cells(cds_male)

#Now plot

#No legend, Figure 3j
plot_cells(cds_male, color_cells_by = "pseudotime", cell_size = 2, trajectory_graph_color = "red", alpha = 1, trajectory_graph_segment_size = 1.2, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(-3, 4)+ylim(-4, 7) + theme(axis.text=element_text(size = 16)) +
  viridis::scale_color_viridis(begin = 0.05, end = 1, option = "D") +
  theme(axis.line.y = element_line(size=0.5, color="black")) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  coord_fixed(ratio = 1) 

#With legend
plot_cells(cds_male, color_cells_by = "pseudotime", cell_size = 2, trajectory_graph_color = "red", alpha = 1, trajectory_graph_segment_size = 1.2, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE) + FontSize(x.title = 20, y.title = 20) + xlim(-3, 4)+ylim(-4, 7) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16), legend.title = element_blank()) +
  viridis::scale_color_viridis(begin = 0.05, end = 1, option = "D") +
  theme(axis.line.y = element_line(size=0.5, color="black")) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  coord_fixed(ratio = 1) 

#Do female analysis
cds_female <- as.cell_data_set(female_for_monocle)

#Pick a trajectory
cds_female <- cluster_cells(cds_female, reduction_method = c("UMAP"), k = 40, resolution = 5e-4)

plot_cells(cds_female, show_trajectory_graph = FALSE)

cds_female <- learn_graph(cds_female, 
                          close_loop = T,
                          learn_graph_control=list(ncenter=130, minimal_branch_len=8),
                          use_partition = T)

cds_female <- order_cells(cds_female)

#Now plot

#No legend, Figure 3c
plot_cells(cds_female, color_cells_by = "pseudotime", cell_size = 2, trajectory_graph_color = "red", alpha = 1, trajectory_graph_segment_size = 1.2, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE) + NoLegend() + FontSize(x.title = 20, y.title = 20) + xlim(0, 11.5)+ylim(-14, 8) + theme(axis.text=element_text(size = 16)) +
  viridis::scale_color_viridis(begin = 0.05, end = 1, option = "D") +
  theme(axis.line.y = element_line(size=0.5, color="black")) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  coord_fixed(ratio = 1) 

#With legend
plot_cells(cds_female, color_cells_by = "pseudotime", cell_size = 2, trajectory_graph_color = "red", alpha = 1, trajectory_graph_segment_size = 1.2, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE) + FontSize(x.title = 20, y.title = 20) + xlim(0, 11.5)+ylim(-14, 8) + theme(axis.text=element_text(size = 16), legend.key.size = unit(1.5, 'cm'), legend.text=element_text(size=16), legend.title = element_blank()) +
  viridis::scale_color_viridis(begin = 0.05, end = 1, option = "D") +
  theme(axis.line.y = element_line(size=0.5, color="black")) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  coord_fixed(ratio = 1) 


#Label pseudotime and add column to metadata of Seurat object
#Male
pseudotime_start_male <- cds_male@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
male_for_monocle <- AddMetaData(male_for_monocle, metadata = pseudotime_start_male, col.name = "Pseudotime_Male")

#Female
pseudotime_start_female <- cds_female@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
female_for_monocle <- AddMetaData(female_for_monocle, metadata = pseudotime_start_female, col.name = "Pseudotime_Female")

#Save to CSV
male_data_to_write_out <- as.data.frame(as.matrix(male_for_monocle@meta.data))
fwrite(x = male_data_to_write_out, row.names = TRUE, file = "2023_1_23_Male_with_pseudotime.csv")

female_data_to_write_out <- as.data.frame(as.matrix(female_for_monocle@meta.data))
fwrite(x = female_data_to_write_out, row.names = TRUE, file = "2023_1_23_Female_with_pseudotime.csv")

#Plot female genes across pseudotime
DefaultAssay(female_for_monocle) <- "RNA"

all.genes.crypto.female <- rownames(female_for_monocle)

all.parasites.female <- colnames(female_for_monocle)

#Oocyst wall LOPIT for Figure 3d
list.of.oocyst.wall <- c("cgd1-3550",
                         "cgd1-800",
                         "cgd2-1590",
                         "cgd2-2510",
                         "cgd2-3040",
                         "cgd2-4350",
                         "cgd2-850",
                         "cgd3-1540",
                         "cgd3-1860",
                         "cgd3-190",
                         "cgd3-3700",
                         "cgd4-3090",
                         "cgd4-500",
                         "cgd4-670",
                         "cgd5-3073",
                         "cgd6-1450",
                         "cgd6-200",
                         "cgd6-2090",
                         "cgd6-210",
                         "cgd6-2470",
                         "cgd6-2900",
                         "cgd6-2920",
                         "cgd6-3730",
                         "cgd6-4440",
                         "cgd6-4640",
                         "cgd6-4840",
                         "cgd6-670",
                         "cgd6-710",
                         "cgd7-180",
                         "cgd7-1800",
                         "cgd7-4310",
                         "cgd7-4560",
                         "cgd7-5150",
                         "cgd7-5400",
                         "cgd7-850",
                         "cgd8-2670",
                         "cgd8-3350",
                         "cgd8-3870",
                         "cgd8-4230",
                         "cgd8-4660",
                         "cgd8-4830",
                         "cgd8-5080",
                         "cgd8-5090",
                         "cgd8-620")

oocyst_wall_pseudotime_female_results <- data.frame(matrix(ncol = 2136, nrow = 44))


rownames(oocyst_wall_pseudotime_female_results) <- c("cgd1-3550",
                                                     "cgd1-800",
                                                     "cgd2-1590",
                                                     "cgd2-2510",
                                                     "cgd2-3040",
                                                     "cgd2-4350",
                                                     "cgd2-850",
                                                     "cgd3-1540",
                                                     "cgd3-1860",
                                                     "cgd3-190",
                                                     "cgd3-3700",
                                                     "cgd4-3090",
                                                     "cgd4-500",
                                                     "cgd4-670",
                                                     "cgd5-3073",
                                                     "cgd6-1450",
                                                     "cgd6-200",
                                                     "cgd6-2090",
                                                     "cgd6-210",
                                                     "cgd6-2470",
                                                     "cgd6-2900",
                                                     "cgd6-2920",
                                                     "cgd6-3730",
                                                     "cgd6-4440",
                                                     "cgd6-4640",
                                                     "cgd6-4840",
                                                     "cgd6-670",
                                                     "cgd6-710",
                                                     "cgd7-180",
                                                     "cgd7-1800",
                                                     "cgd7-4310",
                                                     "cgd7-4560",
                                                     "cgd7-5150",
                                                     "cgd7-5400",
                                                     "cgd7-850",
                                                     "cgd8-2670",
                                                     "cgd8-3350",
                                                     "cgd8-3870",
                                                     "cgd8-4230",
                                                     "cgd8-4660",
                                                     "cgd8-4830",
                                                     "cgd8-5080",
                                                     "cgd8-5090",
                                                     "cgd8-620")

#Run loop to get expression data, only need to run pseudotime once
for (current.gene in 1:length(list.of.oocyst.wall)) {
  gene.name <- list.of.oocyst.wall[current.gene]
  FindMyGene <- which(all.genes.crypto.female==gene.name)
  all_expression_data_female <- GetAssayData(female_for_monocle, slot = 'data')[FindMyGene, all.parasites.female]
  all_pseudotime_female <- as.data.frame(as.matrix(female_for_monocle@meta.data[["Pseudotime_Female"]]), row.names = colnames(female_for_monocle))
  colnames(all_pseudotime_female) <- "Pseudotime"
  print(gene.name)
  as.vector(all_expression_data_female)
  oocyst_wall_pseudotime_female_results[current.gene,] <- all_expression_data_female
  
}


#Get barcode cell labels as column names
colnames(oocyst_wall_pseudotime_female_results) <- female_for_monocle@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
oocyst_wall_pseudotime_female_results_transposed_female <- as.data.frame(t(oocyst_wall_pseudotime_female_results))

#Scale the data
oocyst_wall_pseudotime_female_results_transposed_female <- as.data.frame(scale(oocyst_wall_pseudotime_female_results_transposed_female))

#Add pseudotime to column for plotting
oocyst_wall_pseudotime_female_results_transposed_female$Pseudotime <- all_pseudotime_female$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(oocyst_wall_pseudotime_female_results_transposed_female) <- c("cgd1_3550",
                                                                       "cgd1_800",
                                                                       "cgd2_1590",
                                                                       "cgd2_2510",
                                                                       "cgd2_3040",
                                                                       "cgd2_4350",
                                                                       "cgd2_850",
                                                                       "cgd3_1540",
                                                                       "cgd3_1860",
                                                                       "cgd3_190",
                                                                       "cgd3_3700",
                                                                       "cgd4_3090",
                                                                       "cgd4_500",
                                                                       "cgd4_670",
                                                                       "cgd5_3073",
                                                                       "cgd6_1450",
                                                                       "cgd6_200",
                                                                       "cgd6_2090",
                                                                       "cgd6_210",
                                                                       "cgd6_2470",
                                                                       "cgd6_2900",
                                                                       "cgd6_2920",
                                                                       "cgd6_3730",
                                                                       "cgd6_4440",
                                                                       "cgd6_4640",
                                                                       "cgd6_4840",
                                                                       "cgd6_670",
                                                                       "cgd6_710",
                                                                       "cgd7_180",
                                                                       "cgd7_1800",
                                                                       "cgd7_4310",
                                                                       "cgd7_4560",
                                                                       "cgd7_5150",
                                                                       "cgd7_5400",
                                                                       "cgd7_850",
                                                                       "cgd8_2670",
                                                                       "cgd8_3350",
                                                                       "cgd8_3870",
                                                                       "cgd8_4230",
                                                                       "cgd8_4660",
                                                                       "cgd8_4830",
                                                                       "cgd8_5080",
                                                                       "cgd8_5090",
                                                                       "cgd8_620",
                                                                       "Pseudotime")


#Plot in different colors, Figure 3d
ggplot(oocyst_wall_pseudotime_female_results_transposed_female, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_3550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_800), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_4350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_3700), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4560), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2470), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_1590), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2510), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_850), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1860), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3730), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4840), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_670), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_710), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5400), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2670), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4660), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4830), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3090), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_500), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_670), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_200), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2090), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_210), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1800), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5150), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3350), fill=TRUE, colour="hotpink1", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3040), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1540), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_190), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3073), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1450), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2900), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_2920), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4440), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_4640), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_180), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_850), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4230), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5080), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5090), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_620), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 6, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Dense granule LOPIT for Figure 3e
list.of.dense.granules <- c("cgd1-1230",
                            "cgd1-1250",
                            "cgd1-1433",
                            "cgd1-3780",
                            "cgd1-3790",
                            "cgd1-590",
                            "cgd1-620",
                            "cgd2-2530",
                            "cgd2-2580",
                            "cgd2-3030",
                            "cgd2-3210",
                            "cgd2-3360",
                            "cgd2-430",
                            "cgd3-1690",
                            "cgd3-1700",
                            "cgd3-1750",
                            "cgd3-530",
                            "cgd3-560",
                            "cgd3-600",
                            "cgd3-650",
                            "cgd3-693",
                            "cgd4-230",
                            "cgd4-2510",
                            "cgd4-3460",
                            "cgd4-3530",
                            "cgd4-3580",
                            "cgd4-3970",
                            "cgd4-850",
                            "cgd5-1440",
                            "cgd5-1450",
                            "cgd5-1480",
                            "cgd5-1490",
                            "cgd5-1530",
                            "cgd5-3210",
                            "cgd6-1030",
                            "cgd6-1110",
                            "cgd6-1170",
                            "cgd6-1180",
                            "cgd6-3050",
                            "cgd6-3080",
                            "cgd6-5280",
                            "cgd6-5310",
                            "cgd6-5330",
                            "cgd6-5360",
                            "cgd6-5420",
                            "cgd6-5440",
                            "cgd6-60",
                            "cgd7-1270",
                            "cgd7-1300",
                            "cgd7-1340",
                            "cgd7-2070",
                            "Cgd7-2213",
                            "cgd7-2340",
                            "cgd7-2620",
                            "cgd7-3800",
                            "cgd7-3810",
                            "cgd7-3820",
                            "cgd7-3860",
                            "cgd7-3870",
                            "cgd7-4280",
                            "cgd7-4340",
                            "cgd7-4420",
                            "cgd7-4490",
                            "cgd7-4500",
                            "cgd7-4530",
                            "cgd7-4680",
                            "cgd8-1780",
                            "cgd8-1800",
                            "cgd8-2160",
                            "cgd8-3030",
                            "cgd8-3600",
                            "cgd8-3610",
                            "cgd8-5300",
                            "cgd8-5310",
                            "cgd8-5350",
                            "cgd8-5360",
                            "cgd8-5380",
                            "cgd8-680",
                            "cgd8-690")


dense_granule_pseudotime_results_female <- data.frame(matrix(ncol = 2136, nrow = 79))

rownames(dense_granule_pseudotime_results_female) <- c("cgd1-1230",
                                                       "cgd1-1250",
                                                       "cgd1-1433",
                                                       "cgd1-3780",
                                                       "cgd1-3790",
                                                       "cgd1-590",
                                                       "cgd1-620",
                                                       "cgd2-2530",
                                                       "cgd2-2580",
                                                       "cgd2-3030",
                                                       "cgd2-3210",
                                                       "cgd2-3360",
                                                       "cgd2-430",
                                                       "cgd3-1690",
                                                       "cgd3-1700",
                                                       "cgd3-1750",
                                                       "cgd3-530",
                                                       "cgd3-560",
                                                       "cgd3-600",
                                                       "cgd3-650",
                                                       "cgd3-693",
                                                       "cgd4-230",
                                                       "cgd4-2510",
                                                       "cgd4-3460",
                                                       "cgd4-3530",
                                                       "cgd4-3580",
                                                       "cgd4-3970",
                                                       "cgd4-850",
                                                       "cgd5-1440",
                                                       "cgd5-1450",
                                                       "cgd5-1480",
                                                       "cgd5-1490",
                                                       "cgd5-1530",
                                                       "cgd5-3210",
                                                       "cgd6-1030",
                                                       "cgd6-1110",
                                                       "cgd6-1170",
                                                       "cgd6-1180",
                                                       "cgd6-3050",
                                                       "cgd6-3080",
                                                       "cgd6-5280",
                                                       "cgd6-5310",
                                                       "cgd6-5330",
                                                       "cgd6-5360",
                                                       "cgd6-5420",
                                                       "cgd6-5440",
                                                       "cgd6-60",
                                                       "cgd7-1270",
                                                       "cgd7-1300",
                                                       "cgd7-1340",
                                                       "cgd7-2070",
                                                       "Cgd7-2213",
                                                       "cgd7-2340",
                                                       "cgd7-2620",
                                                       "cgd7-3800",
                                                       "cgd7-3810",
                                                       "cgd7-3820",
                                                       "cgd7-3860",
                                                       "cgd7-3870",
                                                       "cgd7-4280",
                                                       "cgd7-4340",
                                                       "cgd7-4420",
                                                       "cgd7-4490",
                                                       "cgd7-4500",
                                                       "cgd7-4530",
                                                       "cgd7-4680",
                                                       "cgd8-1780",
                                                       "cgd8-1800",
                                                       "cgd8-2160",
                                                       "cgd8-3030",
                                                       "cgd8-3600",
                                                       "cgd8-3610",
                                                       "cgd8-5300",
                                                       "cgd8-5310",
                                                       "cgd8-5350",
                                                       "cgd8-5360",
                                                       "cgd8-5380",
                                                       "cgd8-680",
                                                       "cgd8-690")

#Run loop to get expression data
for (current.gene in 1:length(list.of.dense.granules)) {
  gene.name <- list.of.dense.granules[current.gene]
  FindMyGene <- which(all.genes.crypto.female==gene.name)
  all_expression_data_female <- GetAssayData(female_for_monocle, slot = 'data')[FindMyGene, all.parasites.female]
  print(gene.name)
  as.vector(all_expression_data_female)
  dense_granule_pseudotime_results_female[current.gene,] <- all_expression_data_female
  
}

#Get barcode cell labels as column names
colnames(dense_granule_pseudotime_results_female) <- female_for_monocle@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
dense_granule_pseudotime_results_transposed_female <- as.data.frame(t(dense_granule_pseudotime_results_female))

#Scale the data
dense_granule_pseudotime_results_transposed_female <- as.data.frame(scale(dense_granule_pseudotime_results_transposed_female))

#Add pseudotime to column for plotting
dense_granule_pseudotime_results_transposed_female$Pseudotime <- all_pseudotime_female$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(dense_granule_pseudotime_results_transposed_female) <- c("cgd1_1230",
                                                                  "cgd1_1250",
                                                                  "cgd1_1433",
                                                                  "cgd1_3780",
                                                                  "cgd1_3790",
                                                                  "cgd1_590",
                                                                  "cgd1_620",
                                                                  "cgd2_2530",
                                                                  "cgd2_2580",
                                                                  "cgd2_3030",
                                                                  "cgd2_3210",
                                                                  "cgd2_3360",
                                                                  "cgd2_430",
                                                                  "cgd3_1690",
                                                                  "cgd3_1700",
                                                                  "cgd3_1750",
                                                                  "cgd3_530",
                                                                  "cgd3_560",
                                                                  "cgd3_600",
                                                                  "cgd3_650",
                                                                  "cgd3_693",
                                                                  "cgd4_230",
                                                                  "cgd4_2510",
                                                                  "cgd4_3460",
                                                                  "cgd4_3530",
                                                                  "cgd4_3580",
                                                                  "cgd4_3970",
                                                                  "cgd4_850",
                                                                  "cgd5_1440",
                                                                  "cgd5_1450",
                                                                  "cgd5_1480",
                                                                  "cgd5_1490",
                                                                  "cgd5_1530",
                                                                  "cgd5_3210",
                                                                  "cgd6_1030",
                                                                  "cgd6_1110",
                                                                  "cgd6_1170",
                                                                  "cgd6_1180",
                                                                  "cgd6_3050",
                                                                  "cgd6_3080",
                                                                  "cgd6_5280",
                                                                  "cgd6_5310",
                                                                  "cgd6_5330",
                                                                  "cgd6_5360",
                                                                  "cgd6_5420",
                                                                  "cgd6_5440",
                                                                  "cgd6_60",
                                                                  "cgd7_1270",
                                                                  "cgd7_1300",
                                                                  "cgd7_1340",
                                                                  "cgd7_2070",
                                                                  "Cgd7_2213",
                                                                  "cgd7_2340",
                                                                  "cgd7_2620",
                                                                  "cgd7_3800",
                                                                  "cgd7_3810",
                                                                  "cgd7_3820",
                                                                  "cgd7_3860",
                                                                  "cgd7_3870",
                                                                  "cgd7_4280",
                                                                  "cgd7_4340",
                                                                  "cgd7_4420",
                                                                  "cgd7_4490",
                                                                  "cgd7_4500",
                                                                  "cgd7_4530",
                                                                  "cgd7_4680",
                                                                  "cgd8_1780",
                                                                  "cgd8_1800",
                                                                  "cgd8_2160",
                                                                  "cgd8_3030",
                                                                  "cgd8_3600",
                                                                  "cgd8_3610",
                                                                  "cgd8_5300",
                                                                  "cgd8_5310",
                                                                  "cgd8_5350",
                                                                  "cgd8_5360",
                                                                  "cgd8_5380",
                                                                  "cgd8_680",
                                                                  "cgd8_690",
                                                                  "Pseudotime")

#Plot dense granule across pseudotime, Figure 3e
ggplot(dense_granule_pseudotime_results_transposed_female, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1250), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_1433), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_3790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_590), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd1_620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_3360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_430), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1700), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_1750), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_560), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_600), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_650), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd3_693), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_230), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_2510), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3460), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3580), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_3970), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd4_850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1450), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1480), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_1530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd5_3210), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1170), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_1180), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3050), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_3080), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5330), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_5440), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_60), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2070), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=Cgd7_2213), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_2620), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3800), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3810), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3860), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_3870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4280), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4340), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4420), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4490), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4500), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4530), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1780), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_1800), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_2160), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3030), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3600), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_3610), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5350), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5360), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_5380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_690), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 6, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Crystalloid body LOPIT for Figure 3f
list.of.crystalloid.body <- c("cgd1-1510",
                              "cgd2-2100",
                              "cgd2-2110",
                              "cgd2-790",
                              "cgd6-313",
                              "cgd6-820",
                              "cgd7-1730",
                              "cgd7-300",
                              "cgd7-4810",
                              "cgd7-5140",
                              "cgd8-4290",
                              "cgd8-4300",
                              "cgd8-4310")

crystalloid_body_pseudotime_results_female <- data.frame(matrix(ncol = 2136, nrow = 13))

rownames(crystalloid_body_pseudotime_results_female) <- c("cgd1-1510",
                                                          "cgd2-2100",
                                                          "cgd2-2110",
                                                          "cgd2-790",
                                                          "cgd6-313",
                                                          "cgd6-820",
                                                          "cgd7-1730",
                                                          "cgd7-300",
                                                          "cgd7-4810",
                                                          "cgd7-5140",
                                                          "cgd8-4290",
                                                          "cgd8-4300",
                                                          "cgd8-4310")


#Run loop to get expression data
for (current.gene in 1:length(list.of.crystalloid.body)) {
  gene.name <- list.of.crystalloid.body[current.gene]
  FindMyGene <- which(all.genes.crypto.female==gene.name)
  all_expression_data_female <- GetAssayData(female_for_monocle, slot = 'data')[FindMyGene, all.parasites.female]
  print(gene.name)
  as.vector(all_expression_data_female)
  crystalloid_body_pseudotime_results_female[current.gene,] <- all_expression_data_female
  
}


#Get barcode cell labels as column names
colnames(crystalloid_body_pseudotime_results_female) <- female_for_monocle@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
crystalloid_body_pseudotime_results_transposed_female <- as.data.frame(t(crystalloid_body_pseudotime_results_female))

#Scale the data
crystalloid_body_pseudotime_results_transposed_female <- as.data.frame(scale(crystalloid_body_pseudotime_results_transposed_female))

#Add pseudotime to column for plotting
crystalloid_body_pseudotime_results_transposed_female$Pseudotime <- all_pseudotime_female$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(crystalloid_body_pseudotime_results_transposed_female) <- c("cgd1_1510",
                                                                     "cgd2_2100",
                                                                     "cgd2_2110",
                                                                     "cgd2_790",
                                                                     "cgd6_313",
                                                                     "cgd6_820",
                                                                     "cgd7_1730",
                                                                     "cgd7_300",
                                                                     "cgd7_4810",
                                                                     "cgd7_5140",
                                                                     "cgd8_4290",
                                                                     "cgd8_4300",
                                                                     "cgd8_4310",
                                                                     "Pseudotime")

#Plot crystalloid body across pseudotime, Figure 3f
ggplot(crystalloid_body_pseudotime_results_transposed_female, aes(Pseudotime)) + 
  stat_smooth(aes(y=cgd1_1510), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2100), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_2110), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd2_790), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_313), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd6_820), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_1730), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_4810), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd7_5140), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4290), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4300), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  stat_smooth(aes(y=cgd8_4310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 6, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Secretory Female Exclusive for Figure 3g
list.of.secretory.female <- c("cgd1-670",
                              "cgd1-3830",
                              "cgd4-1300",
                              "cgd4-2430",
                              "cgd5-1423",
                              "cgd6-5340",
                              "cgd7-1290",
                              "cgd8-2680",
                              "cgd8-1160",
                              "cgd5-743",
                              "cgd1-260",
                              "cgd4-4270",
                              "cgd5-1220",
                              "cgd5-3930",
                              "cgd7-5450",
                              "cgd8-3670",
                              "cgd7-4880",
                              "cgd2-3913",
                              "cgd4-1133",
                              "cgd4-2553",
                              "cgd7-400",
                              "cgd8-3133")

secretory_female_pseudotime_results_female <- data.frame(matrix(ncol = 2136, nrow = 22))

rownames(secretory_female_pseudotime_results_female) <- c("cgd1-670",
                                                          "cgd1-3830",
                                                          "cgd4-1300",
                                                          "cgd4-2430",
                                                          "cgd5-1423",
                                                          "cgd6-5340",
                                                          "cgd7-1290",
                                                          "cgd8-2680",
                                                          "cgd8-1160",
                                                          "cgd5-743",
                                                          "cgd1-260",
                                                          "cgd4-4270",
                                                          "cgd5-1220",
                                                          "cgd5-3930",
                                                          "cgd7-5450",
                                                          "cgd8-3670",
                                                          "cgd7-4880",
                                                          "cgd2-3913",
                                                          "cgd4-1133",
                                                          "cgd4-2553",
                                                          "cgd7-400",
                                                          "cgd8-3133")

#Run loop to get expression data
for (current.gene in 1:length(list.of.secretory.female)) {
  gene.name <- list.of.secretory.female[current.gene]
  FindMyGene <- which(all.genes.crypto.female==gene.name)
  all_expression_data_female <- GetAssayData(female_for_monocle, slot = 'data')[FindMyGene, all.parasites.female]
  print(gene.name)
  as.vector(all_expression_data_female)
  secretory_female_pseudotime_results_female[current.gene,] <- all_expression_data_female
  
}

#Get barcode cell labels as column names
colnames(secretory_female_pseudotime_results_female) <- female_for_monocle@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
secretory_female_pseudotime_results_transposed_female <- as.data.frame(t(secretory_female_pseudotime_results_female))

#Scale the data
secretory_female_pseudotime_results_transposed_female <- as.data.frame(scale(secretory_female_pseudotime_results_transposed_female))

#Add pseudotime to column for plotting
secretory_female_pseudotime_results_transposed_female$Pseudotime <- all_pseudotime_female$Pseudotime

#Rename columns because ggplot does not like dashes
colnames(secretory_female_pseudotime_results_transposed_female) <- c("cgd1_670",
                                                                     "cgd1_3830",
                                                                     "cgd4_1300",
                                                                     "cgd4_2430",
                                                                     "cgd5_1423",
                                                                     "cgd6_5340",
                                                                     "cgd7_1290",
                                                                     "cgd8_2680",
                                                                     "cgd8_1160",
                                                                     "cgd5_743",
                                                                     "cgd1_260",
                                                                     "cgd4_4270",
                                                                     "cgd5_1220",
                                                                     "cgd5_3930",
                                                                     "cgd7_5450",
                                                                     "cgd8_3670",
                                                                     "cgd7_4880",
                                                                     "cgd2_3913",
                                                                     "cgd4_1133",
                                                                     "cgd4_2553",
                                                                     "cgd7_400",
                                                                     "cgd8_3133",
                                                                     "Pseudotime")

#Female secretory fertilization receptor candidates (two waves), Figure 3g
ggplot(secretory_female_pseudotime_results_transposed_female, aes(Pseudotime)) + 
  geom_smooth(aes(y=cgd1_670), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd1_3830), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_1300), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_2430), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_1423), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_5340), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_1290), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_2680), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_1160), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_743), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd1_260), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_4270), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_1220), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_3930), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_5450), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_3670), fill=TRUE, colour="snow4", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 6, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Now plot male genes across pseudotime
DefaultAssay(male_for_monocle) <- "RNA"

all.genes.crypto.male <- rownames(male_for_monocle)

all.parasites.male <- colnames(male_for_monocle)


#Secretory and Membrane Male Exclusive for Figure 3l
list.of.secretory.membrane.genes.male <- c("cgd8-2220",
                                     "cgd3-720",
                                     "cgd2-2610",
                                     "cgd5-3570",
                                     "cgd7-5500",
                                     "cgd8-1740",
                                     "cgd3-1130",
                                     "cgd3-1830",
                                     "cgd3-3550",
                                     "cgd3-90",
                                     "cgd6-4010",
                                     "cgd4-3870",
                                     "cgd4-4380",
                                     "cgd2-2850",
                                     "cgd1-1260",
                                     "cgd5-2680",
                                     "cgd7-1240",
                                     "cgd8-2310",
                                     "cgd6-3060",
                                     "cgd6-380")

secretory_membrane_genes_pseudotime_results_male <- data.frame(matrix(ncol = 731, nrow = 20))

rownames(secretory_membrane_genes_pseudotime_results_male) <- c("cgd8-2220",
                                                          "cgd3-720",
                                                          "cgd2-2610",
                                                          "cgd5-3570",
                                                          "cgd7-5500",
                                                          "cgd8-1740",
                                                          "cgd3-1130",
                                                          "cgd3-1830",
                                                          "cgd3-3550",
                                                          "cgd3-90",
                                                          "cgd6-4010",
                                                          "cgd4-3870",
                                                          "cgd4-4380",
                                                          "cgd2-2850",
                                                          "cgd1-1260",
                                                          "cgd5-2680",
                                                          "cgd7-1240",
                                                          "cgd8-2310",
                                                          "cgd6-3060",
                                                          "cgd6-380")

#Run loop to get expression data, I only need to run pseudotime once though
for (current.gene in 1:length(list.of.secretory.membrane.genes.male)) {
  gene.name <- list.of.secretory.membrane.genes.male[current.gene]
  FindMyGene <- which(all.genes.crypto.male==gene.name)
  all_expression_data_male <- GetAssayData(male_for_monocle, slot = 'data')[FindMyGene, all.parasites.male]
  all_pseudotime_male <- as.data.frame(as.matrix(male_for_monocle@meta.data[["Pseudotime_Male"]]), row.names = colnames(male_for_monocle))
  colnames(all_pseudotime_male) <- "Pseudotime"
  print(gene.name)
  as.vector(all_expression_data_male)
  secretory_membrane_genes_pseudotime_results_male[current.gene,] <- all_expression_data_male
  
}

#Get barcode cell labels as column names
colnames(secretory_membrane_genes_pseudotime_results_male) <- male_for_monocle@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
secretory_membrane_genes_pseudotime_results_transposed_male <- as.data.frame(t(secretory_membrane_genes_pseudotime_results_male))

#Scale the data
secretory_membrane_genes_pseudotime_results_transposed_male <- as.data.frame(scale(secretory_membrane_genes_pseudotime_results_transposed_male))

#Add pseudotime to column after transposing 
secretory_membrane_genes_pseudotime_results_transposed_male$Pseudotime <- all_pseudotime_male$Pseudotime   

#Rename columns because ggplot does not like dashes
colnames(secretory_membrane_genes_pseudotime_results_transposed_male) <- c("cgd8_2220",
                                                                     "cgd3_720",
                                                                     "cgd2_2610",
                                                                     "cgd5_3570",
                                                                     "cgd7_5500",
                                                                     "cgd8_1740",
                                                                     "cgd3_1130",
                                                                     "cgd3_1830",
                                                                     "cgd3_3550",
                                                                     "cgd3_90",
                                                                     "cgd6_4010",
                                                                     "cgd4_3870",
                                                                     "cgd4_4380",
                                                                     "cgd2_2850",
                                                                     "cgd1_1260",
                                                                     "cgd5_2680",
                                                                     "cgd7_1240",
                                                                     "cgd8_2310",
                                                                     "cgd6_3060",
                                                                     "cgd6_380",
                                                                     "Pseudotime")

#Plot all together
ggplot(secretory_membrane_genes_pseudotime_results_transposed_male, aes(Pseudotime)) + 
  geom_smooth(aes(y=cgd8_2220), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd2_2610), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_3570), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_5500), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_1740), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_1130), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_1830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_3550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_90), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_4010), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_3870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_4380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd2_2850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd1_1260), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_2680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_1240), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_2310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_3060), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 


#Change size of graph to tall and thin, Figure 3l
#First wave
ggplot(secretory_membrane_genes_pseudotime_results_transposed_male, aes(Pseudotime)) + 
  geom_smooth(aes(y=cgd3_720), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_1130), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_4010), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_3870), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_4380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd2_2610), fill=TRUE, colour="#3FA9F5", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_5500), fill=TRUE, colour="#3FA9F5", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 15, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16))  +
  coord_cartesian(ylim=c(-1.5, 1.5))

#Second wave
ggplot(secretory_membrane_genes_pseudotime_results_transposed_male, aes(Pseudotime)) + 
  geom_smooth(aes(y=cgd2_2850), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd1_1260), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_1830), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_3550), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd3_90), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_2680), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_1240), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_2310), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_3060), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_380), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_2220), fill=TRUE, colour="black", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd5_3570), fill=TRUE, colour="#3FA9F5", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_1740), fill=TRUE, colour="#3FA9F5", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 15, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16))  +
  coord_cartesian(ylim=c(-1.5, 1.5))


#Signaling Male Exclusive for Figure 3n
list.of.signaling.male <- c("cgd2-1270",
                            "cgd6-4020",
                            "cgd1-2980",
                            "cgd7-2470",
                            "cgd8-5323",
                            "cgd4-3240")

signaling_pseudotime_results_male <- data.frame(matrix(ncol = 731, nrow = 6))

rownames(signaling_pseudotime_results_male) <- c("cgd2-1270",
                                                 "cgd6-4020",
                                                 "cgd1-2980",
                                                 "cgd7-2470",
                                                 "cgd8-5323",
                                                 "cgd4-3240")

#Run loop to get expression data
for (current.gene in 1:length(list.of.signaling.male)) {
  gene.name <- list.of.signaling.male[current.gene]
  FindMyGene <- which(all.genes.crypto.male==gene.name)
  all_expression_data_male <- GetAssayData(male_for_monocle, slot = 'data')[FindMyGene, all.parasites.male]
  print(gene.name)
  as.vector(all_expression_data_male)
  signaling_pseudotime_results_male[current.gene,] <- all_expression_data_male
  
}

#Get barcode cell labels as column names
colnames(signaling_pseudotime_results_male) <- male_for_monocle@assays[["RNA"]]@data@Dimnames[[2]]

#Transpose the extracted expression data, make sure it is a data frame
signaling_pseudotime_results_transposed_male <- as.data.frame(t(signaling_pseudotime_results_male))

#Scale the data
signaling_pseudotime_results_transposed_male <- as.data.frame(scale(signaling_pseudotime_results_transposed_male))

#Add pseudotime to column after transposing, for plotting
signaling_pseudotime_results_transposed_male$Pseudotime <- all_pseudotime_male$Pseudotime   

#Rename columns because ggplot does not like dashes
colnames(signaling_pseudotime_results_transposed_male) <- c("cgd2_1270",
                                                            "cgd6_4020",
                                                            "cgd1_2980",
                                                            "cgd7_2470",
                                                            "cgd8_5323",
                                                            "cgd4_3240",
                                                            "Pseudotime")

#Now plot, change size, Figure 3n, deleted cgd8_5323 and cgd4_3240 in figure to highlight other genes
ggplot(signaling_pseudotime_results_transposed_male, aes(Pseudotime)) + 
  geom_smooth(aes(y=cgd2_1270), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd6_4020), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd1_2980), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd7_2470), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd8_5323), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  geom_smooth(aes(y=cgd4_3240), fill=TRUE, colour="lightgray", se=FALSE, show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(0.5, 15, 0.5, 0.5, "cm")) + 
  scale_x_continuous(expand=c(0,0)) +
  ylab("Scaled Gene Expression") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="black")) +
  theme(axis.text.y=element_text(colour="black")) +
  theme(axis.text=element_text(size = 28),
        axis.title = element_text(size=32),
        plot.title=element_text(face = "bold", size = 20),
        legend.text = element_text(size=16)) 



