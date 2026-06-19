#### This script is to pre-process individule sample

args = commandArgs(trailingOnly=TRUE)

v_MT=as.numeric(args[1]) # 100: no filter
v_minCell=as.numeric(args[2])  # 0: no filter
v_minFeatures=as.numeric(args[3])  # 0: no filter
v_log10GenesPerUMI_RNA=as.numeric(args[4])  # 1: no filter
v_nfeatures=as.numeric(args[5])
v_runpca=as.numeric(args[6])
v_UMAPPCA=as.numeric(args[7])
v_resolution=as.numeric(args[8])
sample_name=args[9]
dir = args[10]  # path to project folder
species=args[11] # human or mouse
v_ADT=as.logical(args[12]) # TRUE: ADT assay; FALSE: no ADT assay 
v_VDJ=as.logical(args[13]) # TRUE: VDJ assay; FALSE: no VDJ assay 
v_multiome=as.logical(args[14]) # TRUE: Multiome assay


v_MT
v_minCell
v_minFeatures
v_log10GenesPerUMI_RNA
v_nfeatures
v_runpca
v_UMAPPCA
v_resolution
sample_name
dir
species
v_ADT
v_VDJ
v_multiome

library(Polychrome) # need to install 

.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))
library(SingleR)
library(RColorBrewer)
library(reshape2)
library(Seurat) 
library(dplyr)
library(readxl)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(tidyr)
library(reticulate)

# allow png
# R/4.0.2 on old cluster may have issue with generating png 
options(bitmapType='cairo')

# 0. define color ------------------------------------------------------------------
set.seed(86757)

cell_type_color <- createPalette(50, c("#2a6ebb", "#de3831", "#007367"), range = c(30, 80))

doublet_method_color <- c("darkgrey", "#440154", "#30678D", "#36B677")
names(doublet_method_color) <- c(as.character(seq(0, 3)))

QC_colors_list <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid","orange", "gold", "gray")

# 1. set up working directory and folder name ------------------------------------------------------------------
work_dir <- paste0(dir, "/work/individual/")
#folder <- paste0(sample_name,"-Count" )
folder <- sample_name
date <- Sys.Date()

wd_name <- paste0(work_dir,"/",sample_name,"/","iCel",v_minCell,"_iFtu",v_minFeatures,"_xMT",v_MT,"log10GPU",v_log10GenesPerUMI_RNA)
if (!dir.exists(wd_name)) dir.create(wd_name)
setwd(wd_name)

# 2. create seurat object ------------------------------------------------------------------
data_dir <- paste(dir, "/data/", folder, "/outs/filtered_feature_bc_matrix", sep = "")
data_dir
data <- Read10X(data.dir = data_dir)
  
if (v_ADT){
  # Initialize the Seurat object with the raw (non-normalized data)
  sample.obj <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = v_minCell, project = sample_name)
  #add Antidoby Capture data inot Seurat object
  sample.obj[["ADT"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
    
  #ADT_feature=rownames(sample.obj[["ADT"]])
  raw_cell_number = dim(data$`Gene Expression`)[2]
    
}else if (v_multiome==TRUE){
    sample.obj <- CreateSeuratObject(counts = data$'Gene Expression', min.cells = v_minCell, project = sample_name)
    raw_cell_number = dim(data$'Gene Expression')[2]
}else{
    sample.obj <- CreateSeuratObject(counts = data, min.cells = v_minCell, project = sample_name)
    raw_cell_number = dim(data)[2]
}

raw_cell_number

# 3. plot QC metrics ------------------------------------------------------------------
## add cell cycle scores and MT, and log10GenesPerUMI_RNA
sample.obj[["log10GenesPerUMI_RNA"]] <-log10(sample.obj$nFeature_RNA) / log10(sample.obj$nCount_RNA)

## normalize before assigning cell cycle to aviod error in certain Seurat vesion
sample.obj <- NormalizeData(object = sample.obj, verbose = FALSE)

if (tolower(species)=="human"){
  load(paste0(dir,"/scripts/scRNA_beta/regev_lab_cell_cycle_genes_human.RData"))
  sample.obj[["percent.mt"]] <- PercentageFeatureSet(sample.obj, pattern = "^MT-")  
}else{
  load(paste0(dir,"/scripts/scRNA_beta/regev_lab_cell_cycle_genes_mouse.RData"))
  sample.obj[["percent.mt"]] <- PercentageFeatureSet(sample.obj, pattern = "^mt-")
}
sample.obj <- CellCycleScoring(object = sample.obj, s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)

## QC vlnplot

p3.1<-VlnPlot(sample.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI_RNA","S.Score"),
              ncol = 5, pt.size = 0,
        group.by = "orig.ident",cols = "gray") 
png(paste(sample_name, "_raw_QC_vln_", date, ".png",sep = ""), width = 12, height =12,res=300,units = "in")
print(p3.1)
dev.off()

if(v_ADT){
  p3.2<-VlnPlot(sample.obj, features = c("nFeature_ADT", "nCount_ADT"), ncol = 2, pt.size = 0,
              group.by = "orig.ident",cols = "gray") 
  png(paste(sample_name, "_raw_ADT_QC_vln_", date, ".png",sep = ""), width = 4, height =4,res=300,units = "in")
  print(p3.2)
  dev.off()
  
}

if(v_VDJ){
  folder.vdj=paste0(sample_name,"-VDJ-Count" )
  vdj <- read.csv(paste(dir, "/data/", folder.vdj, "/outs/filtered_contig_annotations.csv", sep = ""), header = T)
  bc.vdj <- gsub("-1","-1",unique(as.character(vdj[,"barcode"])))
  sample.obj@meta.data$VDJ <- "NA"
  sample.obj@meta.data$VDJ[match(bc.vdj, names(sample.obj@active.ident))] <- "Y"

  ## map clonotypes ID to cells
  clonotype <- read.csv(paste(dir, "/data/", folder.vdj, "/outs/clonotypes.csv", sep = ""), header = T)
  clonotype <- clonotype[order(clonotype$frequency, decreasing = T), ]
  #clonotype_ids <- as.character(clonotype$clonotype_id)
  
  bc <- list()
  sample.obj@meta.data$VDJ_clonotype <- "NA"
  sample.obj@meta.data$VDJ_cdr3a_aa <- "NA"
  
  for (n in 1:length(clonotype$clonotype_id)) {
    bc[[n]] <- as.character(vdj$barcode)[as.character(vdj$raw_clonotype_id)==clonotype$clonotype_id[n]]
    bc[[n]] <- gsub("-1", "-1", bc[[n]]) %>% unique()
    sample.obj@meta.data$VDJ_clonotype[match(bc[[n]], names(sample.obj@active.ident))] <- clonotype$clonotype_id[n]
    sample.obj@meta.data$VDJ_cdr3a_aa[match(bc[[n]], names(sample.obj@active.ident))] <- clonotype$cdr3s_aa[n]
  }
  
}

##  nFeature_RNA vs. nCount_RNA scatter plot
meta<-sample.obj@meta.data
p3.3<-meta %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  stat_smooth(method = lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  scale_colour_gradient(low = "gray90", high = "black") 

#p3<-QC_Plot_UMIvsGene(seurat_object = sample.obj, meta_gradient_name = "percent.mt")
png(paste(sample_name, "_raw_QC_scatter1_", date, ".png",sep = ""), width =7, height =5,res=300,units = "in")
print(p3.3)
dev.off()

p3.4<-meta %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = log10GenesPerUMI_RNA)) + 
  geom_point() + 
  stat_smooth(method = lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  scale_colour_gradient(low = "gray90", high = "black") 

#p4<-QC_Plot_UMIvsGene(seurat_object = sample.obj, meta_gradient_name = "log10GenesPerUMI_RNA")
png(paste(sample_name, "_raw_QC_scatter2_", date, ".png",sep = ""), width = 7, height =5,res=300,units = "in")
print(p3.4)
dev.off()

# 4. add doublet results ------------------------------------------------------------------ 

filename = paste0(dir,"/work/dblet/",sample_name, ".metadat.doublet.RData")
if (file.exists(filename)){
  load(filename)
  sample.obj <- AddMetaData(object=sample.obj, metadata=sample.list.metadat[,-c(1:3)])
  
}else{
  sample.obj$scDbl=sample.obj$scDblFinder.score=sample.obj$scrublet.score=sample.obj$scrublet.mask=
    sample.obj$doubletCell.score=sample.obj$doubletCell.mask=sample.obj$DoubletFinder.score=sample.obj$DoubletFinder.class<-NA
}

sample.obj@meta.data$doublet_num <-apply(sample.obj@meta.data[,c("scDbl","scrublet.mask","DoubletFinder.class")], 1, function(x) sum(x=="Doublet") + sum(x=="doublet"))
table(sample.obj@meta.data$doublet_num)


#5. singleR
source(paste0(dir,"/scripts/scRNA_beta/singleR.R"))
if (tolower(species)=="human"){
 singler.preds <- do.singler(sample.obj,single.ref.rds,ref.by.species, species="Human")
 preds.labels <- preds.to.labels(singler.preds[,ref.by.species[["Human"]]])
 sample.obj <- AddMetaData(object=sample.obj, metadata=preds.labels)

 }else{
  singler.preds <- do.singler(sample.obj,single.ref.rds,ref.by.species, species="Mouse")
  preds.labels <- preds.to.labels(singler.preds[,ref.by.species[["Mouse"]]])
  sample.obj <- AddMetaData(object=sample.obj, metadata=preds.labels)

}


# save raw metadata
metadata_single<-sample.obj@meta.data
save(metadata_single,file=paste0("../",sample_name,"_raw_metadata.RData"))


# 6. add filter ------------------------------------------------------------------
if( !(v_minCell==0 && v_minFeatures==0 && v_MT==100 && v_log10GenesPerUMI_RNA==1)){ 
  
  # save a copy of unfiltered
  #saveRDS(sample.obj, paste0(sample_name, "_unfiltered.rds"))
  
 #feature_counts <- rowSums(as.matrix(sample.obj@assays$RNA@counts)>0)
  feature_counts <- rowSums(as.matrix(GetAssayData(sample.obj, assay = "RNA", layer = "counts")) > 0)
  if(v_ADT){
    sample.obj <- subset(sample.obj, subset = nFeature_RNA > v_minFeatures & percent.mt < v_MT & log10GenesPerUMI_RNA > v_log10GenesPerUMI_RNA, 
                         features=c(names(feature_counts)[feature_counts>v_minCell], rownames(sample.obj[["ADT"]])))
  }else{
    sample.obj <- subset(sample.obj, subset = nFeature_RNA > v_minFeatures & percent.mt < v_MT & log10GenesPerUMI_RNA > v_log10GenesPerUMI_RNA, features=names(feature_counts)[feature_counts>v_minCell])
  }
  
  filter_cell_number = length(sample.obj@active.ident)
  tmpdf <- data.frame( sample_name = sample_name, raw_cell_count=raw_cell_number,filter_cell_number = filter_cell_number)
  write.table(tmpdf, file=paste0(sample_name,"_filter_cell_number.txt"),quote=FALSE, sep="\t", row.names = FALSE, col.names= FALSE)
  
  # QC plot of filtered data
  
  p6.1<-VlnPlot(sample.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI_RNA","S.Score"), ncol = 5, pt.size = 0,
          group.by = "orig.ident",cols = "gray") 
  png(paste(sample_name, "_filter_QC_vln_", date, ".png",sep = ""), width = 16, height =8,res=300,units = "in")
  print(p6.1)
  dev.off()
  
  meta<-sample.obj@meta.data
  p6.2<-meta %>% 
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
    geom_point() + 
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    scale_colour_gradient(low = "gray90", high = "black") 
  
  #p3<-QC_Plot_UMIvsGene(seurat_object = sample.obj, meta_gradient_name = "percent.mt")
  png(paste(sample_name, "_filter_QC_scatter1_", date, ".png",sep = ""), width =7, height =5,res=300,units = "in")
  print(p6.2)
  dev.off()
  
  p6.3<-meta %>% 
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = log10GenesPerUMI_RNA)) + 
    geom_point() + 
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    scale_colour_gradient(low = "gray90", high = "black") 
  
  #p4<-QC_Plot_UMIvsGene(seurat_object = sample.obj, meta_gradient_name = "log10GenesPerUMI_RNA")
  png(paste(sample_name, "_filter_QC_scatter2_", date, ".png",sep = ""), width = 7, height =5,res=300,units = "in")
  print(p6.3)
  dev.off()
  
  if(v_ADT){
    p6.4<-VlnPlot(sample.obj, features = c("nFeature_ADT", "nCount_ADT"), ncol = 2, pt.size = 0,
                  group.by = "orig.ident",cols = "gray") 
    png(paste(sample_name, "_filter_ADT_QC_vln_", date, ".png",sep = ""), width = 4, height =4,res=300,units = "in")
    print(p6.4)
    dev.off()
    
  }
  
}else{
  sample.obj<-sample.obj # no filter
}


# 7. Normalization, find variable features------------------------------------------------------------------
sample.obj <- FindVariableFeatures(object = sample.obj, selection.method = "vst", nfeatures = v_nfeatures, verbose = FALSE)
top10 <- head(VariableFeatures(sample.obj), 10)


## plot variable features with labels
p7.1 <- VariableFeaturePlot(sample.obj)
p7.2 <- LabelPoints(plot = p7.1, points = top10, repel = TRUE)
png(paste(sample_name, "_QC_variable_features_", date, ".png",sep = ""), width = 6, height = 5, res=300, units = "in")
print(p7.2)
dev.off()

## remove TCR related gene and IG genes from variableFeatures
Vgene<-VariableFeatures(sample.obj)
if (tolower(species)=="human"){
  Vgene.remove<-Vgene[!grepl("^IG",Vgene) & !grepl("^TR",Vgene) & !grepl("^RP",Vgene)]
}else{
  Vgene.remove<-Vgene[!grepl("^Ig",Vgene) & !grepl("Tra|Trb|Tcr",Vgene) & !grepl("^Rp",Vgene)]
}
VariableFeatures(sample.obj)<-Vgene.remove # put variable genes back

# 8. data scale, dimention reduction ------------------------------------------------------------------
sample.obj <- ScaleData(sample.obj, vars.to.regress = c("S.Score","G2M.Score","percent.mt", "nCount_RNA")) 
sample.obj <- RunPCA(sample.obj, npcs = v_runpca, features = VariableFeatures(object = sample.obj))

png(paste(sample_name, "_PCA_", date, ".png",sep = ""), width =7, height = 9, res=300, units = "in")
p8.1 <- VizDimLoadings(sample.obj, dims = 1:4, reduction = "pca")
print(p8.1)
dev.off()

png(paste(sample_name, "_PCA_heatmap_", date, ".png",sep = ""), width = 7, height = 9, res=300, units = "in")
p8.2<-DimHeatmap(sample.obj, dims = 1:4, cells = 500, balanced = TRUE, ncol = 2)
print(p8.2)
dev.off()

## UMAP
sample.obj <- RunUMAP(sample.obj, reduction = "pca", dims = 1:v_UMAPPCA)
#sample.obj@meta.data$umap_1 <- sample.obj@reductions$umap@cell.embeddings[, "UMAP_1"]
#sample.obj@meta.data$umap_2 <- sample.obj@reductions$umap@cell.embeddings[, "UMAP_2"]
umap_coords <- Embeddings(sample.obj, reduction = "umap")
sample.obj@meta.data$umap_1 <- umap_coords[, 1]
sample.obj@meta.data$umap_2 <- umap_coords[, 2]


# 9. plot QC on UMAP ------------------------------------------------------------------
## 9.1 basic QC
p9.1<-FeaturePlot(sample.obj,features = c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score","log10GenesPerUMI_RNA"))
#                  cols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
png(paste(sample_name, "_QC_UMAP_", date, ".png",sep = ""), width = 10, height = 10, res=300, units = "in")
print(p9.1)
dev.off()

## 9.2 doublet
sample.obj@meta.data$doublet_num <-apply(sample.obj@meta.data[,c("scDbl","scrublet.mask","DoubletFinder.class")], 1, function(x) sum(x=="Doublet") + sum(x=="doublet"))
table(sample.obj@meta.data$doublet_num)

png(paste(sample_name, "_filtered_doublet_overlap_", date, ".png",sep = ""), width = 6, height = 5, res=300, units = "in")
p9.2<- DimPlot(sample.obj, reduction = "umap",group.by="doublet_num",pt.size=1,cols=doublet_method_color)
p9.2[[1]]$layers[[1]]$aes_params$alpha = ifelse ( sample.obj@meta.data$doublet_num==0, 0.2, 1 )
print(p9.2)
dev.off()

## 9.3 singleR

# assign singleR estimation to cells

# assign singleR estimation to cells
for (m in 1:length(preds.labels)){
  set.seed(123)
  singleR_color= createPalette(length(unique(preds.labels[m][[1]])), c("#2a6ebb", "#de3831", "#007367"), range = c(30, 80))
  names(singleR_color)=unique(preds.labels[m][[1]])
  if(length(singleR_color)<18){
    pngwidth=6
  }else if(length(singleR_color)<=40){
    pngwidth=12
  }else if(length(singleR_color)<=80) {
    pngwidth=20
  }else{
    pngwidth=35
  }

  png(paste("cell_types_",names(preds.labels)[m],"_",date,".png", sep=""), width = pngwidth, height = pngwidth,res=300,units = "in")
  p1<-DimPlot(sample.obj, group.by = names(preds.labels)[m],pt.size=1,  cols=singleR_color)
  print(p1)
  dev.off()
}



# 10. clustering analysis ------------------------------------------------------------------
sample.obj <- FindNeighbors(object = sample.obj, reduction = "pca", dims = 1:v_UMAPPCA, force.recalc=TRUE)
sample.obj <- FindClusters(object = sample.obj, resolution = v_resolution)
saveRDS(sample.obj, paste0(sample_name, "_filtered_scaled_reduced_clustered.rds"))

## 10.1. plot QC by cluster
cluster_color<-createPalette(length(levels(Idents(sample.obj))), c("#2a6ebb", "#de3831", "#007367"), range = c(30, 80))
names(cluster_color)<-levels(Idents(sample.obj))

p10.1<-VlnPlot(sample.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI_RNA","S.Score","doublet_num"),
                    cols = cluster_color,pt.size = 0, ncol = 1)+theme(axis.title.x = element_blank())
png(paste(sample_name, "_vln_cluster_", date, ".png",sep = ""), width = length(levels(Idents(sample.obj)))/2, height = 15, res=300, units = "in")
print(p10.1)
dev.off()

## 10.2. cluster UMAP
p10.2<-DimPlot(sample.obj, cols= cluster_color,label = T)
png(paste(sample_name, "_UMAP_cluster_", date, ".png",sep = ""), width = 6, height = 5, res=300, units = "in")
print(p10.2)
dev.off()


# 11. ADT  (optional)  ------------------------------------------------------------------
if(v_ADT){
  ADT_feature=rownames(sample.obj[["ADT"]])
  ADT_cmt<-length(ADT_feature)
  
  ## normalize ADT daa
  sample.obj <- NormalizeData(sample.obj, assay = "ADT", normalization.method = "CLR")
  sample.obj <- ScaleData(sample.obj, assay = "ADT")
  
  DefaultAssay(sample.obj)<-"ADT"
  
  # choose singleR for ADT
  if (tolower(species)=="human"){
    singleR.type<-"monaco.main"
  }else{
    singleR.type<-"immGen.main"
  }
  
  ## plot ADT on umap
  

  for(i in 1:ceiling(ADT_cmt/8)){
    p11.1<-FeaturePlot(sample.obj, features = ADT_feature[c((8*(i-1)+1):(8*i))], ncol = 4,max.cutoff = "q98")
    p11.2<-RidgePlot(sample.obj, features = ADT_feature[c((8*(i-1)+1):(8*i))], cols = cluster_color , ncol = 4)
   #for(j in 1:length(p11.2)){
   #   p11.2[[j]]<-p11.2[[j]]+theme(plot.title = element_text(size=10))
   #}
    singleR_color=cluster_color[1:length(unique(sample.obj@meta.data[,singleR.type]))]
    names(singleR_color)<-unique(unique(sample.obj@meta.data[,singleR.type]))
    p11.3<-RidgePlot(sample.obj, features = ADT_feature[c((8*(i-1)+1):(8*i))], group.by = singleR.type,cols = singleR_color,ncol = 4 )
    
    png(paste(sample_name, "_UMAP_ADT_",i,"_" ,date, ".png",sep = ""), width = 12, height = 6, res=300, units = "in")
    print(p11.1)
    dev.off()
    
    png(paste(sample_name, "_ADT_RidgePlot_cluster_",i,"_" ,date, ".png",sep = ""), width = 12, height =length(levels(Idents(sample.obj)))/2, res=300, units = "in")
    print(p11.2)
    dev.off()
    
    png(paste(sample_name, "_ADT_RidgePlot_",singleR.type,"_",i,"_" ,date, ".png",sep = ""), width = 12, height =7, res=300, units = "in")
    print(p11.3)
    dev.off()
 
  }

}



# 12. VDJ (optional)  ------------------------------------------------------------------
if(v_VDJ){
  TCR_cell_number = sum(sample.obj@meta.data$VDJ=="Y")
  
  ## plot cells with vdj data
  png(paste0(sample_name, "_umap_VDJ_", date, ".png"), width = 6, height = 5, res=300, units = "in")
  p12.1 <- ggplot(sample.obj@meta.data, aes(x = umap_1, y = umap_2, color = VDJ)) + geom_point(size = 2) +
    scale_color_manual(breaks = c("Y", "NA"),values = c("grey", "red")) +
    theme_cowplot()
  print(p12.1)
  dev.off()
  
  
  ## plot top 4 clonotypes
clonotype_ids<-clonotype$clonotype_id
  meta<-sample.obj@meta.data
  meta$top<- "Others"
  meta$top[sample.obj$VDJ_clonotype%in%clonotype_ids[1:4]]<-meta$VDJ_clonotype[meta$VDJ_clonotype%in%clonotype_ids[1:4]]
  meta$top<-factor(meta$top,levels=c(clonotype_ids[1:4],"Others"))
  cbPalette <- c(QC_colors_list[1:4],"grey80")
  p12.2<- ggplot(meta %>% arrange(desc(top)), aes(x = umap_1, y = umap_2, color = top)) + geom_point(size = 2) +
    scale_color_manual(values = cbPalette) +
    theme_light()
  mytable <- table(meta$top)
  mytable <- data.frame(mytable)[1:4,]; colnames(mytable) <- c("VDJ_clonotype", "Freq")
  p12.2<-p12.2 + annotation_custom(tableGrob(mytable, theme = ttheme_default(base_size = 10), rows = NULL), xmin = max(meta$umap_1)+(max(meta$umap_1)-min(meta$umap_1))/4, xmax = max(meta$umap_1), ymin =  max(meta$umap_2)-(max(meta$umap_2)-min(meta$umap_2))/4, ymax =  max(meta$umap_2))

  
  png(paste0(sample_name, "_umap_top_TCR_clonotype_", date, ".png"), width = 7, height = 6, res=300, units = "in")
  print(p12.2)
  dev.off()
}


saveRDS(sample.obj, paste0(sample_name, "_filtered_scaled_reduced_clustered.rds"))

