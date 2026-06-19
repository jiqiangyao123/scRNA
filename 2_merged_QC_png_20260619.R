#### This script is to provide a merged QC including all samples
args = commandArgs(trailingOnly=TRUE)
dir = args[1] 
v_ADT=as.logical(args[2]) # TRUE: ADT assay; FALSE: no ADT assay 
v_VDJ=as.logical(args[3]) # TRUE: VDJ assay; FALSE: no VDJ assay 

.libPaths("/usr/local/lib/R/library")
library(Polychrome) 
library(Seurat)
library(tibble)
library(ggplot2)
library(dplyr)
library(cowplot)
#library(UpSetR)
#library(tidyverse)


# setup dir, path to smaple list, name of singleR database to plot
if (Sys.info()['sysname'] == "Windows"){
  dir <- 'M:/dept/Dept_BBSR/Projects/Monteiro_Alvaro/3900_Monteiro_10x_2022/' # project dir
  
}


work_dir <- paste0(dir, "/work") # work dir
output_folder<- paste(work_dir,'QC/',sep="/") # define output dir
#input_folder<-"/iCel3_iFtu200_xMT10log10GPU0.8/" # define folder names of input .rds data
sample_list<-read.delim(paste(work_dir,"sample.list",sep="/"),header=F) # read in sample name
sample_list
singleR_type<-"hpca.main" # define name of singleR database
large_sample_cutoff<-12000 # define samples with large number of cells

tmp_size =ifelse(dim(sample_list)[1] > 10, 12 + dim(sample_list)[1]*0.5, 12)

# allow png
# R/4.0.2 on old cluster may have issue with generating png 
options(bitmapType='cairo')



# 1 merge raw metadata --------------------------------------------------------

sample_list <- sample_list[,1]
metadata_sum <- tibble()

for(i in sample_list){
  #seurat_obj <- readRDS(paste0(work_dir,"/individual/", i,input_folder, i, '_filtered_scaled_reduced_clustered.rds'))
  load(paste0(work_dir,"/individual/",i,"/",i,"_raw_metadata.RData"))
  metadata_sum <- rbind(metadata_sum, metadata_single)
}

## calculate log10GenesPerUMI_RNA
metadata_sum <- metadata_sum %>%
  mutate(log10GenesPerUMI_RNA = log10(nFeature_RNA) / log10(nCount_RNA))
save(metadata_sum, file = paste0(output_folder,'/metadata_sum.Rdata'))


# 2 define color ------------------------------------------------------------------
set.seed(86757)

sample_color <- createPalette(length(sample_list), c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))  
names(sample_color) <- sample_list
swatch(sample_color)

cell_type_color <- createPalette(length(unique(metadata_sum[,singleR_type])), c("#2a6ebb", "#de3831", "#007367"), range = c(30, 80))
names(cell_type_color) <- unique(metadata_sum[,singleR_type])
swatch(cell_type_color)

doublet_method_color <- c("darkgrey", "#440154", "#30678D", "#36B677","#FDE725")
names(doublet_method_color) <- c(as.character(seq(0, 4)))

save(sample_color, cell_type_color, doublet_method_color, file =  paste0(output_folder,'/color_book.Rdata'))


# 3 plot ------------------------------------------------------------------

load(paste0(output_folder,'/metadata_sum.Rdata'))
load( paste0(output_folder,'/color_book.Rdata'))
source(paste0(dir,'/scripts/scRNA_beta/S_QC_plots.R'))



# 3.1 scQC plots ----------------------------------------------------------

## sc QC plots
sc_QC_plots <- plot_sc_QC(metadata_sum, 
                       myUMIcutoff = 200,
                       myGenecutoff = 200,
                       mymitoratiocut = 15,
                       mylog10GenesPerUMIcut = 0.8,
                       sample_color,
                       cell_type_var = singleR_type)

for(i in 1:length(sc_QC_plots)){
  png(paste0(output_folder,names(sc_QC_plots)[i],".png"), width = tmp_size, height = tmp_size,res=300,units = "in")
  print(sc_QC_plots[[i]])
  dev.off()
}




# 3.2 ADT plot (optional) ------------------------------------------------------------


## ADT QC plots
if(v_ADT){
  ADT_QC_plots <- plot_ADT_QC(metadata_sum, 
                              sample_color)
  
  for(i in 1:length(ADT_QC_plots)){
    png(paste0(output_folder,names(ADT_QC_plots)[i],".png"), width = tmp_size, height = tmp_size,res=300,units = "in")
    print(ADT_QC_plots[[i]])
    dev.off()
  }
}


# 3.3 VDJ plot (optional) ------------------------------------------------------------


## VDJQC plots
if (v_VDJ){
  VDJ_QC_plots <- plot_VDJ_QC(metadata_sum, 
                              sample_color)
  
  for(i in 1:length(VDJ_QC_plots)){
    png(paste0(output_folder,names(VDJ_QC_plots)[i],".png"), width = tmp_size, height = tmp_size,res=300,units = "in")
    print(VDJ_QC_plots[[i]])
    dev.off()
  }
  
  
}

## doublet QC plots

# 3.4 Doublet QC ----------------------------------------------------------

#### count doublet methods

colnames(metadata_sum)

## for each cell, calculate the number of methods that detect doublet (scDbl, scrublet,DoubletFinder )
metadata_sum <- metadata_sum %>%
  dplyr::mutate(scDbl = as.integer(scDbl == "doublet"),   ## change doublet to 1
         scrublet.mask = as.integer(scrublet.mask == "doublet"),
         DoubletFinder.class = as.integer(DoubletFinder.class == "doublet")) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(doublet_pos_methods = sum(scDbl + scrublet.mask  + DoubletFinder.class))  ## calculate the methods that detect a doublet



save(metadata_sum, file =paste0(output_folder,'/metadata_sum.Rdata'))

## plot all samples

double_method_count_plots <- plot_doublet_method_count(doublet_df = metadata_sum,
                                                       doublet_color = doublet_method_color,
                                                       cell_type_var = singleR_type, 
                                                       cluster_var = NULL, 
                                                       method_n_var = 'doublet_pos_methods',
                                                       myUMIcutoff=200,
                                                       myGenecutoff=200,
                                                       my_ncol = 3)


# the width and height need to be modified based on the sample size and my_ncol
for(i in 1:length(double_method_count_plots)){
  png(paste0(output_folder,names(double_method_count_plots)[i],".png"), width = tmp_size, height = tmp_size,res=300,units = "in")
  print(double_method_count_plots[[i]])
  dev.off()
}

## plot samples with a high cell count (optional)
selected_samples <- metadata_sum %>%
  dplyr::select(orig.ident) %>%
  group_by(orig.ident) %>%
  summarise(cell_count = n()) %>%
  filter(cell_count > large_sample_cutoff) %>%
  select(orig.ident) %>%
  unlist()

if(length(selected_samples)>0){
  
 double_method_count_plots_selected <- plot_doublet_method_count(  
 doublet_df = metadata_sum %>% filter(orig.ident %in% selected_samples),
  doublet_color = doublet_method_color,
  cell_type_var = singleR_type, 
  cluster_var = NULL, 
  method_n_var = 'doublet_pos_methods',
  myUMIcutoff = 200,
  myGenecutoff = 200,
  my_ncol = 1)

  # the width and height need to be modified based on the sample size and my_ncol
  for(i in 1:length(double_method_count_plots_selected)){
    png(paste0(output_folder,names(double_method_count_plots_selected)[i],".png"), width = tmp_size, height = tmp_size*2,res=300,units = "in")
    print(double_method_count_plots_selected[[i]])
    dev.off()
  }
  
}


##upsetR

if (FALSE){
  metadata_sum_UpSetR <- metadata_sum %>%
    select(scDbl, scrublet.mask, DoubletFinder.class) %>%
    rownames_to_column(var = "cell_id") %>%
    mutate(cell_id = as.character(cell_id)) %>%
    as.data.frame()
  
  
  png(paste0(output_folder,"/double_method_overlapping_count_upsetR.png"), width = tmp_size, height = tmp_size,res=300,units = "in")
  upset(metadata_sum_UpSetR,
        sets = c("scDbl", "scrublet.mask",  "DoubletFinder.class"), 
        order.by="degree", 
        matrix.color="blue", 
        point.size=5)
  dev.off()
  
  png(paste0(output_folder,"/double_method_overlapping_count_upsetR.2.png"), width = tmp_size, height = tmp_size,res=300,units = "in")
  upset(metadata_sum_UpSetR,
        sets = c("scDbl", "scrublet.mask", "DoubletFinder.class"), 
        order.by="freq", 
        matrix.color="blue", 
        point.size=5)
  dev.off()
}








# 4 apply additional filters and retest QC plots  ----------------------------



load(paste0(output_folder,'/metadata_sum.Rdata'))
load(paste0(output_folder,'/color_book.Rdata'))
source(paste0(dir,'/scripts/scRNA_beta/S_QC_plots.R'))

# customize the filter criteria
# doublet_pos_methods < 2: remove called by >=2 callers
# log10GenesPerUMI_RNA > 0.8: remove log10GenesPerUMI_RNA<=0.8
# log10GenesPerUMI_RNA > 0.8: remove log10GenesPerUMI_RNA<=0.8
# nFeature_RNA >200: remove nFeature_RNA<=200
# percent.mt<15: remove percent.mt>=15

metadata_sum_new <- metadata_sum %>%
  filter(doublet_pos_methods < 2 & log10GenesPerUMI_RNA > 0.8 & nFeature_RNA >200 & percent.mt<15)


##scRNA qc
sc_QC_plots <- plot_sc_QC(metadata_sum_new, 
                          myUMIcutoff = 200,
                          myGenecutoff = 200,
                          mymitoratiocut = 15,
                          mylog10GenesPerUMIcut = 0.8,
                          sample_color,
                          cell_type_var = singleR_type)

for(i in 1:length(sc_QC_plots)){
  png(paste0(output_folder,names(sc_QC_plots)[i],".new_filtered.png"), width = tmp_size, height = tmp_size,res=300,units = "in")
  print(sc_QC_plots[[i]])
  dev.off()
}


##scADT qc

## ADT QC plots
if (v_ADT){
  ADT_QC_plots <- plot_ADT_QC(metadata_sum_new, 
                              sample_color)
  
  for(i in 1:length(ADT_QC_plots)){
    png(paste0(output_folder,names(ADT_QC_plots)[i],".new_filtered.png"), width = tmp_size, height = tmp_size,res=300,units = "in")
    print(ADT_QC_plots[[i]])
    dev.off()
  }
  
}


## VDJQC plots

if (v_VDJ){
  VDJ_QC_plots <- plot_VDJ_QC(metadata_sum_new, 
                              sample_color)
  
  for(i in 1:length(VDJ_QC_plots)){
    png(paste0(output_folder,names(VDJ_QC_plots)[i],".new_filtered.png"), width = tmp_size, height = tmp_size,res=300,units = "in")
    print(VDJ_QC_plots[[i]])
    dev.off()
  }
}



## doublet
double_method_count_plots <- plot_doublet_method_count(doublet_df = metadata_sum_new,
                                                       doublet_color = doublet_method_color,
                                                       myUMIcutoff = 200,
                                                       myGenecutoff = 200,
                                                       cell_type_var =singleR_type, 
                                                       cluster_var = NULL, 
                                                       method_n_var = 'doublet_pos_methods',
                                                       my_ncol = 3)


# the width and height need to be modified based on the sample size and my_ncol
for(i in 1:length(double_method_count_plots)){
  png(paste0(output_folder,names(double_method_count_plots)[i],".new_filtered.png"), width = tmp_size, height = tmp_size,res=300,units = "in")
  print(double_method_count_plots[[i]])
  dev.off()
}

