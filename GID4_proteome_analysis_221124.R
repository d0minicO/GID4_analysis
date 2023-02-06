## script for analysis of GID4 proteomics data Owens et al 2022

## Fully working 221124


library(tidyverse)
library(magrittr)
library(DEP)
library(data.table)
library(biomaRt)
library(patchwork)
library(RColorBrewer)
library(SummarizedExperiment)
library(pheatmap)
library(ggrepel)
library(ggrastr)
library(Biostrings)
library(ggtext)



#################
#### OPTIONS ####
#################

options(max.print=200)



###################
#### Functions ####
###################

# source custom functions
devtools::source_url("https://github.com/d0minicO/phosphoR/blob/main/customFunctions.R?raw=TRUE")

devtools::source_url("https://github.com/d0minicO/unimachR/blob/main/unimachR.R?raw=TRUE")



################
#### inputs ####
################


base = "C:/Users/dowens/OneDrive/GID4_analysis/"

write.table(sess_dat,paste0(base,"SessionInfo.txt"),
            quote=F,
            sep="\t")

in_file = paste0(base,"RAWproteome.txt")


## full list of deleted uniprot IDs already cleaned up
## downloaded 12th October 2022 https://www.uniprot.org/help/deleted_accessions
del_data=paste0(base,"Uniprot_Deleted_accessions_all.Rds")


## load new HGNC map to correct gene symbols to be HGNC compliant
# downloaded April 2022 https://github.com/waldronlab/HGNChelper/blob/master/data/hgnc.table.rda
load(file=paste0(base,"hgnc.table.rda"))
"https://github.com/d0minicO/unimachR/blob/main/hgnc.table.rda?raw=TRUE"
hgnc.table %<>%
  dplyr::select(1,2)


## CTLH complex members
CTLH_symbols =
  c("MKLN1",
    "RANBP9",
    "RMND5A",
    "RMND5B",
    "MAEA",
    "GID8",
    "WDR26",
    "ARMC8",
    "GID4")

CTLH_ids =
  c(
    "Q9UL63",
    "Q96S59",
    "Q9H871",
    "Q96G75",
    "Q7L5Y9",
    "Q9NWU2",
    "Q9H7D7",
    "Q8IUR7",
    "Q8IVV7"
  )

CTLH = 
  tibble(
    CTLH_symbols,
    CTLH_ids
  )


## use custom unimachR function to map IDs to gene symbols
CTLH_out = unimachR(CTLH_ids,hgnc.table)



## the uniprot human proteome
proteome_file = paste0(base,"2021_05_14_Homo_sapiens_UP000005640.fasta")



  
#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")
dir.create(plot_dir,showWarnings = F)

tables = paste0(base,"tables/")
dir.create(tables,showWarnings = F)

###################
#### DATA LOAD ####
###################


## read the overall counts data
dat = 
  read_tsv(in_file) %>%
  dplyr::slice(-1)

## count the rev decoy peptides and contaminant proteins
dat %>%
  group_by(`Potential contaminant`,Reverse) %>%
  summarise(count=n())


## remove rev decoy peptides and contaminant proteins
dat %<>%
  filter(is.na(`Potential contaminant`)) %>%
  filter(is.na(Reverse))



## clean up table to just intensity cols
dat %<>%
  dplyr::select(`Protein IDs`, starts_with("LFQ")) %>%
  dplyr::rename(Uniprot = `Protein IDs`)



####################
#### ID TO GENE ####
####################

## now get vector of all uniprot IDs found
ids =
  dat %>%
  dplyr::select(Uniprot) %>%
  makeChar() %>%
  unique()


ids = unlist(str_split(ids, ";"))

## remove any remaining IDs with the "CON__" flag from Maxquant
## or the "REV__" flag
ids= ids[-grep("CON_",ids)]
ids= ids[-grep("REV_",ids)]

## use custom unimachR function to map IDs to gene symbols with biomart
mapping_out = unimachR(ids,hgnc.table)

#saveRDS(mapping_out,file="IDS_matched_genes.Rds")

#mapping_out = readRDS(file="IDS_matched_genes.Rds")

# get the output as a df
mapping_df =
  as_tibble(mapping_out[[1]]) %>%
  dplyr::rename(ID=uniprot_id)
  


########################
#### DATA WRANGLING ####
########################

## fix the column names
colnames(dat) = gsub("LFQ intensity ","",colnames(dat))


## get into long format for later use
## parsing the sample names to get meaningful info from them and standardise
long = 
  dat %>%
  gather(Sample,Intensity,BIOID2_n1:BIOID2GID4_TREATED_MG_n4) %>%
  mutate(Intensity = as.numeric(Intensity)) %>%
  mutate(Log2Int = log2(Intensity)) %>%
  mutate(Log2Int = if_else(is.infinite(Log2Int),0,Log2Int)) %>%
  mutate(Sample2=Sample) %>%
  separate(Sample,into=c("sample","rep"),sep="_n") %>%
  separate(sample, into=c("Sample","treat1","treat2"),sep="_") %>%
  mutate(treat1 = if_else(
    is.na(treat1),
    "DMSO",
    treat1
  )) %>%
  mutate(treat3 = if_else(
    treat1 == "TREATED",
    treat1,
    treat2
  )) %>%
  mutate(treat1 = if_else(
    !is.na(treat3),
    treat2,
    treat1
  )) %>%
  mutate(treat3 = if_else(
    is.na(treat3),
    "DMSO",
    "PFI7"
  )) %>%
  mutate(treat1 = if_else(
    treat3 == "PFI7" & is.na(treat1),
    "DMSO",
    treat1
  )) %>%
  mutate(treat2 = treat3) %>%
  dplyr::select(-treat3)


##  tidy the data so each uniprot gets one row
long %<>%
  make_unique(name="Uniprot",ids="Uniprot") %>%
  dplyr::select(ID,Sample,treat1,treat2,rep,Intensity, Log2Int,Sample2)




###########################
#### DETECTED PROTEINS ####
###########################

sums =
  long %>%
  filter(Intensity>0) %>%
  mutate(treat1 = gsub("MG","MG132",treat1)) %>%
  mutate(group = paste(Sample,treat1,treat2,sep="_")) %>%
  group_by(group,treat1,treat2,rep) %>%
  summarise(count=n())


# total 6068
length(unique(long$ID))

# mean 5427
mean(sums$count)

# median 5426
median(sums$count)

# sd 63
sd(sums$count)

# min 5331
min(sums$count)

## plot sample numbers
ggplot(sums, aes(group,count,colour=treat2))+
  stat_summary(fun = mean, geom = "bar",position=position_dodge(.6),width=.6,fill="white",size=.1) +
  geom_point(size=.1,position=position_dodge2(.6),shape=21,alpha=.7)+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "uperrorbar",
               position=position_dodge2(.6),width=.2,size=.1,colour="black") +
  scale_y_continuous(limits=c(0,6000))+
  scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
  ylab("Detected proteins count")+
  xlab(NULL)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill=NA,size=.1),
    axis.ticks=element_line(size=.1),
    axis.ticks.length=unit(.3,"mm"),
    panel.background = element_blank(),
    legend.key.size = unit(.5,"cm"),
    axis.text.x = element_text(size=8,angle=45,hjust=1,vjust=1),   
    plot.title = element_text(hjust = 0.5,vjust=-1.5),
    plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
    panel.border=element_rect(size=.1),
    panel.grid=element_blank(),
    text=element_text(size=6))


ggsave(filename=paste0(plot_dir,"QC_detected_proteins.pdf"),
       width=1.7,
       height=1.8)





#########################
#### PREPARE FOR DEP ####
#########################

## adapted from DEP package documentation here https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("BIOID", colnames(dat)) # get LFQ column numbers (area / intensity)


## use long format table to specify experimental design
experimental_design =
  long %>%
  dplyr::select(Sample2,Sample,treat1,treat2,rep) %>%
  mutate(condition=paste(Sample,treat1,treat2,sep="_")) %>%
  mutate(condition = gsub("MG_","MG132_",condition)) %>%
  dplyr::rename(label=Sample2,replicate=rep,) %>%
  dplyr::select(label,condition,replicate) %>%
  distinct() %>%
  ungroup()




## make cols numeric
dat %<>% mutate(across(starts_with("BIOID"),makeNum))

## run make_unique to get each protein on its own row
dat %<>% make_unique(names = "Uniprot",ids="Uniprot")
  
## make the summarized experiment object
data_se <- make_se(dat, LFQ_columns, experimental_design)






#################################
#### FILTERING / NORMALIZING ####
#################################

## following suggestion Goeminne et al 2020, default DEP filtering does not give best results
# https://doi.org/10.1021/acs.analchem.9b04375

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 4 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 2)

# VST normalize the data
data_filt <- normalize_vsn(data_filt)




####################
#### IMPUTATION ####
####################

# set seed to have imputation be reproducible
set.seed(123)
## performing mixed imputation
proteins_MNAR <- get_df_long(data_filt) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(data_filt) %in% proteins_MNAR

# Perform a mixed imputation as done in DEP package documentation
set.seed(123)
mixed_imputation <- impute(
  data_filt, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "MinProb") # imputation function for MNAR


## compare imputation options with a plot
plot_imputation(data_filt,mixed_imputation)

ggsave(paste0(plot_dir,"QC_imputations.pdf"),
       width=6,
       height=6)



#############################
#### MAIN DEP STATS TEST ####
#############################


## proceed with DEP stats test on mixed imp data
# Differential enrichment analysis  based on linear models and empirical Bayes statistics
data_diff_manual <- test_diff(mixed_imputation, type = "all", 
                              test = comp)


# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(1.5))

# Generate a results table and export it
data_results =
  as_tibble(get_results(dep))

sig_prots =
  data_results %>%
  filter(significant) %>%
  pull(ID) %>%
  unique()


## adjusted p values
pvals =
  data_results %>%
  dplyr::select(ID,contains("p.adj"))


keycol <- "condition"
valuecol <- "p.adj"
gathercols <- colnames(dplyr::select(pvals,-ID))

pvals %<>%
  gather_(keycol, valuecol, gathercols) %>%
  mutate(condition=gsub("_p.adj","",condition))


### log2 fold changes
fcs =
  data_results %>%
  dplyr::select(ID,contains("ratio"))


keycol <- "condition"
valuecol <- "ratio"
gathercols <- colnames(dplyr::select(fcs,-ID))

fcs %<>%
  gather_(keycol, valuecol, gathercols) %>%
  mutate(condition=gsub("_ratio","",condition))


# combine them
res = full_join(pvals,fcs,by=c("ID","condition"))


## filter on just the significant genes
sig_res =
  res %>%
  filter(p.adj<=.05 & abs(ratio)>log2(1.5))


## number of sig proteins (should be 427)
length(unique(sig_res$ID))



## summary of siginificant genes in each comparison
sum_stats =
  sig_res %>%
  group_by(condition) %>%
  summarise(count = n())

write.table(sum_stats,
            file=paste0(tables,"DEP_mixedImp_summary_stats.txt"),
            sep="\t",
            col.names=T,
            row.names = F,
            quote=F)





####################
#### CLUSTERING ####
####################

## perform heirarhcical clustering of imputed normalized data to inspect differences between conditions

## get imputed data as numeric matrix
mat = assays(mixed_imputation)[[1]]

## subset on just significant genes
## save the rownames for later to extract which proteins are in each cluster
row_names =
  mat %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  filter(ID %in% sig_prots) %>%
  pull(ID)


## now just get the numeric matrix
mat %<>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  filter(ID %in% sig_prots) %>%
  dplyr::select(-ID) %>%
  mutate(across(.cols=everything(),makeNum))





## give better colnames based on simplified naming defined above
new_colnames =
  c(paste0(unique(experimental_design$condition)[1],"_",1:4),
    paste0(unique(experimental_design$condition)[2],"_",1:4),
    paste0(unique(experimental_design$condition)[3],"_",1:4),
    paste0(unique(experimental_design$condition)[4],"_",1:4),
    paste0(unique(experimental_design$condition)[5],"_",1:4),
    paste0(unique(experimental_design$condition)[6],"_",1:4)
  )

# just line up to check
tibble(old_cols=colnames(mat),
       new_cols=new_colnames)

colnames(mat) = new_colnames


#create data frame for column annotations on heatmap
anno =
  data.frame(sample=as.character(colnames(mat))) %>%
  mutate(sample2 = sample) %>%
  column_to_rownames("sample2") %>%
  separate(sample, into=c("Protein","MG132","PFI7","Rep")) %>%
  dplyr::select(-Rep)



## determine number of clusters for rows based on an elbow plot
## using custom function 
elbow_plot(mat,plot_dir,"protein_clusters")



pdf(paste0(plot_dir,"QC_DEPs_heatmap.pdf"),5,6)
set.seed(123)
hmap.res = pheatmap(mat,scale="row",
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cutree_cols=2,
         cutree_rows=4,
         annotation_col = anno,
         show_rownames=F,
         clustering_method="complete")
print(hmap.res)
dev.off()

closeAllConnections()



#################################################
#### PLOT AVG PROTEIN LEVELS ACROSS CLUSTERS ####
#################################################

# use results of the heatmap above to link the proteins to each cluster
res.clust = 
  cbind(mat,cluster = cutree(hmap.res$tree_row, k = 4)) %>%
  as_tibble() %>%
  mutate(ID = row_names) %>%
  dplyr::select(ID,cluster)


## find out how the cluster numbers here relate to top to bottom on the heatmap
## use the number of proteins in each cluster as a proxy
res.clust %>%
  group_by(cluster) %>%
  dplyr::count()

# cluster mapper df
clust_map =
  tibble(heatmap_clust = c(1,2,3,4),
         cluster = c(2,4,3,1))


## join to change the old clusters to the new heatmap ones
res.clust %<>%
  left_join(clust_map,by="cluster") %>%
  dplyr::select(-cluster) %>%
  dplyr::rename(cluster=heatmap_clust)

## count the proteins in clusters again
## to help identify them!
res.clust %>%
  group_by(cluster) %>%
  dplyr::count()


## center each protein expression by subtracting median (median-centered Log2 intensities)
## have to transpose as scale/sweep operates on columns by default
new_mat = t(mat)

# get protein medians
med.att <- apply(new_mat, 2, median)

# perform centering, join to protein IDs
cent_dat = 
  sweep(data.matrix(new_mat), 2, med.att) %>%
  t() %>%
  as_tibble() %>%
  mutate(ID = row_names) %>%
  left_join(res.clust, by="ID")

## get into long format
keycol <- "condition"
valuecol <- "centd"
gathercols <- colnames(dplyr::select(cent_dat,-ID,-cluster))

cent_dat %<>%
  gather_(keycol, valuecol, gathercols) %>%
  separate(condition, into=c("Sample","treat1","treat2","rep"),sep="_") %>%
  mutate(group = paste(Sample,treat1,treat2,sep="_"))




# check all IDs to check how close to zero the protein medians are
chk = 
  cent_dat %>%
  group_by(ID) %>%
  summarise(median_id = median(centd)) %>%
  arrange(median_id)

## tiny error, fine to proceed
min(chk$median_id)
max(chk$median_id)


## now plot median centered expression across groups by cluster
ggplot(cent_dat, aes(group,centd,colour=treat2))+
  geom_hline(yintercept=0,linetype="dashed",size=.1)+
  geom_boxplot(outlier.shape = NA,size=.1,fill="white")+
  scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
  facet_wrap(~cluster,ncol=4)+
  coord_cartesian(ylim=c(-2.8,2.8))+
  theme_bw()+
  theme(
    strip.background = element_rect(fill=NA,size=.1),
    strip.text = element_text(margin = margin(.5,.5,.5,.5, "mm")),
    axis.ticks=element_line(size=.1),
    axis.ticks.length=unit(.3,"mm"),
    panel.background = element_blank(),
    legend.key.size = unit(.5,"cm"),
    axis.text.x = element_text(size=8,angle=45,hjust=1,vjust=1),   
    plot.title = element_text(hjust = 0.5,vjust=-1.5),
    plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
    panel.border=element_rect(size=.1),
    panel.grid=element_blank(),
    text=element_text(size=6))


ggsave(filename=paste0(plot_dir,"Clusters_expression.pdf"),
       width=3,
       height=2)



#############
#### PCA ####
#############


## do the PCA
pc = prcomp(t(mat), scale. = F, center=T)

## set up values and sample labels to plot using raw ggplot
pc_vals =
  pc$x %>%
  data.frame() %>%
  rownames_to_column("Sample") %>%
  separate(Sample,into=c("protein","treat1","treat2","rep"))



## get the variance explained for each component
var_explained <- pc$sdev^2/sum(pc$sdev^2)

# plot PC1 and PC2
pc_vals %>%
  ggplot(aes(x=PC1,y=PC2))+
  geom_point(data=subset(pc_vals, protein == "BIOID2GID4"),colour = "black", size = 2,aes(shape=treat1)) +
  geom_point(size=1,aes(shape=treat1,colour=treat2))+
  scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
  ggtitle("PCA plot of proteomics","just on DEPs")+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border = element_rect(size=.1),
        axis.ticks = element_line(size=.1),
        text=element_text(size=6),
        legend.key.size = unit(.5,"cm"))


ggsave(paste0(plot_dir,"QC_PCA.pdf"),
       width=2.2,
       height=1.6)




## now PCA just on BioID2GID4 samples

## get imputed data as numeric matrix
mat2 = assays(mixed_imputation)[[1]]

## subset on just significant genes
## remove ID and make all columns numeric
mat2 %<>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  filter(ID %in% sig_prots) %>%
  dplyr::select(contains("GID4")) %>%
  mutate(across(.cols=everything(),makeNum))



## do the PCA
pc = prcomp(t(mat2), scale. = F, center=T)

## set up values and sample labels to plot using raw ggplot
pc_vals =
  pc$x %>%
  data.frame() %>%
  rownames_to_column("Sample") %>%
  separate(Sample,into=c("protein","treat1","treat2","rep"))



## get the variance explained for each component
var_explained <- pc$sdev^2/sum(pc$sdev^2)

# plot PC1 and PC2
pc_vals %>%
  ggplot(aes(x=PC1,y=PC2))+
  geom_point(size=1,aes(shape=treat1,colour=treat2))+
  scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
  ggtitle("PCA plot of proteomics","just on DEPs")+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border = element_rect(size=.1),
        axis.ticks = element_line(size=.1),
        text=element_text(size=6),
        legend.key.size = unit(.5,"cm"))


ggsave(paste0(plot_dir,"QC_PCA_just_BIOID2GID4.pdf"),
       width=2.2,
       height=1.6)



###############################
#### EXPLORE PROTEIN LISTS ####
###############################


## combine the sig proteins with hgnc_symbols
mapped =
  res %>%
  left_join(mapping_df, by="ID")


## have some that are not unique, now they are duplicated over rows which is okay
test =
  mapped %>%
  group_by(ID,condition) %>%
  summarise(count=n()) %>%
  arrange(desc(count))


## 3 IDs not mapped to hgncy symbols (checked and they were all deleted April 2021)
mapped %>% filter(is.na(hgnc_symbol)) %>% pull(ID) %>% unique()
#" P0DN76" "P0DN79" "Q6ZSR9"

## remove these proteins
mapped %<>%
  filter(!is.na(hgnc_symbol))


## add better column info for the comparisons
mapped %<>%
  separate(condition, into=c("cond1","cond2"),sep="_vs_") %>%
  separate(cond1, into=c("cells_1","treat_MG_1","treat_PFI_1")) %>%
  separate(cond2, into=c("cells_2","treat_MG_2","treat_PFI_2")) %>%
  mutate(abs_ratio = abs(ratio))

## save the entire table
write.table(mapped,
            file=paste0(tables,"Proteome_stats_ALL_Proteins.txt"),
            row.names=F,
            col.names = T,
            quote=F,
            sep="\t")


# save just the significant proteins table, and get absolute ratio for filtering big changes
sig_prots_dat = 
  mapped %>%
  filter(p.adj<=.05 & abs_ratio>log2(1.5)) %>%
  arrange(desc(abs_ratio))


## combine with which cluster each protein is in
sig_prots_dat %<>%
  left_join(res.clust,by="ID")



## label the genes that were significant in different biologically meaningful comparisons
sig_prots_dat %<>%
  mutate(comparison = case_when(
    cells_1==cells_2 &
      treat_MG_1 == "DMSO" &
      treat_MG_2 == "DMSO" ~ "PFI7_dependent_noMG",
    cells_1==cells_2 &
      treat_MG_1 == "MG132" &
      treat_MG_2 == "MG132" ~ "PFI7_dependent_MG",
    treat_PFI_1==treat_PFI_2 &
      treat_MG_1 == "DMSO" &
      treat_MG_2 == "DMSO" ~ "GID4_dependent_noMG",
    treat_PFI_1==treat_PFI_2 &
      treat_MG_1 == "MG132" &
      treat_MG_2 == "MG132" ~ "GID4_dependent_MG",
    treat_PFI_1==treat_PFI_2 &
      cells_1==cells_2 ~ "MG132_dependent"
  ))


## export table of significant proteins for supplementary tables
write.table(sig_prots_dat,
            file=paste0(tables,"Proteome_stats_SIG_Proteins_P_FC.txt"),
            row.names=F,
            col.names = T,
            quote=F,
            sep="\t")





#######################
#### DEPS PLOTTING ####
#######################
# DEP plotter

## get imputed data as numeric matrix
imp_dat = assays(mixed_imputation)[[1]]

## give better colnames based on simplified naming defined above
new_colnames =
  c(paste0(unique(experimental_design$condition)[1],"_",1:4),
    paste0(unique(experimental_design$condition)[2],"_",1:4),
    paste0(unique(experimental_design$condition)[3],"_",1:4),
    paste0(unique(experimental_design$condition)[4],"_",1:4),
    paste0(unique(experimental_design$condition)[5],"_",1:4),
    paste0(unique(experimental_design$condition)[6],"_",1:4)
  )

# just line up to check
tibble(old_cols=colnames(imp_dat),
       new_cols=new_colnames)

colnames(imp_dat) = new_colnames


row_names = rownames(imp_dat)

## now do median centering for each protein
## center each protein expression by subtracting median (median-centered Log2 intensities)
## have to transpose as scale/sweep operates on columns by default
new_mat = t(imp_dat)

# get protein medians
med.att <- apply(new_mat, 2, median)

# perform centering, join to protein IDs
cent_dat = 
  sweep(data.matrix(new_mat), 2, med.att) %>%
  t() %>%
  as_tibble() %>%
  mutate(ID = row_names)

## get into long format
keycol <- "condition"
valuecol <- "centd"
gathercols <- colnames(dplyr::select(cent_dat,-ID))

cent_dat %<>%
  gather_(keycol, valuecol, gathercols) %>%
  separate(condition, into=c("Sample","treat1","treat2","rep"),sep="_",remove = F) %>%
  mutate(group = paste(Sample,treat1,treat2,sep="_"))


# check all IDs to check how close to zero the protein medians are
chk = 
  cent_dat %>%
  group_by(ID) %>%
  summarise(median_id = median(centd)) %>%
  arrange(median_id)

## tiny errors, fine to proceed
min(chk$median_id)
max(chk$median_id)


## join protein IDs to hgnc symbols
cent_dat %<>%
  left_join(mapping_df, by="ID")



## plot CTLH complex member expression
cent_dat %>%
  filter(ID %in% CTLH_ids) %>%
  ggplot(aes(group,centd))+
  #stat_summary(fun = mean, geom = "line",position=position_dodge(.6),width=.6) +
  geom_boxplot(size=.1,outlier.shape = NA,aes(colour=Sample))+
  geom_point(shape=21,size=.5,stroke=.3,position=position_dodge(.6),aes(colour=treat2),fill="white")+
  scale_fill_manual(values=c("#377eb8ff","#e4211cff"))+
  scale_colour_manual(values=c("#ccccccff","#000000ff","#377eb8ff","#e4211cff"))+
  facet_wrap(~hgnc_symbol,ncol=4)+
  scale_x_discrete(labels=c('-', '+', '-',"-","+","+"))+
  xlab("MG132")+
  ylab("Median centered Log2(Intensity)")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill=NA,size=.1),
    strip.text = element_text(margin = margin(.5,.5,.5,.5, "mm")),
    axis.ticks=element_line(size=.1),
    axis.ticks.length = unit(.3,"mm"),
    panel.background = element_blank(),
    legend.key.size = unit(.5,"cm"),
    plot.title = element_text(hjust = 0.5,vjust=-1.5),
    plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
    panel.border=element_rect(size=.1),
    panel.grid=element_blank(),
    text=element_text(size=6))


ggsave(paste0(plot_dir,"CTLH_expression.pdf"),
       width=6,
       height=2.5)






## add which cluster each gene is in
## NA for the non sig genes
cent_dat %<>%
  left_join(res.clust,by="ID")




# plot a specific list of proteins
prot_list_2 = c(
  "DDX21",
  "DDX50",
  "IFI16",
  "DDX39A",
  "EIF4A2",
  "LMNB2",
  "DHX40",
  "KIN",
  "DICER1",
  "CHAF1A"
)



## get data on proteins to plot
temp =
  cent_dat %>%
  filter(hgnc_symbol %in% prot_list_2) %>%
  distinct()

## add degron labels to the plot
## first join to table
temp %<>%
  left_join(proteome,by="ID") %>%
  dplyr::select(-seq)



## plot median centered log 2
p1 =
  temp %>%
  ggplot(aes(group,centd))+
  geom_boxplot(size=.1,outlier.shape = NA,aes(colour=Sample))+
  geom_point(shape=1,size=.07,stroke=.5,position=position_dodge(.6),aes(colour=treat2))+
  scale_colour_manual(values=c("#ccccccff","#000000ff","#377eb8ff","#e4211cff"))+
  ggtitle("selected proteins")+
  facet_wrap(hgnc_symbol.x~first_ten,ncol=length(prot_list_2),scales="free_y")+
  scale_x_discrete(labels=c('-', '+', '-',"-","+","+"))+
  xlab("MG132")+
  ylab("Median centered Log2(Intensity)")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill=NA,size=.1),
    strip.text = element_text(margin = margin(.1,.1,.1,.1, "mm"),size=4),
    axis.ticks=element_line(size=.1),
    axis.ticks.length = unit(.3,"mm"),
    panel.background = element_blank(),
    legend.key.size = unit(2,"mm"),
    legend.title=element_blank(),
    plot.title = element_text(hjust = 0.5,vjust=-1.5),
    plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
    panel.border=element_rect(size=.1),
    panel.grid=element_blank(),
    text=element_text(size=5))

ggsave(p1,
       filename=paste0(plot_dir,"Selected_prots_median_centered",".pdf"),
       device="pdf",
       width=7.4,
       height=1.1)




#######################
#### VOLCANO PLOTS ####
#######################


## prepare the df for volcano analysis
volc =
  mapped %>%
  mutate(minuslogp = -log10(p.adj)) %>%
  arrange(desc(minuslogp))


## set up a column for significance
volc %<>%
  mutate(sig = if_else(
    (p.adj<=.05 & abs_ratio > log2(1.5)),
    "sig",
    "notsig"
  ))

## set up a grouping column again
volc %<>%
  mutate(condition=paste(cells_1,treat_MG_1,treat_PFI_1,"vs", cells_2,treat_MG_2,treat_PFI_2,sep="_"))


## swap the direction of the ratio to make more intuitive sense
## points will be to the right if they were
## higher in the sample listed on the right of the comparison
volc %<>%
  mutate(ratio = -ratio)

#comps = unique(volc$condition)
pfi_comp = "BIOID2GID4_DMSO_DMSO_vs_BIOID2GID4_DMSO_PFI7"#comps[1]
mg_comp = "BIOID2_DMSO_DMSO_vs_BIOID2_MG132_DMSO"
gid4_comp = "BIOID2_DMSO_DMSO_vs_BIOID2GID4_DMSO_DMSO"

comps = list(gid4_comp,
             pfi_comp,
          mg_comp
          )


this_comp = pfi_comp



## chose selected proteins to label on the volcano plots

to_label_pfi =
  c(
    "DHX40",
  "DICER1",
  "KIN",
  "CHAF1A",
  "DDX21",
  "DDX50")


to_label_gid =
c(
"IFI16",
"GID4",
"DDX39A",
"EIF4A2",
"LMNB2"
  )




# set up universal features of the plots
boxpad = .5
pointpad = .5
minlength = .01
legpos = "bottom"



unique(volc$condition)


comps2 = c(
  "BIOID2_DMSO_DMSO_vs_BIOID2GID4_DMSO_DMSO",
  "BIOID2_MG132_DMSO_vs_BIOID2GID4_MG132_DMSO",
  "BIOID2GID4_DMSO_DMSO_vs_BIOID2GID4_DMSO_PFI7",
  "BIOID2GID4_MG132_DMSO_vs_BIOID2GID4_MG132_PFI7"
)

gid4_comps = comps2[1:2]

pfi_comps = comps2[3:4]



plot_list = list()
for(i in 1:length(comps2)){
  
  
  this_comp = comps2[[i]]
  
  temp =
    volc %>%
    filter(condition==!!this_comp)
  
  ## label desired proteins
  
  if(this_comp %in% gid4_comps){
    to_label = to_label_gid
  } else if (this_comp %in% pfi_comps){
    to_label = to_label_pfi
    
  }
  
  to_label = 
    temp %>%
    filter(hgnc_symbol %in% to_label) %>%
    pull(ID)
  
  ## set up label column to use in geom_text_repel
  temp %<>%
    mutate(label = if_else(
      ID %in% to_label,
      hgnc_symbol,
      "nolabel"
    ))
  
  
  
  ## only plot a y axis for the first plot
  if(i==1){
    ylab_val = "-log10(p.adjusted)"
  } else {
    ylab_val = NULL
  }
  
  
  ## create dummy variable to get legend to look nice
  #temp %<>%
  #  mutate(allhits = "All proteins")
  
  
  ## get sig prots labels to annotate to upper RH corner of plot
  prot_nums =
    temp %>%
    filter(sig=="sig") %>%
    dplyr::select(ratio) %>%
    mutate(pos_or_neg = if_else(
      ratio>0,
      "down",
      "up"
    )) %>%
    group_by(pos_or_neg) %>%
    dplyr::count()
  
  
  ## find x and y values to position proteins labelling
  label_y = max(temp$minuslogp)*.9
  label_x_left = min(temp$ratio)*.9
  label_x_right = max(temp$ratio)*.9
  
  # custom volcano plot
  plot_list[[i]]=
    ggplot(temp,aes(ratio,minuslogp, label = label))+
    geom_hline(yintercept = -log10(.05),linetype="dashed",size=.1)+
    geom_vline(xintercept = log2(1.5),linetype="dashed",size=.1)+
    geom_vline(xintercept = -log2(1.5),linetype="dashed",size=.1)+
    geom_point_rast(data=subset(temp, sig!="sig"),alpha=.4,raster.dpi=900,shape=".",colour="gray60")+ # rasterise as many points are slow
    geom_point_rast(data=subset(temp, sig=="sig"&label=="nolabel"),alpha=.8,raster.dpi=900,shape=".")+
    geom_point_rast(data=subset(temp, label!="nolabel"),alpha=.95,raster.dpi=900,shape=20,size=.1,colour="red")+
    # decreased proteins number label
    annotate(geom="text", x=label_x_left, y=label_y, label=prot_nums[1,2],color="black",size=1)+
    # increased proteins number label
    annotate(geom="text", x=label_x_right, y=label_y, label=prot_nums[2,2],color="black",size=1)+
    ggtitle(this_comp)+
    ylab(ylab_val)+
    xlab("Log2FC")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(size=.1),
      axis.ticks = element_line(size=.1),
      text=element_text(size=5),
      legend.key.size = unit(5,"mm"),
      title=element_text(size=2.5),
      legend.position = "bottom"
    ) +
    geom_text_repel(data          = subset(temp, label!="nolabel"),
                    colour="black",
                    size          = 1,
                    box.padding   = boxpad,
                    point.padding = pointpad,
                    max.overlaps = 40,
                    force         = 100,
                    segment.size  = 0.1,
                    min.segment.length = minlength,
                    segment.color = "grey50",
                    #direction     = "x")
    )
  
}



wrap_plots(plot_list) +
  plot_layout(guides="collect",ncol=4) +
  plot_annotation(title = "Protein changes dependent on GID4 or PFI-7",
                  theme = 
                    theme(
                      plot.title = element_text(size = 5),
                      legend.position = legpos
                    ))

ggsave(paste0(plot_dir,"Volcano_all.pdf"),
       width=3.8,
       height=1.3)





########################################
#### EXPORT ALL QUANTIFIED PROTEINS ####
########################################

## export the table for supplementary tables
to_export =
  dat %>%
  left_join(mapping_df) %>%
  dplyr::select(ID,hgnc_symbol,starts_with("BIOID")) %>%
  arrange(hgnc_symbol)

## fix col naming
colnames(to_export) = gsub("TREATED","PFI7",colnames(to_export))


write.table(to_export,
            file=paste0(tables,"Proteome_LFQ_Intensities.txt"),
            col.names = T,
            row.names = F,
            quote=F,
            sep="\t")



##########################################
#### ANNOTATING PROTEOME WITH DEGRONS ####
##########################################



proteome = readAAStringSet(proteome_file)

proteome %<>%
  aass2tib() %>% # custom function to make AA String set into a tibble
  separate(names, into=c(NA,"uniprot_id","rest"),sep="\\|") %>%
  separate(rest, into=c(NA,"temp"),sep="GN=") %>%
  separate(temp, into=c("hgnc_symbol"))

# 20380 proteins
nrow(proteome)

## get latest hgnc_symbols, important as it must match our symbols when joining the tables
proteome$hgnc_symbol =
  HGNChelper::checkGeneSymbols(proteome$hgnc_symbol,map=hgnc.table) %>%
  dplyr::select(Suggested.Symbol) %>%
  makeChar

# 595 failed to map
ids =
  proteome %>%
  filter(is.na(hgnc_symbol)) %>%
  pull(uniprot_id) %>%
  makeChar()

## find these using custom function unimachR
ids_mapped = unimachR(ids,hgnc.table)

## keep just the unique ones
## multi mappers here look weird and not informative protein / gene pairs
uniq_ids =
  ids_mapped[[1]] %>%
  group_by(hgnc_symbol) %>%
  dplyr::count() %>%
  filter(n==1) %>%
  pull(hgnc_symbol)

## get table to join with human proteome
to_join =
  ids_mapped[[1]] %>%
  as_tibble() %>%
  filter(hgnc_symbol %in% uniq_ids)

## join the recovered IDs to the table
## make a single column that combines the best HGNC_Symbol that can be found
## discard 199 proteins not mapped to symbol
proteome %<>%
  left_join(to_join,by="uniprot_id") %>%
  mutate(hgnc_symbol = case_when(
    !is.na(hgnc_symbol.x) ~ hgnc_symbol.x,
    !is.na(hgnc_symbol.x) ~ hgnc_symbol.y,
    is.na(hgnc_symbol.x) & is.na(hgnc_symbol.y) ~ "LOST"
  )) %>%
  filter(hgnc_symbol!="LOST")


## tidy up
proteome %<>%
  dplyr::select(-hgnc_symbol.x,-hgnc_symbol.y)


## now mark each protein's degron based on defined criteria Dong 2018 (perfect) & Dong 2020 (flexible)

## need a column for first 10 AAs
## this is for could be degron (if trimmed)
proteome %<>%
  mutate(first_ten = substr(seq,1,10))


## choose the difference degron options
## PERFECT DONG 2018
## SECOND POSITION TOLERATES PA or PS or PT or PH
perfect = c("PG", "PS", "PH", "PT", "PA")

# FLEXIBLE DONG 2020
# set up the different options for first and second position of flexible degron
flexible_first = c("P","I","L","F","V")
flexible_second = c("S","T","G","V","A")

# use expand.grid to combine in all possible permutations
## now paste together to get vector of all possible flexible degrons,
## remember to exclude the perfect degron PG
flexible_all =
  expand.grid(flexible_first,
            flexible_second) %>%
  as_tibble() %>%
  mutate(comb = paste0(Var1,Var2)) %>%
  filter(comb!=perfect)  %>%
  pull(comb)
  

## now mark the proteome based on these criteria
proteome %<>%
  mutate(Degron = case_when(
    grepl(paste0("^M",perfect,collapse="|"),seq) ~ "Perfect",
    grepl(paste0(perfect,collapse="|"), first_ten) ~ "Perfect:if_trimmed",
    grepl(paste0("^M",flexible_all,collapse="|"), seq) ~ "Flexible",
    grepl(paste0(flexible_all,collapse="|"), first_ten) ~ "Flexible:if_trimmed"
  )) %>%
  replace(is.na(.),"None")


## check some genes to make sure it worked
proteome %>%
  filter(hgnc_symbol=="DDX21") %>%
  dplyr::select(-seq)

## better naming of ID
proteome %<>%
  dplyr::rename(ID=uniprot_id)


## save this proteome table
write.table(proteome,
            file=paste0(proteome_file,"_degron_marked.txt"),
            sep="\t",
            quote=F,
            row.names = F,
            col.names = T)

## now save as RDS
saveRDS(proteome,paste0(base,"2021_05_14_Homo_sapiens_UP000005640_degron_marked.Rds"))

proteome = readRDS(paste0(base,"2021_05_14_Homo_sapiens_UP000005640_degron_marked.Rds"))

