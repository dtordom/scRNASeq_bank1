
################################################################################
#' SingleCell Analysis
#' @Rversion R version 4.0.4
#'
################################################################################

##----------------------------------------------------------············· STEP 0
## Setting environment
#' @mainPath: main folder
#' @dataPath: cellranger data path
#' @metadataPath: metadata.csv full path
#' @resourcesPath: path to file contains R functions
#' @specie: specie
#' @assay: data type
#' @metadataColumns: key columns in Metadata. (dataset: name of samples)
#' @varstoPlot: vars to plots
#' @selected_cells: to remove specific cells (file.csv contains cell barcodes, 
#' i,e AAACCCAGTATCTCTT_blw6xa69_3suwslra)
#' @normFactor: To normalize cell expression counts
#' @pct_mito_range: range of allowed percentage of mitocondrial genes in cells
#' @pct_ribo_range: range of allowed percentage of ribosomal genes in cells
#' @keep_genes: genes that avoid filters, or NULL
#' @remove_non_coding: remove non-coding genes: TRUE / FALSE
#' @remove_gene_family: remove  different groups of genes (mt- : Mitochondrial 
#' genes, rpl-:Ribosomal genes (rps prl))
#' @min_gene_count: filter genes by counts
#' @min_gene_per_cell: filter cells by number of measures genes
#' @vars_to_regress: Variables to correct cell expression (mitocondrial genes,
#' cell cycle)

set.seed(12345678)

opt<-list(mainPath="/home/daniel/Desktop/WORK/scRNASeq/MorellMouse", 
          dataPath="/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/data",
          metadataPath ="/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/metadata.csv",
          resourcesPath ="/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/scripts/functions.R",
          specie="mmusculus",
          assay="RNA",
          metadataColumns= c("dataset","assay","genotype","replicate"),
          varstoPlot = c("assay","genotype"),
          selected_cells ="none",
          pct_mito_range = c(0,25),
          normFactor = 10000,
          pct_ribo_range = c(0,50),
          keep_genes =c("Ighd", "Ighm", "Ighg1", "Ighg2c", "Ighg2b", "Ighg3", "Igha", "Ighe"),
          remove_non_coding =T,
          remove_gene_family= c("mt","rpl","rps","hb","malat1"),
          min_gene_count = 5,
          min_gene_per_cell =200,
          vars_to_regress = c("percent.mito","S.Score", "G2M.Score"))

source(opt$resourcesPath)
invisible(lapply(c('Seurat','dplyr','rafalib','Matrix','parallel','biomaRt',
                   'optparse','utils','matrixStats','patchwork','scCATCH','ggplot2',
                   'SingleCellExperiment','scales','RColorBrewer','vegan','ineq',
                   'igraph','sva','scran','scater','batchelor','clustree','optparse',
                   'scales','fields','data.table','scDblFinder','harmony','ggsci',
                   'tidyr','tibble','reshape','stringr'),require,character.only = TRUE))

## Create Subfolders
if (!file.exists(paste0(opt$mainPath,"/Rds"))){
  dir.create(paste0(opt$mainPath,"/Rds"))} 
opt$rdsPath<-paste0(opt$mainPath,"/Rds")

if (!file.exists(paste0(opt$mainPath,"/QC"))){
  dir.create(paste0(opt$mainPath,"/QC"))} 
opt$QCpath<-paste0(opt$mainPath,"/QC")

if (!file.exists(paste0(opt$mainPath,"/Results"))){
  dir.create(paste0(opt$mainPath,"/Results"))} 
opt$results<-paste0(opt$mainPath,"/Results")

## biomart objects
mart = useMart("ensembl", dataset = paste0(casefold(opt$specie),"_gene_ensembl"),host="https://jul2019.archive.ensembl.org")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://jul2019.archive.ensembl.org") 


##----------------------------------------------------------············· STEP 1
## Load data

dataset_metadata <- as.data.frame(read.csv(opt$metadataPath))
#dataset_metadata<-dataset_metadata[,c("dataset","assay","library_id","genotype","replicate")]
rownames(dataset_metadata)<-dataset_metadata$dataset
#dataset_metadata<-dataset_metadata[dataset_metadata$genotype!="Sidt1",]

datasets <- dataset_metadata$dataset
#datasets <- sort(datasets[datasets %in% as.character(dataset_metadata$dataset)])

# paths: paste0(opt$dataPath,"/",datasets[i],"/outs/filtered_feature_bc_matrix")

if(length(datasets) > 1){
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, varlist = c("datasets","opt") )
  
  data <- parLapplyLB(cl, datasets, function(i){
    if( sum(grepl(".mtx", list.files( paste0(opt$dataPath,"/",i,"/outs/filtered_feature_bc_matrix") ))) >= 1 ){
      a <- Seurat::Read10X( paste0(opt$dataPath,"/",i,"/outs/filtered_feature_bc_matrix") )
    } else{print(paste0("Not .mtx file found for ",sep="",i))}
    colnames(a) <- paste0(sub("-.*","",colnames(a)),"_",as.character(i))
    return(a)
  })
  names(data) <- datasets
  
  all_genes <- unique(unlist(parLapplyLB(cl,data,function(x) return(rownames(x)))))
  clusterExport(cl, varlist = c("all_genes") )
  data <- parLapplyLB(cl, data, all_genes=all_genes,function(x,all_genes) {
    m <- Matrix::Matrix(0,nrow = length(all_genes), ncol = ncol(x),sparse = T,
                        dimnames = list(all_genes,colnames(x)))
    m[rownames(x),] <- x
    return(m)
  })
  print(as.data.frame(lapply(data,dim),row.names = c("genes","cells")))
  DATA <- do.call(cbind,data)
  
  DATA <- CreateSeuratObject(DATA,min.cells = 1,min.features = 1,assay = opt$assay)
  DATA$orig.ident <- setNames(sub("(.*?)_","",colnames(DATA)) , colnames(DATA) )
  rm(data)
  
} else {
  a <- Read10X( paste0(opt$dataPath,"/",i,"/outs/filtered_feature_bc_matrix") )
  colnames(a) <- paste0(colnames(a),"_",as.character(i))
  DATA <- CreateSeuratObject(a,min.cells = 1,min.features = 1)
}
cat("\nThe total dimensions of your merged raw dataset is: ",dim(DATA),"\n")
## 24536 genes 13598 cells

cat("\nAdding metadata ...\n")
for(i in opt$metadataColumns){
  DATA <- AddMetaData(object = DATA, 
                      metadata = setNames(as.character(dataset_metadata[match(as.character(DATA$orig.ident), 
                                                                              as.character(dataset_metadata[,1]) ),i]),rownames(DATA@meta.data)), col.name = i)
}

for(i in 1:length(datasets)){
  tmp<-DATA@meta.data[DATA@meta.data$orig.ident==datasets[i],"genotype"][1]
  print(paste0(datasets[i],": ",tmp))
}
print(dataset_metadata[,c("dataset","genotype")])

##--------------------------------- OPTIONAL
## QC of samples
seeExpression(dataObj=DATA,
              metadata=dataset_metadata,
              splitBy="orig.ident",
              gen=c("genotype","replicate"),
              expBase=0,
              genes=c("cd19"))

#RemoveSamples<-c("odqmndmf_j6skyqh9","y6u8cqm4_vn5l2woi")
#for(i in 1:length(RemoveSamples)){
#  sel<-!grepl(RemoveSamples[i],colnames(DATA))
#  DATA<-DATA[,sel]  
#}
##--------------------------------- OPTIONAL

## Select subset of cells (for subclustering of specific previously annotated cell lines)
if (casefold(opt$selected_cells) != "none"){
  cat("\nRemoving non-selected cells from the Seurat object\n")
  cells_use<-readRDS(opt$selected_cells)
  selected_cells <- colnames(DATA) %in% cells_use
  DATA <- subset(DATA, cells=colnames(DATA)[selected_cells])
}
stopCluster(cl);   print(dim(DATA))
rm(all_genes,datasets,i,cl,tmp)

infoSamples<-data.frame("id"=names(table(DATA$orig.ident)),
                        "readCells"=as.numeric(table(DATA$orig.ident)))

save.image(file = paste0(opt$rdsPath,"/raw_seurat_object.RData") )
## 2 Muestras eliminadas (porcentaje de B cells muy por debajo de lo esperado)

##----------------------------------------------------------············· STEP 2
## Quality control

## 2.1. QC by cells
## Diversity Indexes
cat("\nCalculating data diveristy indexes ...\n")
indexes <- t(apply(DATA@assays[[opt$assay]]@counts,2,function(x) {
  c(vegan::diversity(x,index = "simpson"),
    vegan::diversity(x,index = "invsimpson"),
    vegan::diversity(x,index = "shannon"),
    Gini(x)) }))
DATA$simp_index <- indexes[,1]
DATA$invsimp_index <- indexes[,2]
DATA$shan_index <- indexes[,3]
DATA$gini_index <- indexes[,4]

## Percentage of Gene families
cat("\nCalculating percentage of mitocondrial/ribosomal genes ...\n")
Gene.groups <- substring(rownames(x = DATA@assays[[opt$assay]]@counts),1,3)
seq_depth <- Matrix::colSums(DATA@assays[[opt$assay]]@counts)
temp <- rowsum(as.matrix(DATA@assays[[opt$assay]]@counts),Gene.groups)
perc <- sort(apply( t(temp) / seq_depth,2,median) ,decreasing = T)*100
tot <- sort(rowSums(temp)/sum(temp),decreasing = T)*100

#Compute the relative expression of each gene per cell
rel_expression <- Matrix::t( Matrix::t(DATA@assays[[opt$assay]]@counts) / 
                               (Matrix::colSums(DATA@assays[[opt$assay]]@counts)) ) * 100
#memory.limit(size = 1000000000)
most_expressed <- sort(apply(rel_expression,1,mean),T)[1:100] / ncol(DATA)

png(filename = paste0(opt$QCpath,"/Gene_familty proportions.png"),width = 1400*3,height = 4*1400,res = 300)
mypar(4,1,mar=c(5,5,2,1))
boxplot( as.matrix(Matrix::t(rel_expression[names(most_expressed),])),cex=.1,outline=T,las=2,main="% total count per cell",col=hue_pal()(100))
boxplot( (t(temp)/seq_depth) [,names(perc)[1:100]]*100,outline=T,las=2,main="% reads per cell",col=hue_pal()(100))
boxplot(t(temp)[,names(perc)[1:100]], outline=T,las=2,main="reads per cell",col=hue_pal()(100) )
barplot(tot[names(tot)[1:100]],las=2,xaxs="i",main="Total % reads (all cells)",col=hue_pal()(100))
invisible(dev.off())

parameters<-c("rpl","rps","hb[ab]","mito","hb")
for(i in parameters){
  cat(i,"\t")
  family.genes <- rownames(DATA@assays[[opt$assay]]@counts)[grep(pattern = paste0("^",ifelse(i=="mito","mt-",i)), x = casefold(rownames(DATA@assays[[opt$assay]]@counts)), value = F)]
  if(length(family.genes)>1){DATA <- PercentageFeatureSet(DATA,features = family.genes,assay = opt$assay,col.name = paste0("perc_",ifelse(i=="mt-","mito",i)) )}
}
rm("temp","perc","tot","Gene.groups","i","indexes")

## Calculating gene biotype percentages 
cat("\nCalculating gene biotype percentages ...\n")
annot <- getBM(c("external_gene_name","gene_biotype","transcript_biotype","chromosome_name"),mart = mart)
annot[,"chromosome_name"] <- paste0("Chr_",annot[,"chromosome_name"])
annot[ !grepl("^Chr_[123456789XYMT]",annot[,"chromosome_name"]) ,"chromosome_name"] <- "other"

for(z in c("gene_biotype","transcript_biotype","chromosome_name")){
  item <- annot[match(rownames(DATA@assays[[opt$assay]]@counts) , annot[,1]),z]
  item[is.na(item)] <- "unknown"
  
  png(filename = paste0(opt$QCpath,"/",z,"_proportions.png"),width = 1200*3,height = 1200,res = 300)
  mypar(1,3,mar=c(4,2,2,1))
  pie(sort(table(item),decreasing = T), clockwise = T,col = hue_pal()(length(unique(item))))
  title("before filtering")
  par(mar=c(10,2,2,1))
  barplot(sort(table(item),decreasing = T),las=2,xaxs="i",main="Total reads (all cells)",col=hue_pal()(100))
  
  temp <- rowsum(as.matrix(DATA@assays[[opt$assay]]@counts),group=item)
  o <- order(apply(temp,1,median),decreasing = T)
  boxplot( (t(temp)/Matrix::colSums(DATA@assays[[opt$assay]]@counts))[,o]*100,outline=F,las=2,main="% reads per cell",col=hue_pal()(100))
  invisible(dev.off())
  
  aaa <- setNames(as.data.frame(((t(temp)/Matrix::colSums(DATA@assays[[opt$assay]]@counts))[,o]*100)[,names(sort(table(item),decreasing = T))]),paste0("perc_",names(sort(table(item),decreasing = T))))
  DATA@meta.data <- DATA@meta.data[,!(colnames(DATA@meta.data) %in% colnames(aaa))]
  DATA@meta.data <- cbind(DATA@meta.data,aaa)
}

## Cellcycle scoring
cat("\nLog Normalizing counts for cell cycle scoring...\n")
DATA <- NormalizeData(object = DATA, scale.factor = opt$normFactor)

cat("\nPredicting cell cycle scores with Seurat ...\n")
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

if(casefold(opt$specie) != "hsapiens"){
  s.genes = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = s.genes , mart = mart, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol")[,1]
  g2m.genes = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = g2m.genes , mart = mart, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol")[,1]
}
s.genes<-s.genes[order(s.genes)]
g2m.genes<-g2m.genes[order(g2m.genes)]

DATA <- CellCycleScoring(object = DATA, s.features = s.genes, g2m.features = g2m.genes,set.ident = T)
DATA$G1.Score <- 1 - ( DATA$S.Score + DATA$G2M.Score )
DATA$CC.Diff <- DATA$S.Score - DATA$G2M.Score

for(i in opt$varstoPlot){
  feats <- colnames(DATA@meta.data) [ grepl("nFeature|nCount|_index|[.]Score",colnames(DATA@meta.data) ) ]
  feats <- c(feats,"perc_mito" ,"perc_rps","perc_rpl","perc_hb", "perc_protein_coding" ,"perc_lincRNA","perc_snRNA","perc_miRNA","perc_processed_pseudogene",
             "perc_unknown","perc_Chr_1","perc_Chr_X","perc_Chr_Y","perc_Chr_MT")
  feats <- feats[feats %in% colnames(DATA@meta.data)]
  
  png(filename = paste0(opt$QCpath,"/QC_",i,"_ALL.png"),width = 1600*(length(unique(DATA@meta.data[,i]))/2+1),height = 1200*ceiling(length(feats)/5),res = 300)
  print(VlnPlot(object = DATA, features  = feats, ncol = 5,group.by = i,pt.size = 0,assay = opt$assay))
  invisible(dev.off())
}

## Cell filtering identification 
cat("\nIdentify low quality cells ...\n")
NF <-  DATA@meta.data [ grepl("nFeature",colnames(DATA@meta.data)) ][,1]
NC <-  DATA@meta.data [ grepl("nCount",colnames(DATA@meta.data)) ][,1]

mito_range <- opt$pct_mito_range
ribo_range <- opt$pct_ribo_range

## Remove doublets
dta<-as.SingleCellExperiment(DATA)
dta.cnt <- logNormCounts(dta)
dec <- modelGeneVar(dta.cnt, block = dta.cnt$orig.ident)
hvgs = getTopHVGs(dec, n = 2000)
dta.cnt <- runPCA(dta.cnt, subset_row = hvgs)
dta.cnt <- runUMAP(dta.cnt, pca = 10)
dbl.dens <- computeDoubletDensity(dta.cnt,d=ncol(reducedDim(dta.cnt)))
summary(dbl.dens)
#dta.cnt$DoubletScore <- dbl.dens
#plotUMAP(dta.cnt, colour_by="DoubletScore",point_size=0.01)
dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call")
summary(dbl.calls)
DATA@meta.data$doublets<-dbl.calls

Ts <- data.frame(
  Doublets = DATA$doublets=="singlet",
  MitoT = between(DATA$perc_mito, mito_range[1], mito_range[2]),
  RpsT = between(DATA$perc_rps, ribo_range[1], ribo_range[2]),
  RplT = between(DATA$perc_rpl, ribo_range[1], ribo_range[2]),
  nUMIT = between(NF,quantile(NF,probs = c(0.005)),quantile(NF,probs = c(0.995))),
  nCountT = between(NC,quantile(NC,probs = c(0.005)),quantile(NC,probs = c(0.995))),
  GiniT = between(DATA$gini_index,0.8,1),
  SimpT = between(DATA$simp_index,0.8,1),
  protein_codingT = between(DATA$perc_protein_coding,50,100),
  row.names = rownames(DATA@meta.data) )
print(head(Ts,5))

cell_use <- rownames(Ts)[ rowSums(!Ts) == 0 ] ## Cell_use contains cell ids to keep

cat("\nDimentions of the raw.data objects BEFORE filtering ...\n")
print( dim(DATA@assays[[opt$assay]]@counts) ) # 27567 61777

#save.image(file = paste0(opt$rdsPath,"/PREfilt_seurat_object.RData"))
#load("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds/PREfilt_seurat_object.RData")

DATA<-DATA[,cell_use]

seeExpression(dataObj=DATA,
              metadata=dataset_metadata,
              splitBy="orig.ident",
              gen=c("genotype","assay"),
              expBase=0,
              genes=c("cd79a","cd79b"))

## Select the same number of cells per sample (to avoid problems with heterogeneity in mice and the number of cells per mice)
# if(opt$sampleCount){
#   metadata<-DATA@meta.data[cell_use,]
#   ncells<-min(table(metadata$orig.ident))
#   samples<-unique(metadata$orig.ident)
#   selectedCells<-NULL
#   
#   for(smpl in samples){
#     tmp<-metadata[metadata$orig.ident==smpl,]
#     tmp<-rownames(tmp)[sample(1:nrow(tmp),ncells,replace = F)]
#     selectedCells<-c(selectedCells,tmp)
#   }
#   
#   cell_use<-selectedCells
#   rm(metadata,ncells,samples,selectedCells,tmp,smpl)
# }

rm(aaa,human,mart,rel_expression,Ts,family.genes,parameters,o,item,temp,hvgs,
   mito_range,ribo_range,most_expressed,NC,NF,z,seq_depth,dec,dta,dta.cnt,
   dbl.calls,dbl.dens,i)


## 2.2. QC by genes
## To keep specific genes 
if( is.null(opt$keep_genes) ==FALSE ){
  genes_keep <- trimws(unlist(strsplit(casefold(opt$keep_genes), ',')))
  genes_notfound <- setdiff(genes_keep, casefold(rownames(DATA@assays[[opt$assay]]@counts)))
  genes_keep <- rownames(DATA@assays[[opt$assay]]@counts)[casefold(rownames(DATA@assays[[opt$assay]]@counts)) %in% genes_keep]
  if(length(genes_notfound) > 0){
    cat("\nWARNING: The following requested genes were NOT FOUND in the data:\n")
    cat(genes_notfound, "\n")
  }
} else {
  genes_keep <- NULL ## genes_keep
}

## Select only protein-coding genes
if(opt$remove_non_coding){
  sel <- annot[match(rownames(DATA@assays[[opt$assay]]@counts) , annot[,1]),2] == "protein_coding"
  genes_use <- rownames(DATA@assays[[opt$assay]]@counts)[sel]
  genes_use <- union(as.character(na.omit(genes_use)), genes_keep)
  DATA@assays[[opt$assay]]@counts <- DATA@assays[[opt$assay]]@counts[genes_use,]
}

## Removing some family genes from the data (mt)
if(!is.null(opt$remove_gene_family)){
  for(gn in opt$remove_gene_family){
    genes_use<-rownames(DATA@assays[[opt$assay]]@counts)[!grepl(casefold(paste0("^",gn)), 
                                                                casefold(rownames(DATA@assays[[opt$assay]]@counts)))]
    genes_use <- union(as.character(na.omit(genes_use)), genes_keep)
    DATA@assays[[opt$assay]]@counts <- DATA@assays[[opt$assay]]@counts[genes_use,]
  }
}

## Filtering cells and re-normalizing filtered data
DATA <- CreateSeuratObject(counts = DATA@assays[[opt$assay]]@counts , assay = opt$assay, meta.data = DATA@meta.data, min.cells = as.numeric(opt$min_gene_count),min.features = as.numeric(opt$min_gene_per_cell))

cat("\nDimentions of the raw.data objects AFTER filtering ...\n")
print( dim(DATA@assays[[opt$assay]]@counts) )
# 16852 46923

## Normalize counts
cat("\nNormalizing counts ...\n") ## SCTransform
DATA <- NormalizeData(object = DATA,scale.factor = opt$normFactor, normalization.method = "LogNormalize")

for(i in opt$varstoPlot){
  png(filename = paste0(opt$QCpath,"/QC_",i,"_FILTERED.png"),width = 1600*(length(unique(DATA@meta.data[,i]))/2+1),height = 1200*ceiling(length(feats)/5),res = 300)
  print(VlnPlot(object = DATA, features  = feats, ncol = 5,group.by = i,pt.size = 0,assay = opt$assay))
  invisible(dev.off())}

infoSamples$filterCells<-as.numeric(table(DATA$orig.ident))

cat("\nSaving filtered Seurat object ...\n")
rm(annot,cell_use,feats,genes_use,i,genes_keep,genes_notfound,gn,sel)
save.image(file = paste0(opt$rdsPath,"/filt_seurat_object.RData") )
# load("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds/filt_seurat_object.RData")

##----------------------------------------------------------············· STEP 3
## Integration

DATA <- FindVariableFeatures(DATA, selection.method = "vst", nfeatures = 2500)
DATA <- ScaleData(DATA, features = rownames(DATA))
# Cell cycle scores and plots.
DATA <- CellCycleScoring(object = DATA, s.features = s.genes, g2m.features = g2m.genes,set.ident = T)
# PCA previous to cell cycle scoring.
DATA <- RunPCA(DATA, features = c(s.genes, g2m.genes),npcs = 50)

## Plot before cellcycle correction
png(filename = paste0(opt$QCpath,"/cellCycle_BEFORE.png"),width = 2000,height = 1600,res = 300)
plot(DimPlot(DATA, reduction = "pca"))
invisible(dev.off())

## Scale with regress
DATA <- ScaleData(object = DATA, vars.to.regress = opt$vars_to_regress)
DATA <- RunPCA(DATA, features = c(s.genes, g2m.genes))

png(filename = paste0(opt$QCpath,"/cellCycle_AFTER.png"),width = 2000,height = 1600,res = 300)
p1<-DimPlot(DATA, reduction = "pca")
plot(p1)
invisible(dev.off())

## Set Ident (experiment / genotype)
DATA <- SetIdent(DATA, value = DATA@meta.data$genotype) 
DATA <- RunPCA(DATA, features = VariableFeatures(object = DATA))
DATA.i<-RunHarmony(object = DATA,group.by.vars = "orig.ident")
# matrix: DATA.i@assays$RNA@scale.data

## PCA and Visualize Dimensional Reduction genes.
DATA.i <- RunPCA(DATA.i, features = rownames(DATA.i), ncps = 100, verbose = FALSE)
png(filename = paste0(opt$results,"/1_viz_dim_loadings.png"),width = 2000,height = 1600,res = 300)
plot(VizDimLoadings(DATA.i, dims = 1:2, reduction = "pca") + theme(legend.position="bottom"))
invisible(dev.off())

## PCA projection and integration visualization plot.
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))
png(filename = paste0(opt$results,"/2_dimplot_PCA.png"),width = 2000,height = 1600,res = 300)
plot(DimPlot(DATA.i, reduction = "pca", group.by = "genotype", cols=getPalette(length(levels(as.factor(DATA.i$genotype))))))
invisible(dev.off())

## Principal component study using Elbow plot.
png(filename = paste0(opt$results,"/3_elbowplot.png"),width = 2000,height = 1600,res = 300)
plot(ElbowPlot(DATA.i, ndims = 50) + theme(legend.position="bottom"))
invisible(dev.off())

##  Save expression matrix (scaled, corrected and integrated)
cat("\nSaving Robject with integrated data")
rm(g2m.genes,s.genes,getPalette,p1)
save.image(file = paste0(opt$rdsPath,"/integrated_seurat_object.RData") )
# load("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds/integrated_seurat_object.RData")


##----------------------------------------------------------············· STEP 4
## Clustering

seeExpression(dataObj=DATA.i,
              metadata=dataset_metadata,
              splitBy="orig.ident",
              gen=c("genotype","assay"),
              expBase=0,
              genes=c("cd79a","cd79b"))

subDir <- paste0(opt$results,"/Maincells")
if (!file.exists(subDir)){dir.create(subDir)} 

All.PostSCT <- FindNeighbors(object = DATA.i,reduction = "pca",assay = "RNA",dims = 1:40)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.01)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.05)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.1)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.2)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.3)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.4)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.5)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.6)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.7)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.8)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 0.9)
All.PostSCT <- FindClusters(object = All.PostSCT, resolution = 1)

unsup.clust.colors <- pal_igv("defaul", alpha = 1)(30)
clustree(All.PostSCT, prefix = "RNA_snn_res.")

# 0.1 ...
cat("\nPerforming Clustering")
DATA.i <- FindNeighbors(object = DATA.i,reduction = "pca",dims = 1:40)
DATA.i <- FindClusters(DATA.i, resolution = 0.1, algorithm = 1) ## celltype
DATA.i <- RunUMAP(DATA.i,dims = 1:40, n.components = 2, verbose = FALSE, future.seed = NULL)

print("Clusters found (and cells per cluster)...")
print(table(DATA.i@meta.data$seurat_clusters)) ## 12 clusters

## First plot
plot(DimPlot(DATA.i, reduction = "umap", label = TRUE, pt.size = 0.00001,label.size = 3.5))

#clust.order<-as.character(c(8,7,11,5,10,9,12,0,6,3,1,4,2))
#clust.order<-as.character(c(7,3,6,9,8,0,4,5,1,2))
#ClustersColors<-PutColor(name.variables = clust.order)

# ClustersColors<-c("0"="#f5776f",
#                   "1"="#f2a733",
#                   "2"="#c4a411",
#                   "3"="#92ae0e",
#                   "4"="#48ac30",
#                   "5"="#3acb92",
#                   "6"="#02c1ab",
#                   "7"="#0cbbda",
#                   "8"="#01acfc",
#                   "9"="#8186dc",
#                   "10"="#d77cca",
#                   "11"="#f862dc",
#                   "12"="#fd66a5")

invisible(dev.off())
tiff(filename = paste0(subDir,"/UMAP_initial2.tiff"),width = 1350,height = 1150,res = 300)
plot(DimPlot(DATA.i, reduction = "umap", label = TRUE, pt.size = 0.0001,label.size = 2.5)) #,cols=ClustersColors
invisible(dev.off())

cat("\nSaving Robject with integrated data")
rm(clust.order,All.PostSCT,unsup.clust.colors)

save.image(file = paste0(opt$rdsPath,"/clustered_seurat_object.RData"))
# load("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds/clustered_seurat_object.RData")

genes<-c("cd79a","cd79b")
genes<-c("Tbx21","Itgax","Itgam")
genes<-c("cd3d","cd3e")
genes<-"cd19"

m<-seeExpression(dataObj=DATA.i,
              metadata=dataset_metadata,
              splitBy="seurat_clusters",
              gen=c("genotype","replicate"),
              expBase=0,
              genes=genes)
rownames(m)<-m$id
m[order(m$propCells,decreasing = T),1:5]


## Con esta info nos podemos saltar la anotacion, seleccionamos el cluster 0

##----------------------------------------------------------············· STEP 5
## Annotation 

bycell<-readRDS(paste0(opt$mainPath,sep="","/maincells_CELL.rds"))
bybulk<-readRDS(paste0(opt$mainPath,sep="","/maincells_BULK.rds"))

CellColors<-c("Neutrophil"="#CC3333",
              "Mono"="#CC9900","FCGR3A_mono"="#996666","CD14_mono"="#CC9999",
              "Macrophage"="#FF9966",
              "RBC"="#990000",
              "Platelet"="#FF6633",
              "DC"="#CC99CC",
              "NK_cell"="#66CC33",
              "T_cell"="#339999",
              "B_cell"="#336699", #336699
              "Myocyte"="#FFCCCC",
              "Adipocyte"="#99CCFF",
              "Epithelial"="#FFFFCC",
              "Endothelial"="#99CC99",
              "Fibroblast"="#99FFCC")

for(i in 1:length(bycell)){
  bycell[[i]]<-str_to_title(bycell[[i]])}
for(i in 1:length(bybulk)){
  bybulk[[i]]<-str_to_title(bybulk[[i]])}

mergedAnn<-list()
mergedAnn[["B_cell"]]<-unique(c(bybulk$B_cell,bycell$B_cell))
mergedAnn[["DC"]]<-unique(c(bybulk$DC,bycell$DC)) #bulk
mergedAnn[["Mono"]]<-unique(c(bybulk$Mono))
mergedAnn[["Macrophage"]]<-unique(c(bycell$Macrophage,bycell$Macrophage)) #cell
mergedAnn[["Neutrophil"]]<-unique(c(bybulk$Neutrophil,bycell$Neutrophil)) 
mergedAnn[["NK_cell"]]<-unique(c(bybulk$NK_cell,bycell$NK_cell))
mergedAnn[["T_cell"]]<-unique(c(bybulk$T_cell,bycell$T_cell))
mergedAnn[["RBC"]]<-unique(c(bycell$RBC))
mergedAnn[["Fibroblast"]]<-unique(c(bycell$Fibroblast))
mergedAnn[["Myocyte"]]<-unique(c(bycell$Myocyte))
CellColors<-CellColors[names(mergedAnn)]

DATA.i<-AnnotateCells(seurat=DATA.i,
                      geneset=mergedAnn,
                      wd=1425,
                      ht=995,
                      cor.diff=0.15,
                      cor.min=0,
                      metadata.name="cor_cell",
                      path=subDir,
                      ClustersColors=CellColors)

markerDots(seurat=DATA.i,
           genes = c("Gzma","Nkg7","Klrd1",
                     "Cd3g","Cd3d","Cd28",
                     "Cd19","Cd79a","Ms4a1",
                     "Flt3","Clec10a","Tbc1d4",
                     "Clec4a1","Clec4a3","Ear2",
                     "Ctss","Ms4a6d","Lyz2",
                     "S100a8","Cxcr2","Mmp8",
                     "Alas2","Snca"),
           path=paste0(subDir,"/Dotplot.tiff"))


## Save data
rm(list=setdiff(ls(),c("DATA.i","DATA","opt","dataset_metadata",
                       "CellColors","ClustersColors","infoSamples")))
save.image(file = paste0(opt$rdsPath,"/maincells_seurat_object.RData") )
#load("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds/maincells_seurat_object.RData")

## Save specific celltypes (B cells)
selCells<-ifelse(#DATA.i@meta.data$cor_cell=="B_cell" & 
  (DATA.i@meta.data$seurat_clusters==0  
     | DATA.i@meta.data$seurat_clusters==5 
   #  | DATA.i@meta.data$seurat_clusters==12
  ),T,F)
table(selCells)
cells<-unique(rownames(DATA.i@meta.data[selCells,]))
saveRDS(cells,paste0(opt$rdsPath,"/Bcells_05.rds"))

table(DATA.i@meta.data[cells,"genotype"])
# BANK1   TLR    WT 
# 7090  5001  7149 

selCells<-ifelse(#DATA.i@meta.data$cor_cell=="B_cell" & 
  (DATA.i@meta.data$seurat_clusters==0  
   #  | DATA.i@meta.data$seurat_clusters==5 
   #  | DATA.i@meta.data$seurat_clusters==12
  ),T,F)
table(selCells)
cells<-unique(rownames(DATA.i@meta.data[selCells,]))
saveRDS(cells,paste0(opt$rdsPath,"/Bcells_0.rds"))

table(DATA.i@meta.data[cells,"genotype"])
# BANK1   TLR    WT 
# 6101  4244  5979 


saveRDS(infoSamples,paste0(opt$rdsPath,"/infoSamples.rds"))
