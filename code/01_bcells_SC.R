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
          selected_cells ="/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds/Bcells_05.rds",
          normFactor = 10000,
          pct_mito_range = c(0,25),
          pct_ribo_range = c(0,50),
          keep_genes =c("Ighd", "Ighm", "Ighg1", "Ighg2c", "Ighg2b", "Ighg3", "Igha", "Ighe"),
          remove_non_coding =T,
          remove_gene_family= c("mt","rpl","rps","hb","malat1"),
          min_gene_count = 5,
          min_gene_per_cell =200,
          vars_to_regress = c("percent.mito","S.Score", "G2M.Score"))

source(opt$resourcesPath)
invisible(lapply(c('Seurat','dplyr','rafalib','Matrix','parallel','biomaRt','pheatmap',
                   'optparse','utils','matrixStats','patchwork','scCATCH','ggplot2',
                   'SingleCellExperiment','scales','RColorBrewer','vegan','ineq',
                   'igraph','sva','scran','scater','batchelor','clustree','optparse',
                   'scales','fields','data.table','scDblFinder','harmony','ggsci',
                   'tidyr','tibble','reshape','stringr'),require,character.only = TRUE))

## Create Subfolders
if (!file.exists(paste0(opt$mainPath,"/Rds_bcells"))){
  dir.create(paste0(opt$mainPath,"/Rds_bcells"))} 
opt$rdsPath<-paste0(opt$mainPath,"/Rds_bcells")

if (!file.exists(paste0(opt$mainPath,"/QC_bcells"))){
  dir.create(paste0(opt$mainPath,"/QC_bcells"))} 
opt$QCpath<-paste0(opt$mainPath,"/QC_bcells")

if (!file.exists(paste0(opt$mainPath,"/Results_bcells"))){
  dir.create(paste0(opt$mainPath,"/Results_bcells"))} 
opt$results<-paste0(opt$mainPath,"/Results_bcells")

## biomart objects
mart = useMart("ensembl", dataset = paste0(casefold(opt$specie),"_gene_ensembl"),host="https://jul2019.archive.ensembl.org")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://jul2019.archive.ensembl.org") 

#infoSamples<-readRDS("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds/infoSamples.rds")

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
## 27567 61777 

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
              gen=c("genotype","assay"),
              expBase=0,
              genes=c("cd79a","cd79b"))

# RemoveSamples<-c("odqmndmf_j6skyqh9","y6u8cqm4_vn5l2woi")
# for(i in 1:length(RemoveSamples)){
#   sel<-!grepl(RemoveSamples[i],colnames(DATA))
#   DATA<-DATA[,sel]  
# }

##--------------------------------- OPTIONAL


## Select subset of cells (for subclustering of specific previously annotated cell lines)
if (casefold(opt$selected_cells) != "none"){
  cat("\nRemoving non-selected cells from the Seurat object\n")
  cells_use<-readRDS(opt$selected_cells)
  selected_cells <- colnames(DATA) %in% cells_use
  DATA <- subset(DATA, cells=colnames(DATA)[selected_cells])
}
stopCluster(cl);   print(dim(DATA))
rm(all_genes,datasets,i,cl)

save.image(file = paste0(opt$rdsPath,"/raw_seurat_object.RData") )
# load("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds_bcells/raw_seurat_object.RData")


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

## Los dobletes ya los quitamos en el primer filtrado, no volver a quitar
Ts <- data.frame(
  # Doublets = DATA$doublets=="singlet",
  MitoT = between(DATA$perc_mito, mito_range[1], mito_range[2]),
  RpsT = between(DATA$perc_rps, ribo_range[1], ribo_range[2]),
  RplT = between(DATA$perc_rpl, ribo_range[1], ribo_range[2]),
  nUMIT = between(NF,quantile(NF,probs = c(0.005)),quantile(NF,probs = c(0.995))),
  nCountT = between(NC,quantile(NC,probs = c(0.005)),quantile(NC,probs = c(0.995))),
  GiniT = between(DATA$gini_index,0.8,1),
  SimpT = between(DATA$simp_index,0.8,1),
  protein_codingT = between(DATA$perc_protein_coding,50,100),
  row.names = rownames(DATA@meta.data) )
print(head(Ts,10))

cell_use <- rownames(Ts)[ rowSums(!Ts) == 0 ] ## Cell_use contains cell ids to keep

cat("\nDimentions of the raw.data objects BEFORE filtering ...\n")
print( dim(DATA@assays[[opt$assay]]@counts) ) # 27567 19240

rm(aaa,human,mart,rel_expression,Ts,family.genes,parameters,o,item,temp,hvgs,
   mito_range,ribo_range,most_expressed,NC,NF,z,seq_depth,dec,dta,dta.cnt,
   dbl.calls,dbl.dens,i,tmp)


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
DATA <- CreateSeuratObject(counts = DATA@assays[[opt$assay]]@counts[,cell_use] , assay = opt$assay, meta.data = DATA@meta.data[cell_use,], min.cells = as.numeric(opt$min_gene_count),min.features = as.numeric(opt$min_gene_per_cell))

cat("\nDimentions of the raw.data objects AFTER filtering ...\n")
print( dim(DATA@assays[[opt$assay]]@counts) )
# 14921 19017

## Normalize counts
cat("\nNormalizing counts ...\n")
DATA <- NormalizeData(object = DATA,scale.factor = opt$normFactor, normalization.method = "LogNormalize")

for(i in opt$varstoPlot){
  png(filename = paste0(opt$QCpath,"/QC_",i,"_FILTERED.png"),width = 1600*(length(unique(DATA@meta.data[,i]))/2+1),height = 1200*ceiling(length(feats)/5),res = 300)
  print(VlnPlot(object = DATA, features  = feats, ncol = 5,group.by = i,pt.size = 0,assay = opt$assay))
  invisible(dev.off())}

cat("\nSaving filtered Seurat object ...\n")
rm(annot,cell_use,feats,genes_use,i,genes_keep,genes_notfound,gn)
save.image(file = paste0(opt$rdsPath,"/filt_seurat_object.RData") )

seeExpression(dataObj=DATA,
              metadata=dataset_metadata,
              splitBy="orig.ident",
              gen=c("genotype","assay"),
              expBase=0,
              genes=c("cd79a"))

seeExpression(dataObj=DATA,
              metadata=dataset_metadata,
              splitBy="orig.ident",
              gen=c("genotype","assay"),
              expBase=0,
              genes=c("Tbx21","Itgax","Itgam"))



m$orden<-1:nrow(m)

rownames(infoSamples)<-infoSamples$id
totalcells<-NULL
for(gen in unique(m$genotype)){
  tmp<-m[m$genotype==gen,]
  tmp<-sum(infoSamples[rownames(tmp),]$filterCells)
  totalcells<-c(totalcells,tmp)
}
names(totalcells)<-unique(m$genotype)

propABCTotal<-NULL
for(i in 1:nrow(m)){
  x<-m$nCells[i]/as.numeric(totalcells[m$genotype[i]])
  propABCTotal<-c(propABCTotal,x)
}
m$propABCTotal<-propABCTotal

ggplot(m,aes(x=orden,y=propABCTotal,fill=genotype))+geom_bar(stat="identity")

# res<-seeExpression(dataObj=DATA,
#                    metadata=dataset_metadata,
#                    splitBy="orig.ident",
#                    gen=c("genotype","replicate"),
#                    expBase=0,
#                    genes=c("Tbx21","Itgax","Itgam"))
# res<-res[,-1]
# res<-res[order(res$genotype,res$replicate),]
# 
# res$propABCs<-res$nCells/(11498)
# 
# tiff(filename = "/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Results_bcells/Bcells/ABCproportion.tiff",res = 300,width = 700,height = 700)
# ggplot(res,aes(y=propABCs,x=genotype, fill=genotype))+geom_boxplot(outlier.alpha = 0)+
#   theme_classic()+scale_fill_manual(values=c("WT"="#ed8f7b","Bank1"="#47bf71","TLR7"="#3b7878"))+
#   theme(axis.text.x=element_text(angle=90))+
#   xlab("") + ylab("%ABC")
# invisible(dev.off())


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
rm(g2m.genes,s.genes,getPalette)
save.image(file = paste0(opt$rdsPath,"/integrated_seurat_object.RData") )


##----------------------------------------------------------············· STEP 4
## Clustering

subDir <- paste0(opt$results,"/Bcells")
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

# 0.2 ...
cat("\nPerforming Clustering")
DATA.i <- FindNeighbors(object = DATA.i,reduction = "pca",dims = 1:40)
DATA.i <- FindClusters(DATA.i, resolution = 0.2, algorithm = 1) ## celltype
DATA.i <- RunUMAP(DATA.i,dims = 1:40, n.components = 2, verbose = FALSE, future.seed = NULL)

print("Clusters found (and cells per cluster)...")
print(table(DATA.i@meta.data$seurat_clusters)) ## 9 clusters

## First plot
plot(DimPlot(DATA.i, reduction = "umap", label = TRUE, pt.size = 0.00001,label.size = 3.5))

## Provisional
#clust.order<-as.character(c(0:11))
#clust.order<-as.character(c(3,2,1,0,4))

ClustersColors<-c("0"="#fa8881",
                  "1"="#e89614",
                  "2"="#c69a0a",
                  "3"="#9dba1c",
                  "4"="#27b806",
                  "5"="#1fbc7c",
                  "6"="#02c0ab",
                  "7"="#07bcda",
                  "8"="#2cb4f5",
                  "9"="#8b93ff",
                  "10"="#d980ff",
                  "11"="#f36bdb",
                  "12"="#ce657d")

#ClustersColors<-PutColor(name.variables = clust.order)

invisible(dev.off())
tiff(filename = paste0(subDir,"/UMAP_initial.tiff"),width = 1200,height = 950,res = 300)
plot(DimPlot(DATA.i, reduction = "umap", label = TRUE, pt.size = 0.0001,label.size = 3,cols=ClustersColors)) 
invisible(dev.off())

cat("\nSaving Robject with integrated data")
rm(clust.order,All.PostSCT,unsup.clust.colors,cells_use,p1,selected_cells,sel)

save.image(file = paste0(opt$rdsPath,"/clustered_seurat_object.RData"))
# load("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds_bcells/clustered_seurat_object.RData")

## Filter small clusters (VER CUANDO TENGAMOS TODAS LAS MUESTRAS)
# smallCl<-names(table(DATA.i@meta.data$seurat_clusters))[table(DATA.i@meta.data$seurat_clusters)>=100]
# sel<-ifelse(DATA.i@meta.data$seurat_clusters %in% smallCl,T,F)
# 
# DATA.i <- subset(DATA.i, cells=colnames(DATA.i)[sel])
# DATA.i@meta.data$seurat_clusters<-factor(DATA.i@meta.data$seurat_clusters,
#                                             labels=smallCl)
# 
# clust.order<-as.character(c(3,0,5,7,6,2,1,4))
# ClustersColors<-PutColor(name.variables = clust.order)

ClustersColors2<-c("0"="lightgrey",
                  "1"="lightgrey",
                  "2"="lightgrey",
                  "3"="lightgrey",
                  "4"="lightgrey",
                  "5"="lightgrey",
                  "6"="lightgrey",
                  "7"="#07bcda",
                  "8"="lightgrey",
                  "9"="lightgrey",
                  "10"="lightgrey")

tiff(filename = paste0(subDir,"/figure5_parte1.tiff"),width = 1200,height = 900,res = 300)
plot(DimPlot(DATA.i, reduction = "umap", label = FALSE, pt.size = 0.0001,label.size = 3,
             cols=ClustersColors2)) 
invisible(dev.off())

## FIGURE 2A ###################################################################
tiff(filename = paste0(subDir,"/figure2A.tiff"),width = 1200,height = 900,res = 300)
plot(DimPlot(DATA.i, reduction = "umap", label = TRUE, pt.size = 0.0001,label.size = 3,
             cols=ClustersColors)) 
invisible(dev.off())

## Split by Genotype

data1 <- subset(x = DATA.i, subset = genotype == "WT")
data2 <- subset(x = DATA.i, subset = genotype == "TLR")
data3 <- subset(x = DATA.i, subset = genotype == "BANK1")

p1<-plot(DimPlot(data1, reduction = "umap", label = FALSE, pt.size = 0.0001,label.size = 3,
             cols=ClustersColors)) 
p2<-plot(DimPlot(data2, reduction = "umap", label = FALSE, pt.size = 0.0001,label.size = 3,
                 cols=ClustersColors)) 
p3<-plot(DimPlot(data3, reduction = "umap", label = FALSE, pt.size = 0.0001,label.size = 3,
                 cols=ClustersColors)) 

library(ggpubr)


tiff(filename = paste0(subDir,"/figure1e.tiff"),width = 24,height = 8,res = 300,units="cm")
plot(ggarrange(p1,p2,p3,nrow=1,ncol=3,legend = "none"))
invisible(dev.off())

## Split by Genotype




m<-seeExpression(dataObj=DATA.i,
                 metadata=dataset_metadata,
                 splitBy="seurat_clusters",
                 gen=NULL,
                 expBase=0,
                 genes=c("Tbx21","Itgax","Itgam"))

#m$id<-paste0("cl",sep="",m$id)
#m[order(m$propCells,decreasing=T),]

invisible(dev.off())
tiff(filename = paste0(subDir,"/ABCprop.tiff"),width = 600,height = 400,res = 300)
ggplot(m,aes(x=id,y=propCells,fill=id))+ geom_bar(stat = "identity",color="black",linewidth=0.2)+
  theme_bw()+ylab("%ABCs") + xlab("Cluster")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 5.8,angle = 0,hjust = 1),
        axis.text.y = element_text(size = 5.8,angle = 0,hjust = 1),
        axis.title = element_text(size=6),
        axis.ticks.x = element_line(colour = 'black', size = 0.3),
        axis.ticks.y = element_line(colour = 'black', size = 0.3),
        axis.line = element_line(colour = 'black', size = 0.3)) +
  scale_fill_manual(values=ClustersColors)
invisible(dev.off())



#ggplot(m,aes(x=id,y=propCells,fill=fill))+geom_bar(stat="identity")

#ifn<-c("Ifitm3","Ifitm2","Ifitm1","Isg15") #,

# seeExpression(dataObj=DATA.i,
#               metadata=dataset_metadata,
#               splitBy="seurat_clusters",
#               gen=NULL,
#               expBase=0,
#               genes=ifn)

# markerDots(seurat=DATA.i,
#            genes = ifn,
#            wd=1000,
#            ht=600,
#            path=paste0(subDir,"/IFN.tiff"))



##----------------------------------------------------------············· STEP 5
## Annotation
#https://www.nature.com/articles/s41577-020-00446-2

# - *Naive B cell*: Ighd, Ptprc, Cd24a, Cd22, Btla, Cd5
# - *Marginal zone B cell*:  Ptprc, Cd1d, Cd1d1, Cr2, Cd9, Ighm [-Fcer2a, -Ighd]           # Cr1, "Cd1d"
# - *Follicular B cell*: Fcer2a, Ptprc, Ighd, Fcr2, Fcer1g, Cxcr5 [-Cr1, -Cd1d, -Ighm]
# - *Germinal center B cell*: Bcl6, S1pr2, Aicda, Cd69, Mki67, Foxo1                         # Ly77
# - *Pre-memory B cell*: Ccr6, Bmpr1a, Cd59a "Ncoa1", "Efnb1", "Cd86",
# - *Memory B cell*: Cd80, Cd38, Nt5e [-Sdc1, -Prdm1, -Xbp1] Pdcd1lg2, Cd27, Zeb2, Mndal, Tle3
# - *Plasma cell*: Sdc1, Prdm1, Xbp1 [-Cd19, -Ptprc] Irf4, Cd28, Tnfrsf17, Bhlha15
# - *ABCs: "Tbx21","Itgax","Itgam"

# markers<-c("Ighd", "Ptprc", "Cd24a", "Cd22","Btla","Cd5","Ms4a1","Tcl1","Cd53",
#            "Cr1", "Cd1d", "Cd1d1", "Cr2", "Cd9","Ighm","S1pr3","Tlr3","Adam28","Ebf1","Slc22a2",
#            "Fcer2a", "Fcr2", "Fcer1g", "Cxcr5","Lrrk2","Gdf11","Icos1","Pax5",
#            "Bcl6", "S1pr2", "Cd69", "Mki67", "Foxo1", "Aicda", "Ly77","Fas","Neil1","Plxnb2","Tnfsf9","Rgs13","Tnfrsf13c","Egr1","Pif1",
#            "Ccr6", "Bmpr1a", "Cd59a", "Ncoa1", "Efnb1", "Cd86","Klf2","Pbx3","Etv6",
#            "Cd80", "Cd38", "Pdcd1lg2", "Cd27", "Zeb2", "Mndal", "Tle3", "Nt5e", "Sdc1", "Prdm1", "Xbp1","Itga4","Cxcr3","Cd47","Slamf1","Spib","Pecam1","Pml", 
#            "Cd19","Irf4", "Cd28", "Tnfrsf17", "Bhlha15","Il6r","Cd93","Ighg1","Ighg","Ighg2c","Atf6","Lgals1","Itgb5","Selpig",
#            "Tbx21","Itgax","Itgam",)

markersB<-read.csv(file = "/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/bcell_markers3.csv")

## Correlation

markerCor<-markersB[markersB$Expression=="pos",1:2]
cellsid<-unique(markerCor$CellType)
markersCor.list<-list()
for(i in 1:length(cellsid)){
  markersCor.list[[cellsid[i]]]<-as.character(markerCor[markerCor$CellType==cellsid[i],"Gene"])
}

DATA.i<-AnnotateCells(seurat=DATA.i,
                      geneset=markersCor.list,
                      wd=1425,
                      ht=995,
                      cor.diff=0.15,
                      cor.min=0,
                      metadata.name="cor_cell",
                      path=subDir,
                      ClustersColors=ClustersColors)

## Pensar como automatizar esto
## Elinar genes que no se expresen apenas y los que no aparecen en los datos
## Ver que hacer con los duplicados, lo ideal sería ir uniendo las matrices con los nombres de los genes cambiados
## Sacar Figura heatmap, con anotación positiva y negativa
## Con gaps


m<-markerDots(seurat=DATA.i,
           genes = markersB$Gene,
           wd=1000,
           ht=600,
           path=paste0(subDir,"/Bsubsets.tiff"))

markersB<-markersB[markersB$Gene %in% rownames(DATA.i@assays$RNA),]

notexp<-c("Cd5","Tlr3","Adam28","Aicda","Mki67",
          "S1pr2","Neil1","Plxnb2","Tnfsf9", #,"Fas"
          "Rgs13","Efnb1","Bmpr1a","Cd59a","Sdc1",
          "Prdm1","Bhlha15","Cd28","Tnfrsf17","Cd93",
          "Ighg1","Ighg2c","Itgb5","Pif1","Cd80","Cxcr3")


markersB<-markersB[!markersB$Gene %in% notexp,]

subB<-unique(markersB$CellType)

library('pheatmap')
cells<-NULL
ann<-NULL
for(subp in 1:length(subB)){
  
  tmp<-markersB[markersB$CellType == subB[subp],]
  tmp$geneplot<-paste0(tmp$CellType,sep=" ",tmp$Gene)
  genes = tmp$Gene
  
  ##
  metadata<-DATA.i@meta.data
  metadata$cell_id<-rownames(metadata)
  
  genes<-str_to_title(genes)
  genes<-intersect(genes,rownames(DATA.i@assays$RNA@data))
  expression_matrix <- DATA.i@assays$RNA@data[genes,,drop=F]
  
  cell_markers <- expression_matrix %>%
    Matrix::t() %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    pivot_longer(cols = genes, names_to = "gene", values_to = "expression") %>%
    left_join(metadata[,c("cell_id","seurat_clusters")])
  
  count_per_clusters <- unique(cell_markers[,c("cell_id","seurat_clusters")]) %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(cells = n())
  
  percentage_expression <- cell_markers %>%
    mutate(bool = expression > 0) %>%
    filter(bool) %>%
    group_by(seurat_clusters,gene) %>%
    dplyr::summarise(count = n()) %>%
    left_join(count_per_clusters) %>%
    mutate(percentage = (count / cells)*100)
  
  # Escalar la expresión por gen
  expression_vals <- rbindlist(mclapply(unique(cell_markers$gene), function(gene_to_use){
    expression_vals <- cell_markers[cell_markers$gene == gene_to_use,]
    expression_vals$scale_vals <- scale(expression_vals$expression)
    expression_vals
  })) %>%
    # Calcular la mediana de los valores escalados
    group_by(seurat_clusters,gene) %>%
    dplyr::summarise(median_expression = median(scale_vals)) %>%
    left_join(percentage_expression)
  
  expression_vals$gene <- factor(expression_vals$gene, levels = unique(cell_markers$gene))
  m<-expression_vals[,c("seurat_clusters","percentage","gene")]
  m[is.na(m)]<-0; m<-as.data.frame(m)
  
  clusters<-as.numeric(names(table(m$seurat_clusters)))
  m.tmp<-matrix(ncol=length(clusters),nrow=length(genes))
  colnames(m.tmp)<-clusters; rownames(m.tmp)<-genes
  
  for(i in 1:length(genes)){
    for(j in clusters){
      punt<-j+1
      m.tmp[i,punt]<-m[m$seurat_clusters==clusters[punt] & m$gene==genes[i],]$percentage
    }
  }

  m<-m.tmp[tmp$Gene,]
  rownames(m)<-tmp$geneplot
  ann<-rbind(ann,tmp)
  
  cells<-rbind(cells,m)
  
}

ann<-ann[!duplicated(ann$geneplot),]
cells<-cells[!duplicated(rownames(cells)),]

rownames(ann)<-ann$geneplot
ann<-ann[,c("CellType","Expression")]
colnames(cells)<-paste0("Cl",sep="",colnames(cells))

ann2<-data.frame("celltype"=ann$CellType)
rownames(ann2)<-rownames(ann)

dev.off()
tiff(filename="/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells/Bcells/SubsetMarkers.tiff",res=300, height = 845,width = 1650)
pheatmap(t(cells),scale="column",fontsize = 6.5, 
         border_color = "black",
         gaps_col = c(11,18,22,28,32,44,51),
         #gaps_col = c(11,18,22,28,32,44,50),
         #gaps_col = c(12,21,25,40,47,63,79),
         cluster_cols = F, cluster_rows = T,
         breaks = seq(-1.5,1.5,length.out = 100),
         annotation_col = ann2,
         color = colorRampPalette(c("#4393c3","white", "#d6604d"))(100))
dev.off()

## Plot mas reducido

cells2<-cells[c(1,2,3,),]


#
# 
# 
# abcs<-c("Tbx21","Itgax","Itgam")
# markerDots(seurat=DATA.i,
#            genes = abcs,
#            wd=1000,
#            ht=600,
#            path=paste0(subDir,"/ABC.tiff"))
# 
# bsubsets<-c("ighd",
#             "cd27","cd80",
#             "cd1d1",
#             "Tbx21","Itgax","Itgam")
# 
# ## Naive b cell
# markerDots(seurat=DATA.i,
#            genes = c("cd27","cd80"),
#            wd=1000,
#            ht=600,
#            path=paste0(subDir,"/Bsubsets1.tiff"))
# 
# 
# bycell<-readRDS("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Bcells_CELL.rds")
# bycell<-bycell[-c(6,7)]
# bycell<-bycell[c("NaiveBcell","GC_LightZone","ABC","MemoryCells")] #,"GC_DarkZone","MarginalZone"
# bycell$ABC<-unique(c(bycell$ABC,"Tbx21","Itgax","Itgam"))
# 
# 
# CellColors<-c("NaiveBcell"="#df7656",
#               "MarginalZone"="#f29a42",
#               "PlasmaCells"="#7c67f8",
#               "GC_DarkZone"="#45659e",
#               "GC_LightZone"="#f6c4a7",
#               "MemoryCells"="#a6b8e4",
#               "PreMemCells"="#ff8fad",
#               "ABC"="#0099b2")
# 
# GenotypeColors<-c("WT"="#ed8f7b","BANK1"="#47bf71","TLR"="#3b7878")
# 
# DATA.i<-AnnotateCells(seurat=DATA.i,
#                       geneset=bycell[c(1,3)],
#                       wd=1700,
#                       ht=1300,
#                       cor.diff=0.2,
#                       cor.min=0.25,
#                       pt.size=0.0001,
#                       metadata.name="subB",
#                       path=subDir,
#                       ClustersColors=CellColors)

ProportionPlot(seuratObj=DATA.i,
               wd = 1000,
               ht=500,
               idColors = CellColors,
               path=subDir,
               column_split="seurat_clusters", ## clusters, genotypes...
               column_tags="subB") ## cell assigment column


# sel<-ifelse(DATA.i@meta.data$genotype=="TLR",T,F)
# samples<-unique(DATA.i@meta.data[sel,"orig.ident"])
# for(smpl in samples){
#   print(table(DATA.i@meta.data[ifelse(DATA.i@meta.data$orig.ident==smpl,T,F),"seurat_clusters"]))
# }
# 
# length(DATA.i@meta.data[ifelse(DATA.i@meta.data$genotype=="TLR7" & DATA.i$subB=="ABC",T,F),"subB"])
# 
# length(DATA.i@meta.data[ifelse(DATA.i@meta.data$genotype=="Bank1" & DATA.i$subB=="ABC",T,F),"subB"])
# 
# length(DATA.i@meta.data[ifelse(DATA.i@meta.data$genotype=="WT" & DATA.i$subB=="ABC",T,F),"subB"])


##----------------------------------------------------------············· STEP 7
## Figures Proportion

# ProportionPlot(seuratObj=DATA.i,
#                wd = 1000,
#                ht=500,
#                idColors = GenotypeColors,
#                path=subDir,
#                column_split="seurat_clusters", ## clusters, genotypes...
#                column_tags="genotype") ## cell assigment column

## Figura
ProportionPlot(seuratObj=DATA.i,
               wd = 600,
               ht=1550,
               idColors = ClustersColors,
               path=subDir,
               column_split="genotype", ## clusters, genotypes...
               column_tags="seurat_clusters") ## cell assigment column

GenotypeColors<-c("WT"="#ed8f7b","BANK1"="#47bf71","TLR"="#3b7878")

ProportionPlot(seuratObj=DATA.i,
               wd = 1000,
               ht=1200,
               idColors = GenotypeColors,
               path=subDir,
               column_split="seurat_clusters", ## clusters, genotypes...
               column_tags="genotype") ## cell assigment column

## ProportionPlot
# metadata <- DATA.i@meta.data # Coger metadata 
# metadata.barr<-metadata[,c("seurat_clusters","genotype")]
# metadata.barr$seurat_clusters<-paste0("Bcell_",metadata.barr$seurat_clusters)
# 
# metadata.barr$seurat_clusters<-as.factor(metadata.barr$seurat_clusters)
# #metadata.barr$seurat_clusters<-factor(x=metadata.barr$seurat_clusters,levels=)
# 
# metadata.barr$genotype<-factor(x = metadata.barr$genotype,
#                                levels = c("WT","Bank1","TLR7"))
# 
# p2<-ggplot(metadata.barr, aes(x = seurat_clusters))+
#   geom_bar(aes(fill = genotype), position = "fill", color = "black", lwd = 0.2)+
#   theme_classic()+
#   theme(axis.title = element_text(size = 5),
#         axis.ticks = element_blank(),
#         axis.text.x= element_text(angle = 90),
#         axis.text = element_text(size = 4.5),
#         axis.line = element_blank(),
#         legend.position="none")+
#   coord_flip()+
#   scale_y_continuous(expand = c(0,0))+
#   ylab("Proportion of cells")+xlab("Clusters")+
#   scale_fill_manual(values = list("WT"="#ed8f7b","Bank1"="#47bf71","TLR7"="#3b7878"))
# tiff(filename = paste0(subDir,"/figure2B.tiff"),width = 320,height = 400,res = 300)
# plot(p2)
# invisible(dev.off())


##----------------------------------------------------------············· STEP 7
## Significance BOXpLOTS

metadata <- DATA.i@meta.data # Coger metadata 

sum_cells <- metadata %>%
  group_by(orig.ident) %>%
  dplyr::summarise(n_cells = n())

proportion_dataset <- metadata %>%
  group_by(orig.ident, seurat_clusters, genotype) %>%
  dplyr::summarise(cell = n()) %>%
  left_join(sum_cells) %>%
  mutate(proportion_clust = cell / n_cells)

proportion_dataset <- proportion_dataset[order(proportion_dataset$genotype),]
proportion_dataset$orig.ident <- factor(proportion_dataset$orig.ident, levels = unique(proportion_dataset$orig.ident))

p_values_genotype <- rbindlist(lapply(unique(proportion_dataset$seurat_clusters), function(cluster){ # Para hacer test de wilcoxon
  #print(cluster)
  tmp_proportion_dataset <- proportion_dataset[proportion_dataset$seurat_clusters == cluster,]
  total<-max(tmp_proportion_dataset$proportion_clust)
  
  prop_1 <- tmp_proportion_dataset[tmp_proportion_dataset$genotype == "TLR",]$proportion_clust # Cambiar por condición 1
  prop_2 <- tmp_proportion_dataset[tmp_proportion_dataset$genotype == "BANK1",]$proportion_clust # Cambiar por condición 2
  prop_3 <- tmp_proportion_dataset[tmp_proportion_dataset$genotype == "WT",]$proportion_clust # Cambiar por condición 2
  
  tlr7_wt<-c(mean(prop_1*100),mean(prop_3*100))
  tlr7_bank1<-c(mean(prop_1*100),mean(prop_2*100))
  bank1_wt<-c(mean(prop_2*100),mean(prop_3*100))
             
  # p1<-wilcox.test(c(prop_1), c(prop_3))$p.value
  # p2<-wilcox.test(c(prop_1), c(prop_2))$p.value
  # p3<-wilcox.test(c(prop_2), c(prop_3))$p.value
  
  prop_1<-(prop_1*100)/total
  prop_2<-(prop_2*100)/total
  prop_3<-(prop_3*100)/total
  
  p1<-prop.test(c(mean(prop_1),mean(prop_3)),n = rep(100,2))$p.value
  p2<-prop.test(c(mean(prop_1),mean(prop_2)),n = rep(100,2))$p.value
  p3<-prop.test(c(mean(prop_2),mean(prop_3)),n = rep(100,2))$p.value
  
  pval <- data.frame(seurat_clusters = cluster, comp = c("tlr7_wt","tlr7_bank1","bank1_wt"), 
                     pval = c(p1,p2,p3),prop = c(tlr7_wt[1]-tlr7_wt[2],
                                                 tlr7_bank1[1]-tlr7_bank1[2],
                                                 bank1_wt[1]-bank1_wt[2]),
                     max_val=max((tlr7_wt/100) + 0.1,(tlr7_bank1/100) + 0.1,(bank1_wt/100) + 0.1))
  
  #pval <- data.frame(seurat_clusters = cluster, max_val = max(c(prop_tlr7, prop_wt)), pval = wilcox.test(prop_tlr7, prop_wt)$p.value)
  
}))


p_values_genotype$pval <- ifelse(p_values_genotype$pval < 0.0001, "****",   # Esto es para evitar que salga el número
                                 ifelse(p_values_genotype$pval < 0.001, "***",
                                        ifelse(p_values_genotype$pval < 0.01, "**",
                                               ifelse(p_values_genotype$pval < 0.05,"*",""))))

## PLOT

library("ggpubr")

plotList<-list()
clusters<-unique(proportion_dataset$seurat_clusters)

maxH<-c(0.7,0.67,0.20,0.60,0.2,0.1,0.12,0.08,0.01,0.01)

res<-lapply(1:10,function(i){
  
  pval<-NULL
  
  tmp<-proportion_dataset[proportion_dataset$seurat_clusters==clusters[i],]
  pval<-p_values_genotype[p_values_genotype$seurat_clusters==clusters[i],]
  pval<-cbind(pval,do.call("rbind",strsplit(pval$comp,"_")))
  pval$max_val<-maxH[i]+(maxH[i]*0.1)
  pval$height<-c(pval$max_val[1] * 0.94,pval$max_val[1] * 0.88,pval$max_val[1] * 0.99)
  pval$alfa = ifelse(pval$pval=="",0,1)
  
  tmp$genotype<-ifelse(tmp$genotype == "BANK1","T7.B1KO",
                       ifelse(tmp$genotype=="TLR","T7","WT"))
  
  tmp<-rbind(as.data.frame(tmp[tmp$genotype=="WT",]),
             as.data.frame(tmp[tmp$genotype=="T7",]),
             as.data.frame(tmp[tmp$genotype=="T7.B1KO",]))
  
  tmp$genotype<-factor(tmp$genotype,levels = unique(tmp$genotype))
  
  p1<-ggplot(data = tmp,
         aes(x = genotype, y = proportion_clust))+
    geom_violin(aes(fill = genotype),lwd=0.3)+
    geom_boxplot(fill = "white", outlier.shape = NA, width = 0.4,lwd=0.3)+
    geom_text(data = pval, aes(y = height, x = c(1.5,2.5,2), label = formatC(pval,format = "e",digits = 1)), size = 3)+ 
    geom_segment(aes(x="T7", xend = "WT",y = pval$height[1]-(max(pval$height)*0.01),yend = pval$height[1]-(max(pval$height)*0.01)),
                 linewidth=.1,alpha=pval$alfa[1])+
    geom_segment(aes(x="T7", xend = "T7.B1KO",y = pval$height[2]-(max(pval$height)*0.01),yend = pval$height[2]-(max(pval$height)*0.01)),
                 linewidth=.1,alpha=pval$alfa[2])+
    geom_segment(aes(x="T7.B1KO", xend = "WT",y = pval$height[3]-(max(pval$height)*0.01),yend = pval$height[3]-(max(pval$height)*0.01)),
                 linewidth=.1,alpha=pval$alfa[3])+
    scale_fill_manual(values = c("T7" = "#f7df7f", "WT" = "#c2c2c2","T7.B1KO" = "#7fb9e0"))+ #Cambiar
    theme_classic()+
    ylab("Proportion of cells per sample")+
    theme(title = element_text(size=8),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 7),
          axis.text.x = element_blank(),
          strip.text = element_blank(), # element_text(size = 7)
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(colour = 'black', size = 0.3,),
          axis.ticks.length=unit(.05, "cm"),
          axis.line = element_line(colour = 'black', size = 0.3)
    )+ ylim(0,maxH[i]+(maxH[i]*0.1)) + ggtitle(paste0("Cluster ",sep="",clusters[i]))+
    xlab("genotype")
  p1
  print(pval)
  print(i)
  #plotList[[i]]<-p1
  return(p1)
  
})

dev.off()

p1<-ggarrange(plotlist = res,ncol = 5,nrow = 2,legend = "bottom",common.legend = TRUE)
plot(p1)




# gg3 <- ggplot(data = proportion_dataset,
#               aes(x = genotype, y = proportion_clust))+
#   facet_wrap(~seurat_clusters, scales = "free",ncol = 5)+
#   geom_violin(aes(fill = genotype),lwd=0.3)+
#   geom_boxplot(fill = "white", outlier.shape = NA, width = 0.4,lwd=0.3)+
#   #geom_text(data = p_values_genotype, aes(y = max_val * 0.9, x = 1.5, label = pval), size = 3)+
#   # Usar la línea siguiente si prefieres que salga el número
#   #geom_text(data = p_values_genotype, aes(y = max_val * 0.9, x = 1.5, label = formatC(pval,format = "e",digits = 1)), size = 3)+ 
#   scale_fill_manual(values = c("TLR" = "#f0c35e", "WT" = "#43a2ca","BANK1" = "#00cc66"))+ #Cambiar
#   theme_classic()+
#   ylab("Proportion of cells per sample")+
#   theme(legend.position = "none",
#         axis.title.y = element_text(size = 9),
#         axis.title.x = element_blank(),
#         axis.text = element_text(size = 5),
#         axis.text.x = element_text(size = 5),
#         strip.text = element_blank(), # element_text(size = 7)
#         strip.background = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.ticks.y = element_line(colour = 'black', size = 0.3,),
#         axis.ticks.length=unit(.05, "cm"),
#         axis.line = element_line(colour = 'black', size = 0.3)
#   )+
#   xlab("genotype")
# 
# plot(gg3)
# invisible(dev.off()

##----------------------------------------------------------············· STEP 8
## Trajectory


invisible(lapply(c('Seurat','dplyr','rafalib','Matrix','parallel','biomaRt','pheatmap',
                   'optparse','utils','matrixStats','patchwork','scCATCH','ggplot2',
                   'SingleCellExperiment','scales','RColorBrewer','vegan','ineq',
                   'igraph','sva','scran','scater','batchelor','clustree','optparse',
                   'scales','fields','data.table','scDblFinder','harmony','ggsci',
                   'tidyr','tibble','reshape','stringr'),require,character.only = TRUE))

library("scater")
library("TSCAN")

load("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds_bcells/clustered_seurat_object.RData")
source(opt$resourcesPath)

dta<-as.SingleCellExperiment(x=DATA.i, assay = "RNA")
colLabels(dta)<-dta$seurat_clusters

by.cluster <- aggregateAcrossCells(dta, ids=colLabels(dta))
centroids <- reducedDim(by.cluster, "UMAP")

mst <- createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")

# plotUMAP(dta, colour_by="seurat_clusters") + 
#   geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
# 
map.tscan <- mapCellsToEdges(dta, mst=mst, use.dimred="UMAP")

tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

common.pseudo <- averagePseudotime(tscan.pseudo) 

plotUMAP(dta,colour_by=I(common.pseudo),
         text_by="label", text_colour="red",point_size=0.0001) +
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))


tiff(filename = paste0(subDir,"/SeudoTimes.tiff"),width = 1200,height = 900,res = 300)
plotUMAP(dta,colour_by=I(common.pseudo),point_size=0.0001) +
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
invisible(dev.off())


tiff(filename = paste0(subDir,"/SeudoTimesLAB.tiff"),width = 1200,height = 900,res = 300)
plotUMAP(dta,colour_by=I(common.pseudo),
         text_by="label", text_colour="red",point_size=0.0001) +
  geom_line(data=line.data, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))
invisible(dev.off())

tiff(filename = paste0(subDir,"/SeudoTimesPoints.tiff"),width = 1400,height = 1400,res = 300)
pseudo.og <- quickPseudotime(dta, use.dimred="UMAP", outgroup=TRUE)
set.seed(10101)
plot(pseudo.og$mst)
invisible(dev.off())

pseudo.mnn <- quickPseudotime(dta, use.dimred="UMAP", with.mnn=TRUE)
mnn.pseudo <- averagePseudotime(pseudo.mnn$ordering)
plotUMAP(dta, colour_by=I(mnn.pseudo), text_by="label", text_colour="red") +
  geom_line(data=pseudo.mnn$connected$UMAP, mapping=aes(x=UMAP_1, y=UMAP_2, group=edge))

##

sce.nest <- runUMAP(dta, dimred="PCA")
reducedDim(sce.sling2, "UMAP") <- reducedDim(sce.nest, "UMAP")

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=Dim.1, y=Dim.2), size=1.2)
}

gg


#####################


#save.image("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds_bcells/trajectory_seurat_object.RData")
#load("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK3/Rds_bcells/trajectory_seurat_object.RData")


library("Seurat")
library("destiny")
library("scater")
#library("clusterExperiment")
library("SingleCellExperiment")
library("ggplot2")
library("ggbeeswarm")

DATA.sub<-DATA.i
dta<-as.SingleCellExperiment(x=DATA.sub, assay = "RNA")
deng <- logcounts(dta)

cellLabels <- dta$seurat_clusters
colnames(deng) <- cellLabels

# seleccionamos aquellos clusters que cambian entre condiciones
#selected<-ifelse(dta$seurat_clusters %in% c(0,1,3,4,6,7,8),T,F)
#deng<-deng[,selected]

# Make a diffusion map.
dm <- DiffusionMap(t(as.matrix(deng)))

save.image("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds_bcells/TrayectoryPlots.RData")

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  ClustersColors = cellLabels, #[selected],
                  genotype = dta$genotype) #[selected])

x<-tmp[order(tmp$ClustersColors,decreasing = F),]



invisible(dev.off())
tiff(filename = paste0(opt$results,"/trayectoryAll.tiff"),width = 1400,height = 1000,res = 300)
pAll<-ggplot(x, aes(x = DC1, y = DC2, fill = ClustersColors)) +
  geom_point(size=1.2, shape = 21, colour = "black",stroke = 0.2)  + scale_fill_manual(values=ClustersColors)+
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()+ xlim(-0.077,0.015) + ylim(-0.07,0.02)
plot(pAll)
invisible(dev.off())





# x2<-tmp[tmp$ClustersColors==8,]
# ggplot(x2, aes(x = DC1, y = DC2, fill = ClustersColors)) +
#   geom_point(size=2, shape = 21, colour = "black")  + scale_fill_manual(values=ClustersColors)+
#   xlab("Diffusion component 1") + 
#   ylab("Diffusion component 2") +
#   theme_classic()+ xlim(-0.077,0.015) + ylim(-0.07,0.02)

# x<-tmp[tmp$genotype=="WT",]
# p0<-ggplot(x, aes(x = DC1, y = DC2, colour = ClustersColors)) +
#   geom_point()  + scale_color_manual(values=ClustersColors)+
#   xlab("Diffusion component 1") + 
#   ylab("Diffusion component 2") +
#   theme_classic() #+ xlim(-0.08,0.015) + ylim(-0.07,0.04)
# p0
# 
# x<-tmp[tmp$genotype=="TLR",]
# p1<-ggplot(x, aes(x = DC1, y = DC2, colour = ClustersColors)) +
#   geom_point()  + scale_color_manual(values=ClustersColors)+
#   xlab("Diffusion component 1") + 
#   ylab("Diffusion component 2") +
#   theme_classic() #+ xlim(-0.08,0.015) + ylim(-0.07,0.04)
# p1
# 
# x<-tmp[tmp$genotype=="BANK1",]  # & tmp$ClustersColors==0
# p2<-ggplot(x, aes(x = DC1, y = DC2, colour = ClustersColors)) +
#   geom_point()  + scale_color_manual(values=ClustersColors)+
#   xlab("Diffusion component 1") + 
#   ylab("Diffusion component 2") +
#   theme_classic() #+ xlim(-0.08,0.015) + ylim(-0.075,0.04)
# p2



 

##----------------------------------------------------------············· STEP 11
## Mas plots

genes<-c("Itgax","Tbx21","Itgam",
         "Ahnak","Itgb1","Cd72","Hck","Zeb2","Bhlhe40","Fcrl5","Sox5","Tnfrsf1b",
         "Nt5e","Cd80","Cd44","Pdcd1lg2",
         "Zbtb32","Fas","Foxo1","Ighg1","Ighg2c")
dev.off()

i=i+1
  tiff(paste0("GENES_",genes[i],sep="",".tiff"),res=300,units = "cm",width = 7,height = 5)
    FeaturePlot(DATA.i,features = genes[i], pt.size =  0.0000001, combine = FALSE,cols = c("#e7e7e7","red"))
  dev.off()




## hEATMAP

tabl<-read.csv("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells/All_comparations.csv",sep="\t")
tabl<-tabl[tabl$cluster==7,]
rownames(tabl)<-NULL

sel<-ifelse((tabl$avg_log2FC < -1 | tabl$avg_log2FC > 1),T,F)

tabl.2<-tabl[sel,]

library(pheatmap)

metadata <- DATA.i@meta.data # Coger metadata 

x<-DATA.i@assays$RNA
x<-x[tabl$gene,]

EXP<-NULL
for(i in 0:9){
  selected<-rownames(metadata[metadata$seurat_clusters==i,])
  res<-apply(x[,selected],1,mean)
  EXP<-cbind(EXP,res)
}
colnames(EXP)<-paste0("Cl",sep="",0:9)

pheatmap(t(EXP),scale="column",show_colnames = F)

EXP<-EXP[genes,]
pheatmap(t(EXP),scale="column",show_colnames = T,border_color = "black")


#####################################
## Mas plots

gnes<-read.csv("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells/Bcells/GenotypeDEG.txt",sep="\t")

gnes<-gnes[gnes$Cluster=="cluster7",]
gnes<-gnes$Gene

gnes<-c(gnes,c("Itgax","Tbx21","Itgam",
               "Ahnak","Itgb1","Cd72","Hck","Zeb2","Bhlhe40","Fcrl5","Sox5","Tnfrsf1b",
               "Nt5e","Cd80","Cd44","Pdcd1lg2",
               "Zbtb32","Fas","Foxo1","Ighg1","Ighg2c"))

gnes<-gnes[!gnes %in% c("Bank1")]

metadata <- DATA.i@meta.data # Coger metadata 

x<-DATA.i@assays$RNA
x<-x[gnes,]

samples<-unique(metadata$orig.ident)#[4:9]

EXP<-NULL
for(i in 1:length(samples)){
  selected<-rownames(metadata[metadata$orig.ident==samples[i],])
  res<-apply(x[,selected],1,mean)
  EXP<-cbind(EXP,res)
}
colnames(EXP)<-samples

pheatmap(EXP,scale="row",breaks = seq(-1.5,1.5,length.out = 100),cluster_cols = F,gaps_col = c(3,6),border_color = "black")

####

library(ggplot2)
library("ggrepel")

go<-read.csv("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells/Bcells/functionsGEN.csv",sep=",")
go<-go[go$cluster=="cluster7",c("cluster","pval_adj","relative_enrichment","fold","description")]
go<-go[order(go$pval_adj,decreasing = T),]
paths<-go[1:15,]

go<-go[go$fold<0,]
go<-go[order(go$pval_adj,decreasing = T),]
go$pval_adj<-go$pval_adj* -1
paths<-rbind(paths,go[1:15,])
go<-paths

#go<-go[c(1:15,175:189),]
go$description<-factor(go$description,levels=unique(go$description))
go$Genotype<-ifelse(go$relative_enrichment>0,"T7","T7.B1KO")

texto<-as.character(go$description)

ggplot(go,aes(x=pval_adj,y=description,fill=Genotype))+ theme_classic()+
  geom_bar(stat="identity",colour="black",size=0.5,width=0.8)+
#  geom_vline(xintercept = c(-4,-3,-2,-1,0,1,2,3,4),linetype="dashed",color="black",size=0.1)+
  theme(axis.line.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size=8),
        legend.position="top",
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8))+
  xlab("Enrichment score")+
  annotate("text", x= -0.05,y=1, label= texto[1],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=2, label= texto[2],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=3, label= texto[3],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=4, label= texto[4],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=5, label= texto[5],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=6, label= texto[6],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=7, label= texto[7],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=8, label= texto[8],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=9, label= texto[9],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=10, label= texto[10],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=11, label= texto[11],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=12, label= texto[12],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=13, label= texto[13],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=14, label= texto[14],size = 3,hjust = 1)+
  annotate("text", x= -0.05,y=15, label= texto[15],size = 3,hjust = 1)+
  
  annotate("text", x= 0.05,y=16, label= texto[16],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=17, label= texto[17],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=18, label= texto[18],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=19, label= texto[19],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=20, label= texto[20],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=21, label= texto[21],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=22, label= texto[22],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=23, label= texto[23],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=24, label= texto[24],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=25, label= texto[25],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=26, label= texto[26],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=27, label= texto[27],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=28, label= texto[28],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=29, label= texto[29],size = 3,hjust = 0)+
  annotate("text", x= 0.05,y=30, label= texto[30],size = 3,hjust = 0)+
  

  scale_fill_manual(values=c("T7"="#ff9900","T7.B1KO"="#66cccc"))




## Terminar Figura
## Inscripcion Pedro








