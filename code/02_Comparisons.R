##############################
## SingleCell Analysis
## R version R version 4.0.4
##
##############################
## 06_Differential expression analysis


load("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds_bcells/clustered_seurat_object.RData")
source(opt$resourcesPath)
library("Seurat")
library("pheatmap")
library(RColorBrewer)
library("dplyr")
library("stringr")

setwd("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells")

##------------------------------------------------------ STEP 1
## Finding differentially expressed genes (cluster biomarkers)

##----------
## Comparacion uno contra uno

#clusters<-unique(as.character(DATA.i@meta.data$seurat_clusters))
clusters<-as.character(0:9)

DATA_markers <-getDEG(DATA.i,filterAdj = T,
                      selectClusters=clusters,
                      only.pos = FALSE,
                      nameFile=paste0(opt$results,"/All_comparations.csv"))


cat("\nPlotting heatmap of Cluster Marker genes ...\n")
top<-10
DATA_markers %>% group_by(cluster) %>% top_n(top,avg_log2FC) -> top10


#DATA.tmp<-subset(DATA.i,subset = seurat_clusters %in% as.character(c(7,8,9,16,19)))
DATA.tmp<-DATA.i
genes<-top10$gene[!duplicated(top10$gene)]
#clusters<-as.character(c(7,8,9,16,19))

m <- DATA.tmp@assays$RNA@counts
m<-m[genes,rownames(DATA.tmp@meta.data)[DATA.tmp@meta.data$seurat_clusters %in% clusters]]

## Heatmap by cluster
clusters<-as.character(0:9)
samples<-unique(as.character(DATA.i$orig.ident))

M<-NULL
nmes<-NULL
cls<-NULL
ids<-NULL

for(i in 1:length(clusters)){
  print(i)
  for(j in 1:length(samples)){
    sel<-ifelse(DATA.i$seurat_clusters==clusters[i] & DATA.i$orig.ident == samples[j],T,F)   
    
    tmp<-m[,sel]
    if(is.null(nrow(tmp))){
      res<-tmp
    }else{
      res<-apply(tmp,1,function(x){
        mean(x,na.rm=T)
      })
    }
    cls<-c(cls,clusters[i])
    x<-ifelse(str_detect(samples[j], c("WT")),"WT",
              ifelse(str_detect(samples[j],c("TLR")),"TLR7","BANK1"))
    ids<-c(ids,x)
    M<-cbind(M,res)
    nmes<-c(nmes,paste0(clusters[i],sep="_",samples[j]))
  }
  
}
rownames(M)<-rownames(m)
colnames(M)<-nmes

ann<-data.frame("cluster"=cls,
                "genotype"=ids)
rownames(ann)<-nmes

sel<-!colSums(is.na(M)) > 0
M<-M[,sel]
ann<-ann[sel,]



invisible(dev.off())
tiff(filename = paste0(opt$results,"/Heatmap2.tiff"),width = 14,height = 6.8,res = 300,units = "cm")
pheatmap(t(M),cluster_rows = F,cluster_cols = F,
         show_rownames = F,show_colnames = T,
         breaks = seq(-2.5,2.5,length.out = 100),
         #gaps_row = c(10,20,30,40,50,60,70,80,90),
         annotation_row = ann,
         border_color = NA, fontsize = 6,
         color = colorRampPalette(c("#4393c3","white", "#d6604d"))(100),
         scale = "column")
dev.off()

write.table(rownames(M)[nrow(M):1],file=paste0(opt$results,"/","genes.txt"),sep="\t")
# 
# pheatmap(t(cells),scale="column",fontsize = 6.5, 
#          border_color = "black",
#          gaps_col = c(9,18,27,32,36,49,55),
#          cluster_cols = F, cluster_rows = T,
#          breaks = seq(-1.5,1.5,length.out = 100),
#          #annotation_col = ann,
#          color = colorRampPalette(c("#4393c3","white", "#d6604d"))(100))

##


clusterpos<-as.character(DATA.tmp@meta.data[colnames(m),"seurat_clusters"])

mplot<-matrix(data=0,ncol=length(clusters),nrow=length(genes))
colnames(mplot)<-clusters
rownames(mplot)<-genes

for(i in 1:length(genes)){
  for(j in 1:length(clusters)){
    
    tmp<-m[i,clusterpos==clusters[j]]
    mplot[i,j]<-mean(as.numeric(tmp),na.r=T)
  }
}

colors<-brewer.pal(n = 9, name = 'Blues')

invisible(dev.off())
tiff(filename = paste0(opt$results,"/top.tiff"),width = 600,height = 2300,res = 300)
p1<-pheatmap(mplot,scale="row",cluster_rows = F,cluster_cols = T,
             color=c("white",colors), breaks = seq(-1.6,1.6,length.out = 9),
             clustering_method = "complete",fontsize=6,border_color="#303030")
print(p1)
invisible(dev.off())


##------------------------------------------------------ STEP 2
## Genes by Clusters

##--------
## IFN

#genes<-c("Ighg1","Ighg2","Ighg3","Igh-2a","Igh-1b","Ighg2c","Ighd","Ighm","Igha","Ighe")

DATA.tmp<-DATA.i
genes<-genes<-c("Ifng","Ifna1","Ifnar1","Ifi44l","Isg15","Ifi44","Ifi6","Ifitm3","Ifitm2","Ifitm1",
                "Oas2","Eif2ak2","Oas3","Mx1","Oasl","Irf7","Plscr1","Fcgr1","Irf9")

clusters<-unique(DATA.tmp@meta.data$seurat_clusters)

m <- DATA.tmp@assays$RNA@counts
genes<-intersect(genes,rownames(m))

m<-m[genes,]

clusterpos<-as.character(DATA.tmp@meta.data[colnames(m),"seurat_clusters"])

mplot2<-matrix(data=0,ncol=length(clusters),nrow=length(genes))
colnames(mplot2)<-clusters
rownames(mplot2)<-genes

for(i in 1:length(genes)){
  for(j in 1:length(clusters)){
    
    tmp<-m[i,clusterpos==clusters[j]]
    mplot2[i,j]<-mean(as.numeric(tmp),na.r=T)
  }
}

colors<-brewer.pal(n = 9, name = 'Blues')

tiff(filename = paste0(opt$results,"/IFN.tiff"),width = 800,height = 608,res = 300)
p2<-pheatmap(mplot2,scale="none",cluster_rows = T,cluster_cols = T,
             color=c("white",colors), breaks = seq(0,0.8,length.out = 9),
             clustering_method = "complete",fontsize=6,border_color="#303030")
print(p2)
invisible(dev.off())

#p2<-pheatmap(mplot2,scale="row",cluster_rows = T,cluster_cols = T,
#             color=c("white",colors), #breaks = seq(0,0.5,length.out = 9),
#             clustering_method = "complete",fontsize=6,border_color="#303030")
#print(p2)


##--------
## Atibodies

DATA.tmp<-DATA.i
genes<-c("Ighg1","Ighg2","Ighg3","Igh-2a","Igh-1b","Ighg2c","Ighd","Ighm","Igha","Ighe")

#genes<-genes<-c("Ifng","Ifna1","Ifnar1","Ifi44l","Isg15","Ifi44","Ifi6",
#                "Oas2","Eif2ak2","Oas3","Mx1","Oasl","Irf7","Plscr1","Fcgr1","Irf9")

clusters<-unique(DATA.tmp@meta.data$seurat_clusters)

m <- DATA.tmp@assays$RNA@counts
genes<-intersect(genes,rownames(m))

m<-m[genes,]

clusterpos<-as.character(DATA.tmp@meta.data[colnames(m),"seurat_clusters"])

mplot3<-matrix(data=0,ncol=length(clusters),nrow=length(genes))
colnames(mplot3)<-clusters
rownames(mplot3)<-genes

for(i in 1:length(genes)){
  for(j in 1:length(clusters)){
    
    tmp<-m[i,clusterpos==clusters[j]]
    mplot3[i,j]<-mean(as.numeric(tmp),na.r=T)
  }
}

colors<-brewer.pal(n = 9, name = 'Blues')

hits<-c(0,1,5,10,100)
colorseq<-c("white","#5b9ecc","#083e89","#6e3e89","#9300e3")

bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)


tiff(filename = paste0(opt$results,"/Autoantibodies.tiff"),width = 1000,height = 450,res = 300)
p3<-pheatmap(mplot3,scale="none",cluster_rows = T,cluster_cols = T,
             #color=c("white",colors), breaks = seq(0,1,length.out = 9),
             breaks=bks$breaks,color=bks$colors,
             clustering_method = "complete",fontsize=6,border_color="#303030")
invisible(dev.off())

#p2<-pheatmap(mplot2,scale="row",cluster_rows = T,cluster_cols = T,
#             color=c("white",colors), #breaks = seq(0,0.5,length.out = 9),
#             clustering_method = "complete",fontsize=6,border_color="#303030")
#print(p2)

p3<-pheatmap(mplot3,scale="row",cluster_rows = T,cluster_cols = T,
             color=c("white",colors), #breaks = seq(0,1,length.out = 9),
             clustering_method = "complete",fontsize=6,border_color="#303030")












##

































##----------
## Comparacion ABCs resto: ...................
genesSel<-NULL
comp<-0:19
comp<-comp[!comp %in% c(7,8,9,16,19)] # 7

DATA.tmp <- SetIdent(DATA.i, value = DATA.i@meta.data$seurat_clusters) 
DefaultAssay(DATA.tmp)<-"RNA"

## 7----------
DATA_markers <- FindMarkers(object = DATA.tmp,ident.1 = 7,ident.2 = comp,only.pos = F)
## Guardar matriz (Falta)
DATA_markers$p_val_adj<-p.adjust(DATA_markers$p_val,method = "bonferroni")

genesSel<-c(genesSel,
            rownames(DATA_markers[order(DATA_markers$p_val_adj,-DATA_markers$avg_log2FC),])[1:10],
            rownames(DATA_markers[order(DATA_markers$p_val_adj,DATA_markers$avg_log2FC),])[1:10])

## 8----------
DATA_markers <- FindMarkers(object = DATA.tmp,ident.1 = 8,ident.2 = comp,only.pos = F)
## Guardar matriz (Falta)
DATA_markers$p_val_adj<-p.adjust(DATA_markers$p_val)

genesSel<-c(genesSel,
            rownames(DATA_markers[order(DATA_markers$p_val_adj,-DATA_markers$avg_log2FC),])[1:10],
            rownames(DATA_markers[order(DATA_markers$p_val_adj,DATA_markers$avg_log2FC),])[1:10])

## 9----------
DATA_markers <- FindMarkers(object = DATA.tmp,ident.1 = 9,ident.2 = comp,only.pos = F)
## Guardar matriz (Falta)
DATA_markers$p_val_adj<-p.adjust(DATA_markers$p_val)

genesSel<-c(genesSel,
            rownames(DATA_markers[order(DATA_markers$p_val_adj,-DATA_markers$avg_log2FC),])[1:10],
            rownames(DATA_markers[order(DATA_markers$p_val_adj,DATA_markers$avg_log2FC),])[1:10])

## 16----------
DATA_markers <- FindMarkers(object = DATA.tmp,ident.1 = 16,ident.2 = comp,only.pos = F)
## Guardar matriz (Falta)
DATA_markers$p_val_adj<-p.adjust(DATA_markers$p_val)

genesSel<-c(genesSel,
            rownames(DATA_markers[order(DATA_markers$p_val_adj,-DATA_markers$avg_log2FC),])[1:10],
            rownames(DATA_markers[order(DATA_markers$p_val_adj,DATA_markers$avg_log2FC),])[1:10])

## 19----------
DATA_markers <- FindMarkers(object = DATA.tmp,ident.1 = 19,ident.2 = comp,only.pos = F)
## Guardar matriz (Falta)
DATA_markers$p_val_adj<-p.adjust(DATA_markers$p_val)

genesSel<-c(genesSel,
            rownames(DATA_markers[order(DATA_markers$p_val_adj,-DATA_markers$avg_log2FC),])[1:10],
            rownames(DATA_markers[order(DATA_markers$p_val_adj,DATA_markers$avg_log2FC),])[1:10])

## all abc----------
DATA_markers <- FindMarkers(object = DATA.tmp,ident.1 = c(7,8,9,16,19),ident.2 = comp,only.pos = F)
## Guardar matriz (Falta)
DATA_markers$p_val_adj<-p.adjust(DATA_markers$p_val)

genesSel<-c(genesSel,
            rownames(DATA_markers[order(DATA_markers$p_val_adj,-DATA_markers$avg_log2FC),])[1:10],
            rownames(DATA_markers[order(DATA_markers$p_val_adj,DATA_markers$avg_log2FC),])[1:10])







cat("\nPlotting heatmap of Cluster Marker genes ...\n")
top<-10
DATA_markers %>% group_by(cluster) %>% top_n(top,avg_log2FC) -> top10

DATA.tmp<-subset(DATA.i,subset = seurat_clusters %in% as.character(c(7,8,9,16,19)))
genes<-top10$gene
clusters<-as.character(c(7,8,9,16,19))

m <- DATA.tmp@assays$RNA@counts
m<-m[genes,rownames(DATA.tmp@meta.data)[DATA.tmp@meta.data$seurat_clusters %in% clusters]]

clusterpos<-as.character(DATA.tmp@meta.data[colnames(m),"seurat_clusters"])

mplot<-matrix(data=0,ncol=length(clusters),nrow=length(genes))
colnames(mplot)<-clusters
rownames(mplot)<-genes

for(i in 1:length(genes)){
  for(j in 1:length(clusters)){
    
    tmp<-m[i,clusterpos==clusters[j]]
    mplot[i,j]<-mean(as.numeric(tmp),na.r=T)
  }
}

colors<-brewer.pal(n = 9, name = 'Blues')
pheatmap(mplot,scale="row",cluster_rows = F,cluster_cols = T,
         color=c("white",colors),
         clustering_method = "complete",fontsize=6,border_color="#303030")

#invisible(dev.off())

## plot
m <- DATA.tmp@assays$RNA@counts
genesall<-unique(genesSel)
m<-m[genesall,]

clusterpos<-as.character(DATA.tmp@meta.data[colnames(m),"seurat_clusters"])

mplotall<-matrix(data=0,ncol=length(clusters),nrow=length(genesall))
colnames(mplotall)<-clusters
rownames(mplotall)<-genesall

for(i in 1:length(genesall)){
  for(j in 1:length(clusters)){
    
    tmp<-m[i,clusterpos==clusters[j]]
    mplotall[i,j]<-mean(as.numeric(tmp),na.r=T)
  }
}

colors<-brewer.pal(n = 9, name = 'Blues')
pheatmap(mplotall,scale="row",cluster_rows = T,cluster_cols = T,
         color=c("white",colors),
         clustering_method = "complete",fontsize=6,border_color="#303030")











##--------
DATA_markers<- getDEG(DATA.i,
                      selectClusters=c("8","9","16","13","19"),
                      nameFile=paste0(opt$Results,"/Cluster_marker_genes_ABC2.csv"))

cat("\nPlotting heatmap of Cluster Marker genes ...\n")
DATA_markers %>% group_by(cluster) %>% top_n(top,avg_log2FC) -> top10
DATA.tmp<-subset(DATA.i,subset = seurat_clusters %in% c("8","9","16","13","19"))

#png(filename = paste0("heatmap_top_",sep="",top),width = 1500,height = 1500,res = 300)
DoHeatmap(object = DATA.tmp, features = as.character(unique(top10$gene)))
#invisible(dev.off())

###################





PlotMarkers(ObjetoSC=DATA.i,     ## Single cell R Object
            name="abc",         ## name for output plots
            POS_markers=c("Tbx21","Itgax","Itgam","Ccl5"),  ## vector of positive markers
            NEG_markers=NULL,  ## vector of negative markers
            path="/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK1/Bcells/Results/BCells",         ## Path to save the plots
            sizelabel=2,  ## Size for cluster labels
            gate=NULL)   ## vector with -x,x,-y,y limits

PlotMarkers(ObjetoSC=DATA.i,     ## Single cell R Object
            name="Igg",         ## name for output plots
            POS_markers=c("Ifitm1"),  ## vector of positive markers
            NEG_markers=NULL,  ## vector of negative markers
            path="/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK1/Bcells/Results/BCells",         ## Path to save the plots
            sizelabel=2,  ## Size for cluster labels
            gate=NULL)   ## vector with -x,x,-y,y limits

## Mediar,
## 





