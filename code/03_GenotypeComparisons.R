##############################
## SingleCell Analysis
## R version R version 4.0.4
##
##############################
## 10_GenotypeComparisons

load("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds_bcells/clustered_seurat_object.RData")
source(opt$resourcesPath)
library("Seurat")
library("pheatmap")
library("RColorBrewer")
library("dplyr")

##------------------------------------------------------ STEP 1
## Finding differentially expressed genes (between phenotypes)

## POR AQUI

clusters<-as.character(names(table(DATA.i@meta.data$seurat_clusters)))


comb<-data.frame("comb1"=c("TLR","BANK1"))

byGEN<-list()
for(g in 1:length(clusters)){

  res<-getDEGGen(DATA.tmp=DATA.i,
                 adj = F, comb = comb,
                 cluster=clusters[g])
  
  res<-res[ifelse(res$p_val_adj<=0.05,T,F),]
  byGEN[[g]]<-res
  names(byGEN)[length(byGEN)]<-paste0("cluster",clusters[g],sep="")
}

# res<-getDEGGen(DATA.tmp=DATA.i,
#               cluster="7",adj = F)


save.image("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Rds_bcells/byGEN.RData")

## Save file
setwd("/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells/Bcells/compbyGen")

tabla<-NULL
for(i in 1:length(byGEN)){
  tmp<-byGEN[[i]]
  if(nrow(tmp)>0){
    tmp$cluster<-rep(names(byGEN)[i],nrow(tmp))
    tmp<-tmp[,c("cluster","gene","p_val","p_val_adj","avg_log2FC","pct.1","pct.2")]
    colnames(tmp)<-c("Cluster","Gene","Pvalue","Pvalue_Adj","avg_log2FC(tlr7_vs_bank1)","Pct_tlr7","Pct_bank1")
    tabla<-rbind(tabla,tmp)
  }
}

write.table(tabla,"/home/daniel/Desktop/WORK/scRNASeq/MorellMouse/Results_bcells/Bcells/GenotypeDEG.txt",quote = F,row.names = F,sep="\t")


##------------------------
## Sacar Figuras

genes<-unique(tabla$Gene)

FC<-data.frame(matrix(data=0,ncol=length(unique(tabla$Cluster)),nrow=length(genes)))
colnames(FC)<-unique(tabla$Cluster)
rownames(FC)<-genes

for(i in 1:length(genes)){
  for(j in 1:nrow(tabla)){
    if(tabla$Gene[j]==genes[i]){
      FC[genes[i],tabla$Cluster[j]]<-tabla$`avg_log2FC(tlr7_vs_bank1)`[j]
    }
  }
}


pheatmap(t(FC),cluster_rows = T,cluster_cols = T,
         show_rownames = F,show_colnames = T,
         breaks = seq(-1,1,length.out = 100),
         gaps_col = c(1,2,3,4,5,6,7),
         #annotation_row = ann,
         fontsize = 8,
         color = colorRampPalette(c("#4393c3","white", "#d6604d"))(100),
         border_color = "black",
         scale = "none")


## POR AQUIII

clusters<-unique(tabla$Cluster)

m<-data.frame("Cluster"=clusters,
              "Value"=as.numeric(table(tabla$Cluster)[clusters]))
m<-m[order(m$Value,decreasing=T),]
m$Cluster<-factor(x = m$Cluster,levels = m$Cluster)

p1<-ggplot(m,aes(x=Cluster,y=Value))+geom_bar(stat="identity",color="black")+
  theme_bw()+theme(panel.grid = element_blank())

CommonGenes<-NULL
for(i in 1:nrow(m)){
  
  genes<-tabla[tabla$Cluster==m$Cluster[i],"Gene"]
  NOgenes<-tabla[tabla$Cluster!=m$Cluster[i],"Gene"]
  uniqGenes<-length(genes)-length(intersect(genes,NOgenes))
  
  res<-NULL
  for(j in 1:nrow(m)){
    
    if(j == i){
      res<-c(res,uniqGenes)
    }else{
      genes2<-tabla[tabla$Cluster==m$Cluster[j],"Gene"]
      res<-c(res,length(intersect(genes,genes2)))
    }

  }
  CommonGenes<-rbind(CommonGenes,res)
}
rownames(CommonGenes)<-m$Cluster
colnames(CommonGenes)<-m$Cluster

Prop<-CommonGenes
for(i in 1:nrow(CommonGenes)){
  Prop[i,]<-(CommonGenes[i,] / apply(CommonGenes,2,sum))*100
}
Prop

m2<-melt(Prop)
colnames(m2)<-c("groups","splitby","prop")
m2$groups<-gsub("cluster","",m2$groups)
m2$splitby<-gsub("cluster","Cluster ",m2$splitby)

m2$splitby<-factor(m2$splitby,levels = as.character(unique(m2$splitby)))

m2$groups<-paste0("Cluster ",sep="",m2$groups)

p2 <- ggplot(data=m2, aes(x=splitby, y=prop, fill=groups)) + geom_col(colour="black",size=0.3)+
  scale_fill_manual(values=list("Cluster 0"="#fa8881","Cluster 1"="#e89614","Cluster 2"="#c69a0a","Cluster 3"="#9dba1c",
                                "Cluster 4"="#27b806","Cluster 5"="#1fbc7c","Cluster 6"="#02c0ab","Cluster 7"="#07bcda",
                                "Cluster 8"="#2cb4f5","Cluster 9"="#8b93ff","Cluster 10"="#d980ff","Cluster 11"="#f36bdb",
                                "Cluster 12"="#ce657d"))+theme(axis.text.x = element_text(angle=90))
plot(p2)
dev.off()

library(ggpubr)

tiff(filename = paste0(opt$results,"/GenDEG.tiff"),width = 800,height = 1000,res = 300)
plot(ggarrange(p1,p2,nrow = 2,ncol=1,legend = "none",heights = c(0.5, 1)))
invisible(dev.off())


## FALTA: UMAP diferenciando TLR7 / BANK1 / WT





############################################# old


##
library(stringr)

## Significant genes
genes<-NULL
# tmp = res

for(i in 1:length(byGEN)){
  tmp<-byGEN[[i]]
  tmp<-tmp[ifelse(tmp$comparation=="TLR_vs_BANK1",T,F),]
  tmp<-tmp[tmp$p_val_adj<=0.05,]
  
  ## Remove Ribosomal genes
  tmp<-tmp[!str_detect(string = tmp$gene,pattern = "Rps"),]
  tmp<-tmp[!str_detect(string = tmp$gene,pattern = "Rpl"),]
  
  tmp<-tmp[order(tmp$avg_log2FC,decreasing = T),]
  
  # if(nrow(tmp)>0){
  #   tmpP<-tmp[ifelse(tmp$avg_log2FC>0,T,F),]
  #   if(nrow(tmpP)>0){
  #     tmpP<-tmpP[order(tmpP$avg_log2FC,decreasing=T),]
  #     tmpP<-tmpP$gene[1:15]; tmpP<-as.character(na.omit(tmpP))
  #     genes<-c(genes,tmpP)
  #   }
  #   tmpN<-tmp[ifelse(tmp$avg_log2FC<0,T,F),]
  #   if(nrow(tmpN)>0){
  #     tmpN<-tmpN[order(tmpN$avg_log2FC,decreasing=F),]
  #     tmpN<-tmpN$gene[1:15]; tmpN<-as.character(na.omit(tmpN))
  #     genes<-c(genes,tmpN)
  #   }
  # }
}
# genes<-unique(genes)


##-------------
## Get FoldChanges

comb<-data.frame("comb1"=c("TLR7","TLR7_Bank1-KO"))

EXP<-list()
for(g in 1:length(clusters)){
  
  res<-getDEGGen(DATA.tmp=DATA.i,
                 cluster=clusters[g],
                 adj = F,
                 comb=comb)
  EXP[[g]]<-res
  names(EXP)[length(EXP)]<-paste0("cluster",clusters[g],sep="")
}

M<-NULL
for(i in 1:length(EXP)){
  
  tmp<-EXP[[i]]
  tmp<-tmp[ifelse(tmp$comparation=="TLR7_vs_TLR7_Bank1-KO",T,F),]
  rownames(tmp)<-tmp$gene
  
  vect<-NULL
  for(geni in 1:length(genes)){
    x<-as.numeric(tmp[genes[geni],"avg_log2FC"])
    if(is.na(x)){
      x<-0
    }
    vect<-c(vect,x)
  }
  M<-cbind(M,vect)
  
}
colnames(M)<-names(byGEN)
rownames(M)<-genes

library("pheatmap")

hits<-c(-2.6,0,2.6)
colorseq<-c("#003366","white","#FFCC33")
bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)



setwd("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK1/Bcells/Results/BCells")
dev.off()

tiff(filename = paste0(opt$Results,"/BCells/TOPgenes.tiff"),width = 2000,height = 800,res = 300)
pheatmap(t(M),border_color = "black",cluster_rows = T,cluster_cols = T,
         breaks=bks$breaks,color=bks$colors,scale = "column",
         clustering_method = "complete",fontsize=5.5)
invisible(dev.off())



##------------------------------------------------------ STEP 2
## Functional analysis

setwd("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK1/Bcells/Results")
source("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK1/scripts/GC4programmatic.R")

PATHS<-list()
count<-1
for(i in 1:length(byGEN)){
  
  tmp<-byGEN[[i]]
  tmp<-tmp[ifelse(tmp$comparation=="TLR7_vs_TLR7_Bank1-KO",T,F),]
  
  ## Remove Ribosomal genes
  tmp<-tmp[!str_detect(string = tmp$gene,pattern = "Rps"),]
  tmp<-tmp[!str_detect(string = tmp$gene,pattern = "Rpl"),]
  
  up<-tmp[ifelse(tmp$avg_log2FC>0,T,F),]
  down<-tmp[ifelse(tmp$avg_log2FC<0,T,F),]
  UP<-NULL
  DOWN<-NULL
  
  if(length(up$gene)>1){
    UP <- launchAnalysis(organism = "Mus musculus",
                         inputType = "genes",
                         inputQuery = as.character(up$gene),
                         annotationsDBs = c("GO_BP","GO_CC"),
                         outputType = "dataframe",
                         inputCoannotation = "no",
                         inputName1 = "abcs",
                         universeScope = "annotated")
  }
  
  if(length(down$gene)>1){
    DOWN <- launchAnalysis(organism = "Mus musculus",
                           inputType = "genes",
                           inputQuery = as.character(down$gene),
                           annotationsDBs = c("GO_BP","GO_CC"),
                           outputType = "dataframe",
                           inputCoannotation = "no",
                           inputName1 = "abcs",
                           universeScope = "annotated")
  }

  if(is.null(UP)==FALSE){
    UP<-UP$`abcs-GO_BP`[,c("description","annotation_id","pval_adj")]
    UP<-UP[ifelse(UP$pval_adj<=0.05,T,F),]
    UP$dir<-rep("up",nrow(UP))
  }
  if(is.null(DOWN)==FALSE){
    DOWN<-DOWN$`abcs-GO_BP`[,c("description","annotation_id","pval_adj")]
    DOWN<-DOWN[ifelse(DOWN$pval_adj<=0.05,T,F),]
    DOWN$dir<-rep("down",nrow(DOWN))
  }
  
  res<-rbind(UP,DOWN)
  if(is.null(res)==FALSE){
    res$Cluster<-rep(names(byGEN)[i],nrow(res))
    PATHS[[count]]<-res
    count<-count+1
    names(PATHS)[length(PATHS)]<-names(byGEN)[i]
  }

}

## Save Table
TAB<-NULL
for(i in 1:length(PATHS)){
  TAB<-rbind(TAB,PATHS[[i]])
}
TAB<-TAB[,c("Cluster","dir","pval_adj","annotation_id","description")]
colnames(TAB)<-c("ClusterReference","Direction","PvalAdj","AnnotationID","Description")

write.table(TAB,"Functional_TLRvsKO.txt",
            quote = F,row.names = F,
            sep="\t")


##---------
## Get top paths

## Significant genes
paths<-NULL

for(i in 1:length(PATHS)){
  tmp<-PATHS[[i]]
  
  tmpup<-tmp[ifelse(tmp$dir=="up",T,F),]
  tmpdown<-tmp[ifelse(tmp$dir!="up",T,F),]
  
    if(nrow(tmpup)>0){
      tmpup<-tmpup[order(tmpup$pval_adj,decreasing=F),]
      tmpup<-tmpup$annotation_id[1:5]; tmpup<-as.character(na.omit(tmpup))
      paths<-c(paths,tmpup)
    }

    if(nrow(tmpdown)>0){
      tmpdown<-tmpdown[order(tmpdown$pval_adj,decreasing=F),]
      tmpdown<-tmpdown$annotation_id[1:5]; tmpdown<-as.character(na.omit(tmpdown))
      paths<-c(paths,tmpdown)
    }
}
paths<-unique(paths)


##-------------
## Pvals

M<-NULL
for(i in 1:length(PATHS)){
  
  tmp<-PATHS[[i]]
  tmp<-tmp[!duplicated(tmp$annotation_id),]
  rownames(tmp)<-tmp$annotation_id
  
  vect<-NULL
  for(geni in 1:length(paths)){
    x<-as.numeric(tmp[paths[geni],"pval_adj"])
    
    if(length(x)==0 | is.na(x)){
      x<-0
    }else{
      x<--log10(x)
      if(tmp[paths[geni],"dir"]!="up"){
        x<-x*-1
      }
    }
    vect<-c(vect,x)
  }
  M<-cbind(M,vect)
  
}
colnames(M)<-names(PATHS)
rownames(M)<-paths

library("pheatmap")

setwd("/home/daniel/Desktop/WORK/WORK/scRNASeq/BANK1/Bcells/Results/BCells")

hits<-c(-2.8,-2,-1.3,0,1.3,3,10.5)
colorseq<-c("#003366","#0066CC","lightblue","white","#fefdbc","#f6824c","#9e0142") #
bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)

tab<-TAB[!duplicated(TAB$AnnotationID),]
rownames(tab)<-tab$AnnotationID
tab<-tab[rownames(M),]

M2<-t(M)
colnames(M2)<-tab$Description

tiff(filename = paste0(opt$Results,"/BCells/TOPpaths.tiff"),width = 1750,height = 2100,res = 300)
pheatmap(t(M2),border_color = "black",cluster_rows = T,cluster_cols = T,
         breaks=bks$breaks,color=bks$colors,scale = "none",
         clustering_method = "complete",fontsize=5.5)
invisible(dev.off())




