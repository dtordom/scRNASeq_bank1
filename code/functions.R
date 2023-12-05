################################################################################
#' SingleCell Analysis functions
#' @Rversion R version 4.0.4
#'
################################################################################
## functions.R

##··············································································
## Ordered color palette to clusters
PutColor<-function(name.variables){
  
  palette.colors<-c("#a9182d","#f01515","#f1562a","#df7656","#fc7d50",
                    "#f29a42","#f3b465","#f6c4a7","#fba996","#ff8fad",
                    "#d98293","#b486cb","#8e40b4","#1a1aaa","#43639d",
                    "#7c67f8","#7c8bf8","#0099b2","#00bff9","#a6b8e4",
                    "#bbbdaa","#bace9d","#c0b878","#92ad66","#769544",
                    "#4fbb4b","#009b2e","#009000","#3a5d3a","#717171",
                    "#a4a4a4","#b89357","#c88d07","#c85807")
  
  if(length(name.variables)>length(palette.colors)){
    print("Error: length(name.variables)>length(palette.colors)")
  }else{
    if(length(name.variables)<=(length(palette.colors)/2)){
      palette.colors<-palette.colors[(1:length(palette.colors))%%2 !=0]
    }
    palette.colors<-palette.colors[1:length(name.variables)]
    
    names(palette.colors)<-name.variables
    
    return(palette.colors)
  }
  
}


##··············································································
## AnnotateCells
AnnotateCells<-function(seurat,
                        geneset,
                        wd=1300,
                        ht=900,
                        cor.diff=0.15,
                        cor.min=0,
                        pt.size=0.00001,
                        metadata.name,
                        path=paste0(opt$Results,sep="","/MainCells"),
                        ClustersColors){

  if (!file.exists(path)){dir.create(path)} 
  for(i in 1:length(geneset)){geneset[[i]]<-str_to_title(geneset[[i]])}
  
  cat("\nCreating cell identity matrix\n")
  cell_ident <- unique(unlist(geneset))
  cell_ident <- lapply(geneset,function(x){
    ifelse(cell_ident %in% x,1,ifelse(paste0("-",cell_ident) %in% x,-1,0))})
  cell_ident <- as.data.frame(cell_ident)
  rownames(cell_ident) <- casefold(unique(unlist(geneset)))
  
  cat("\nSelecting detected genes ...\n")
  sel <- rownames(seurat@assays[[opt$assay]]@data)[casefold(rownames(seurat@assays[[opt$assay]]@data)) %in% rownames(cell_ident)]
  cell_ident <- cell_ident[ casefold(sel),]
  
  cat("\nPredicting cell type by correlation method. Computing correlations ... ...\n")
  cors <- apply(seurat@assays[[opt$assay]]@data[sel,],2,function(x) cor(x,cell_ident))
  cors[is.na(cors)] <- -1
  rownames(cors) <- colnames(cell_ident)
  print(cors[,1:5])
  
  cat("\nPredicting cell types ...\n")
  pred <- unlist(apply(cors,2,function(x){
    tmp<-x[order(x,decreasing = T)]
    if(tmp[1]>=cor.min && (tmp[1]-cor.diff)>=tmp[2]){
      return(names(tmp[1]))
    }else{
      return(NA)
    }
    })
  )
  print(table(pred))
  seurat <- AddMetaData(seurat,metadata = pred,col.name = metadata.name)
  
  tiff(filename = paste0(path,"/",sep="",metadata.name,"_UMAP_tags.tiff"),width = wd,height = ht,res = 300)
  p.temp <- DimPlot(seurat,dims = 1:2,reduction = "umap",group.by = metadata.name,
                    pt.size = pt.size,label = F,cols = ClustersColors)
  plot(p.temp)
  invisible(dev.off())
  
  ##--------
  ## Plots of each cell type
  col_scale <- c("grey85","navy")
  myCorGrad <- colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(9)
  
  tiff(filename = paste0(path,"/",metadata.name,"_UMAP_cell_preds.tiff"),width = 450*4,height = 400*ceiling(nrow(cors) / 4),res = 150)
  par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(nrow(cors) / 4),4))
  myCorGrad <- paste0(colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(10))
  for(j in rownames(cors)){
    lim <- max(as.numeric(cors[j,]),na.rm = T)
    temp <- ((cors[j,]) - 0) / ( max(lim,0.6) + 0)
    temp[is.na(temp)] <- min(temp,na.rm = T)
    temp <- round((temp)*9)+1
    temp[temp <= 1] <- 1
    o <- order(temp)
    plot(seurat@reductions[["umap"]]@cell.embeddings[o,],pch=20,cex=pt.size, line=0.5, col=myCorGrad[ temp[o] ], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to ",j))
    image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
  }
  dev.off()
  
  return(seurat)
  
}


##··············································································
markerDots<-function(seurat,
                     genes,
                     wd=1000,
                     ht=680,
                     path){
  
  require('pheatmap')
  metadata<-seurat@meta.data
  metadata$cell_id<-rownames(metadata)
  
  genes<-str_to_title(genes)
  genes<-intersect(genes,rownames(seurat@assays$RNA@data))
  expression_matrix <- seurat@assays$RNA@data[genes,,drop=F]
  
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
  tmp<-matrix(ncol=length(clusters),nrow=length(genes))
  colnames(tmp)<-clusters; rownames(tmp)<-genes
  
  for(i in 1:length(genes)){
    for(j in clusters){
      punt<-j+1
      tmp[i,punt]<-m[m$seurat_clusters==clusters[punt] & m$gene==genes[i],]$percentage
    }
  }
  cl<-pheatmap(tmp,cluster_rows = F,scale="none")
  cl<-cl$tree_col$order

  expression_vals$seurat_clusters<-factor(x = expression_vals$seurat_clusters,
                                          levels=colnames(tmp)[cl])
  
  tiff(filename = path,width = wd,height = ht,res = 300)
  gg4 <- ggplot(expression_vals, aes(x = seurat_clusters, y = gene))+
    geom_point(aes(fill = median_expression, size = percentage),shape = 21,stroke=0.3)+
    theme_classic()+coord_flip()+
    theme(axis.title = element_blank(),
          legend.key.size = unit(0.2,"cm"),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5),
          axis.text.x = element_text(size = 5.8,angle = 90,hjust = 1),
          axis.text.y = element_text(size = 5.8,angle = 0,hjust = 1),
          axis.ticks.x = element_line(colour = 'black', size = 0.3),
          axis.ticks.y = element_line(colour = 'black', size = 0.3),
          axis.line = element_line(colour = 'black', size = 0.3))+
    scale_fill_gradient(low = "white", high = "darkred")+
    scale_size(range=c(0.5,3.3))+
    labs(fill = "Expression", size = "Percentage of\npositive cells")
  plot(gg4)
  invisible(dev.off())
  
  return(tmp)
  
}

##··············································································
## Plot proportion of cells in clusters, genotypes...
ProportionPlot<-function(seuratObj,
                         wd=2200,
                         ht=1300,
                         path=setwd(),
                         idColors,
                         column_split, ## clusters, genotypes...
                         column_tags){ ## cell assigment column
  
  seuratObj <- SetIdent(seuratObj, value = factor(as.character(seuratObj@meta.data[,column_split])))
  #x<-colnames(seuratObj@meta.data)[grepl(pattern = column_tags, x = colnames(seuratObj@meta.data))]
  pred_factor <- as.factor(seuratObj@meta.data[,column_tags])
  cell_type_levels <- levels(pred_factor)
  
  proportion <- as.data.frame(lapply(levels(seuratObj@active.ident),function(x){c(unname(table(pred_factor[seuratObj@active.ident==x]))) [1:length(cell_type_levels)]} ))
  proportion[is.na(proportion)] <- 0
  rownames(proportion) <- cell_type_levels
  colnames(proportion) <- levels(seuratObj@active.ident)
  
  cl_order <- order(as.numeric(levels(seuratObj@active.ident)))
  sa <- cbind(stack(as.data.frame(t(t(proportion[,cl_order])/colSums(proportion[,cl_order])) )), rep(rownames(proportion),ncol(proportion)) )
  colnames(sa) <- c("prop","splitby","groups")
  
  tiff(filename = paste0(path,"/",column_split,"_",column_tags,"_cells_barplot.tiff"),width = wd,height = ht,res = 300)
  plot1 <- ggplot(data=sa, aes(x=splitby, y=prop, fill=groups)) + geom_col(colour="black",size=0.3)+
    scale_fill_manual(values=idColors)+theme(axis.text.x = element_text(angle=90))
  plot(plot1)
  invisible(dev.off())
  
  return()
}

##··············································································
## Ver expresion y ncells que expresan por muestra
seeExpression<-function(dataObj,
                        metadata,
                        splitBy="orig.ident",
                        gen="genotype",
                        expBase=0,
                        genes="cd19"){
  
  genes<-str_to_title(genes)
  genes<-intersect(genes,rownames(dataObj@assays$RNA@counts))
  
  vars<-unique(dataObj@meta.data[,splitBy])
  res<-data.frame("id"=vars,
                  "exp"=rep(0,length(vars)),
                  "nCells"=rep(0,length(vars)),
                  "nTot"=rep(0,length(vars)),
                  "propCells"=rep(0,length(vars)))
  
  for(i in 1:length(vars)){
    cellSel<-unique(rownames(dataObj@meta.data)[dataObj@meta.data[,splitBy]==vars[i]])
    tmp<-dataObj@assays$RNA@counts[genes,cellSel]
    if(is.null(nrow(tmp))){
      res[i,"exp"]<-mean(tmp)
    }else{
      tmp<-apply(tmp,2,mean)
      res[i,"exp"]<-mean(tmp)
    }
    res[i,"nCells"]<-length(tmp[which(tmp>expBase)])
    res[i,"nTot"]<-length(tmp)
    res[i,"propCells"]<-length(tmp[which(tmp>expBase)])/length(tmp)
    
  }
  if(!is.null(gen)){
    addData<-metadata[res$id,gen]
    colnames(addData)<-gen
    res<-cbind(res,addData)
  }
  #res$class<-metadata[res$id,gen]
  #res<-res[order(res$class),]
  
  return(res)
  
}



##-------
## DEG by cluster
getDEG<-function(DATA.tmp,selectClusters,filterAdj=T,nameFile,only.pos =FALSE){
  
  DefaultAssay(DATA.tmp)<-"RNA"
  DATA.tmp <- SetIdent(DATA.tmp, value = DATA.tmp@meta.data$seurat_clusters) 
  
  DATA.tmp<-subset(DATA.tmp,subset = seurat_clusters %in% selectClusters)
  DATA.tmp@meta.data <- DATA.tmp@meta.data[colnames(DATA.tmp)[DATA.tmp@meta.data$seurat_clusters %in% selectClusters ],]
  
  DATA_markers <- FindAllMarkers(object = DATA.tmp, assay = opt$assay, only.pos =only.pos,
                                 min.pct = 0.1, min.diff.pct = 0.05,
                                 max.cells.per.ident = Inf,print.bar = T,
                                 do.print = T,return.thresh = 0.05,
                                 test.use = "t")
  
  if(filterAdj){
    sel<-ifelse(DATA_markers$p_val_adj<=0.05,T,F)
    DATA_markers<-DATA_markers[sel,]
  }
  write.table(DATA_markers,file = nameFile,row.names = F,sep="\t")
  return(DATA_markers)
}


##-------
## Create breaks and custom color ranges
Custom.colorRange<-function(hits,colorseq,size=1000){
  
  breaks=seq(min(hits),max(hits),length.out = size)
  cols<-NULL
  i=1
  while(i<length(hits)){
    sel<-ifelse(breaks>=hits[i] & breaks<hits[i+1],T,F)
    n<-table(sel)["TRUE"]
    if((i+1)==(length(hits))){
      n<-n+1
    }
    tmp<-colorRampPalette(colors=c(colorseq[i],colorseq[i+1]))(n = n)
    cols<-c(cols,tmp)
    i<-i+1
  }
  vals<-list(breaks,cols); names(vals)<-c("breaks","colors")
  return(vals)
}


##-------
## DEG by genotype in one cluster
getDEGGen<-function(DATA.tmp,cluster,adj=TRUE,comb=NULL){
  
  DefaultAssay(DATA.tmp)<-"RNA"
  DATA.tmp <- SetIdent(DATA.tmp, value = DATA.tmp@meta.data$seurat_clusters) 
  
  DATA.tmp<-subset(DATA.tmp,subset = seurat_clusters %in% cluster)
  DATA.tmp@meta.data <- DATA.tmp@meta.data[colnames(DATA.tmp)[DATA.tmp@meta.data$seurat_clusters %in% cluster],]
  
  DATA.tmp <- SetIdent(DATA.tmp, value = DATA.tmp@meta.data$genotype) 
  
  if(is.null(comb)==T){
    comb<-data.frame("comb1"=c("BANK1","WT"),
                     "comb2"=c("TLR","BANK1"),
                     "comb3"=c("TLR","WT"))
  }
  
  RESULTS<-NULL
  for(i in 1:ncol(comb)){
    
    DATA.tmp.i<-subset(DATA.tmp,subset = genotype %in% as.character(comb[,i]))
    DATA.tmp.i@meta.data <- DATA.tmp.i@meta.data[colnames(DATA.tmp.i)[DATA.tmp.i@meta.data$genotype %in% as.character(comb[,i])],]
    
    if(adj){
      
      markers <- FindAllMarkers(object = DATA.tmp.i, assay = opt$assay, only.pos =FALSE,
                                max.cells.per.ident = Inf,print.bar = T,
                                do.print = T,return.thresh = 0.05,
                                test.use = "t")
      
      if(nrow(markers)>0){
        markers<-markers[ifelse(markers$cluster==comb[1,i],T,F),]
        markers<-markers[ifelse(markers$p_val_adj<=0.05,T,F),c("gene","avg_log2FC","p_val_adj")]
        if(nrow(markers)>0){
          markers<-markers[order(markers$avg_log2FC,decreasing=T),]
          markers$comparation<-rep(paste0(comb[1,i],sep="","_vs_",comb[2,i]),nrow(markers))
          rownames(markers)<-NULL
          RESULTS<-rbind(RESULTS,markers)
        }
      }
      RESULTS<-RESULTS[ifelse(abs(RESULTS$avg_log2FC)>=1,T,F),]
    }else{
      
      markers <- FindAllMarkers(object = DATA.tmp.i, assay = opt$assay, only.pos =FALSE,
                                max.cells.per.ident = Inf,print.bar = T,
                                do.print = T,return.thresh = 1,
                                test.use = "t")
      
      if(nrow(markers)>0){
        markers<-markers[ifelse(markers$cluster==comb[1,i],T,F),]
        markers<-markers[order(markers$avg_log2FC,decreasing=T),]
        markers$comparation<-rep(paste0(comb[1,i],sep="","_vs_",comb[2,i]),nrow(markers))
        rownames(markers)<-NULL
        RESULTS<-rbind(RESULTS,markers)
      }
      
      #RESULTS<-RESULTS[ifelse(abs(RESULTS$avg_log2FC)>=1,T,F),]
    }
    
  }
  
  #write.table(DATA_markers,file = nameFile,row.names = F,sep="\t")
  return(RESULTS)
}





