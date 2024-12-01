
rawdataqc <- function(rawdatadir,projectname,tissues) {
    library(Seurat)
    library(SeuratData)
    library(DoubletFinder)
    library(patchwork)
    library(ggplot2)
    library(harmony)
    if (!dir.exists("../output/")) dir.create("../output/")
    if (!dir.exists(paste0("../output/01.QC"))) dir.create(paste0("../output/01.QC"))
    if (!dir.exists(paste0("../output/02.cluster"))) dir.create(paste0("../output/02.cluster"))
    samples=list.files(paste0(rawdatadir,tissues))
    ####依次读入样本
    sceList = lapply(samples,function(pro){
        folder=file.path(paste0(rawdatadir,tissues),pro)
        CreateSeuratObject(counts = Read10X(folder),
                           project = pro )
        })
    seuratdata <- merge(sceList[[1]],
                        y=c(sceList[[2]],sceList[[3]],sceList[[4]],
                            sceList[[5]],sceList[[6]],sceList[[7]],
                            sceList[[8]],sceList[[9]],sceList[[10]],
                            sceList[[11]]),
                        add.cell.ids=c(paste0(rep("DLE",4),seq(1,4),"_B"),
                                       paste0(rep("HC",3),seq(1,3),"_B"),
                                       paste0(rep("SLE",4),seq(1,4),"_B")),
                        project="B1_BLOOD")
    #不整合
    # run standard anlaysis workflow
    seuratdata <- NormalizeData(seuratdata)
    seuratdata <- FindVariableFeatures(seuratdata)
    seuratdata <- ScaleData(seuratdata)
    seuratdata <- RunPCA(seuratdata)
    seuratdata <- RunHarmony(seuratdata,"orig.ident")#
    
    seuratdata <- RunUMAP(seuratdata, dims = 1:30, reduction = "harmony", reduction.name = "umap.unintegrated")
    #scmeta <- seuratdata@meta.data
    umap_ui <- DimPlot(seuratdata, reduction = "umap.unintegrated", group.by = c("orig.ident"))
    
    #执行整合
    seuratdata <- IntegrateLayers(object = seuratdata,assay = ,
                                  method = CCAIntegration, 
                                  orig.reduction = "harmony", 
                                  new.reduction = "integrated.cca",
                                  verbose = FALSE)
    
    # re-join layers after integration
    seuratdata[["RNA"]] <- JoinLayers(seuratdata[["RNA"]])
    
    seuratdata <- FindNeighbors(seuratdata, reduction = "integrated.cca", dims = 1:30)
    seuratdata <- FindClusters(seuratdata, resolution = 0.5)
    seuratdata <- RunUMAP(seuratdata, dims = 1:30, reduction = "integrated.cca")
    
    umap_i <- DimPlot(seuratdata, reduction = "umap", group.by = c("orig.ident"))
    p <- umap_ui+umap_i  
    p  
    ggsave("../output/01.QC/integrated.pdf",plot = p,width = 12,height = 5)
    # 去除双细胞
    ## pK Identification (no ground-truth)
    sweep.res.list<- paramSweep(seuratdata, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    sweep.stats[order(sweep.stats$BCreal),]
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- seuratdata@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(seuratdata@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seuratdata <- doubletFinder(seuratdata, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, 
                                reuse.pANN = FALSE, sct = FALSE)
    colnames(seuratdata@meta.data)[7] <- "DF.classifications"
    seuratdata <- subset(seuratdata,subset = `DF.classifications`=="Singlet")
    
    seuratdata$condition <- gsub("\\d_\\w$","",  seuratdata$orig.ident)
    #rm("counts")
    ####生成seurat对象
    save(seuratdata,file=paste0("../output/01.QC/2.",projectname,".raw.Rdata") )
}