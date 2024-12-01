#species="mouse"
cellchat_single <- function(projectname,species="mouse") {
  library(tidyverse)
  library(ggplot2)
  library(CellChat)
  library(patchwork)
  library(Seurat)
  #library(future)
  options(stringsAsFactors = FALSE)
  
if (!dir.exists(paste0("../output/08_cell_com/"))) dir.create(paste0("../output/08_cell_com/"))

  my.sapply <- if (future::nbrOfWorkers() == 1) {
      pbapply::pbsapply
  } else {
      future.apply::future_sapply
  }
  
####Part I: Data input & processing and initialization of CellChat object####

# Load data
loadmerge=load(paste0("../output/01.QC/8.",projectname,".celltype.sample.Rdata"))
seuratdata <- downsampled_obj
Idents(seuratdata) <- seuratdata$Tissue

#tissue <- "GZ01-14dpi"
for (tissue in unique(seuratdata$Tissue[-1])) {
  #### Part I: Data input & processing and initialization of CellChat object
  if (!dir.exists(paste0("../output/08_cell_com/",tissue))) dir.create(paste0("../output/08_cell_com/",tissue))
  savedir <- paste0("../output/08_cell_com/",tissue,"/")
  subtissue <- subset(seuratdata, idents = tissue)
  # Prepare input data for CelChat analysis
  data.input = subtissue[["RNA"]]$counts # counts data matrix
  index <- duplicated(rownames(data.input))
  data.input <- data.input[!index,]
  data.input <- normalizeData(data.raw = data.input) # normalized data matrix
  meta = subtissue@meta.data # a dataframe with rownames containing cell mata data
  rm(subtissue)
  gc()
  # Create a CellChat object
  print("Create a CellChat object")
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "manual_L3")
  cellchat <- setIdent(cellchat, ident.use = "manual_L3") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  # Set the ligand-receptor interaction database
  if(species=="mouse") CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  if(species=="human") CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)
  print("CellChatDB")
  
  # use a subset of CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use

  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  #plan(multisession) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  #cellchat <- ProjectData(cellchat, PPI.mouse)

  #### Part II: Inference of cell-cell communication network
  cellchat@idents <- factor(cellchat@idents,levels = unique(cellchat@idents))
  cellchat <- computeCommunProb(cellchat)
  
  #> triMean is used for calculating the average gene expression per cell group. 
  #> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
  #> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat)
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  # visualize the aggregated cell-cell communication network.
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  
  pdf(paste0(savedir,"0_cc_com_network.pdf"),width = 12,height = 12)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE, title.name = "Interaction weights/strength")
  dev.off()
  pdf(paste0(savedir,"1_cc_com_network_single.pdf"),width = 30,height = 50)
  mat <- cellchat@net$weight
  par(mfrow = c(8,5), xpd=TRUE)
  #all.n <- rownames(mat)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
   # sr <- rownames(mat)[1]
   # tr <- 
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE,
                     #sources.use = ,targets.use = ,
                     edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  ####Part IV: Systems analysis of cell-cell communication network
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  
  
  # Part III: Visualization of cell-cell communication network
  if (!dir.exists(paste0(savedir,"pathway"))) dir.create(paste0(savedir,"pathway"))
  savedir1 <- paste0(savedir,"pathway/")
  pathways.show <- cellchat@netP$pathways[2]
  for (pathways.show in cellchat@netP$pathways) {
    pdf(paste0(savedir1,pathways.show,"_network.pdf"),width =10,height = 10)
    # Hierarchy plot
    # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
    
    vertex.receiver = seq(1,length(cellchat@net$prob %>% colnames)) # a numeric vector. 
    netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
    # Circle plot
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    #> Heatmap
    nh <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    plot(nh)
    #> Do heatmap based on a single object
    dev.off() 
    vertex.receiver <- seq(1,20)

    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
      gg <- netAnalysis_contribution(cellchat, signaling = pathways.show)
      ggsave(filename=paste0(savedir1,pathways.show,"_L-R_contribution.pdf"), plot=gg,
             width = 3, height = 2, units = 'in', dpi = 300)
      
      #### Plot the signaling gene expression distribution using violin/dot plot
      pge <- plotGeneExpression(cellchat, signaling = "CXCL")
      ggsave(filename=paste0(savedir1,pathways.show,"_GeneExpression.pdf"), plot=pge,
             width = 12, height = 8, units = 'in', dpi = 300)
      # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
      pdf(paste0(savedir1,pathways.show,"_signalingRole_network.pdf"),width = 18,height = 8)
      netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, 
                                        height = 2.5, font.size = 10)
      dev.off()
      # Visualize the dominant senders (sources) and receivers (targets) in a 2D space
      # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      gg1 <- netAnalysis_signalingRole_scatter(cellchat)
      #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      # Signaling role analysis on the cell-cell communication networks of interest
      gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
      #> Signaling role analysis on the cell-cell communication network from user's input
      gg <- gg1 + gg2
      ggsave(filename=paste0(savedir1,pathways.show,"_sender_rece_2D.pdf"), plot=gg,
             width = 16, height = 8, units = 'in', dpi = 300)
      # Identify signals contributing most to outgoing or incoming signaling of certain cell groups
      
      # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 24,width = 20)
      ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 24,width = 20)
      hh <- ht1 + ht2
      
      pdf(paste0(savedir1,pathways.show,"_signalingRole_heatmap.pdf"),width = 24,height = 16)
      plot(hh)
      dev.off()
      
   
  }
  # Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
  library(NMF)
  
  saveRDS(cellchat, file = paste0(savedir,"cellchat_LS.rds"))
}
  
}

