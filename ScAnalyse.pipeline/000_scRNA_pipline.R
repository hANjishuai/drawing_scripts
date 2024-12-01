#### scRNA analysis pipline
rm(list = ls())
gc()
write( "times", file = "../output/times.txt")
sink(file = "../output/times.txt",append = T);
Sys.time();sink()
options(future.globals.maxSize = 50000 * 1024^2)  # 设置为 1000 MB

#### load constant for scRNA analysis ####
.libPaths("~/R/4.4.1/library/")
source("./env_setting/000_constant_scRNA.R")
constant_scRNA()
load("../output/constant.RData")
#source("./scRNA_utils.R")
####处理rawdata
source("./stepsOfpipeline/00_rawdatamerge.R")
rawdataqc(rawdatadir <- rawdatadir,
          projectname <- projectname)
####质控
source("./stepsOfpipeline/01_scRNA.QC.R")
scRNA.qc(projectname = projectname, species =  species)
####聚类
source("./stepsOfpipeline/02_0_cluster.R")
cluster_marker(projectname =  projectname)
source("./stepsOfpipeline/03_0_cellmarkerplot.R")#10.2
manual_marker_plot(projectname =projectname )
#### manual check outfiles identity final celltype ####
#### 细胞注释 
   # singleR 注释
   source("./stepsOfpipeline/04_0_cellanno.R")
   cellanno(projectname = projectname,specie = "mouse" )
   cellanno.check(projectname=projectname,
                  newcluster = "../output/03.celltype/0.singleR_cellplot.csv")
   # 保存cluster marker 
   source("./stepsOfpipeline/04_1_clustermarkers_table.R")
   clustermaker_table()
   source("./stepsOfpipeline/04_2_totalcell_statistics.R")
   totalcell_statistics()
   source("./stepsOfpipeline/04_3_totalcell_statistics_ratio.R")
   totalcell_statistics_ratio()
#### get sub cellcluster
source("./stepsOfpipeline/05_0_subsetcell.R")
subcell(projectname = projectname) 

#### identification of DEGs
source("./06_0_celldeg.R")
celldeg(cell = NULL,compares = compares,specie = "mouse")

#### 下采样，获取小数据集
source("./07_0_sample_seuratdata.R")
sample_seuratdata()

#### 轨迹分析 ####
source("./08_0_monocle2.R")
#monocle2_self(cell = "Microglia",T0name = "Microglia_Ramified_Ecscr")
#monocle2_self(cell = "Bcell",T0name = "Bcell_Myl4")

#### 通讯分析 ####
source("./09_0_cellchat_single_subcell.R")
cellchat_single(projectname=projectname, species="mouse")

sink(file = "../output/times.txt",append = T);
Sys.time();sink()

#### 亚群分析 ####
source("./stepsOfpipeline/10_0_subset.R")
subcell.qc(species = "human")

#### 亚群注释 ####
source("./stepsOfpipeline/10_1_subcellanno.R")
subcellanno("../output/04.subcell/03.annotation/manual_level2.csv")


