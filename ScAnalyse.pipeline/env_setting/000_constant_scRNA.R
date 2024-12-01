#按照自己的文件结构及目录名，构建环境变量
constant_scRNA <- function(){
  rawdatadir = "../data/"
  projectname = "B1"
  species = "human"
  #cell="CD8"
  tissues="Blood"
  samples = list.files(paste0("../data/",tissues))
  #compares = c("LE VS HC")
  level_treat = c("HC","DLE","SLE")
  save(rawdatadir,projectname,species,samples,#compares,
       level_treat,
       file = "../output/constant.RData")
  }
