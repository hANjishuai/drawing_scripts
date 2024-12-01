convert_human_to_mosue <- function(gene_list){
    library(dplyr)
    if (file.exists("../db/HOM_MouseHumanSequence.rpt.txt")) {
        mouse_human_genes <- read.delim("../db/HOM_MouseHumanSequence.rpt.txt",
                                        sep = "\t")
    } else {
        mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")}
    
    humangenes <- mouse_human_genes[mouse_human_genes$Symbol %in% gene_list,]
    allgenes <- mouse_human_genes[mouse_human_genes$DB.Class.Key %in% humangenes$DB.Class.Key  ,]
    mousegenes <- allgenes[allgenes$Common.Organism.Name %in% "mouse, laboratory", ]
}

calTissueDist <- function(dat.tb,byPatient=F,colname.cluster="majorCluster",simulate.p.value=FALSE,
                          colname.patient="patient",colname.tissue="loc")
{
    if(byPatient==F){
        N.o <- table(dat.tb[[colname.cluster]],dat.tb[[colname.tissue]])
        res.chisq <- chisq.test(N.o,simulate.p.value=simulate.p.value)
        R.oe <- (res.chisq$observed)/(res.chisq$expected)
    }else{
        N.o.byPatient <- table(dat.tb[[colname.patient]],
                               dat.tb[[cluster.colname]], dat.tb[[colname.tissue]])
        R.oe <- apply(N.o.byPatient,1,function(x){
            res.chisq <- chisq.test(x)
            return((res.chisq$observed)/(res.chisq$expected))
        })
    }
    return(R.oe)
}
