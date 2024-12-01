#! /bin/bash
#=======================================================================================
#                                   pySCENIC数据分析
#=======================================================================================
reference_path='/home/jifanghan/Documents/R_pipeline/reference/pyscenic/'
TF_list=$reference_path"hs_hgnc_tfs.txt"
motifs=$reference_path"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
feather=$reference_path"hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"

#分析第一步：GRN---运行完得到sce.adj.csv文件
##这一步的目的
#推断转录因子与提供的表达矩阵基因的共表达模块，基于grnboost2，R中时GENIE3
pyscenic grn --num_workers 10 \
--sparse \
--method grnboost2 \
--output sce.adj.csv \
sce.loom \
$TF_list

#分析第二步：RcisTarget---运行完得到sce.regulons.csv文件
#这一步的目的
#进行TF-motif富集分析，识别直接靶标
#得到转录因子(TF)与其对应的直接作用的靶点,称为regulon(每一个regulon是1个TF和其调控的靶基因)
pyscenic ctx --num_workers 10 \
--output sce.regulons.csv \
--expression_mtx_fname sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname $motifs \
sce.adj.csv \
$feather


#分析第三步：AUCell---运行完得到sce_SCENIC.loom文件，即分析结果
#这一步的目的
#使用AUCell对每个细胞的每个regulon活性进行评分。
pyscenic aucell --num_workers 3 \
--output sce_SCENIC.loom \
sce.loom \
sce.regulons.csv
