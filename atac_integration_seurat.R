#!/usr/bin/env Rscript

#Based on the Seurat tutorial 'PBMC scATAC-seq Vignette'
#Reference: https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html
#Ignore lines starting w/ '##' as they merely explain R grammer
#USAGE: Rscript atac_integration_seurat.R -a /path/to/atac_data -b /path/to/annot_data -c path/to/csv -d nCount_num -r path/to/rds -o /path/to/output -n 'custom_project_name'

library(dplyr)
library(Seurat)
library(optparse)
library(cowplot)
library(grid)
library(ggpubr)
library(ggplot2)
library(hdf5r)

#this scripts receives working directories as arguments
#directory names must end with '/'; i.e. dir/
validate.folder <- function(s){
    if (substring(s,nchar(s),nchar(s))=='/'){
        return(s)
    }else{
        return(paste(s,'/',sep=''))
    }
}

option_list <- list(
    make_option(c('-a','--atac_data'),type='character',default='./',help='[Required] Input peak data (peak_bc_matrix.h5) [default ./]', metavar='character'),
    make_option(c('-b','--annot_data'),type='character',default='./',help='[Required] Input annotation data (.gtf) [default ./]', metavar='character'),
    make_option(c('-c','--atac_csv'),type='character',default='./',help='[Required] Input csv data (singlecell.csv) [default ./]', metavar='character'),
    make_option(c('-d','--atac_count'),type='integer',default=1000,help='[Setup; Filtering data] Minimum of nCount_ATAC [default 1000]', metavar='number'),
    make_option(c('-r','--rds_data'),type='character',default='./',help='[Required] Input rds data (10x_v3.rds) [default ./]', metavar='character'),
    make_option(c('-o','--output_dir'),type='character',default='./',help='[Recommended] Output plot/text data directory [default ./]', metavar='character'),
    make_option(c('-n','--project_name'),type='character',default='sample',help='[Recommended] name of project. output csv files will be given this name [default sample]', metavar='character')
)

opt <- parse_args(OptionParser(option_list=option_list))
#data1.directory = validate.folder(opt$input_dir1)
#data2.directory = validate.folder(opt$input_dir2)
out.dir = validate.folder(opt$output_dir)

Setup <- function(peaks, activity.matrix){
    sample.atac <- CreateSeuratObject(counts=peaks, assay="ATAC", project="10x_ATAC")
    sample.atac[["ACTIVITY"]] <- CreateAssayObject(counts=activity.matrix)
    meta <- read.table(opt$atac_csv, sep=",", header=TRUE, row.names=1,stringsAsFactors=FALSE)
    meta <- meta[colnames(sample.atac),]
    sample.atac <- AddMetaData(sample.atac, metadata=meta)
    sample.atac <- subset(sample.atac, subset=nCount_ATAC > opt$atac_count)
    sample.atac$tech <- "atac"
    return(sample.atac)
}

Preprocess_coembed <- function(sample.atac){
    DefaultAssay(sample.atac) <- "ACTIVITY"
    sample.atac <- FindVariableFeatures(sample.atac)
    sample.atac <- NormalizeData(sample.atac)
    sample.atac <- ScaleData(sample.atac)
    
    DefaultAssay(sample.atac) <- "ATAC"
    VariableFeatures(sample.atac) <- names(which(Matrix::rowSums(sample.atac) > 100))
    sample.atac <- RunLSI(sample.atac, n=50, scale.max=NULL)
    sample.atac <- RunUMAP(sample.atac, reduction="lsi", dims=1:50)
    
    sample.rna <- readRDS(opt$rds_data)
    sample.rna$tech <- "rna"
    p1 <- DimPlot(sample.atac, reduction="umap") + NoLegend() + ggtitle("scATAC-seq")
    p2 <- DimPlot(sample.rna, group.by = "celltype", label=TRUE, repel=TRUE) + NoLegend() + ggtitle("scRNA-seq")
    png(filename=paste(out.dir,"cc_combineplot.png",sep=''))
    print(CombinePlots(plots=list(p1,p2)))
    dev.off()
    
    transfer.anchors <- FindTransferAnchors(reference=sample.rna, query=sample.atac, features=VariableFeatures(object=sample.rna), reference.assay="RNA", query.assay="ACTIVITY", reduction="cca")
    celltype.predictions <- TransferData(anchorset=transfer.anchors, refdata=sample.rna$celltype, weight.reduction=sample.atac[["lsi"]])
    sample.atac <- AddMetaData(sample.atac, metadata=celltype.predictions)
    png(filename=paste(out.dir,"histogram.png",sep=''))
    print(hist(sample.atac$prediction.score.max))
    print(abline(v=0.5, col="red"))
    dev.off()

    sink(paste(out.dir,'table.txt',sep=''))
    print(table(sample.atac$prediction.score.max > 0.5))
    sink()

    sample.atac.filtered <- subset(sample.atac, subset=prediction.score.max>0.5)
    sample.atac.filtered$predicted.id <- factor(sample.atac.filtered$predicted.id, levels=levels(sample.rna))
    p1 <- DimPlot(sample.atac.filtered, group.by="predicted.id", label=TRUE, repel=TRUE) + ggtitle("scATAC-seq cells") + NoLegend() + scale_colour_hue(drop=FALSE)
    p2 <- DimPlot(sample.rna, group.by="celltype", label=TRUE, repel=TRUE) + ggtitle("scRNA-seq cells") + NoLegend()
    png(filename=paste(out.dir,"filtered_combineplots.png",sep=''))
    print(CombinePlots(plots = list(p1,p2)))
    dev.off()

    genes.use <- VariableFeatures(sample.rna)
    refdata <- GetAssayData(sample.rna, assay="RNA", slot="data")[genes.use,]
    imputation <- TransferData(anchorset=transfer.anchors, refdata=refdata, weight.reduction=sample.atac[["lsi"]])
    sample.atac[["RNA"]] <- imputation
    coembed <- merged(x=sample.rna, y=sample.atac)
    coembed <- ScaleData(coembed, features=genes.use, do.sclae=FALSE)
    coembed <- RunPCA(coembed, features=genes.use, verbose=FALSE)
    coembed <- RunUMAP(coembed, dims=1:30)
    coembed$celltype <- ifelse(!is.na(coembed$celltype),coembed$celltype, coembed$predicted.id)
    p1 <- DimPlot(coembed, group.by="tech")
    p2 <- DimPlot(coembed, group.by="celltype", label=TRUE, repel=TRUE)
    png(filename=paste(out.dir,"coembed_combineplots.png",sep=''))
    print(CombinePlots(plots = list(p1,p2)))
    dev.off()

    return(coembed)
}

main <- function(){
    #Loading dataset
    peaks <- Read10X_h5(opt$atac_data)
    activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file=opt$annot_data, verbose=TRUE)
    print('Loading - DONE')
    #Setting up objects
    sample.atac <- Setup(peaks, activity.matrix)
    print('Setup - DONE')
    #Data preprocessing and co-embedding
    coembed <- Preprocess_coembed(sample.atac)
    print('Preprocessing and Co-embedding - DONE')
    #Saving data
    saveRDS(coembed, file=paste(out.dir,opt$projoect,'_final.rds',sep=''))
    print('ANALYSIS COMPLETE!')
}

main()
