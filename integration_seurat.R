#!/usr/bin/env Rscript

#Based on the Seurat tutorial 'Integrating stimulated vs. control PBMC datasets to learn cell-type specific responses'
#Reference: https://satijalab.org/seurat/v3.0/immune_alignment.html
#Ignore lines starting w/ '##' as they merely explain R grammer
#USAGE: Rscript integration_seurat.R -1 /path/to/dir1/ -2 /path/to/dir2/ -o /path/to/output -n1 normal -n2 tumor -n 'custom_project_name'

library(dplyr)
library(Seurat)
library(optparse)
library(cowplot)
library(grid)
library(ggpubr)

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
    make_option(c('-1','--input_dir1'),type='character',default='./',help='[Required] Input data directory1 [default ./]', metavar='character'),
    make_option(c('-2','--input_dir2'),type='character',default='./',help='[Required] Input data directory2 [default ./]', metavar='character'),
    make_option(c('-3','--name_dir1'),type='character',default='sample1',help='[Recommended] Name for directory1 [default sample1]', metavar='character'),
    make_option(c('-4','--name_dir2'),type='character',default='sample2',help='[Recommended] Name for directory2 [default sample2]', metavar='character'),
    make_option(c('-o','--output_dir'),type='character',default='./',help='[Recommended] Output plot/text data directory [default ./]', metavar='character'),
    make_option(c('-n','--project_name'),type='character',default='sample',help='[Recommended] name of project. output csv files will be given this name [default sample]', metavar='character'),
    make_option(c('-c','--cluster_id'),type='character',default='',help='[Assigning cell type identity] path/name of textfile that contains custom cluster IDs', metavar='character'),
    make_option(c('-g','--marker_genes'),type='character',default='',help='[Finding DEGs] path/name of textfile that contains custom marker genes', metavar='character'),
    make_option(c('-f','--find_marker'),action='store_true',default=FALSE,help='[Finding biomarkers] Plot all Marker candidates in seperate scatter plot [default FALSE]')
)

opt <- parse_args(OptionParser(option_list=option_list))
data1.directory = validate.folder(opt$input_dir1)
data2.directory = validate.folder(opt$input_dir2)
out.dir = validate.folder(opt$output_dir)

Setup <- function(sample.data, name){
    sample <- CreateSeuratObject(counts = sample.data, min.cells=5)
    sample$stim <- name
    sample <- subset(sample, subset=nFeature_RNA > 300)
    sample <- NormalizeData(sample, verbose=FALSE)
    sample <- FindVariableFeatures(sample, selection.method="vst",nfeatures=3000)
    return(sample)
}

Integrate <- function(sample1, sample2){
    #Performing Integration
    sample.anchors <- FindIntegrationAnchors(object.list=list(sample1, sample2), dims=1:20)
    sample.combined <- IntegrateData(anchorset = sample.anchors, dims=1:20)
    #Performingan integarted analysis
    DefaultAssay(sample.combined) <- "integrated"
    sample.combined <- ScaleData(sample.combined, verbose=FALSE)
    sample.combined <- RunPCA(sample.combined, npcs=30, verbose=FALSE)
    sample.combined <- RunUMAP(sample.combined, reduction="pca", dims=1:20)
    sample.combined <- FindNeighbors(sample.combined, reduction="pca", dims=1:20)
    sample.combined <- FindClusters(sample.combined, resolution=0.5)
    p1 <- DimPlot(sample.combined, reduction="umap", group.by="stim")
    p2 <- DimPlot(sample.combined, reduction="umap", label=TRUE)
    png(filename=paste(out.dir,"dimplot_combined.png",sep=''))
    print(plot_grid(p1,p2))
    dev.off()
    png(filename=paste(out.dir,"dimplot_split.png",sep=''))
    print(DimPlot(sample.combined, reduction="umap",split.by="stim"))
    dev.off()
    return(sample.combined)
}

Identify <- function(sample.combined){
    #Finding conserved markers
    DefaultAssay(sample.combined) <- "RNA"
    
    sample.markers <- FindAllMarkers(sample.combined, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
    markers10 <- sample.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
    write.table(sample.markers,file=paste(out.dir,opt$project_name,'_all_markers.tsv',sep=''),sep='\t',na='',row.names=T,col.names=NA,quote=F)
    write.table(markers10,file=paste(out.dir,opt$project_name,'_top10_markers.tsv',sep=''),sep='\t',na='',row.names=T,col.names=NA,quote=F)

    if(opt$cluster_id != ''){
        new.cluster.ids <- scan(opt$cluster_id,charater())
        names(new.cluster.ids) <- levels(sample.combined)
        sample.combined <- RenameIdents(sample.combined, new.cluster.ids)
        png(filename=paste(out.dir,"celltype_scatter.png",sep=''))
        print(DimPlot(sample.combined, label=TRUE))
        dev.off()
    }
    
    if(opt$find_marker == TRUE){
        png(filename=paste(out.dir,"deg_featureplot.png",sep=''),width=6000,height=9000)
        print(FeaturePlot(sample.combined, features=markers10$gene))
        dev.off()
    }
    return(sample.combined)
}

main <- function(){
    #Loading dataset
    sample1.data <- Read10X(data.dir = data1.directory)
    sample2.data <- Read10X(data.dir = data2.directory)
    #Setting up objects
    sample1 <- Setup(sample1.data, opt$name_dir1)
    sample2 <- Setup(sample2.data, opt$name_dir2)
    print('Setup - DONE')
    #Performing integartion
    sample.combined <- Integrate(sample1, sample2)
    print('Integration - DONE')
    #Identifying conserved cell type markers
    sample.combined <- Identify(sample.combined)
    print('Identifying markers - DONE')
    saveRDS(sample.combined, file=paste(out.dir,opt$projoect,'_final.rds',sep=''))
    print('ANALYSIS COMPLETE!')
}

main()
