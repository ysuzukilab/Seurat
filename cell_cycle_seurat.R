#!/usr/bin/env Rscript

#Based on the Seurat tutorial 'Cell-Cycle Scoring and Regression'
#Reference: https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html
#Ignore lines starting w/ '##' as they merely explain R grammer
#USAGE: Rscript cell_cycle_seurat.R -i /path/to/filtered_feature_bc_matrix -o /path/to/output/ 

library(Seurat)
library(optparse)
library(Matrix)

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
    make_option(c('-i','--input_dir'),type='character',default='./',help='[Recommended] Input data directory [default ./]', metavar='character'),
    make_option(c('-o','--output_dir'),type='character',default='./',help='[Recommended] Output plot/text data directory [default ./]', metavar='character'),
    make_option(c('-n','--project_name'),type='character',default='sample',help='[Recommended] name of project. output csv files will be given this name [default sample]', metavar='character')
)


opt <- parse_args(OptionParser(option_list=option_list))
data.directory = validate.folder(opt$input_dir)
out.dir = validate.folder(opt$output_dir)


Initialize <- function(sample){
    #Performing Initialization
    sample <- NormalizeData(sample)
    sample <- FindVariableFeatures(sample, selection.method="vst")
    sample <- ScaleData(sample, features=rownames(sample))
    return(sample)
}

Run_pca <- function(sample){
    #Perform linear dimensional reduction
    sample <- RunPCA(sample, features = VariableFeatures(sample))
    sink(paste(out.dir,'pca_results.txt',sep=''))
    print(sample[["pca"]], dims = 1:10, nfeatures = 10)
    sink()

    #Visualizing PCA results in a few different ways	
    png(filename=paste(out.dir,"pca_genes.png",sep=''))
    print(VizDimLoadings(sample, dims = 1:4, reduction = "pca"))
    dev.off()

    png(filename=paste(out.dir,"pca_scatter.png",sep=''))
    print(DimPlot(sample, reduction = "pca"))
    dev.off()
	
    png(filename=paste(out.dir,"pca_heatmap.png",sep=''),width=4000,height=6000,pointsize=75)
    print(DimHeatmap(sample, dims = 1:15, cells=500, balanced = TRUE))
    dev.off()
    
    return(sample)
}

Score <- function(sample){
    #Cell cycle scoring
    sample <- CellCycleScoring(sample, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
    #View cell cycle scores and phase assignments
    write.table(sample[[]],file=paste(out.dir,opt$project_name,'_cellcycle.tsv',sep=''),sep='\t',row.names=T,col.names=NA,quote=F)
    return(sample)
}

Visualize <- function(sample){
    #Visualizing the distribution of cell cycle markers across sample
    png(filename=paste(out.dir,"ridgeplot.png",sep=''))
    print(RidgePlot(sample, features=c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol=2))
    dev.off()
    #Looking at PC genes
    sample <- RunPCA(sample, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
    sink(paste(out.dir,'pca_results_ccgenes.txt',sep=''))
    print(sample[["pca"]], dims = 1:5)
    sink()
    #Plotting sample by cell cycle
    png(filename=paste(out.dir,"cc_dimplot.png",sep=''))
    print(DimPlot(sample))
    dev.off()

}

main <- function(){
    # Load dataset
    sample.data <- Read10X(data.dir = data.directory)
    #A list of cell cycle markers
    sample <- CreateSeuratObject(counts = sample.data)
    #Initialization
    sample <- Initialize(sample)
    print('Initialization - DONE')
    #Running PCA
    sample <- Run_pca(sample)
    print('Running PCA - DONE')
    #Assigning cell-cycle scores
    sample <- Score(sample)
    print('Cell-cycle Scoring - DONE')
    #Visualization
    sample <- Visualize(sample)
    print('Visualization - DONE')
    print('ANALYSIS COMPLETE!')
}

main()
