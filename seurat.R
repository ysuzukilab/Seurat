#!/usr/bin/env Rscript

#This script runs Seurat and generates plots and input matrices for cibersort
#Based on the Seurat tutorial 'Guided Clustering Tutorial (2,700 PBMCs)'
#Reference: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html

#USAGE: Rscript  seurat.R -i input_dir/ -o output_dir/ 

library(dplyr)
library(Seurat)
library(optparse)
library(grid)
library(ggpubr)

validate.folder <- function(s){
    #convert directory paths so that they end with '/'; i.e. dir/
    if (substring(s,nchar(s),nchar(s))=='/'){
        return(s)
    }else{
        return(paste(s,'/',sep=''))
    }
}

option_list <- list(
    make_option(c('-i','--input_dir'),type='character',default='./',help='[Required] Input data directory [default ./]', metavar='character'),
    make_option(c('-o','--output_dir'),type='character',default='./',help='[Recommended] Output plot/text data directory [default ./]', metavar='character'),
    make_option(c('-n','--project_name'),type='character',default='sample',help='[Recommended] Name of project. Output files will be given this name [default sample]', metavar='character'),
    make_option(c('--visualise_rawdata'),action='store_true',default=FALSE,help='[Recommended/Preprocessing] Plot raw data to determine preprocessing cutoff values. [default FALSE]'),
    make_option(c('-s','--nFeature_RNA_min'),type='integer',default=200,help='[Preprocessing; Filtering data] Minimum of nFeature_RNA [default 200]', metavar='number'),
    make_option(c('-l','--nFeature_RNA_max'),type='integer',default=0,help='[Preprocessing; Filtering data] Maximum of nFeature_RNA [default FALSE]', metavar='number'),
    make_option(c('-m','--percent.mt'),type='integer',default=5,help='[Preprocessing; Filtering data] Maximum percentage of mitochondria genome. Higher percent.mt indicates dead cell [default 5]', metavar='number'),
    make_option(c('-d','--pca_dims'),type='integer',default=10,help='[Clustering/Dimensionality reduction] Number of principal components to use. Cf. pca_jackstraw.png, pca_elbowPlot.png [default 10]', metavar='number'),
    make_option(c('-j','--jackstrawed'),type='character',default='',help='[Resume] Set RDS file name (and the path to that file) if you want to use previously calculated JackStraw results (e.g. *_jackstrawed.rds) [default FALSE]',metavar='character'),
    make_option(c('--deg'),type='character',default='',help='[Resume] Set RDS file name (and the path to that file) if you want to use previously calculated results and resume from finding DEGs(e.g. *_final.rds) [default FALSE]',metavar='character'),
    make_option(c('-g','--marker_genes'),type='character',default='',help='[Finding DEGs] path/name of textfile that contains custom marker genes', metavar='character'),
    make_option(c('-c','--cluster_id'),type='character',default='',help='[Assigning cell type identity] path/name of textfile that contains custom cluster IDs', metavar='character'),
    make_option(c('--marker_threshold'),type='integer',default=0,help='[Finding biomarkers] avg_logFC threshold  [default FALSE]',metavar='number'),
    make_option(c('--marker_numbers'),type='integer',default=0,help='[Finding biomarkers] Top n genes to retreive as marker candidate   [default 10]',metavar='number'),
    make_option(c('-f','--find_marker'),action='store_true',default=FALSE,help='[Finding biomarkers] Plot all Marker candidates in seperate scatter plot [default FALSE]')
)

opt <- parse_args(OptionParser(option_list=option_list))
data.directory = validate.folder(opt$input_dir)
out.dir = validate.folder(opt$output_dir)

pre_processing <- function(sample){
    #Data Preprocessing: remove low quality data (cells).
    #Below example calculates percentage of reads mapped to mitochondria genome, since high mito => dying cell
    #This preprocessing part should ideally be adjusted according to project. (use --visualise_rawdata !)
    sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
    
    #Visualise data quality
    png(filename=paste(out.dir,"preprocessing_violin.png",sep=''))
    plot(VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    dev.off()

    #Visualise feature-feature relationships
    plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    png(filename=paste(out.dir,"preprocessing_scatter.png",sep=''),width=2000,height=1000)
    print(CombinePlots(plots = list(plot1, plot2)))
    dev.off()

    if (opt$visualise_rawdata==TRUE){
        return()
    }
    #print('hoo')
    #Apply filter; remove low quality data
    #Adjust --nFeature_RNA_min, --percent.mt, --nFeature_RNA_max according to how data looks in preprocessing_violin.png & preprocessing_scatter.png
    if (opt$nFeature_RNA_max == 0){
        sample <- subset(sample, subset = nFeature_RNA > opt$nFeature_RNA_min & percent.mt < opt$percent.mt)
    }else{
        sample <- subset(sample, subset = nFeature_RNA > opt$nFeature_RNA_min & nFeature_RNA < opt$nFeature_RNA_max & percent.mt < opt$percent.mt)
    }

    #Data normalisation; normalised data stored in sample[["RNA"]]@data
    sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)

    #Indentify highly variable features (feature selection)
    #Calculate a subset of features that exhibit high cell-to-cell variation in the dataset
    sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)

    # Plot variable features with and without labels
    png(filename=paste(out.dir,'preprocessing_variable_features.png',sep=''), width = 2000, height = 1000)
    top10 <- head(VariableFeatures(sample), 10) # Identify the 10 most highly variable genes
    plot1 <- VariableFeaturePlot(object=sample)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge=0,ynudge=0)
    plot(CombinePlots(plots = list(plot1, plot2)))
    dev.off()

    #Scaling data: Standard pre-processing step prior to dimensional reduction
    all.genes <- rownames(sample)
    sample <- ScaleData(sample, features = all.genes)

    print('Data preprocessing - DONE')
    return(sample)
}

run_pca <- function(sample){
    #Perform linear dimensional reduction
    sample <- RunPCA(sample, features = VariableFeatures(object = sample))
    sink(paste(out.dir,'pca_results.txt',sep=''))
    print(sample[["pca"]])
    sink()

    # Examine and visualize PCA results in a few different ways	
    png(filename=paste(out.dir,"pca_genes.png",sep=''))
    print(VizDimLoadings(sample, dims = 1:2, reduction = "pca"))
    dev.off()

    png(filename=paste(out.dir,"pca_scatter.png",sep=''))
    print(DimPlot(sample, reduction = "pca"))
    dev.off()
	
    png(filename=paste(out.dir,"pca_heatmap.png",sep=''),width=4000,height=6000,pointsize=75)
    DimHeatmap(sample, dims = 1:15, cells = 500, balanced = TRUE)
    dev.off()

    #Determine the 'dimensionality' of the dataset
    sample <- JackStraw(sample, num.replicate = 100) #This may take a while
    print('JackStraw - DONE')
    sample <- ScoreJackStraw(sample, dims = 1:20)

    png(filename=paste(out.dir,"pca_jackstraw.png",sep=''))
    plot(JackStrawPlot(sample, dims = 1:15)) #Jackstraw/Jackstrawplot may be replaced by Elbowplot 
    dev.off()
	
    png(filename=paste(out.dir,"pca_elbowPlot.png",sep=''))
    plot(ElbowPlot(sample))
    dev.off()

    saveRDS(sample,file=paste(out.dir,opt$project_name,'_jackstrawed.rds',sep='')) #Save this to avoid running jackstraw when possible
    print('PCA - DONE')
    return(sample)
}

cluster <- function(sample){
    #Cluster the cells
    sample <- FindNeighbors(sample, dims = 1:opt$pca_dims)
    sample <- FindClusters(sample, resolution = 0.5)
    print('Clustering - DONE')
    return(sample)
}

Dimensionality_reduction <- function(sample,method){
    ##Run non-linear dimensional reduction (UMAP/tSNE)
    sample <- RunUMAP(sample, dims = 1:opt$pca_dims)
    png(filename=paste(out.dir,"cluster_scatter.png",sep=''))
    print(DimPlot(sample, reduction = "umap",label=TRUE))
    dev.off()
    print('Dimensionality reduction - DONE')
    return(sample)
}

find_DEG <- function(sample){
    #Find differentially expressed features (Identifying cluster biomarkers)
    sample.markers <- FindAllMarkers(sample, min.pct = 0.25, logfc.threshold = 0.25)
    sample.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    markers10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    write.table(sample.markers,file=paste(out.dir,opt$project_name,'_all_markers.tsv',sep=''),sep='\t',na='',row.names=T,col.names=NA,quote=F)
    write.table(markers10,file=paste(out.dir,opt$project_name,'_top10_markers.tsv',sep=''),sep='\t',na='',row.names=T,col.names=NA,quote=F)

    #Visualise user specified marker gene expressions
    if (opt$marker_genes!=''){
        feature_genes <- scan(opt$marker_genes,character())
        png(filename=paste(out.dir,"deg_feature.png",sep=''))
        print(FeaturePlot(sample, features = feature_genes))
        dev.off()
    }

    #Plot highly expressed genes (biomarker candidate)
    if (opt$marker_numbers != 0){
        subset <- sample.markers %>% group_by(cluster) %>% top_n(n = opt$marker_numbers, wt = avg_logFC)
    } else if (opt$marker_threshold != 0){
        subset <- sample.markers[sample.markers$avg_logFC>opt$marker_threshold,]
    }else{
        subset <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    }
    if (opt$find_marker == TRUE){
        png(filename=paste(out.dir,"deg_featureplot.png",sep=''),width=6000,height=9000)
        print(FeaturePlot(sample, features =subset$gene))
        dev.off()
    }
    #Plot top 10 highly expressed genes vs cluster/Cell phenotype
    png(filename=paste(out.dir,"deg_heatmap.png",sep=''),width=1000,height=1000)
    print(DoHeatmap(sample, features = markers10$gene) + NoLegend(),width=1000,height=1000)
    dev.off()

    print('Finding DEGs - DONE')
    return(sample)
}

assign_cellType_to_cluster <- function(sample){
    #Assign cell type identity to clusters	
    new.cluster.ids <- scan(opt$cluster_id,character())
    names(new.cluster.ids) <- levels(sample)
    sample <- RenameIdents(sample, new.cluster.ids)

    #Scatter plot coloured by cell phenotype/cluster
    png(filename=paste(out.dir,"celltype_scatter.png",sep=''))
    plot(DimPlot(sample, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
    dev.off()

    print('Assigning cell type to cluster - DONE')
    return(sample)
}

write_to_output <- function(sample){
    #Barcode - Phenotype list
    write.table(Idents(sample),file=paste(out.dir,opt$project_name,'_cell_ID_phenotype.tsv',sep=''),sep='\t',na='',row.names=T,col.names=NA,quote=F)

    #Gene - Cell Phenotype matrix
    oldnames = colnames(sample@assays$RNA@counts)
    newnames <- c()
    #Review the following rows!
    for(barcode in oldnames){
        newnames <- append(newnames,as.character(sample@active.ident[barcode]))#factor (dct-like) of barcode-phenotype: sample@active.ident
    }
    df = as.data.frame(sample@assays$RNA@counts)
    names(df) <- newnames
    #sample@assays$RNA@counts %>% rename_at(vars(oldnames), ~ newnames)
    write.table(df,file=paste(out.dir,opt$project_name,'_cell_expression.tsv',sep=''),sep='\t',na='',row.names=T,col.names=NA,quote=F)
}

main <- function(){
    if (opt$deg != ''){
        #Resume from finding Differentially Expressed Genes (DEGs)
        print('resuming...')
        sample <- readRDS(file=opt$deg)
    }else if (opt$jackstrawed == ''){
        #Start a fresh analysis from Step one!
        # Load the PBMC dataset
        sample.data <- Read10X(data.dir = data.directory)
        # Initialize the Seurat object with the raw (non-normalized data).
        sample <- CreateSeuratObject(counts = sample.data, project = opt$project_name, min.cells = 3, min.features = 200)
        sample = pre_processing(sample)
        if (opt$visualise_rawdata == TRUE){
            print('DONE! Raw data has been visualised.')
            print('Next: Determine cutoff values.Then run the same script without --visualise_rawdata')
            return()
    }
        sample = run_pca(sample)
        sample = cluster(sample)
        sample = Dimensionality_reduction(sample,'umap')
    }else{
        #Resume from after Jackstraw calculation
        print('resuming...')
        sample <- readRDS(file=opt$jackstrawed)
        sample = cluster(sample)
        sample = Dimensionality_reduction(sample,'umap')
    }
    sample = find_DEG(sample)
    if (opt$cluster_id != ''){
        sample = assign_cellType_to_cluster(sample)
    }
    saveRDS(sample,file=paste(out.dir,opt$project_name,'_final.rds',sep=''))
    write_to_output(sample)
    print('ANALYSIS COMPLETE!')

    print('*****Ignore error message about as.character() being deprecated*****')
}

main()
