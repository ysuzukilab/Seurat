# Seurat
Repository for Seurat analysis

## Overview
### Input  
10x chromium output.  
./input_dir  
	L--- barcodes.tsv.gz  
	L--- features.tsv.gz  
	L--- matrix.mtx.gz  
### Output  
- .png format plots  
- .txt files  
- .tsv files -> may be used as Cibersort input  
- .rds files  

## Usage  
### Basic usage      
```Bash
Rscript  seurat.R -i ./input_dir/ -o ./output_dir/ --project_name 'custom_project_name'     
```  
### Step 1: Initial run  
Visualise raw data and determine cutoff values to use for dropping low quality  data/cells  
```Bash
Rscript  seurat.R -i ./input_dir/ -o ./output_dir/ --project_name 'custom_project_name' --visualise_rawdata  
```  
### Step 2:   
```Bash
Rscript seurat.R  -i ./input_dir/ -o ./output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name'   
```  
### Step 3:  
```Bash
Rscript seurat.R  -i ./input_dir/ -o ./output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name' --marker_genes ./path_to/marker_gene_file.txt --cluster_id ./path_to/cluster_id.txt  
```


