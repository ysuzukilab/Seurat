# Seurat
Repository for Seurat analysis

## Usage  
### Basic usage      
```R
Rscript  seurat.R -i input_dir/ -o output_dir/ --project_name 'custom_project_name'     
```  
### Step 1: Initial run  
Visualise raw data and determine cutoff values to use for dropping low quality  data/cells  
```R
Rscript  seurat.R -i input_dir/ -o output_dir/ --project_name 'custom_project_name' --visualise_rawdata  
```  
### Step 2:   
```R
Rscript seurat.R  -i input_dir/ -o output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name'   
```  
### Step 3:  



