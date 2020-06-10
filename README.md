
Table of Contents  
=================  
<!--ts-->
      * [Prerequisites - Install R packages](#prerequisites---install-r-packages)  
      * [10x chromium single cell analysis](#10x-chromium-single-cell-analysis)
         * [Input](#input)
         * [Output](#output)
      * [Usage](#usage)
         * [Basic usage](#basic-usage)
         * [Step 1: Initial run](#step-1-initial-run)
         * [Step 2: Determine dimensions to use](#step-2-determine-dimensions-to-use)
         * [Step 3: Determine cluster phenotypes](#step-3-determine-cluster-phenotypes)
         * [Step 4: Final/Complete run](#step-4-finalcomplete-run)
         * [Step 5: Optional: Identifying novel biomarkers](#step-5-optional-identifying-novel-biomarkers)
      * [Cell Cycle Analysis](#cell-cycle-analysis)
         * [Input](#input-1)
         * [Output](#output-1)
      * [Usage](#usage-1)
         * [Basic usage](#basic-usage-1)
      * [Integrating datasets](#integrating-datasets)
         * [Input](#input-2)
         * [Output](#output-2)
      * [Usage](#usage-2)
         * [Basic usage](#basic-usage-2)
      * [Integrating scATAC and scRNA](#integrating-scatac-and-scrna)
         * [Input](#input-3)
         * [Output](#output-3)
      * [Usage](#usage-3)
         * [Basic usage](#basic-usage-3)
<!--te-->

## Prerequisites - Install R packages
```R
options(repos="http://cran.ism.ac.jp")
install.packages('dplyr')			
install.packages('Seurat')
install.packages('optparse')
install.packages('ggpubr')
install.packages('cowplot')
install.packages('grid')
install.packages('ggplot2')
install.packages('hdf5r')
```
Downloading the package "hdf5r" may not work on ssh.


## 10x chromium single cell analysis
10x_seurat.R automates the analysis shown in the [Seurat tutorial](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html 'Seurat pbmc tutorial').  
### Input  
10x chromium output files.    
./input_dir  
	├--- barcodes.tsv.gz  
	├--- features.tsv.gz  
	└--- matrix.mtx.gz  
### Output  
- .png format plots  
- .txt files  
- .tsv files -> may be used as Cibersortx input  
- .rds files  

## Usage  
### Basic usage      
```Bash
$ Rscript  seurat.R -i ./input_dir/ -o ./output_dir/ --project_name 'custom_project_name'     
```  
### Step 1: Initial run 
Visualise raw data and determine cutoff values to use for dropping low quality  data/cells.  
```Bash
$ Rscript  seurat.R -i ./input_dir/ -o ./output_dir/ --project_name 'custom_project_name' --visualise_rawdata  
```  
### Step 2: Determine dimensions to use
According to your plots generated in step 1, assign integers to --nFeature_RNA_min, --nFeature_RNA_max and --percent.mt custom_value.
```Bash
$ Rscript seurat.R  -i ./input_dir/ -o ./output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name'  
```  
### Step 3: Determine cluster phenotypes  
According to pca_*.png plots from Step2, reassign --pca_dims with the number of dimensions suitable for your dataset.  
In addition, visualise the expression of previously identified genes by making a .txt file with such genes and specifying it under --marker_genes.  
Note: If you do not have any a priori knowledge concerning markers, you may run without --marker_genes option.  
Example of marker gene file:
```
"MS4A1"
"GNLY"
```
```Bash
$ Rscript seurat.R  -i ./input_dir/ -o ./output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name' --jackstrawed ./output_dir/custom_project_name_jackstrawed.rds --marker_genes ./path_to/marker_gene_file.txt --find_marker
```
### Step 4: Final/Complete run
Use deg_heatmap.png, deg_features.png and cluster_scatter.png to assign each cluster a phenotype. Assign Phenotype by making a .txt file where each row contains the phenotype of a cluster. Specify this file using --cluster_id.  
Below is an example of text file. In this example, the cluster 1 represents Platelets, and cluster 2 represents CD14+ Mono cells.
```
Platelet
CD14+ Mono
```

```Bash
$ Rscript seurat.R  -i ./input_dir/ -o ./output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name'  --jackstrawed ./output_dir/custom_project_name_jackstrawed.rds --marker_genes ./path_to/marker_gene_file.txt --cluster_id ./path_to/cluster_id.txt  
```
### Step 5: Optional: Identifying novel biomarkers
Use deg_heatmap.png, deg_features.png, *_all_markers.tsv and *_top10_markers.tsv to identify novel biomarkers OR to visualise the expression of known biomarkers. List the genes in interest in a .txt file and specify under --marker_genes.  
```Bash
$ Rscript seurat.R  -i ./input_dir/ -o ./output_dir/ --nFeature_RNA_min custom_value --nFeature_RNA_max custom_value --percent.mt custom_value --project_name 'custom_project_name' --jackstrawed ./output_dir/custom_project_name_jackstrawed.rds --marker_genes ./path_to/marker_gene_file.txt --cluster_id ./path_to/cluster_id.txt
```

```Bash
Options:
    -i CHARACTER, --input_dir=CHARACTER
        [Required] Input data directory [default ./]
    -o CHARACTER, --output_dir=CHARACTER
        [Recommended] Output plot/text data directory [default ./]
    -n CHARACTER, --project_name=CHARACTER
        [Recommended] Name of project. Output files will be given this name [default sample]
    --visualise_rawdata
        [Recommended/Preprocessing] Plot raw data to determine preprocessing cutoff values. [default FALSE]
    -s NUMBER, --nFeature_RNA_min=NUMBER
        [Preprocessing; Filtering data] Minimum of nFeature_RNA [default 200]
    -l NUMBER, --nFeature_RNA_max=NUMBER
        [Preprocessing; Filtering data] Maximum of nFeature_RNA [default FALSE]
    -m NUMBER, --percent.mt=NUMBER
        [Preprocessing; Filtering data] Maximum percentage of mitochondria genome. Higher percent.mt indicates dead cell [default 5]
    -d NUMBER, --pca_dims=NUMBER
        [Clustering/Dimensionality reduction] Number of principal components to use. Cf. pca_jackstraw.png, pca_elbowPlot.png [default 10]
    -j CHARACTER, --jackstrawed=CHARACTER
        [Resume] Set RDS file name (and the path to that file) if you want to use previously calculated JackStraw results (e.g. *_jackstrawed.rds) [default FALSE]
    --deg=CHARACTER
        [Resume] Set RDS file name (and the path to that file) if you want to use previously calculated results and resume from finding DEGs (e.g. *_final.rds) [default FALSE]
    -g CHARACTER, --marker_genes=CHARACTER
       [Finding DEGs] path/name of textfile that contains custom marker genes
    -c CHARACTER, --cluster_id=CHARACTER
        [Assigning cell type identity] path/name of textfile that contains custom cluster IDs
    --marker_threshold=NUMBER
        [Finding biomarkers] avg_logFC threshold  [default FALSE]
	--marker_numbers=NUMBER
        [Finding biomarkers] Top n genes to retreive as marker candidate   [default 10]
    -f, --find_marker
        [Finding biomarkers] Plot all Marker candidates in seperate scatter plot [default FALSE]
    -h, --help
        Shows this help message and exit
```

## Cell Cycle Analysis
cell_cycle_seurat.R automates the analysis shown in the [Seurat tutorial](https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html 'Cell-Cycle Scoring and Regression'). However this does not regress out scores but only assigns cell-cycle scores.

### Input  
10x chromium output files.    
./input_dir  
	├--- barcodes.tsv.gz  
	├--- features.tsv.gz  
	└--- matrix.mtx.gz  
### Output  
- .png format plots  
- .txt files  
- .tsv files 

## Usage  
### Basic usage      
```Bash
$ Rscript cell_cycle_seurat.R -i /path/to/input_dir/ -o /path/to/output_dir/
```  

## Integrating datasets
integration_seurat.R automates the analysis shown in the [Seurat tutorial](https://satijalab.org/seurat/v3.0/immune_alignment.html 'Tutorial: Integrating stimulated vs. control PBMC datasets to learn cell-type specific responses').

### Input  
Two 10x chromium output files.    
./input_dir1  
	├--- barcodes.tsv.gz  
	├--- features.tsv.gz  
	└--- matrix.mtx.gz  
./input_dir2  
	├--- barcodes.tsv.gz  
	├--- features.tsv.gz  
	└--- matrix.mtx.gz  

### Output  
- .png format plots  
- .tsv files 
- .rds files

## Usage  
### Basic usage      
```Bash
$ Rscript integration_seurat.R --input_dir1 /path/to/dir1/ --input_dir2 /path/to/dir2/ -o /path/to/output -n1 name1 -n2 name2 -n 'custom_project_name'
```  

## Integrating scATAC and scRNA
atac_integration_seurat.R automates the analysis shown in the [Seurat tutorial](https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html 'PBMC scATAC-seq Vignette').

### Input  
- peak_matrix.h5 file
- .gtf file
- singlecell.csv file
- scRNA.rds file

### Output  
- .png format plots  
- .txt files 
- .rds files

## Usage  
### Basic usage      
```Bash
$ Rscript atac_integration_seurat.R -a path/to/peak_bc_matrix.h5 -b path/to/gtf -c path/to/singlecell.csv -r path/to/scRNA.rds -o 'custom_ouput_name' -n 'custom_project_name'
```  












































