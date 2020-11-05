# Separating the Matrix into Sample Objects

The two matrices that I have need to be subset and turned into Seurat objects. 

First, set working directory and load the libraries needed (Seurat and dplyr)

    setwd("/media/sjblair/Ext4Scratch/sc_Spine")
    library(Seurat)
    library(dplyr)
    
The two matrices should be loaded into R.  The dot before the folder tells the system to look relative to the working directory.  Otherwise the system would look for a folder called data in root. 

    sc1.matrix <- read.table('./data/sc1.counts.annotated.matrix.repGene')
    sc2.matrix <- read.table('./data/sc2.counts.annotated.matrix.repGene')

I would like to make the visualizations that will come more readable.  So tidyverse and dplyr will be used to parse the gene name from the contig name string.  Meanwhile we want to remove the ambiguous contigs from the dataset.  It is possible that there may be something of interest in these data but it is not part of this project to study them more.

First we should make a copy of the data.

    sc1.anno <- sc1.matrix
    sc2.anno <- sc2.matrix
Using tidyverse and dplyr, we can parse the contig string to remove the identifying marker.  The contig name is separated by the pipe character.  So, we will split the name by searching for the pipe character.

    library(tidyverse)
    library(dplyr)
    sc1.anno <- separate(sc1.anno, 
                         col = X, 
                         into = c("contig", "uniprot.entry", "gene"), 
                         sep = "\\|")

We can clean this up even more by separating the contig name by the caret to make it look cleaner.  It is important to keep the contig information for comparing later downstream if we want tocidentify contig by "gene" name.

    sc1.anno <- separate(sc1.anno, 
                         col = contig, 
                         into = c("contig"), 
                         sep = "\\^")
    sc1.anno <- separate(sc1.anno, 
                         col = gene, 
                         into = c("gene"), 
                         sep = "\\^")
We can now omit the ambigous data.  Further study of transcripts that have no homologs are beyond the scope of this study.  Finally, we can tidy up the gene names further by cutting the names further.  Later if we want to add them back we can use the original matrix and contig to do so.

    sc1.anno <-na.omit(sc1.anno)
    sc1.anno <- separate(sc1.anno, 
                         col = gene, 
                         into = c("gene"), 
                         sep = "\\_")

You could change the gene names to lower case.  I decided to leave them uppercase though.

    #sc1.anno[[3]] <- tolower(sc1.anno[[3]])

Now we change the rownames to the "gene" name.  This makes it simple to name pca graphs later.

    row.names(sc1.anno) <- make.names(sc1.anno$gene, 
                                      unique = T)

Save the matrix for later.

    saveRDS(sc1.anno, file  =  "./filtered/sc1.anno")


Now do the same with the sc2 dataset.

    sc2.anno <- separate(sc2.anno, 
                         col = X, 
                         into = c("contig", "uniprot.entry", "gene"), 
                         sep = "\\|")
    sc2.anno <- separate(sc2.anno, 
                         col = contig, 
                         into = c("contig"), 
                         sep = "\\^")
    sc2.anno <- separate(sc2.anno, 
                         col = gene, 
                         into = c("gene"), 
                         sep = "\\^")
    sc2.anno <-na.omit(sc2.anno)
    sc2.anno <- separate(sc2.anno, 
                         col = gene, 
                         into = c("gene"), 
                         sep = "\\_")
    #sc2.anno[[3]] <- tolower(sc2.anno[[3]])
    row.names(sc2.anno) <- make.names(sc2.anno$gene, 
                                      unique = T)
    saveRDS(sc2.anno, file  =  "./filtered/sc2.anno")

 

Hani gave me an excel sheet that has the following sample info in it:
||Intact.Limb|Limb 3dpa|Limb 14dpa|Intact Spine|Spine 3dpa|Spine 14dpa|
| :-------------: | :----------: | :-----------: |:-------------: | :-------------: | :-------------: | :-------------: | 
|sc1|||S4, S7|||S9, S12, S14|
|sc2|S3, S7|S1, S2, S26||S4, S5 |S8, S9, S11||
The other data collected were for internal and adjacent models and will not be used here.
## Creating Seurat Objects 
Using R we can easily  create objects directly from the matrices as follows

    limb.intact <- 
      CreateSeuratObject(
        sc2.anno[, grep('^S(3|7)_', colnames(sc2.anno)),], 
        project = "homeostasisLimb",
        assay = "RNA",
        min.cells = 3, 
        min.features = 100,
        names.field = 2,
        names.delim = "_")
    saveRDS(limb.intact, file  =  "./seurat_objects/limb.intact.rds")

I will break down these arguments before moving on and creating the additional objects
|  | Description |
|--|--|
| limb.intact | This is the name of the Seurat object.|
|CreateSeuratObject|The Seurat function that creates objects.|
|counts|unormalized data in the form of a matrix, more below.|
|project|creates a project name for the object that is visible in the metadata.|
|assay|name for initial input data.|
|min.cells|this is a cutoff to only include features detected in at least this many cells. **I choose to filter out features that have a total frequency of 3 or less.** |
|min.features|this cutoff will include cells that have at least this many features.|
|names.field|you can give the cells identification based on their sample name in the matrix. **I choose to include all cells that have at least 100 features in it.** cells that have less than 100 features are prone to be unreliable.|
|names.delim|this is the delim that sepeterates the name.  in our case the name is SX_xxxxx so a delim of _ is specifed.|
|saveRDS|this saves the object limb.intact to the filename that file = points to.|
**You may have noticed that I omitted the counts argument** 
sc2.anno[, grep('^S(3|7)_', colnames(sc2.anno)),] is the counts argument but I used a feature instead.  So, instead of telling it that I wanted to use the file *sc2.anno* I instead used a grep command to remove only the samples that include data for the samples of interest as explained here.
||Description|
|--|--|
|sc2.anno|The tabular data that contains the samples we want.|
|[, grep('^S(3\|7)_',|grep is searching for a 3 character string.  The ^ tells grep to find any strings that begin with the character capital S. The parenthesis and logical operator \| tells grep the next character must be 3 OR 7.  The third character must be underscore.  Combined, grep knows it is looking for strings that begin with S3_ or S7_.     |
|colnames(sc2.anno)|This is the location grep is searching.  Just like in terminal, grep follows the syntax grep string location. This is telling grep to use the search string in the column names of the sc2.anno table.|
**Now it is time to run similar commands for each sample**

For limb 3dpa object the S1_, S2_, S26_ samples from s2 will be subset into a Seurat object

    limb.3dpa <- 
          CreateSeuratObject(
            sc2.anno[, grep('^S(1|2|26)_', colnames(sc2.anno)),], 
            project = "limb_3dpa",
            assay = "RNA",
            min.cells = 3, 
            min.features = 100,
            names.field = 2,
            names.delim = "_")
        saveRDS(limb.3dpa, file  =  "./seurat_objects/limb.3dpa.rds")
        
For intact spine the S4_, S5_ samples from sc2 matrix will be subset into a Seurat object

    spine.intact <- 
          CreateSeuratObject(
            sc2.anno[, grep('^S(4|5)_', colnames(sc2.anno)),], 
            project = "homeostasis_Spine",
            assay = "RNA",
            min.cells = 3, 
            min.features = 100,
            names.field = 2,
            names.delim = "_")
        saveRDS(spine.intact, file  =  "./seurat_objects/spine.intact.rds")
        
For spine 3dpa the S8_, S9_, and S11_ samples from sc2 matrix will be subset into a Seurat object

    spine.3dpa <- 
          CreateSeuratObject(
            sc2.anno[, grep('^S(8|9|11)_', colnames(sc2.anno)),], 
            project = "spine_3dpa",
            assay = "RNA",
            min.cells = 3, 
            min.features = 100,
            names.field = 2,
            names.delim = "_")
        saveRDS(spine.3dpa, file  =  "./seurat_objects/spine.3dpa.rds")

For limb 14dpa the S4_ and S7_ samples from sc1 matrix will be subset into a Seurat object

    limb.14dpa <- 
          CreateSeuratObject(
            sc1.anno[, grep('^S(4|7)_', colnames(sc1.anno)),], 
            project = "limb_14dpa",
            assay = "RNA",
            min.cells = 3, 
            min.features = 100,
            names.field = 2,
            names.delim = "_")
        saveRDS(limb.14dpa, file  =  "./seurat_objects/limb.14dpa.rds")
        
For spine 14dpa the S9_, S12_, and S14 samples from sc1 matrix will be subset into a Seurat object

    spine.14dpa <- 
          CreateSeuratObject(
            sc1.anno[, grep('^S(9|23|14)_', colnames(sc1.anno)),], 
            project = "spine_14dpa",
            assay = "RNA",
            min.cells = 3, 
            min.features = 100,
            names.field = 2,
            names.delim = "_")
        saveRDS(spine.14dpa, file  =  "./seurat_objects/spine.14dpa.rds")

>sessionInfo()  
R version 4.0.3 (2020-10-10). 
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS
Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     
other attached packages:
[1] Seurat_3.2.2










