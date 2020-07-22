![Logo](https://github.com/elagralinska/APL/blob/master/img_folder/APL_logo.png)


### Introduction

**APL** is an interactive application that can be used for data exploration as well as finding condition-specific genes in complex data sets. After uploading the input matrix and choosing the parameters, the tool generates the 2D- and 3D-representation of the data, as well as the Association Plot for any selected condition or set of conditions. To identify the condition-specific genes one can zoom in the generated Association Plot and click the mouse over the genes presented in the plot. The tool allows also for enrichment analysis of the condition-specific genes for Gene Ontology terms (based on R package `topGO`). 

Below we explain step-by-step how to use the tool.

### Table of contents:

  + [1.  Installation](#1--installation)
    + [2.  First steps](#2--first-steps)
      - [Data upload](#data-upload)
      - [Number of genes](#number-of-genes)
      - [Number of dimensions](#number-of-dimensions)
    + [3.  2D- and 3D-representation of the data](#3--2d--and-3d-representation-of-the-data)
    + [4.  Association Plots](#4--association-plots)
    + [5.  Gene enrichment analysis](#5--gene-enrichment-analysis)
    + [6.  Final remarks](#6--final-remarks)


### 1.  Installation

APL can be run directly from R. Before running the tool, make sure that following R packages are installed on your machine: `shiny`, `reticulate`, `ggplot2`, `plotly`, `dplyr`, `shinycssloaders`, `shinyjs`, `topGO`, `Rgraphviz`, `org.Hs.eg.db`, `org.Mm.eg.db`. To download and launch APL we recommend using following commands in R session:

 ```R
     # load shiny library
     library(shiny)
   
     # run the tool
     runGitHub("APL","elagralinska")  
  ```

For further alternatives please see the [documentation](https://github.com/elagralinska/associationplots/blob/master/README.md).

After launching the tool you will see the following window:

![Screenshot1](screenshots/Screenshot_1.png)

### 2.  First steps

#### Data upload
To demonstrate you how to use the tool, we applied it to an example data set, which is provided [here](https://github.com/elagralinska/APL/blob/master/TUTORIAL/input_matrix_DEMO.txt). The data constitutes a subset of GTEx data, containing gene expression values across multiple human tissues. The provided file contains in total information on 22 tissues, 10 samples per each. 

To analyze the data using our tool please save the file on your machine and then upload it in the application in the field `"Gene expression data input"`.


#### Number of genes
The factor, which greatly contributes to the APL computation time, is the size of the input matrix. Reducing its size by discarding genes with the lowest variance across conditions (and thus, the ones which are not likely to be highly condition-specific) would significantly shorten the computation time. 

In our program the user can choose the number of genes with the highest variance, which should be retained for further analysis. To do this, after uploading the input matrix a notification with two gene numbers will be displayed. The first value indicates the total number of genes, and the second one - number of genes with a non-zero variance. To choose the number of genes for further analysis, please use the field `"Number of genes"`. The selected number should be smaller or equal to the second displayed value. The default parameter is 1000.

For the analysis of our example data `5000` genes were chosen.

#### Number of dimensions
The step of reducing the size of the input matrix is followed by dimensionality reduction of the data. Reduced number of dimensions not only eliminates noise in the data but can also significantly shorten the computation time. To choose the number of dimensions, which will be further used for generating the Association Plots, please use the parameter `"Number of dimensions"`.

To give the user an impression on how many dimensions should be kept, four methods from the literature have been implemented in the program. The short description of the methods as well as the resulting numbers of dimensions are displayed in the tab `"Suggested number of dimensions"`. Below we explain each of them.

![Screenshot2](screenshots/Screenshot_2.png)



The implemented methods:
1.  *”Average rule“* - retains those dimensions that explain more inertia than one dimension on average. Such average value is calculated by dividing 100% by the total number of columns from the input matrix.
2.  *”80% rule“* - retains the minimum number of dimensions, which in total account for more than 80% of total variance in the data.
3.  *”Elbow rule“* - is based on scree plots of randomized data. The dimension is read off from the point where the actual
scree plot of the data enters the band of randomized eigenvalues.
4.  *"Scree plot of the real data"* - decision can be made based on the shape of the scree plot by reading the number of dimensions corresponding to the jump in the plot.


In the example analysis presented in this tutorial the first three above listed methods gave following results: `31`, `21`, `25` dimensions. At the end `25` dimensions (number obtained using the "elbow rule") were chosen.
Now, to continue with the analysis, please click on a blue button `"Start calculations"`.

### 3.  2D- and 3D-representation of the data
Based on the input matrix provided by the user APL generates 2D- and 3D- representation of the data. The generated plots are interactive and will be displayed in the tabs `"2D plot"` and `"3D plot"`. Below you can see the screenshot of a 2D plot generated for demo data.

![Screenshot3](screenshots/Screenshot_3.png)

As it was shown in the screenshot, below the 2D plot a barplot with expression values across all of the conditions for gene `GH2` is displayed. To generate such barplots, click a mouse over a gene (blue dot) in the 2D plot. 
One can also color the certain genes in all plots generated in APL (including 2D- and 3D-plots). To do this, please type the names of genes to color in the field `Genes to highlight`. (*Important: the gene names have to be in the same form as in the input matrix. When highlighting multiple genes, please use commas to separate the names, e.g. `KLF,GH2,LMN1`.*) The selected genes will be marked in red (example: gene GH2 in the screenshot above). 


### 4.  Association Plots
The next plot, which can by generated by APL, is the so-called Association Plot - a planar representation of gene-condition relationships in high-dimensional data. Association Plots are intended for finding and visualizing genes specific for a selected condition (or a group of conditions) from the input data. To generate such plot we have to first specify on which conditions we would like to focus. This can be done in two different ways:

- **by providing IDs of conditions manually**

To do this, please use the column IDs that are indicated in the main panel in the tab `"Which conditions?"` and type them in the corresponding field located in the side panel (*Important: please use only comma-separated numbers, without spaces between them, e.g. `2,3,4`*). After selecting the IDs of the cluster or of the conditions the Association Plot will be displayed in the tab `"Association Plot"`.

- **by uploading a cluster annotation file**

When working with single-cell transcriptomics data (or data with a very high number of conditions), clustering of cells (conditions) will typically be part of data preprocessing. Thus, to know which genes are characteristic for a given cluster, APL provides a possibility to upload a cluster annotation file instead of providing the condition IDs manually. The uploaded file should consists of two columns - the first one with the condition ID and the second one with the cluster ID. (*Important: condition IDs in the cluster annotation file can be different from IDs in the input expression matrix. However, in both files the condition IDs should appear in the same order, which means that the first row of the cluster annotation file corresponds to the first condition from the input expression matrix, the second - to the second, etc.)*. For instance, when processing single-cell data using the R package `Seurat`, such file can be obtained using the function `Idents`:

  ```R
      # create table with cells and corresponding cluster IDs
      cluster_ids=Idents(nameofyourdata)

      
      # save the table on your machine
      write.table(cluster_ids, "cluster_annotation_file.txt", sep=",",col.names=F)
  ```

For further details please see an [example](https://github.com/elagralinska/APL/blob/master/TUTORIAL/input_clusterAnnotation_FORMAT.txt) of such file. 


In APL, after uploading the cluster annotation file a selection box with all cluster IDs will appear in the side panel. In order to generate the Association Plot for conditions that belong to a given cluster, choose a corresponding cluster IDs in the selection box. The Association Plot will appear in the tab `"Association Plot"`. 



The Association Plots can be used for finding condition-specific genes. The interpretation of such plots is simple: The farther out to the right a genes lies in the plot, the more associated it is to the given set of conditions. For the demonstration purposes the Association Plot in this tutorial was generated for liver samples (condition IDs: `212,213,214,215,216,217,218,219,220,221` or cluster ID: `22` from [this](https://github.com/elagralinska/APL/blob/master/TUTORIAL/input_clusterAnnotation_DEMO.txt) cluster annotation file) and can be seen below in form of red crosses. The green crosses show the location of all other conditions in the plot (use the select box below the Association Plot to show or hide certain conditions). Additionally, genes known from the literature to be highly expressed in liver are colored in red (by typing the gene names `F9,CFHR2,HPX,SPP2,MBL2,FGA,HAO1,APOH,F2,C9` in the field `Genes to highlight`). As expected, all of them are located in the right tail of the plot.

To obtain the bar plot of gene expression values across all conditions, click a mouse over a given gene. The bar plot will show up below the Association Plot. 


![Screenshot4](screenshots/Screenshot_4.png)

Additionally, below the Association Plot a download button is located. Use it in order to download the list of the genes ranked according to their specificity for the given subset of conditions.

### 5.  Gene enrichment analysis
APL allows also for enrichment analysis of the condition-specific genes for Gene Ontology terms. This part of the program is based on `topGO` - an R package that facilitates such type of analyses.

To conduct gene enrichment analysis in APL first make sure that in the field `"Which conditions?"` the correct condition IDs are chosen or that the cluster annotation file is uploaded and the correct cluster IDs is selected. Next, in the field `"Annotation file"` in the side panel choose the correct gene ontology file (mouse or human) or provide the name of one of the Bioconductor AnnotationData Packages for another organism:
- Anopheles `org.Ag.eg.db`
- Bovine `org.Bt.eg.db`
- Worm `org.Ce.eg.db`
- Canine `org.Cf.eg.db`
- Fly `org.Dm.eg.db`
- Zebrafish `org.Dr.eg.db`
- E coli strain K12 `org.EcK12.eg.db`
- Chicken `org.Gg.eg.db`
- Rhesus `org.Mmu.eg.db`
- Malaria `org.Pf.plasmo.db`
- Chimp `org.Pt.eg.db`
- Rat `org.Rn.eg.db`
- Yeast `org.Sc.sgd.db`
- Pig `org.Ss.eg.db`
- Xenopus `org.Xl.eg.db`

Check [here](http://bioconductor.org/packages/3.10/data/annotation/) for the most recent list of existing annotation packages. *(Important: when using a new annotation package, please make sure that it is installed on your machine before launching APL.)* 

Finally, in the selection box `"Gene ontology"` you can choose between the type of the gene ontology (BP - biological processes, CC - cellular components, MF - molecular function).

After choosing the correct parameters the GO enrichment test will be conducted using the Kolmogorov-Smirnov test and the results will be displayed in the tab `"GSEA plot"` in two forms:

- **tabular** - The generated table will contain top ten biological processes enriched among top 1000 condition-specific genes obtained from the Association Plot (or if the total number of genes used for generating the Association Plot was lower than 1000 - 25% of number of genes used for this analysis).

- **graphical** - Another way of understanding the gene enrichment results is to investigate the distribution of the enriched GO terms over the entire GO graph. In APL, we use a function from `topGO` that allows on generating such subgraph. After running the gene enrichment test, the resulting subgraph will be displayed below the table. The nodes of the graph represent different GO terms: the significant ones are presented in form of rectangles, whereas the others in form of ellipses. Moreover, depending on the relative significance of GO terms the nodes have different colors: dark red relates to the highest significance and light yellow - to the lowest. Additionally, the most recent GO subgraph is automatically saved as a pdf file in the working directory under the name `"GSEA_classic_10_def.pdf"`.
For more details please see the documentation of the `topGO` package.


For the demonstration purpose, the enrichment analysis in this tutorial was conducted for liver samples (IDs: `212,213,214,215,216,217,218,219,220,221` or cluster ID: `22` from [this](https://github.com/elagralinska/APL/blob/master/TUTORIAL/input_clusterAnnotation_DEMO.txt) cluster annotation file) using the human annotation file. The obtained results are shown below.

![Screenshot5](screenshots/Screenshot_5.png)


### 6.  Final remarks
- To save the state of APL application and get a URL, which will restore the APL session with that state, please use the button `"Bookmarking..."` at the bottom of the side panel.
- The plots in APL are generated using the `plotly` package. To save a plot, use the camera icon in the modebar.


