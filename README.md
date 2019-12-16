![Logo](https://github.com/elagralinska/APL/blob/master/img_folder/APL_logo.png)

  # APL
  Finding condition-specific genes using the Association Plots.
  
  For more details, please see our [tutorial](https://github.com/elagralinska/associationplots/blob/master/TUTORIAL/Tutorial.md).
  
  
  
  ## Introduction
  APL is an interactive application developed using R package Shiny. After uploading the input matrix and choosing the parameters, the tool generates the **2D- and 3D-representation** of the data, as well as the **Association Plot** for any selected  condition or set of conditions. To identify the **condition-specific genes** one can zoom in the generated Association Plot and click the mouse over the genes presented in the plot. The tool allows also for **enrichment analysis** of the condition-specific genes based on the human or mouse gene ontology.
  
  
  
  ## Installation and usage
  APL can be run directly from R. Before running the tool, make sure that following R packages are installed on your machine:
  - `shiny`
  - `ggplot2`
  - `plotly`
  - `dplyr`
  - `shinycssloaders`
  - `shinyjs`
  - `topGO`
  - `Rgraphviz`
  - `org.Hs.eg.db`
  - `org.Mm.eg.db`
 
  To use APL we recommend first to download our GitHub repository to your working directory and unzip the downloaded file.  Next, to run the tool - open an R session and type following commands:
   
  ```R
      # load shiny library
      library(shiny)
      
      # run the tool (please make sure that the folder "APL-master" is located in your working directory)
      runApp("APL-master")
  ```
  
  
  Alternatively, you can download and launch APL using following commands in R:
 
   ```R
       # load shiny library
       library(shiny)
      
       # run the tool
       runGitHub("APL","elagralinska")  
   ```
 
 
 
  To start the analysis using APL an input matrix in a tab-delimited .txt format has to be provided. The example of such file can be found [here](https://github.com/elagralinska/APL/blob/master/TUTORIAL/input_matrix_FORMAT.txt).
  
  ## Demo
  
  For demonstration purposes an example subset of GTEx data have been provided [here](https://github.com/elagralinska/APL/blob/master/TUTORIAL/input_matrix_DEMO.txt). The provided file contains gene expression values across 22 distinct human tissues, 10 samples per each. To learn more how to analyze such data using APL see our [tutorial](https://github.com/elagralinska/associationplots/blob/master/TUTORIAL/Tutorial.md). 
