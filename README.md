***LOGO***

  # ***NAME OF THE TOOL***
  Finding condition-specific genes using the Association Plots.
  For more details, please see our tutorial.
  
  
  
  ## Introduction
  Blablabla is an interactive application developed using R package Shiny. After uploading the input matrix and choosing the parameters, the tool generates the **2D- and 3D-representation** of the data, as well as the **Association Plot** for any selected  condition or set of conditions. To identify the **condition-specific genes** one can zoom in the generated Association Plot and click the mouse over the genes presented in the plot. The tool allows also for **enrichment analysis** of the condition-specific genes based on the human or mouse gene ontology.
  
  
  
  ## Installation and usage
  Blablabla can be run directly from RStudio. Before running the tool, make sure that following R packages are installed on your machine:
  - `shiny`
  - `ggplot2`
  - `plotly`
  - `dplyr`
  - `shinycssloaders`
  
  
  
  To run Blablabla open **RStudio** and type following commands:
 ```R
      # load shiny library
      library(shiny)
      
      # run the tool
      runGitHub("associationplots","elagralinska")  
 ```
 
  
  Alternatively, you can first download the tool from our GitHub repository and then run the tool in RStudio. To do this, first create a new directory **"Blablabla"** on your machine and save all the files and directories from our GitHub repository in the newly created directory. Next, open RStudio and set **"Blablabla"** as your working directory. Now you can run the tool using following commands:
  
  ```R
      library(shiny)
      runApp("blablabla")
  ```
  
  ## Demo
  
  
