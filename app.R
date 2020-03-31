# Load libraries
library(ggplot2)
library(plotly)
library(shiny)
library(dplyr)
library(shinycssloaders)
library(Rgraphviz)
library(shinyjs)


# -------------------------------------------------------------------------------------------   
# Define UI for app ({fluidPage...} for bookmarking)

ui <- function(request) {fluidPage(
	
	useShinyjs(),
	# App title ----
	titlePanel(""),
	# Sidebar layout with input and output definitions ----
	sidebarLayout(
		# Sidebar panel for inputs ----
		sidebarPanel(
			h4("Settings", style = "color:red"),
			"Please, upload a gene expression data matrix and choose optimal parameters.", 
			br(),br(),
			fileInput("dataIn", label = h5(strong("Gene expression data input:"))),
			uiOutput("number_of_genes"),
			tags$hr(), 
			numericInput("num_genes", h5(strong("Number of genes:")), value = 1000),
			numericInput("num_dims", h5(strong("Number of dimensions:")), value = 4),
			actionButton("do", "Start calculations",icon("paper-plane"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
			tags$hr(), 
			radioButtons("which_sampl_buttons", h5(strong("Choose conditions by:")), choices = list("providing IDs of conditions" = 1, "uploading cluster annotation file" = 2),selected = 1),
			textInput("which_sampl", h5("Which conditions?"), value = "2,3,4"),

			fileInput("clusterFile", label = h5("Upload cluster annotation file:")),
			selectInput("which_cluster", h5("Which cluster?"), choices = list("Cluster 1" = 1, "Cluster 2" = 2, "Cluster 3" = 3,"Cluster 4" = 4), selected = 1),
			br(),br(),
			tags$hr(),
			textInput("which_genes", h5(strong("Genes to highlight:")), value = ""),
			tags$hr(), 
			radioButtons("radio", h5(strong("Annotation file:")), choices = list("Human" = 1, "Mouse" = 2,"Another annotation (please provide the name of the annotation package):" = 3),selected = 3),
			textInput("annotIn", label = NULL),
			selectInput("which_ontology", h5(strong("Gene ontology:")), choices = list("BP (biological processes)" = "BP", "CC (cellular components)" = "CC", "MF (molecular function)" = "MF"), selected = 1),			
			actionButton("run_gsea", "Run enrichment analysis",icon("paper-plane"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
			tags$hr(), 
			bookmarkButton()
		),
		# Main panel for displaying outputs ----
		mainPanel(
			fluidRow(
				column(4, br(), img(src = "img/APL_logo.png", height = 76, width = 173), ""),
				column(1),
				column(7,
					div(h4("APL"), style = "color:blue"),
					p("Explore your data and find condition-specific genes."),
					p("For an introduction", "visit the ", a("Github repository", href = "http://github.com/elagralinska/APL"), " and see our"),
					a("tutorial.", href = "https://github.com/elagralinska/APL/blob/master/TUTORIAL/Tutorial.md")
				)
			),
			br(),
			h4("Results", style = "color:red"),
			textOutput("selected_num_dims"),
			tabsetPanel(type = "tabs",
                tabPanel("Suggested number of dimensions", 
					uiOutput("methods") %>% withSpinner(color="#0066FF"),br(),br(),
					plotlyOutput("method4"),br(),br()
				),
                  
                tabPanel("Which conditions?", br(),
					textOutput("sentence_conditions") %>% withSpinner(color="#0066FF"), br(),
					tableOutput("sample_ids")  
                ),
                  
                tabPanel("2D plot", br(),br(),
					plotlyOutput("plot2D") %>% withSpinner(color="#0066FF"),
					br(),br(),br(),br(),
					div(verbatimTextOutput("selection2D"), align = "center"),
					br(),
					plotlyOutput("barplot2D"),
					br(), br()
				),
                  
                tabPanel("3D plot", br(),br(), 
					plotlyOutput("plot3D") %>% withSpinner(color="#0066FF"),
					br(),br(),br(),br()
				),

                tabPanel("Association Plot", br(),br(), 
					plotlyOutput("plotsincos") %>% withSpinner(color="#0066FF"),
					br(),br(),
					selectInput("which_conditions_to_highlight", h5("Show:"), choices = list("only selected conditions" = 1, "all conditions" = 2,"no conditions" = 3), selected = 1),
					br(),br(),
					downloadButton("downloadRanking_assocPlot", "Download gene ranking"),
					br(),br(),
					div(verbatimTextOutput("selectionsc"), align = "center"),
					br(),
					plotlyOutput("barplotsc"),
					br(), br()
				),

				tabPanel("GSEA", br(),br(), 
					div(verbatimTextOutput("gsea_errors"), align = "center"),
					tableOutput("ranking_GSEA_table") %>% withSpinner(color="#0066FF"),
					plotOutput("ranking_GSEA_graph",width = "1300px", height="1300px")
				)
			)

		)
	)
)}

# -------------------------------------------------------------------------------------------   
# Define server logic required to draw a histogram ----

server <- function(input, output,session) {
	
	# Choosing the condition IDs for the Association Plot - choose between a list of IDs or file with clusters to upload
	observe({
		shinyjs::hide("which_sampl")
		shinyjs::hide("clusterFile")
		shinyjs::hide("which_cluster")

		if(input$which_sampl_buttons == 1)
		shinyjs::show("which_sampl")
		
		if(input$which_sampl_buttons == 2)
		shinyjs::show("clusterFile")
		
		req(input$clusterFile)
		if(input$which_sampl_buttons == 2)
		shinyjs::show("which_cluster")
		
	})
	

  
	# change size limit of input file
	options(shiny.maxRequestSize=5000*1024^2) #5000 - number of Mb,


	# -------------------------------------------------------------------------------------------   
	# conduct first part of CA calculations
	
	# read input matrix	
	in_table1 <- reactive({
		req(input$dataIn)
		in_table <- read.table(input$dataIn$datapath,header=T, sep = "\t", quote="", fill=FALSE)
	})
	
	# calculate number of columns
	col_num <- reactive({
		ncol(in_table1())
	})

	# calculate total number of genes in the input matrix
	total_genes <- reactive({
		nrow(in_table1())
	})
	
	# save names of the columns in form of a vector
	header <- reactive({
		colnames(in_table1())
	})
	
	# create a submatrix - without gene names and only with expression values			  
	data <- reactive({
		data <- unname(as.matrix(in_table1()[,-1]))
	})

	# calculate number of genes in the input matrix, which have non-zero expression in at least one sample
	non_zero_index <- reactive({	  
		  # delete genes (rows) that contain only zeros across all samples
		  non_zero_index <- which(rowSums(data()) > 0)
	})
	
	# calculate the gene ranking according to their variance across samples
	ix_var <- reactive({

		  data_no0 <- data()[non_zero_index(),]
		  expr_data <- data_no0

		  ## calculating the chi-square components ranking (according to this ranking we will choose later the top XXX genes "topvar")
		  # CA calculations (I)
			  
		  P <- expr_data/sum(expr_data)             # proportions matrix
		  p_i_plus <- rowSums(P)                    # row masses
		  p_plus_j <- colSums(P)                    # column masses
			  
		  E <- p_i_plus %o% p_plus_j      			# expected proportions
		  Z <-  (P - E) / sqrt(E)         			# standardized residuals
		  Z[is.nan(Z)] <- 0
			  
		  chi_data <- sum(expr_data) * (Z ^ 2)		# chi-square components matrix

		  # sort chi-square components matrix by variance (descending)
		  variances <- apply(chi_data,1,var)
		  ix_var <- order(-variances)
		  ix_var
	})
	
	# create a subset of expression matrix - only with selected number of top genes with the highest variance	
	expr_data <- reactive({	  

		  expr_data <- data()[non_zero_index(),]
		  ## reduce the expression matrix to the chosen number of genes
		  expr_data <- expr_data[ix_var()[1:input$num_genes],]
		  expr_data
				  
	})

	# create a list of names of the top genes with the highest variance (corresponds to matrix expr_data())		
	used_gene_names <- reactive({
		  req(input$dataIn)

		  ## create a vector with names of investigated genes (sorted by the calculated variances, descending)
		  all_gene_names <- unname(in_table1()[1])
		  used_gene_names <- all_gene_names[non_zero_index(),]      # these are the gene_names actually used (after purging all-0 rows)
		  used_gene_names <- used_gene_names[ix_var()[1:input$num_genes],]
		  used_gene_names
	})
				
	# unlist used_gene_names		
	used_gene_names_list <- reactive({
		  req(input$dataIn)
		  used_gene_names_list <- unlist(as.list(used_gene_names()))
	})

 
 
	# -------------------------------------------------------------------------------------------   
	## CA calculations (I)
	
	# conduct SVD
	SVD <- reactive({
		req(input$dataIn)
		validate(need(length(non_zero_index()) > input$num_genes , "Too high number of genes selected for further analysis! To continue, reduce the number of genes in the side panel!")) 
		
	    P <- expr_data()/sum(expr_data())         # proportions matrix
		p_i_plus <- rowSums(P)                    # row masses
		p_plus_j <- colSums(P)                    # column masses
			  
		E <- p_i_plus %o% p_plus_j      # expected proportions
		Z <-  (P - E) / sqrt(E)         # standardized residuals
		Z[is.nan(Z)] <- 0
		SVD <- svd(Z)
	})

	# based on conducted SVD - calculate eigenvalues and explained inertia (explained by each dimension, sorted in descending order)
	expl_inertia <- reactive({
		req(input$dataIn)

		ev <- SVD()$d^2							# eigenvalues = singularvalues^2
		expl_inertia <- (ev/sum(ev)) *100		# explained inertia (%)
	})
	
	# calculate maximal number of dimensions in this analysis
	max_num_dims <- reactive({
		req(input$dataIn)
		max_num_dims <- length(SVD()$d)			# maximal number of dimensions
	})	  

	# how much inertia is explained on average by one dimension?
	avg_inertia <- reactive({
		req(input$dataIn)

		avg_inertia <- 100/max_num_dims()		# percentage of inertia explained by 1 dimension (on average)
	})



	# Calculate a sufficient number of dimensions using different methods
	# -------------------------------------------------------------------------------------------   
	# Method 1: Dim's > average inertia
	dim_number_m1 <- reactive({	  
		req(input$dataIn)

	    dim_number_m1 <- sum(expl_inertia() > avg_inertia())  # result: number of dimensions, all of which explain more than avg_inertia
		dim_number_m1
	})
	  
	# -------------------------------------------------------------------------------------------   
	# Method 2: Sum of dim's > 80% of the total inertia
	dim_number_m2 <- reactive({
		req(input$dataIn)

		dim_number_m2 <- min(which(cumsum(expl_inertia())>80)) # the first dimension for which the cumulative sum of inertias (from dim1 up to given dimension) is higher than 80%
		dim_number_m2
	})
	   
	# -------------------------------------------------------------------------------------------   
	# Method 3: Graphical representation of explained inertias (scree plot)
	# the user can set the threshold based on the scree plot
	  
	# the scree plot below method 4

	# -------------------------------------------------------------------------------------------   
	# Method 4: Formalization of the elbow rule
	dim_number_m4 <- reactive({  
		req(input$dataIn)

		  
		# 1.Generate an artificial data matrix by randomly permuting the responses (rows)
		#   of each sample (column)  of the original expression matrix
		# 2.Perform the CA calculations (I) for the generated data
		# 3.Repeat the steps 1. to 3. 3 times and save the calculated
		#   explained inertia vector (art_expl_inertia) for each artificial matrix
		#   as a row in an overall inertia matrix (matrix_art_expl_inertia)
			  
		matrix_art_expl_inertia <- matrix(nrow = 3, ncol = max_num_dims())
			  
		for (k in 1:3) {
				
			art_expr_data <- apply(expr_data(), 2, FUN=sample)

			art_P <- art_expr_data/sum(art_expr_data)
			art_p_i_plus <- rowSums(art_P)                    # row masses
			art_p_plus_j <- colSums(art_P)                    # column masses
				
			art_E <- art_p_i_plus %o% art_p_plus_j            # expected proportions
				
			art_Z <-  (art_P - art_E) / sqrt(art_E)           # standardized residuals
			art_Z[is.nan(art_Z)] <- 0
			art_SVD <- svd(art_Z)
				
			art_sv <- art_SVD$d
			art_ev <- art_sv^2
			art_expl_inertia <- (art_ev/sum(art_ev)) *100
				
			matrix_art_expl_inertia[k,] <- art_expl_inertia
		}
			  
		# 4.Superimpose the scree plot of each artificial matrix to the graph
		#   of the original matrix by plotting the rows of matrix_art_expl_inertia
		#   as additional lines in the barplot
	
		ER_scree_plots <- renderPlotly({
			plot1=plot_ly(type='bar') %>%
			add_trace(x = 1:max_num_dims(), y = expl_inertia(), type = 'bar', marker = list(color = 'rgb(49,130,189)'), name ='original data') %>%
			add_segments(x = 0, xend = max_num_dims(), y = avg_inertia(), yend = avg_inertia(), name = 'average inertia', showlegend = TRUE,line=list(color = 'red')) %>%
			layout(title = 'Scree plot of original and simulated data', xaxis = list(title ='Number of dimensions'), yaxis = list(title ='Explained inertia (%)'), legend = list(x = 0.8, y = 0.9)) 
			# add traces - simulated data
			for (k in 1:3) {
				plot1 <- plot1 %>% add_trace(
					x = 1:max_num_dims(), y = matrix_art_expl_inertia[k,], line = list(dash = 'dash'), 
					mode = 'lines', type = 'scatter', name = paste(k,'. simulated data')
				)
			}
			plot1
		})
	  
		output$method4 = ER_scree_plots #***for gtex data window is too small!!!***
			    
		# 5.Identify the largest number of dimension based on intersection of the real data's scree plot
		#   with the average of simulated scree plots
			  
		avg_art_inertia <- colMeans(matrix_art_expl_inertia)		# average artificial inertia vector

		tmp=as.integer(expl_inertia()>avg_art_inertia)						# check for each position (dimension) if expl_inertia is higher than the avg_artificial_inertia. If yes =1, if no =0
		if (sum(tmp)==0 || sum(tmp)==max_num_dims()){
			dim_number_m4 <- max_num_dims()									# result: if the lines do not intersect, choose max_num_dims
		}else{
			dim_number_m4 <- length(tmp[cumsum(tmp == 0)<1 & tmp!=0])		# result: if the lines intersect at at 1 or more than 1 positions, choose the intersection with the lowest (BEFORE: it was with highest) x-coordinate (:= the first intersection!)
		}
		dim_number_m4
	})


	# -------------------------------------------------------------------------------------------   
	# Output - left panel - print total number of genes and number of non-zero genes

	output$number_of_genes <- renderUI({
		tagList(tags$p(
			HTML('The uploaded file contains expression data of'), total_genes(), 
			HTML('genes. '), length(non_zero_index())-1, 
			HTML('of them are expressed in at least one sample.') 
		))
	})

	# -------------------------------------------------------------------------------------------   
	# Output - first tab - with suggested number of dimensions

	output$methods <- renderUI({
		req(input$dataIn)

		tagList(tags$p(
				
			HTML('
				<div class="panel panel-primary">
				  <div class="panel-heading">Here we propose <b>four methods</b> for selecting the number of dimensions that will be kept for further analysis:</div>
				  <div class="panel-body">
					<ol>
					 <li> <i>”Average rule“</i> - retains those dimensions that explain more inertia than one dimension on average.</li>
					 <li> <i>”80% rule“</i> - retains the minimum number of dimensions, which in total account for more than 80% of total variance in the data.</li>
					 <li> <i>”Elbow rule“</i> - is based on scree plots of randomized data. The dimension is read off from the point where the original scree plot of the data enters the band of ranomized eigenvalues.</li>
					 <li> <i>"Scree plot of the real data"</i> - decision can be made based on the shape of the scree plot.</li>
					</ol>
				  </div>
				</div>'
			),
		  
			HTML('
				<div class="panel panel-primary">
					<div class="panel-heading">These are the <b>results</b>:</div>
					<div class="panel-body">
					  <ol>
					    <li>   ', dim_number_m1(), ' dimensions </li>
					    <li>   ', dim_number_m2(), ' dimensions </li>
					    <li>   ', dim_number_m4(), ' dimensions (see the scree plot below) </li>
					    <li>   You can also decide on your own based on the scree plot displayed below </li> 
					  </ol>
					</div>
				</div>'
			)
		))
	})
	
	# -------------------------------------------------------------------------------------------   
	# Output - second tab  - print sample names and their IDs
	output$sentence_conditions <- renderText({
		req(input$dataIn)
		HTML('Below you will find a full list of conditions from the input data. Please, specify the ID numbers of conditions, for which you would like to find condition-specific genes. ')
	})
 
  	output$sample_ids <- renderTable({
		req(input$dataIn)
		sample_names_idx <- (2:length(header()))
		T <- as.data.frame(header()[sample_names_idx])
		rownames(T) <- sample_names_idx
		colnames(T) <- 'Sample_name'
		T
	}, include.rownames=T)
	
	
	# -------------------------------------------------------------------------------------------   
	# Read cluster annotation file for the Association Plot and gene enrichment analysis
  
	# Uploading file with cluster annotation
	in_table_clusters <- reactive({
		req(input$clusterFile)
		in_table_clusters <- read.table(input$clusterFile$datapath, header = F, sep = ",", fill=FALSE)
	})
	
	# Construct list of all cluster IDs
	clusterids_list <- reactive({
		req(input$clusterFile)
		clusterids = sort(unique(in_table_clusters()[,2]))
		clusterids_list=paste(c("Cluster", clusterids[1]), collapse = " ")
		for (i in 2:length(clusterids)){
			clusterids_list=append(clusterids_list, paste(c("Cluster", clusterids[i]), collapse = " "))
		}
		clusterids_list
	})

	# update the "which_cluster" selection input by adding the IDs of the clusters
	observe({
		updateSelectInput(
			session,
			"which_cluster",
			choices = clusterids_list()
		)
	})
	
	# Construct a list of cluster IDs and cell IDs (or sample IDs) belonging to a given cluster
	clusters_composition <- reactive({
		req(input$clusterFile)
		clusterids = sort(unique(in_table_clusters()[,2]))
		cluster_composition = list(which(in_table_clusters()[,2] == clusterids[1])+1)
		for (j in 2:length(clusterids)){
			cluster_composition[[j]] = which(in_table_clusters()[,2] == clusterids[j])+1
		}
		# add also names to each vector in the list (with cluster ids)
		for (j in 1:length(clusterids)){
			names(cluster_composition)[j] <- paste("Cluster", clusterids[j], sep = " ") 
		}
		cluster_composition
	})
	
	
	
	# -----------------------------------------------------------------------------------------------------------------------------			
	# -----------------------------------------------------------------------------------------------------------------------------
	# Start running the analysis after the number of dimensions and number of genes have been chosen and the button "Start calculations!" has been pressed
	observeEvent(input$do, {


	    # -------------------------------------------------------------------------------------------   
	    # CA calculations (II)

		# calculate principal coordinates of genes (depending on chosen number of dimensions)
		genes <- reactive({
			  P <- expr_data()/sum(expr_data())         # proportions matrix
			  p_i_plus <- rowSums(P)                    # row masses
			  
			  standard_coordinates_rows <- sweep(SVD()$u[,1:input$num_dims], 1, sqrt(p_i_plus), "/")
			  principal_coordinates_rows <- sweep(standard_coordinates_rows, 2, SVD()$d[1:input$num_dims], "*")
			  genes <- t(principal_coordinates_rows)
		})
			  
		# calculate standard coordinates of samples (depending on chosen number of dimensions)
		samples <- reactive({
			  P <- expr_data()/sum(expr_data())             # proportions matrix
			  p_plus_j <- colSums(P)                    # column masses
			  
			  standard_coordinates_cols <- sweep(SVD()$v[,1:input$num_dims], 1, sqrt(p_plus_j), "/")
			  samples <- t(standard_coordinates_cols)
		})


		# -------------------------------------------------------------------------------------------  
		## choose samples that should be considered in further analysis for generating Association Plot (based either on sample IDs or the uploaded annotation file)
  		
		which_samples <- reactive({
			# for providing sample IDs manually
			if(input$which_sampl_buttons == 1){
				which_samples <- input$which_sampl
				which_samples <- as.numeric(unlist(strsplit(which_samples, split=","))) #convert input into a numeric vector
			}
			
			# for uploading the cluster annotation file
			if(input$which_sampl_buttons == 2){
				which_element_in_list <- which(names(clusters_composition()) == input$which_cluster)
				if(length(which_element_in_list)>0){
					which_samples <- clusters_composition()[[which_element_in_list]]
				}else{return(NULL)}
			}
			
			which_samples
		})

		# create vector with IDs of all other samples (= not chosen samples)
		which_samples_not_chosen <- reactive({
		  all_samples <- c(2:length(header()))
		  which_samples_not_chosen <- setdiff(all_samples, which_samples())
		  
		  which_samples_not_chosen
		})
		
		# IDs of the samples
		all_sample_names <- reactive({
			sample_names_idx <- (2:length(header()))
			header()[sample_names_idx]
		})
		
		# names of the chosen samples
		used_sample_names <- reactive({ ###first save the input as a vector
			header()[which_samples()]
		})
	 
		# names of the NOT chosen samples
		NOT_used_sample_names <- reactive({ ###first save the input as a vector
		  header()[which_samples_not_chosen()]
		})
		
		# -------------------------------------------------------------------------------------------
		## choose which marker genes to highlight
	  
		markers <- reactive({
			as.character(unlist(strsplit(input$which_genes, split=',')))
		})

		# get indices of highlighted genes of the used_gene_names_list
		hl_genes <- reactive({
			sort(match(markers(), used_gene_names_list()))
		})



		# -------------------------------------------------------------------------------------------   
		## samples - coordinates in n-dimensional space	
	
		# calculate average coordinates of chosen samples
		avg_sample_coord <- reactive({
			if (is.null(which_samples()))
				return(NULL)
			avg_sample_coord <- rep(4,0)
			for(j in 1:input$num_dims){
				avg_sample_coord[j] <- 0
				for(i in 1:length(which_samples())){
					temp_index <- which_samples()[i]
					avg_sample_coord[j] <- avg_sample_coord[j] + samples()[j,temp_index-1]
				}
				avg_sample_coord[j] <- avg_sample_coord[j] / length(which_samples())
			}  
			avg_sample_coord
			
		})

		# calculate length_vector_sample (the distance from the origin to the "average" sample)
		length_vector_sample <- reactive({
			length_vector_sample <- 0
			for(j in 1:input$num_dims){
				length_vector_sample <- length_vector_sample + avg_sample_coord()[j]^2
			}
			length_vector_sample <- sqrt(length_vector_sample)
		})
	  


		# -------------------------------------------------------------------------------------------   
		## genes and samples - coordinates in Association Plot	

		# calculate coordinates of genes in Association Plot (for chosen samples)
		G_coordinates_sincos <- reactive({
			if (is.null(which_samples()))
				return(NULL)
			G_coordinates_sincos <-  matrix(0,3,ncol(genes()))              # contains coordinates for Association Plot - first row is x, second row is y, third row is gene ranking
			for(gene_nr in 1:ncol(genes())){
				tmp <- 0
				length_vector_gene <- 0
				for(j in 1:input$num_dims){
					tmp <- tmp + avg_sample_coord()[j] * genes()[j,gene_nr]                           # multiplication of two vectors - avg_sample and gene
					length_vector_gene <- length_vector_gene + (genes()[j,gene_nr]^2)
				}
				length_vector_gene <- sqrt(length_vector_gene)
				G_coordinates_sincos[1,gene_nr] <- tmp / length_vector_sample()                                         # coordinate x for given gene
				G_coordinates_sincos[2,gene_nr] <- sqrt((length_vector_gene ^2) - (G_coordinates_sincos[1,gene_nr]^2))  # coordinate y for given gene
				G_coordinates_sincos[3,gene_nr] <- G_coordinates_sincos[1,gene_nr]										# gene rank
			}
			G_coordinates_sincos[3,]=rank(-G_coordinates_sincos[3,]) # create ranking of genes grom the most specific gene (1) to the least specific one (highest number)
			G_coordinates_sincos
		})


		# calculate coordinates of samples in Association Plot
		S_coordinates_sincos <- reactive({
			if (is.null(which_samples()))
				return(NULL)
			S_coordinates_sincos <- matrix(0,2,ncol(samples()))                       #contains coordinates for Association Plot - first row is x, second row is y
			for(sample_nr in 1:ncol(samples())){
				tmp1 <- 0
				length_vector_sampleNr <- 0
				for(j in 1:input$num_dims){
					tmp1 <- tmp1 + avg_sample_coord()[j] * samples()[j,sample_nr]                             #multiplication of two vectors - avg_sample and sample
					length_vector_sampleNr <- length_vector_sampleNr + (samples()[j,sample_nr]^2)
				}
				length_vector_sampleNr <- sqrt(length_vector_sampleNr)
				S_coordinates_sincos[1,sample_nr] <- tmp1 / length_vector_sample()                                               #coordinate x for given sample
				S_coordinates_sincos[2,sample_nr] <- sqrt((length_vector_sampleNr ^2) - (S_coordinates_sincos[1,sample_nr]^2))  #coordinate y for given sample
			}
			S_coordinates_sincos
		})

		

		# -------------------------------------------------------------------------------------------   
		## output file 1 - ranking of the genes from the Association Plot	
	
		tmp_ranking <- reactive({
			tmp_ranking = cbind(used_gene_names(),t(G_coordinates_sincos()))
			head(tmp_ranking)
			tmp_ranking=tmp_ranking[order(tmp_ranking[,4]),c(4,1)]
			colnames(tmp_ranking)=c("Ranking","Gene_name")
			
			tmp_ranking
		})
	  
		# Downloadable csv of the ranking
		output$downloadRanking_assocPlot <- downloadHandler(
			filename = function() {
				paste("Gene_ranking_",input$num_genes,"_genes_",input$num_dims,"_dims.csv", sep = "")
			},
			content = function(file) {
				write.csv(tmp_ranking(), file, row.names = FALSE)
			}
		)
	  

		# -------------------------------------------------------------------------------------------
		# -------------------------------------------------------------------------------------------
		## Plots
	  
		## plot 1 - 2D subspace

		output$plot2D <- renderPlotly({
			plot_ly(type='scatter', source='plot2D', mode='markers') %>%
			add_trace(x = samples()[1,], y = samples()[2,], mode = 'markers+text', 
				text = paste(all_sample_names()), textposition = "left", opacity = 1,
				marker = list(color = ' #990000 ', symbol = 'x', size = 5), name = 'samples', hoverinfo = 'text', type = 'scatter'
			) %>%
			add_trace(x = genes()[1,], y = genes()[2,], mode = 'markers', text = paste(used_gene_names_list()), opacity = 0.7,
				marker = list(color ='#0066FF', symbol = 'circle-open', size = 2.5), name = 'genes', hoverinfo = 'text', type = 'scatter'
			) %>%
			add_trace(x = genes()[1,hl_genes()], y = genes()[2,hl_genes()], mode = 'markers+text',
				text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
				marker = list(symbol = 'circle', color ='#FF0000', size = 5),
				name = 'marked gene(s)', hoverinfo = 'text', type = 'scatter'
			) %>%
			layout(autosize=T, title = '2D CA plot',showlegend = FALSE, xaxis = list(title = 'Dim1'),
				yaxis = list(title = 'Dim2')
			)
		})
	 
		## plot 2 - 3D subspace

		output$plot3D <- renderPlotly({
			plot1=plot_ly() %>%
			add_trace(x = samples()[1,], y = samples()[2,], z = samples()[3,], mode = 'markers+text', text = paste(all_sample_names()),
				textposition = "left", marker = list(size = 2,  color = '#990000', symbol = 'x'), name = 'samples', hoverinfo = 'text', type = 'scatter3d'
			) %>%
			add_trace(x = genes()[1,], y = genes()[2,], z = genes()[3,], mode = 'markers', text = paste(used_gene_names_list()), opacity = 0.7,
				marker = list(size = 1, color ='#0066FF', symbol = 'circle-open'), name = 'genes', hoverinfo = 'text', type = 'scatter3d'
			) %>%
			add_trace(x = genes()[1,hl_genes()], y = genes()[2,hl_genes()], z = genes()[3,hl_genes()], mode = 'markers+text',
				text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
				marker = list(size = 2, symbol = 'circle', color ='#FF0000'),
				name = 'marked gene(s)', hoverinfo = 'text', type = 'scatter3d'
			) %>%
			layout(autosize=T, title = '3D CA plot', showlegend = FALSE, scene = list(xaxis = list(title = 'Dim1'), yaxis = list(title = 'Dim2'), zaxis = list(title = 'Dim3')))
		})

		## plot 3 - barplot for 2D plot

		output$selection2D <- renderPrint({
			req(input$dataIn)
			s <- event_data("plotly_click", source = 'plot2D')	
			if (length(s) == 0){
				cat("Click on a gene in the 2D CA plot to display the GE levels for each sample.")
			}else if(s$curveNumber == 1){
				cat("Please try again and click on a gene in the 2D CA plot to display the GE levels for each sample.")
			}
		})
		  
		output$barplot2D <- renderPlotly({
			s <- event_data("plotly_click", source = 'plot2D')
			as.list(s)
			if (length(s)) {
				if (s$curveNumber == 2 | s$curveNumber == 3){
					varsi <- c(s$x, s$y)
					get_gene_idx_x <- grep(varsi[1], genes()[1,])
					get_gene_idx_y <- grep(varsi[2], genes()[2,])
					get_gene_idx <- intersect(get_gene_idx_x, get_gene_idx_y)
					bar_x <- expr_data()[get_gene_idx,]
					if(length(get_gene_idx)>0){
						plot_ly(x = bar_x, y = all_sample_names(), type = 'bar', orientation = 'h', height=length(header())*12) %>%
						layout(title = used_gene_names_list()[get_gene_idx],
						yaxis = list(categoryorder = "array",categoryarray = all_sample_names(), tickmode = 'linear'))
					}
				} else {
					plotly_empty(type = 'scatter', mode = 'none')
				}
			}
		})
		  
		## plot 4 - Association plot
		
		output$plotsincos <- renderPlotly({
			
		  associationplot <- plot_ly(type = 'scatter', source = "plotsincos", mode = 'markers') %>%
#		    add_trace(x =  S_coordinates_sincos()[1,which_samples()-1], y =  S_coordinates_sincos()[2,which_samples()-1], mode = 'markers',
#		              text = paste(used_sample_names()), textposition = "left",
#		              marker = list(color = '#990000', symbol = 'x', size = 5), name = 'samples', hoverinfo = 'text', type = 'scatter'
#		    ) %>%
		    # add empty trace so that clicking on genes and generating of a corresponding bar plot works
		    add_trace(x =  S_coordinates_sincos()[1,1]-S_coordinates_sincos()[1,1], y =  S_coordinates_sincos()[2,1]-S_coordinates_sincos()[2,1],
		              marker = list(size = 1, color = "black")
		    ) %>%
		    # add blue dots - all genes
		    add_trace(x = G_coordinates_sincos()[1,], y = G_coordinates_sincos()[2,], mode = 'markers', text = paste(used_gene_names_list()),
		              opacity = 0.7, marker = list(color ='#0066FF', symbol = 'circle-open', size = 2.5), name = 'genes', hoverinfo = 'text', type = 'scatter'
		    ) %>%
		    # add blue dots - highlighted genes
		    add_trace(x = G_coordinates_sincos()[1,hl_genes()], y = G_coordinates_sincos()[2,hl_genes()], mode = 'markers+text',
		              text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
		              marker = list(symbol = 'circle', color = '#FF0000', size = 5),
		              name = 'marked genes', hoverinfo = 'text', type = 'scatter'
		    ) %>%
		    layout(title = paste('Association Plot \n', input$num_dims,' first dimensions, condition IDs: ', paste(which_samples(), collapse = ','),'\n'),
		           xaxis = list(title = 'Distance from origin (x)', rangemode = "tozero"),
		           yaxis = list(title = 'Distance from gene to its projection (y)', rangemode = "tozero"),showlegend = FALSE
		    )
			
			if(input$which_conditions_to_highlight == 1){
			  
			  # add red crosses - chosen samples
			    associationplot <- add_trace(associationplot, x =  S_coordinates_sincos()[1,which_samples()-1], y =  S_coordinates_sincos()[2,which_samples()-1], mode = 'markers',
			    text = paste(used_sample_names()), textposition = "left",
			    marker = list(color = '#e60000', symbol = 'x', size = 6), name = 'samples', hoverinfo = 'text', type = 'scatter')
			                          
			}else{
			    if(input$which_conditions_to_highlight == 2){
			  
			      # add green crosses - all other samples
			       associationplot <- add_trace(associationplot, x =  S_coordinates_sincos()[1,which_samples_not_chosen()-1], y =  S_coordinates_sincos()[2,which_samples_not_chosen()-1], mode = 'markers',
			       text = paste(NOT_used_sample_names()), textposition = "left",
			       marker = list(color = '#008000', symbol = 'x', size = 6), name = 'samples', hoverinfo = 'text', type = 'scatter')

			       # add red crosses - chosen samples
			       associationplot <- add_trace(associationplot, x =  S_coordinates_sincos()[1,which_samples()-1], y =  S_coordinates_sincos()[2,which_samples()-1], mode = 'markers',
			       text = paste(used_sample_names()), textposition = "left",
			       marker = list(color = '#e60000', symbol = 'x', size = 6), name = 'samples', hoverinfo = 'text', type = 'scatter')
			    
			    }else{
			       # don't show any samples in the Association Plot (option 3)
			       associationplot
			    }
			}
		})


		## plot 5 - barplot for Association plot


		output$selectionsc <- renderPrint({
			t <- event_data("plotly_click", source = "plotsincos")
			if (length(t) == 0){
				cat("Click on a gene in the Association Plot to display the GE levels for each sample.")
			}else if(t$curveNumber == 1){
				cat("Please try again and click on a gene in the Association Plot to display the GE levels for each sample.")
			}
		})
		  
		  
		output$barplotsc <- renderPlotly({
			t <- event_data("plotly_click", source = "plotsincos")
			as.list(t)
			if (length(t)) {
				if (t$curveNumber == 2 | t$curveNumber == 3){
					varsi <- c(t$x, t$y)
					get_gene_idx_x <- grep(varsi[1], G_coordinates_sincos()[1,])
					get_gene_idx_y <- grep(varsi[2], G_coordinates_sincos()[2,])
					get_gene_idx <- intersect(get_gene_idx_x, get_gene_idx_y)
					bar_x <- expr_data()[get_gene_idx,]
					if(length(get_gene_idx)>0){
						plot_ly(x = bar_x, y = all_sample_names(), type = 'bar', orientation = 'h', height=length(header())*12) %>%
						layout(title = used_gene_names_list()[get_gene_idx],
						yaxis = list(categoryorder = "array",categoryarray = all_sample_names(), tickmode = 'linear'))
					}
				} else {
					plotly_empty(type = 'scatter', mode = 'none')
				}
			}
		})
		
 
		# -------------------------------------------------------------------------------------------   
		# -------------------------------------------------------------------------------------------   
		## Gene set enrichment analysis (GSEA) using topGO package

		observeEvent(input$run_gsea, {		  

			# choose annotation file
			if(input$radio == 1){
				output$gsea_errors <- renderPrint({
					paste0("Analysis conducted using human annotation file.")
				})
				chosen_annot_file <- "org.Hs.eg.db"	 # Human
			}else{ 
				if(input$radio == 2){
					output$gsea_errors <- renderPrint({
						paste0("Analysis conducted using mouse annotation file.")
					})			
				chosen_annot_file <- "org.Mm.eg.db"	 # Mouse
				}else{
					if(input$radio == 3){
						chosen_annot_file <- paste(input$annotIn)	# using annotation package for another organism
					}
				}
			}
			
			# check if the chosen annotation is installed and if its name is correct
			output$gsea_errors <- renderPrint({		
				if(!require(chosen_annot_file, character.only = TRUE, quietly=T))
				paste0("The package ", chosen_annot_file, " is either not installed or not available. Please install it or choose another annotation.")
			})							




			if(require(chosen_annot_file, character.only = TRUE, quietly = T)){
				
				# create a ranking of the condition-specific genes - based on the current Association Plot
				tmp_ranking = cbind(isolate(used_gene_names()),t(isolate(G_coordinates_sincos())))
				head(tmp_ranking)
				tmp_ranking=tmp_ranking[order(tmp_ranking[,4]),c(4,1)]
				colnames(tmp_ranking)=c("Ranking","Gene_name")
				# select the number (1000) of top ranked genes for the enrichment analysis
				nr_genes_gsea = 1000
				# if the number of genes used for the analysis is smaller than 1000 genes then take the floor of 25% of selected number of genes
				if (isolate(input$num_genes) < nr_genes_gsea){
					nr_genes_gsea = floor(0.25 * isolate(input$num_genes))
				}
				# create a vector of ranking values (1:1000) and the names of the genes for topGO analysis
				rankedGenes <- c(1:nr_genes_gsea)
				names(rankedGenes) <- tmp_ranking[1:nr_genes_gsea,2]
					
				# check if package topGO is installed and load it
				require(topGO)
				
				# specify gene ontology type (BP, CC, MF)
				selected_ontology <- input$which_ontology
				
				# conduct gene enrichment analysis using topGO functions (Kolmogorov-Smirnov test)
				selection <- function(x) TRUE
				allGO2genes <- annFUN.org(whichOnto=selected_ontology, feasibleGenes=NULL, mapping= paste(chosen_annot_file), ID="symbol")
				GOdata <- new("topGOdata", ontology=selected_ontology, allGenes=rankedGenes, annot=annFUN.GO2genes, GO2genes=allGO2genes, geneSel=selection, nodeSize=5)
				results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
			}

			# generate the output table with top 10 enriched GO terms
			output$ranking_GSEA_table <- renderTable({
				if(require(chosen_annot_file, character.only = TRUE, quietly = T)){
				
					goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=10)

					#goEnrichment <- goEnrichment[goEnrichment$KS < 0.05,]
					goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]

					## print 10 top GO results enriched in the Association Plot	
					goEnrichment		
				}
			})
			
			# generate the output subgraph with top 10 enriched GO terms and neighboring terms
			output$ranking_GSEA_graph <- renderPlot({
				if(require(chosen_annot_file, character.only = TRUE, quietly = T)){

					showSigOfNodes(GOdata, score(results.ks), firstSigNodes = 10, useInfo = 'def')
					printGraph(GOdata, results.ks, firstSigNodes = 10, fn.prefix = "GSEA", useInfo = "def", pdfSW = TRUE)

				}
			})			
			

			
					
		})
	



	}) # end of observe event (button "Start calculations!)
}

# for showing the APL logo
addResourcePath('img','img_folder')

# enable bookmarking
enableBookmarking(store = "url")

shinyApp(ui = ui, server = server)
