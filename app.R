# Load libraries
library(ggplot2)
library(plotly)
library(shiny)
library(dplyr)
library(shinycssloaders)

# -------------------------------------------------------------------------------------------   
# Define UI for app ({fluidPage...} for bookmarking)

ui <- function(request) {fluidPage(

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
			actionButton("do", "Start calculations!",icon("paper-plane"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
			tags$hr(), 
			textInput("which_sampl", h5(strong("Which conditions?")), value = "2,3,4"),
			textInput("which_genes", h5(strong("Genes to highlight:")), value = ""),
			tags$hr(), 
			radioButtons("radio", h5(strong("Annotation file:")), choices = list("Human" = 1, "Mouse" = 2,"User's own file (please upload):" = 3),selected = 3),
			fileInput("annotIn", label = NULL),
			br(),br(),
			bookmarkButton()
		),
		# Main panel for displaying outputs ----
		mainPanel(
			fluidRow(
				column(4, br(), img(src = "img/APL_logo.png", height = 76, width = 173), ""),
				column(1),
				column(7,
					div(h4("APL"), style = "color:blue"),
					p("Finding condition-specific genes in your data."),
					p("For an introduction", "visit the ", a("Github repository", href = "http://github.com/elagralinska/APL"), " and the"),
					a("tutorial.", href = "https://github.com/elagralinska/APL/blob/master/TUTORIAL/Tutorial.md")
				)
			),
			br(),
			h4("Results", style = "color:red"),
			textOutput("selected_num_dims"),
			tabsetPanel(type = "tabs",
                tabPanel("Suggested number of dimensions", 

                uiOutput("methods"),br(),br(),
                plotlyOutput("method4"),br(),br()),
                  
                tabPanel("Which conditions?", br(),
					textOutput("sentence_conditions"), br(),
					tableOutput("sample_ids")  %>% withSpinner(color="#0066FF")
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

                tabPanel("Association plot", br(),br(), 
					plotlyOutput("plotsincos") %>% withSpinner(color="#0066FF"),
					br(),br(),
					downloadButton("downloadRanking_assocPlot", "Download gene ranking"),
					br(),br(),
					div(verbatimTextOutput("selectionsc"), align = "center"),
					br(),
					plotlyOutput("barplotsc"),
					br(), br()
				),

				tabPanel("GSEA plot", br(),br(), 
					plotlyOutput("plotsincos_gsea") ,
					plotlyOutput("plotsincos_gsea_legend") ,
					tableOutput("ranking_GSEA_ten") %>% withSpinner(color="#0066FF")
				)
			)

		)
	)
)}

# -------------------------------------------------------------------------------------------   
# Define server logic required to draw a histogram ----

server <- function(input, output,session) {
	
  
	# change size limit of input file
	options(shiny.maxRequestSize=5000*1024^2) #5000 - number of Mb,


	dist3 <- reactive({
		input$num_dims+7
	})

	# -------------------------------------------------------------------------------------------   
	# conduct first part of CA calculations
		
	in_table1 <- reactive({
		req(input$dataIn)
		in_table <- read.table(input$dataIn$datapath,header=T, sep = "\t", quote="", fill=FALSE)
	})
	

	col_num <- reactive({
		ncol(in_table1())
	})

	total_genes <- reactive({
		nrow(in_table1())
	})
	
	header <- reactive({
		colnames(in_table1())
	})
				  
	data <- reactive({
		data <- unname(as.matrix(in_table1()[,-1]))
	})


	non_zero_index <- reactive({	  
		  # delete genes (rows) that contain only zeros across all samples
		  non_zero_index <- which(rowSums(data()) > 0)
	})

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
	
			
	expr_data <- reactive({	  

		  expr_data <- data()[non_zero_index(),]
		  ## reduce the expression matrix to the chosen number of genes
		  expr_data <- expr_data[ix_var()[1:input$num_genes],]
		  expr_data
				  
	})

			
	used_gene_names <- reactive({
		  req(input$dataIn)

		  ## create a vector with names of investigated genes (sorted by the calculated variances, descending)
		  all_gene_names <- unname(in_table1()[1])
		  used_gene_names <- all_gene_names[non_zero_index(),]      # these are the gene_names actually used (after purging all-0 rows)
		  used_gene_names <- used_gene_names[ix_var()[1:input$num_genes],]
		  used_gene_names
	})
				
				
	used_gene_names_list <- reactive({
		  req(input$dataIn)
		  used_gene_names_list <- unlist(as.list(used_gene_names()))
	})

 
 
	# -------------------------------------------------------------------------------------------   
	## CA calculations (I)
	SVD <- reactive({
		req(input$dataIn)

	    P <- expr_data()/sum(expr_data())         # proportions matrix
		p_i_plus <- rowSums(P)                    # row masses
		p_plus_j <- colSums(P)                    # column masses
			  
		E <- p_i_plus %o% p_plus_j      # expected proportions
		Z <-  (P - E) / sqrt(E)         # standardized residuals
		Z[is.nan(Z)] <- 0
		SVD <- svd(Z)
	})


	expl_inertia <- reactive({
		req(input$dataIn)

		ev <- SVD()$d^2							# eigenvalues = singularvalues^2
		expl_inertia <- (ev/sum(ev)) *100		# explained inertia (%)
	})

	max_num_dims <- reactive({
		req(input$dataIn)
		max_num_dims <- length(SVD()$d)			# maximal number of dimensions
	})	  

	avg_inertia <- reactive({
		req(input$dataIn)

		avg_inertia <- 100/max_num_dims()		# percentage of inertia explained by 1 dimension (on average)
	})



	
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
	# Output - left panel - total number of genes

	output$number_of_genes <- renderUI({
		tagList(tags$p(
			HTML('The uploaded file contains expression data of'), total_genes(), 
			HTML('genes. '), length(non_zero_index()), 
			HTML('of them are expressed in at least one sample.') 
		))
	})

	# -------------------------------------------------------------------------------------------   
	# Output - first tab

	output$methods <- renderUI({
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
	# Output - second tab (sample names and IDs)
	output$sentence_conditions <- renderText({
		req(input$dataIn)
		HTML('Below you will find a full list of conditions from input data. Please, specify the ID numbers of conditions, for which you would like to find condition-specific genes. ')
	})
 
  	output$sample_ids <- renderTable({
		req(input$dataIn)
		sample_names_idx <- (2:length(header()))
		T <- as.data.frame(header()[sample_names_idx])
		rownames(T) <- sample_names_idx
		colnames(T) <- 'Sample_name'
		T
	}, include.rownames=T)
	
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
		## choose samples that should be considered in further analysis for generating Association Plot
  		
		which_samples <- reactive({
			which_samples <- input$which_sampl
			as.numeric(unlist(strsplit(which_samples, split=","))) #convert input into a numeric vector
		})


		all_sample_names <- reactive({
			sample_names_idx <- (2:length(header()))
			header()[sample_names_idx]
		})
	  
		used_sample_names <- reactive({ ###first save the input as a vector
			header()[which_samples()]
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
		## genes and samples - coordinates in association plot	

		# calculate coordinates of genes in association plot (for chosen samples)
		G_coordinates_sincos <- reactive({
			G_coordinates_sincos <-  matrix(0,3,ncol(genes()))              # contains coordinates for association plot - first row is x, second row is y, third row is gene ranking
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



		# calculate coordinates of samples in association plot (for chosen samples)
		S_coordinates_sincos <- reactive({
			S_coordinates_sincos <- matrix(0,2,ncol(samples()))                       #contains coordinates for association plot - first row is x, second row is y
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
		## output file 1 - ranking of the genes from the association plot	
	
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
				marker = list(color = ' #990000 ', symbol = 'x', size = 5), name = 'samples', hoverinfo = 'x+y+text', type = 'scatter'
			) %>%
			add_trace(x = genes()[1,], y = genes()[2,], mode = 'markers', text = paste(used_gene_names_list()), opacity = 0.7,
				marker = list(color ='#0066FF', symbol = 'circle-open', size = 2.5), name = 'genes', hoverinfo = 'x+y+text', type = 'scatter'
			) %>%
			add_trace(x = genes()[1,hl_genes()], y = genes()[2,hl_genes()], mode = 'markers+text',
				text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
				marker = list(symbol = 'circle', color ='#FF0000', size = 5),
				name = 'marked gene(s)', hoverinfo = 'x+y+text', type = 'scatter'
			) %>%
			layout(autosize=T, title = '2D CA plot',showlegend = FALSE, xaxis = list(title = 'Dim1'),
				yaxis = list(title = 'Dim2')
			)
		})
	 
		## plot 2 - 3D subspace

		output$plot3D <- renderPlotly({
			plot1=plot_ly() %>%
			add_trace(x = samples()[1,], y = samples()[2,], z = samples()[3,], mode = 'markers+text', text = paste(all_sample_names()),
				textposition = "left", marker = list(size = 2,  color = '#990000', symbol = 'x'), name = 'samples', hoverinfo = 'x+y+text', type = 'scatter3d'
			) %>%
			add_trace(x = genes()[1,], y = genes()[2,], z = genes()[3,], mode = 'markers', text = paste(used_gene_names_list()), opacity = 0.7,
				marker = list(size = 1, color ='#0066FF', symbol = 'circle-open'), name = 'genes', hoverinfo = 'x+y+text', type = 'scatter3d'
			) %>%
			add_trace(x = genes()[1,hl_genes()], y = genes()[2,hl_genes()], z = genes()[3,hl_genes()], mode = 'markers+text',
				text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
				marker = list(size = 2, symbol = 'circle', color ='#FF0000'),
				name = 'marked gene(s)', hoverinfo = 'x+y+text', type = 'scatter3d'
			) %>%
			layout(autosize=T, title = '3D CA plot', showlegend = FALSE, scene = list(xaxis = list(title = 'Dim1'), yaxis = list(title = 'Dim2'), zaxis = list(title = 'Dim3')))
		})

		## plot 3 - barplot for 2D plot

		output$selection2D <- renderPrint({
			req(input$dataIn)
			s <- event_data("plotly_click", source = 'plot2D')	
			if (length(s) == 0){
				cat("Click on a gene in the 2D CA map to display the GE levels for each sample.")
			}else if(s$curveNumber == 1){
				cat("Please try again and click on a gene in the 2D CA map to display the GE levels for each sample.")
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
					plot_ly(x = bar_x, y = all_sample_names(), type = 'bar', orientation = 'h', height=length(header())*12) %>%
					layout(title = used_gene_names_list()[get_gene_idx],
						yaxis = list(categoryorder = "array",categoryarray = all_sample_names(), tickmode = 'linear')
					)
				} else {
					plotly_empty(type = 'scatter', mode = 'none')
				}
			}
		})
		  

		output$plotsincos <- renderPlotly({
			plot_ly(type = 'scatter', source = "plotsincos", mode = 'markers') %>%
			add_trace(x =  S_coordinates_sincos()[1,which_samples()-1], y =  S_coordinates_sincos()[2,which_samples()-1], mode = 'markers+text',
				#text = paste(used_sample_names()), textposition = "left",
				marker = list(color = '#990000', symbol = 'x', size = 5), name = 'samples', hoverinfo = 'x+y+text', type = 'scatter'
			) %>%
			add_trace(x = G_coordinates_sincos()[1,], y = G_coordinates_sincos()[2,], mode = 'markers', text = paste(used_gene_names_list()),
				opacity = 0.7, marker = list(color ='#0066FF', symbol = 'circle-open', size = 2.5), name = 'genes', hoverinfo = 'x+y+text', type = 'scatter'
			) %>%
			add_trace(x = G_coordinates_sincos()[1,hl_genes()], y = G_coordinates_sincos()[2,hl_genes()], mode = 'markers+text',
				text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
				marker = list(symbol = 'circle', color = '#FF0000', size = 5),
				name = 'marked genes', hoverinfo = 'x+y+text', type = 'scatter'
			) %>%
			layout(title = paste('association plot \n', input$num_dims,' first dimensions, samples: ', paste(which_samples(), collapse = ','),'\n'),
				xaxis = list(title = 'Distance from origin (x)', rangemode = "tozero"),
				yaxis = list(title = 'Distance from gene to its projection (y)', rangemode = "tozero"),showlegend = FALSE
			)
		})
		  

		output$selectionsc <- renderPrint({
			t <- event_data("plotly_click", source = "plotsincos")
			if (length(t) == 0){
				cat("Click on a gene in the association plot to display the GE levels for each sample.")
			}else if(t$curveNumber == 1){
				cat("Please try again and click on a gene in the association plot to display the GE levels for each sample.")
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
					plot_ly(x = bar_x, y = all_sample_names(), type = 'bar', orientation = 'h', height=length(header())*12) %>%
					layout(title = used_gene_names_list()[get_gene_idx],
						yaxis = list(categoryorder = "array",categoryarray = all_sample_names(), tickmode = 'linear')
					)
				} else {
					plotly_empty(type = 'scatter', mode = 'none')
				}
			}
		})
		
 
		# -------------------------------------------------------------------------------------------   
		# -------------------------------------------------------------------------------------------   
		## Gene set enrichment analysis (GSEA)

		## read the annotation file
		GOterms_table <- reactive({
			if(input$radio == 1){
				GOterms_table <- read.table("annotation_files/annotation_human.txt", header = TRUE, fill=FALSE)			
			}else{ 
				if(input$radio == 2){
					GOterms_table <- read.table("annotation_files/annotation_mouse.txt", header = TRUE, fill=FALSE)
				}else{
					if(input$radio == 3){
						req(input$annotIn)
						GOterms_table <- read.table(input$annotIn$datapath, header = TRUE, fill=FALSE)
					}
				}
			}
			colnames(GOterms_table)[1] <- "GO_id"
			GOterms_table <- GOterms_table %>% distinct() # removes duplicated rows 
		})
		
		## choose subset of genes located in the 20-degree-cone in the association plot
		G_coordinates_sincos_GOsubset <- reactive({
			G_coordinates_sincos()[, G_coordinates_sincos()[2,] < (G_coordinates_sincos()[1,] * tan(pi/9))] 
		})

		## extract the names of the selected genes (in the 3rd line used_gene_names changed to used_gene_names_list!)
		used_gene_names_GOsubset_list <- reactive({
			gene_index_GOsubset <- which(G_coordinates_sincos()[2,] < (G_coordinates_sincos()[1,] * tan(pi/9)))
			used_gene_names_GOsubset <- used_gene_names()[gene_index_GOsubset,]      # these are the used_gene_names in association plot in 20-degree-cone
			used_gene_names_GOsubset_list <- unlist(as.list(used_gene_names_GOsubset))
		})

		## match genes from the cone with the GO terms
  
		# create a new table with only 1 row per 1 GO_term
		GOterms_table_unique <- reactive({
			GOterms_table_unique <- GOterms_table()[!duplicated(GOterms_table()[,1]), 1:2]
			colnames(GOterms_table_unique)[1] <- "GO_id"
			GOterms_table_unique
		})
		
		# count the number of genes in each GO_term
		GOterms_size <- reactive({
			GOterms_size <- aggregate(data.frame(count = GOterms_table()[,1]), list(value = GOterms_table()[,1]), length)
			colnames(GOterms_size)[1] <- "GO_id"
			GOterms_size
		})

		# combine those 2 tables by the ID of GO term
		GOterms_summarized <- reactive({
			GOterms_summarized <- merge(GOterms_table_unique(), GOterms_size(), by="GO_id")
		})

		## modify GOterms_table

		# add column with the count number to original table
		GOterms_table2 <- reactive({ 		# (from here "GOterms_table" is called "GOterms_table2")
			merge(GOterms_table(), GOterms_summarized()[,c(1,3)], by="GO_id")
		})
		
		# add 5-th column: with 1 if gene is located within the cone, or -1 otherwise
		GOterms_table3 <- reactive({
			GOterms_table3 <- GOterms_table2()
			output <- numeric(nrow(GOterms_table2()))
			condition <- (GOterms_table2()$gene_name) %in% used_gene_names_GOsubset_list()  # condition check outside the loop
			for (i in 1:nrow(GOterms_table2())) {
				if (condition[i]) {
					output[i] <- 1
				}else{
					output[i] <- -1
				}
			}
			GOterms_table3$is_in_cone <- output
			GOterms_table3
		})

		# add 6-th column: with 1 if gene is within the set of analyzed genes, or 0 otherwise 
		GOterms_table4 <- reactive({
			GOterms_table4 <- GOterms_table3()
			output <- numeric(nrow(GOterms_table3()))
			condition <- (GOterms_table3()$gene_name) %in% used_gene_names_list()  # condition check outside the loop
			for (i in 1:nrow(GOterms_table3())) {
				if (condition[i]) {
					output[i] <- 1
				} else {
					output[i] <- 0
				}
			}
			GOterms_table4$is_in_sincos <- output
			GOterms_table4
		})
 
		# add 7-th column: temporal GO_rank (based on position of the gene in the association ranking) *MAKE FASTER* #maybe filter G_coordinates_sincos for GO genes
		GOterms_table5 <- reactive({
			GOterms_table5 <- GOterms_table4()
			output <- numeric(nrow(GOterms_table4()))
			tmp_gene_id <- match(GOterms_table4()[,3], used_gene_names_list())
			condition1 <- (GOterms_table4()$is_in_sincos) == 1  # condition check outside the loop
			condition2 <- (GOterms_table4()$is_in_cone) == 1  # condition check outside the loop
			for (i in 1:nrow(GOterms_table4())){
				if(condition1[i]){ #if the gene is in the association plot
					if(condition2[i]){ #if the gene is located within the cone
						output[i] <- input$num_genes + 1 - G_coordinates_sincos()[3,tmp_gene_id[i]] #change rank
					}else{ #if the gene is not located within the cone
						output[i] <- G_coordinates_sincos()[3,tmp_gene_id[i]] * (-1) #negative value
					}
				}else{ #if the gene is not in the association plot
					output[i] <- (input$num_genes+1) * (-1)
				}
			}
			GOterms_table5$gene_rank_tmp <- output
			GOterms_table5
		})

		# summing values (number of genes in association plot, ranking value of the GO term) according to GO_id
		tmp2 <- reactive({
			tmp <- GOterms_table5()[,c(1,6,7)]
			tmp2 <- tmp %>% 
			group_by(GO_id) %>% 
			summarise_at(c("is_in_sincos","gene_rank_tmp"),list(sum))
		})	

		# calculate final ranking value
		GOterms_summarized2 <- reactive({
			GOterms_summarized2 <- GOterms_summarized()
			GOterms_summarized2[,4:5] <- tmp2()[2:3] # add 2 columns - number of genes in sincos plot (per GO term) and sum of temporal rank (per GO term)
			GOterms_summarized2[,6] <- GOterms_summarized2[,5] / GOterms_summarized2[,3] # divide the sum of temporal rank (per GO term) by the number of genes in the given GO term
			colnames(GOterms_summarized2)[6] <- "rank_value"
			GOterms_summarized2=GOterms_summarized2[complete.cases(GOterms_summarized2[,6]), ]  #remove rows with NaN as final ranking value (6th column) (== GO terms, which don't have genes in association plot)
			GOterms_summarized2
		})
		
		GO_ranking <- reactive({
			# sort the GO terms descending (from the top enriched GO set to the least enriched GO set)
			GO_ranking <- GOterms_summarized2()[order(-GOterms_summarized2()[,6]),]
			# choose gene terms that contain more than 3 genes
			GO_ranking <- GO_ranking[GO_ranking[,3]>3,]
			GO_ranking
		})
		
		tmp <- reactive({	
			# create table with XXX top GO terms and the names of genes that are in association plot
			GO_ranking_genes <- GOterms_table5()[(GOterms_table5()[,1] %in% GO_ranking()[1:10,1]) & (GOterms_table5()[,6] == 1), c(1:3)]
			# create table: gene name + all GO terms the gene is involved in
			tmp=aggregate(GO_ranking_genes[,2], list(GO_ranking_genes[,3]), function(x) paste0(unique(x)))
		})
		
		GO_ranking_genes2 <- reactive({
			# match the gene name to gene ID
			GO_ranking_genes2 <- match(tmp()[,1], used_gene_names_list())
		})


		## colors for association plot containing enriched GO categories
		# generate separate matrix
		tmpx <- reactive({
			tmpx=cbind(tmp()[,1:2],t(G_coordinates_sincos()[1:2,GO_ranking_genes2()])) 
			colnames(tmpx)=c("gene","GOterm","xc","yc")
			# modify GO terms
			tmpx[,2] <- gsub('c("', "", tmpx[,2], fixed="TRUE")
			tmpx[,2] <- gsub('")', "", tmpx[,2], fixed="TRUE")
			tmpx[,2] <- gsub('", "', " | ", tmpx[,2], fixed="TRUE")
			tmpx
		})

		# -------------------------------------------------------------------------------------------   
		## print 10 top GO results enriched in association-plot	
	
		tmp_ranking_GSEA <- reactive({
			# Suppress warnings related to GSEA plot (Warning message: In RColorBrewer::brewer.pal(N, "Set2") : minimal value for n is 3, returning requested palette with 3 different levels)
			storeWarn<- getOption("warn")
			options(warn = -1)
			# Create list of top 10 enriched terms
			tmp_ranking_GSEA = GO_ranking()[1:10,c(6,1,2)]
			tmp_ranking_GSEA[,1] = c(1:10)
			tmp_ranking_GSEA
		})
	  
		output$ranking_GSEA_ten <- renderTable({
			req(input$dataIn)
			tmp_ranking_GSEA()
		})


		## plot 1 - 2D subspace

		output$plotsincos_gsea <- renderPlotly({
			plot_ly(type='scatter', source='plotsincos_gsea', mode='markers') %>%
			add_trace(x = G_coordinates_sincos()[1,], y = G_coordinates_sincos()[2,], mode = 'markers',text = paste(used_gene_names_list()), opacity = 0.7, textposition = "left",
				marker = list(color ='#0066FF', symbol = 'circle-open', size = 2.5), name = 'genes', hoverinfo = 'x+y+text', type = 'scatter'#, legendgroup = 'group1'
			) %>%
			add_trace(x = G_coordinates_sincos()[1,hl_genes()], y = G_coordinates_sincos()[2,hl_genes()], mode = 'markers+text',
				text = paste(used_gene_names_list()[hl_genes()]), textposition = "left", textfont=list(color='#FF0000'),
				marker = list(symbol = 'circle', color = '#FF0000', size = 5),
				name = 'marked genes', hoverinfo = 'x+y+text', type = 'scatter', showlegend=F
			) %>%
			add_trace(data=tmpx(), x = ~xc, y = ~yc,  type="scatter", color = ~GOterm , mode = 'markers+text',
				text = paste(used_gene_names_list()[GO_ranking_genes2()]), textposition = "left", textfont=list(color='black',size=13), 
				#legendgroup = 'group2',
				marker=list( size=10 , opacity=0.8)
			) %>%
			add_segments(x = 0, xend = max(G_coordinates_sincos()[1,])+0.5, y = 0, yend = (max(G_coordinates_sincos()[1,])+0.5) * tan(pi/9), showlegend=F) %>%
			layout(autosize=T,
				title = paste('association plot with gene set enrichment analysis \n', input$num_dims, ' first dimensions, samples: ', paste(which_samples(), collapse = ','),'\n'),
				xaxis = list(title = 'Distance from origin (x)', rangemode = "tozero"),
				yaxis = list(title = 'Distance from gene to its projection (y)', rangemode = "tozero"),
				showlegend = F
			)
		})


		output$plotsincos_gsea_legend <- renderPlotly({
			plot_ly(type='scatter', source='plotsincos_gsea', mode='markers') %>%
			add_trace(data=tmpx(), x = ~xc, y = ~yc,  type="scatter", color = ~GOterm , 
				legendgroup = 'group2',
				marker=list( size=10 , opacity=0.8)
			) %>%
			layout(autosize=T, xaxis = list(
					  title = "",
					  zeroline = FALSE,
					  showline = FALSE,
					  showticklabels = FALSE,
					  showgrid = FALSE
				),
				yaxis = list(
					  title = "",
					  zeroline = FALSE,
					  showline = FALSE,
					  showticklabels = FALSE,
					  showgrid = FALSE
				),
				showlegend = T,legend = list(x = 0.25, y = 0.8)
			)
		})


	}) # end of observe event (button "Start calculations!)
}

addResourcePath('img','img_folder')

# enable bookmarking
enableBookmarking(store = "url")

shinyApp(ui = ui, server = server)
