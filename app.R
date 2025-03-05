list_of_packages <- c("shiny", "bslib", "Seurat", "ggplot2", "dplyr", "bsicons", "plotly", "presto", "DT")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)) install.packages(new_packages)

invisible(lapply(list_of_packages, library, character.only = TRUE))

# --- Load Data BEFORE ui or server (NON-REACTIVE) ---
load("/Volumes/DMGE$/teamfolder/Shiny_sc/sc_obj.RData")  # Correct path.  This MUST be a valid path!

# Ensure 'identity' is a factor
sc_obj <- SetIdent(sc_obj, value = factor(Idents(sc_obj),
                                          levels = c("C1", "C2", "C3", "M1", "M2", "M3", "M4", "P1", "P2", "P3", "P4")))  # Adapt if needed


#Removing Ribossomal, blood, etc
Rik.genes <- grep(pattern = "\\d+.*Rik$", x = rownames(x = sc_obj@assays$RNA), value = TRUE)
Rps.genes <- grep(pattern = "^Rps", x = rownames(x = sc_obj@assays$RNA), value = TRUE)
Rpl.genes <- grep(pattern = "^Rpl", x = rownames(x = sc_obj@assays$RNA), value = TRUE)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = sc_obj@assays$RNA), value = TRUE)
pseudo.genes <- grep(pattern = "^Gm", x = rownames(x = sc_obj@assays$RNA), value = TRUE)
blood.genes <- c('Hbb-y', 'Hba-a1', 'Hbb-x', 'Hbb-bh1', 'Hba-a2', 'Hbb-bh1', "Hba-x", "Hbb-bs", "Hbb-bt")
removal <- read.csv(file = "/Users/victorio/ShinyApp/sc_app/remove.csv", header = F) #make sure the path is right
ychrom <- removal$V2
ychrom <- ychrom[1:1555]
xchrom <- removal$V3
removing <- c(Rps.genes, Rpl.genes, mito.genes, pseudo.genes, blood.genes, ychrom, xchrom, Rik.genes)
rm(Rps.genes, Rpl.genes, mito.genes, pseudo.genes, blood.genes, ychrom, xchrom) # This line doesn't really do anything useful in a Shiny context.
all_genes <- rownames(sc_obj@assays$RNA) # Get all available genes from the Seurat object
genes_to_test <- setdiff(all_genes, removing) # Filter out the excluded genes

# Define the User Interface ####
ui <- fluidPage(
  titlePanel("Gene Expression Visualization"),
  theme = bs_theme(version = 5, bootswatch = "lumen"),
  # Define the sidebar
  sidebarLayout(
    sidebarPanel(
      tags$img(src = "https://static.wixstatic.com/media/cf7f0d_baf87b2bcb6d47a687a25cad51b2ccd5~mv2.png",
               style = "width: 100%;", alt = "Developmental Genetics"),
      br(), br(), br(),
      downloadButton("downloadDimplot", "Download Dimplot"),
      br(), br(), br(),
      h4("Differential Expression"),
      selectInput("identity_test", "Select Identity", unique(sc_obj$identity), selected = unique(sc_obj$identity)[5]),
      selectInput("genotype1", "Genotype 1", choices = unique(sc_obj$genotype), selected = unique(sc_obj$genotype)[1]),       
      selectInput("genotype2", "Genotype 2", choices = unique(sc_obj$genotype), selected = unique(sc_obj$genotype)[2]),
      actionButton("run_de", "Run Differential Expression"),
      downloadButton("downloadDE", "Download DE Results"),
      downloadButton("downloadDEVolcano", "Download Volcano Plot"),
      br(), br(), br(), br(),
      h4("FeaturePlot & VlnPlot Inputs"),
      textInput("gene", "Enter Gene Name:"),
      actionButton("updatePlot", "Update Plot"), # Keep the action button
      br(),
      downloadButton("downloadFeatureplot", "Download FeaturePlot"),
      downloadButton("downloadVlnplot", "Download VlnPlot"),
      br(), br(), br(),
      h4("Scatter Plot Inputs"),
      textInput("gene_x", "X-axis Gene:"),
      textInput("gene_y", "Y-axis Gene:"),
      textInput("gene_color", "Color by Gene:"),
      downloadButton("downloadScatterplot", "Download Scatterplot"),
      h4("Total Cells per Genotype"),
      tableOutput("total_cell_count_table"),
      h4("Cell Count by Genotype"),
      textInput("gene_count", "Enter Gene Name:"),
      tableOutput("cell_count_table"),
      width = 3
    ),
    # Define the cards
    mainPanel(
      fluidRow(
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("Dimplot"))),
                 plotOutput("dimPlot"),
                 actionButton("showHeatmap", "Show Heatmap"),
                 downloadButton("downloadMarkers", "Download Markers")
               )
        ),
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("Volcano plot"))),
                 plotOutput("de_volcano"),
                 # Add numericInputs *inside* the Volcano Plot card, below the plot
                 hr(),  # Add a horizontal line for visual separation
                 fluidRow(
                   column(4, numericInput("xlim_min", "X-axis Min:", value = -1)),
                   column(4, numericInput("xlim_max", "X-axis Max:", value = 1)),
                   column(4, numericInput("ylim_max", "Y-axis Max:", value = 5))
                 )
               )
        )
      ),
      fluidRow(
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("FeaturePlot"))),
                 plotOutput("featurePlot"),
                 actionButton("showSplitFeaturePlot", "Split view by genotype")
               )
        ),
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("VlnPlot"))),
                 plotOutput("vlnPlot"),
                 actionButton("showSplitVlnPlot", "Split view by genotype")
               )
        )
      ),
      fluidRow(
        column(12,
               card(
                 card_header(div(style = "text-align: center;", strong("Scatterplot"))),
                 plotOutput("scatterPlot")
               )
        )
      ),
      fluidRow(
        column(12,
               card(
                 card_header(div(style = "text-align: center;", strong("Blended FeaturePlot"))),
                 plotOutput("blendedFeaturePlot"),
                 hr(), # Horizontal line for separation
                 fluidRow(
                   column(6, textInput("gene1", "Gene 1:", value = "Msx1")),
                   column(6, textInput("gene2", "Gene 2:", value = "Sox9"))
                 ),
                 actionButton("updateBlendedPlot", "Update Blended Plot") # Button to trigger plot update
               )
        )
      ),
      fluidRow(
        column(12,
               card(
                 card_header(div(style = "text-align: center;", strong("Differential Expression Results"))),
                 conditionalPanel(
                   condition = "input.run_de > 0",
                   DT::dataTableOutput("de_results")
                 )
               )
        )
      ),
      width = 9
    )
  )
)

# Define the server ####
server <- function(input, output, session) {
  
  # Create a reactiveVal to store the markers
  top10_markers <- reactiveVal(NULL)
  
  # --- Update FeaturePlot and VlnPlot on Button Click ---
  observeEvent(input$updatePlot, {
    req(input$gene)
    
    # Input Validation: Check if the gene exists
    validate(
      need(input$gene %in% rownames(sc_obj),
           paste("Gene", input$gene, "not found in Seurat object."))
    )
    
    # Create Feature Plots *INSIDE* observeEvent (and only if the gene is valid)
    output$featurePlot <- renderPlot({
      FeaturePlot(sc_obj, features = input$gene, order = TRUE)
    })
    output$featurePlotSplit <- renderPlot({
      FeaturePlot(sc_obj, features = input$gene, split.by = "genotype", order = TRUE)
    })
    
    
    # Create Violin Plots *INSIDE* observeEvent (and only if the gene is valid)
    output$vlnPlot <- renderPlot({
      VlnPlot(sc_obj, features = input$gene, pt.size = 0.1)
    })
    output$vlnPlotSplit <- renderPlot({
      VlnPlot(sc_obj, features = input$gene, pt.size = 0.1, split.by = 'genotype')
    })
    
  })
  
  # --- Show Modal with Split FeaturePlot ---
  observeEvent(input$showSplitFeaturePlot, {
    req(input$gene) #Ensure gene input
    
    # Input validation
    validate(
      need(input$gene %in% rownames(sc_obj),
           paste("Gene", input$gene, "not found in Seurat Object."))
    )
    # Show modal only if gene is valid
    if (input$gene %in% rownames(sc_obj)){
      showModal(
        modalDialog(
          title = "FeaturePlot Split by Genotype",
          plotOutput("featurePlotSplit", height = "300px", width = "900px"),  #Plot inside modal
          easyClose = TRUE,
          footer = modalButton("Close")
        )
      )
    }
  })
  
  # --- Show Modal with Split VlnPlot ---
  observeEvent(input$showSplitVlnPlot, {
    req(input$gene)
    #Input validation
    validate(
      need(input$gene %in% rownames(sc_obj), "Gene not found")
    )
    #Show modal only if gene is valid
    if (input$gene %in% rownames(sc_obj)){
      showModal(modalDialog(
        title = "Violin Plot - Split by Genotype",
        plotOutput("vlnPlotSplit",  height = "300px", width = "900px"), # Plot inside modal
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    }
  })
  
  # --- Run and Display Differential Expression Results ---
  observeEvent(input$run_de, {
    req(input$identity_test, input$genotype1, input$genotype2)
    
    withProgress(message = "Calculating Differential Expression...", value = 0, {
      # Subset Seurat object
      subset_obj <- subset(sc_obj, idents = c(input$identity_test))
      subset_obj <- subset(subset_obj, subset = genotype %in% c(input$genotype1, input$genotype2))
      incProgress(0.1, detail = "Subsetting data...")
      
      # Normalize, Find Variable Features, and Scale Data
      subset_obj <- NormalizeData(subset_obj) %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "orig.ident"))
      incProgress(0.3, detail = "Normalizing data...")
      # Run DE analysis
      tryCatch({
        de_wilcox <- wilcoxauc(subset_obj, group_by = "genotype") %>%
          filter(group == input$genotype2) %>%
          mutate(DE = abs(logFC) > log(1.1) & padj < 0.01) %>%
          mutate(DEG = ifelse(DE, feature, NA))
        incProgress(0.7, detail = "Calculating DE...")
        
        # Get initial axis limits
        initial_xlim <- c(min(de_wilcox$logFC), max(de_wilcox$logFC))
        initial_ylim <- c(0, max(-log10(de_wilcox$padj)))
        
        # Update numericInputs with initial values
        updateNumericInput(session, "xlim_min", value = round(initial_xlim[1],2))
        updateNumericInput(session, "xlim_max", value = round(initial_xlim[2],2))
        updateNumericInput(session, "ylim_max", value = round(initial_ylim[2],2))
        
        # Render DE results table (inside tryCatch)
        output$de_results <- DT::renderDataTable({
          de_wilcox %>% arrange(padj)
        })
        
        
        # Render Volcano Plot (inside tryCatch, with dynamic limits)
        output$de_volcano <- renderPlot({
          ggplot(de_wilcox, aes(x = logFC, y = -log10(padj), col = DE, label = DEG)) +
            geom_point() +
            ggrepel::geom_text_repel() +
            geom_vline(xintercept = c(-log(1.1), log(1.1), 0), linetype = "dotted") +
            geom_hline(yintercept = -log10(0.01), linetype = "dotted") +
            scale_color_manual(values = c("#909090", "red")) +
            theme_minimal() +
            xlim(input$xlim_min, input$xlim_max) +
            ylim(0, input$ylim_max)
        })
        incProgress(1.0, detail = "Rendering plot...") # Complete
        
      }, error = function(e) {
        # ... (error handling, same as before) ...
        output$de_volcano <- renderPlot({
          ggplot() +
            geom_text(aes(x = 0.5, y = 0.5, label = paste("Error in DE:", e$message)), size = 5) +
            theme_void()
        })
        output$de_results <- DT::renderDataTable(NULL)
      })
    })
  })
  
  # --- Download Handlers ---
  
  output$downloadDEVolcano <- downloadHandler(
    filename = function() {
      paste("VolcanoPlot_", input$identity_test, "_", input$genotype1, "_vs_", input$genotype2, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(input$identity_test, input$genotype1, input$genotype2)
      # Subset Seurat object
      subset_obj <- subset(sc_obj, idents = c(input$identity_test))
      subset_obj <- subset(subset_obj, subset = genotype %in% c(input$genotype1, input$genotype2))
      # Normalize, Find Variable Features, and Scale Data (within the subset)
      subset_obj <- NormalizeData(subset_obj) %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "orig.ident"))
      
      # Differential expression using wilcoxauc
      
      de_wilcox <- wilcoxauc(subset_obj, group_by = "genotype") %>%
        filter(group == input$genotype2) %>%
        mutate(DE = abs(logFC) > log(1.1) & padj < 0.01) %>%
        mutate(DEG = ifelse(DE, feature, NA))
      pdf(file, width = 8, height = 6)  # Create PDF
      print(ggplot(de_wilcox, aes(x = logFC, y = -log10(padj), col = DE, label = DEG)) + #Plot
              geom_point() +
              ggrepel::geom_text_repel() +
              geom_vline(xintercept = c(-log(1.1), log(1.1), 0), col = "#303030", linetype = "dotted") +
              geom_hline(yintercept = -log10(0.01), col = "#303030", linetype = "dotted") +
              scale_color_manual(values = c("#909090", "red")) +
              theme_minimal()+
              xlim(input$xlim_min, input$xlim_max) +  # Apply x-axis limits
              ylim(0, input$ylim_max)
      )
      dev.off() # Close the PDF device
    }
  )
  
  output$downloadDE <- downloadHandler(
    filename = function() {
      paste("DE_results_", input$identity_test, "_", input$genotype1, "_vs_", input$genotype2, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$identity_test, input$genotype1, input$genotype2)
      # Subset Seurat object
      subset_obj <- subset(sc_obj, idents = c(input$identity_test))
      subset_obj <- subset(subset_obj, subset = genotype %in% c(input$genotype1, input$genotype2))
      # Normalize, Find Variable Features, and Scale Data (within the subset)
      subset_obj <- NormalizeData(subset_obj) %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "orig.ident"))
      
      # Differential expression using wilcoxauc
      
      de_wilcox <- wilcoxauc(subset_obj, group_by = "genotype") %>%
        filter(group == input$genotype2) %>%
        mutate(DE = abs(logFC) > log(1.1) & padj < 0.01) %>%
        mutate(DEG = ifelse(DE, feature, NA))
      write.csv(de_wilcox, file, row.names = FALSE) # Write to CSV
    }
  )
  
  output$downloadScatterplot <- downloadHandler(
    filename = function() {
      paste("Scatterplot_", Sys.Date(), "_",
            input$gene_x, "_", input$gene_y, "_", input$gene_color,
            ".pdf", sep = "")
    },
    content = function(file) {
      req(input$gene_x, input$gene_y, input$gene_color) # Ensure genes are selected
      # Gene Validation
      validate(
        need(all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj)),
             "One or more selected genes not found in Seurat object.")
      )
      # Continue if genes exist
      if(all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj))){
        pdf(file, width = 7.5, height = 4.5)
        expr_data <- FetchData(sc_obj, vars = c(input$gene_x, input$gene_y, input$gene_color, "genotype"))
        print(ggplot(expr_data, aes_string(x = input$gene_x, y = input$gene_y)) +
                geom_point(aes_string(color = input$gene_color), size = 2, alpha = 0.8) +
                scale_color_gradient(low = "grey", high = "red") +
                labs(color = paste(input$gene_color, "Expression")) +
                geom_vline(xintercept = 0, linetype = "dashed") +
                geom_hline(yintercept = 0, linetype = "dashed") +
                theme_minimal() +
                facet_wrap(~ genotype)) # Facet by genotype
        dev.off()
      }
    }
  )
  
  output$downloadDimplot <- downloadHandler(
    filename = function() {
      paste("Dimplot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 7.5, height = 4.5)
      print(DimPlot(sc_obj, group.by = 'identity', reduction = "umap", label = TRUE))
      dev.off()
    }
  )
  
  
  output$downloadFeatureplot <- downloadHandler(
    filename = function() {
      paste("Featureplot_", input$gene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(input$gene) # Ensure gene input
      if (input$gene %in% rownames(sc_obj)) {
        pdf(file, width = 7.5, height = 4.5)
        print(FeaturePlot(sc_obj, features = input$gene, order = TRUE))
        dev.off()
      }
    }
  )
  
  output$downloadVlnplot <- downloadHandler(
    filename = function() {
      paste("Vlnplot_", input$gene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(input$gene) # Ensure gene input
      if (input$gene %in% rownames(sc_obj)) { # Gene Validation
        pdf(file, width = 7.5, height = 4.5)
        print(VlnPlot(sc_obj, features = input$gene, pt.size = 0.1))
        dev.off()
      }
    }
  )
  
  
  # Cell count calculation (per genotype, for a given gene)
  output$cell_count_table <- renderTable({
    req(input$gene_count) # Ensure gene is entered
    #Gene validation
    if(input$gene_count %in% rownames(sc_obj)){
      gene_expr <- FetchData(sc_obj, vars = c(input$gene_count, "genotype"))
      gene_sym <- input$gene_count # More robust way to handle gene names
      cell_count <- table(subset(gene_expr, gene_expr[[gene_sym]] > 0)$genotype)
      t(as.matrix(cell_count)) # Transpose for better display
    }
  }, rownames = TRUE, colnames = TRUE)
  
  # Total cell count (per genotype)
  output$total_cell_count_table <- renderTable({
    total_cell_count <- table(sc_obj$genotype)
    t(as.matrix(total_cell_count))  # Transpose
  }, rownames = TRUE, colnames = TRUE)
  
  
  # --- Initial Dimplot (Rendered Once) ---
  output$dimPlot <- renderPlot({
    DimPlot(sc_obj, group.by = 'identity', reduction = "umap", label = TRUE)
  })
  
  # --- Show Heatmap in Modal (and calculate markers) ---
  observeEvent(input$showHeatmap, {
    
    # 1. Show the modal *immediately* with a loading message
    showModal(modalDialog(
      title = "Heatmap of Top 10 Markers per Cluster",
      "Calculating markers and rendering heatmap...",  # Initial message
      plotOutput("heatmap", height = "700px"), # Still include plotOutput
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("downloadHeatmap", "Download Heatmap")
      )
    ))
    
    
    # 2. Calculate Markers *INSIDE* the observeEvent, within withProgress
    withProgress(message = "Calculating Markers...", value = 0, {
      cluster_markers_identity <- FindAllMarkers(sc_obj, features = genes_to_test,
                                                 only.pos = TRUE, min.pct = 0.25,
                                                 logfc.threshold = log(1.2))
      incProgress(0.6, detail = "Grouping Markers...")
      top10 <- cluster_markers_identity %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)
      # 3. Store the results in the reactiveVal
      top10_markers(top10)  # Store the calculated markers
      incProgress(1.0, detail = "Done")
    })
    
    # --- Render Heatmap (using the reactiveVal) ---
    output$heatmap <- renderPlot({
      req(top10_markers())  # Use req() to wait for markers
      
      DoHeatmap(sc_obj, features = top10_markers()$gene) +
        theme(axis.text.y = element_text(size = 8))
    })
  })
  
  # --- Download Handler (using the reactiveVal)
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste("Heatmap_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(top10_markers()) # Make sure markers exist.
      pdf(file, width = 6, height = 8)
      print(DoHeatmap(sc_obj, features = top10_markers()$gene) +
              theme(axis.text.y = element_text(size = 8)))
      dev.off()
    }
  )
  
  # --- Download Markers ---
  output$downloadMarkers <- downloadHandler(
    filename = function() {
      paste("Cluster_Markers_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(top10_markers()) #wait for markers
      write.csv(top10_markers(), file, row.names = FALSE) # Use top10_markers()
    }
  )
  
  
  
  output$scatterPlot <- renderPlot({
    req(input$gene_x, input$gene_y, input$gene_color)
    # Gene Validation
    validate(
      need(all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj)),
           "One or more selected genes not found in Seurat object.")
    )
    # Continue if genes exist
    if(all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj))){
      expr_data <- FetchData(sc_obj, vars = c(input$gene_x, input$gene_y, input$gene_color, "genotype"))
      ggplot(expr_data, aes_string(x = input$gene_x, y = input$gene_y)) +
        geom_point(aes_string(color = input$gene_color), size = 2, alpha = 0.8) +
        scale_color_gradient(low = "grey", high = "red") +
        labs(color = paste(input$gene_color, "Expression")) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        theme_minimal() +
        facet_wrap(~ genotype) # Facet by genotype
    }
  })
  # --- Blended FeaturePlot ---
  observeEvent(input$updateBlendedPlot, {
    req(input$gene1, input$gene2)
    
    # Validate that both genes exist
    validate(
      need(input$gene1 %in% rownames(sc_obj) && input$gene2 %in% rownames(sc_obj),
           "One or both of the specified genes were not found.")
    )
    
    output$blendedFeaturePlot <- renderPlot({
      FeaturePlot(sc_obj, features = c(input$gene1, input$gene2), blend = TRUE, 
                  order = T, blend.threshold = 0.1)
    })
  })
}

shinyApp(ui = ui, server = server)
