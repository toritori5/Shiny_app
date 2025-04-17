Single-Cell RNA-Seq Data Visualization APP

This document guides you through using a Shiny app to explore single-cell RNA sequencing (scRNA-seq) data. The app lets you visualize gene expression, compare different groups of cells, and identify genes that are differentially expressed.

Part 1: Setting Up Your Environment (One-Time Setup)
1.	Install R:
o	Go to https://cran.r-project.org/
o	Click the link for your operating system (Windows, macOS, or Linux).
o	Download and install the latest version of R. Follow the on-screen instructions.
2.	Install RStudio:
o	Go to https://posit.co/download/rstudio-desktop/
o	Download and install the free RStudio Desktop version for your operating system.
3.	Install Xcode:
o	Go to the App Store and download Xcode. Install the app.

Part 2: Running the App
1.	Open RStudio.
2.	Connect to the server. The single cell object is hosted in the dmg server and needs to be available for the app to run. 
3.	Run the App: In the RStudio Console (usually at the bottom left), copy and paste the following code, then press Enter:
install.packages("shiny")
library(shiny)
runGitHub("Shiny_app", "toritori5")

This will download the app from GitHub and launch it in a pop up window. You only need to run “install.packages("shiny")” the first time. 

Part 3: Using the App
The app interface has two main sections: a sidebar on the left and a main panel on the right.
Sidebar (Left Side):
•	Download Dimplot: Download the main "Dimplot" as a PDF.
•	Differential Expression:
o	Select Identity: Choose a cell type or cluster to analyze (e.g., C1, C2, M1).
o	Genotype 1: Select the first group you want to compare (e.g., 1_wt).
o	Genotype 2: Select the second group you want to compare (e.g., 2_E1C8).
o	Run Differential Expression: Click this to find genes that are expressed differently between the two selected groups. A progress bar will appear.
o	Download DE Results: Download a table of the differential expression results.
o	Download Volcano Plot: Download the volcano plot as a PDF.
•	FeaturePlot & VlnPlot Inputs:
o	Enter Gene Name: Type the name of a gene you're interested in (e.g., Tfap2b). Make sure to use the correct capitalization.
o	Update Plot: Click this to display the gene's expression on the FeaturePlot and VlnPlot.
o	Download Featureplot: Download the current FeaturePlot.
o	Download Vlnplot: Download the current VlnPlot.
•	Scatter Plot Inputs:
o	X-axis Gene: Enter a gene for the x-axis of a scatter plot.
o	Y-axis Gene: Enter a gene for the y-axis.
o	Color by Gene: Enter a gene to color the points by.
o	Download Scatterplot: Download the scatter plot.
•	Total Cells per Genotype: Shows the total number of cells in each group.
•	Cell Count by Genotype:
o	Enter Gene Name: Type a gene name.
o	The table will show how many cells in each group express that gene.
Main Panel (Right Side):
•	Dimplot: This is a main overview plot. Each point is a cell, and cells are grouped by their type.
o	Show Heatmap: Click to see a heatmap of the top marker genes for each cluster (this may take some seconds to load). You can also download the markers.
o	Download Markers: Download the list of markers used in the Heatmap.
•	Volcano plot: This appears after you run the "Differential Expression" analysis. It shows which genes are significantly different between the two groups you selected. Genes further to the right are upregulated in Genotype 2, genes to the left are upregulated in Genotype 1. Genes higher up are more statistically significant. You can change the x and Y axis of the volcano plot using the numeric inputs below the plot.
•	FeaturePlot: Shows the expression of the gene you entered in the sidebar on the Dimplot. Brighter colors usually mean higher expression.
o	Show Split FeaturePlot: Click to see the FeaturePlot split by genotype (in a pop-up window). Right click to download.
•	VlnPlot: Another way to visualize gene expression. It shows the distribution of expression levels for the gene you entered.
o	Show Split VlnPlot: Click to see the VlnPlot split by genotype (in a pop-up window). Right click to download.
•	Scatterplot: Shows the relationship between the expression of three genes (X-axis, Y-axis and colored) you entered in the sidebar.
•	Blended FeaturePlot: Combines the expression of two genes on the Dimplot.
o	Gene 1: Enter the first gene.
o	Gene 2: Enter the second gene.
o	Update Blended Plot: Click to create the blended plot.
•	Differential Expression Results: This table appears after you run the "Differential Expression" analysis. It lists the genes that are significantly different between the two groups.
Important Notes:
•	Case Sensitivity: Gene names are case-sensitive. "Sox9" is different from "sox9."
•	Waiting: Some calculations (like "Run Differential Expression" and "Show Heatmap") may take some time. Be patient and wait for the progress bar to complete.
•	Error Messages: If a message like "Gene not found" appears, double check that the gene you entered is correct.

