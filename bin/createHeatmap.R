################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, metavar="path", help="Input sample file"),
  make_option(c("-g", "--geneFunctions_file"), type="character", default=NULL, metavar="path", help="Gene Functions file."),
  make_option(c("-a", "--annoData_file"), type="character", default=NULL, metavar="path", help="Annotation Data file."),
  make_option(c("-p", "--outprefix"), type="character", default='projectID', metavar="string", help="Output prefix.")
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

sampleInput=opt$input_file
geneInput=opt$geneFunctions_file
annoInput=opt$annoData_file
outprefix=opt$outprefix

testing="Y"
if (testing == "Y"){
  sampleInput="sampleData.csv"
  geneInput="geneFunctions.csv"
  annoInput="annoData.csv"
  outprefix="test"
}


if (is.null(sampleInput)){
  print_help(opt_parser)
  stop("Please provide an input file.", call.=FALSE)
}

################################################
################################################
## READ IN FILES##
################################################
################################################
sampleData=read.csv(sampleInput,row.names=1)
annoData=read.csv(annoInput,row.names=1)
geneFunctions=read.csv(geneInput,row.names=1)

################################################
################################################
## Set colors##
################################################
################################################
annoColors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)

################################################
################################################
## Create a basic heatmap##
################################################
################################################
basic_output <- paste0("basic_heatmap_", outprefix, ".pdf")
pdf(basic_output)
pheatmap(
  sampleData,
  clustering_distance_rows = "euclidean",   # Row clustering using Euclidean distance
  clustering_distance_cols = "euclidean",   # Column clustering using Euclidean distance
  clustering_method = "ward.D",             # Use Ward’s clustering method
  main = "Basic Heatmap",                   # Set the title of the heatmap
  fontsize_row = 8,                         # Font size for row labels
  fontsize_col = 8                          # Font size for column labels
)
dev.off()

################################################
################################################
## Create a basic heatmap##
################################################
################################################
complex_output <- paste0("complex_heatmap_", outprefix, ".pdf")
pdf(complex_output)
pheatmap(
  sampleData,
  annotation_col = annoData,                # Add column annotations
  annotation_colors = annoColors,          # Set colors for annotations
  clustering_distance_rows = "euclidean",  # Row clustering using Euclidean distance
  clustering_distance_cols = "euclidean",  # Column clustering using Euclidean distance
  clustering_method = "ward.D",            # Use Ward’s clustering method
  main = "Complex Heatmap",                # Set the title of the heatmap
  fontsize_row = 8,                        # Font size for row labels
  fontsize_col = 8,                        # Font size for column labels
  show_rownames = FALSE,                   # Hide row names
  show_colnames = FALSE,                   # Hide column names
  legend_breaks = c(min(sampleData), mean(sampleData), max(sampleData)), # Legend breaks
  legend_labels = c("Low", "Medium", "High")                              # Legend labels
)
dev.off()
