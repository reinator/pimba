# ---
# title: Manipulate PIMBA output in R
# description: Diversity analysis for metabacoding data
# Take in output from PIMBA and manipulate files in R
# ---

# Setup environment
args<-commandArgs(TRUE)
#Rscript Phyloseq_metabarcoding.vPimba.R <otu_table> <tax_assignment> <metadata> 
####PHYLOSEQ PIPELINE######################################


# Setup environment
library(phyloseq)
library(ggplot2)
library(ape)
library(scales)
library(plyr)
library(Biostrings)
library(vegan)
library(grid)
library(reshape2)
library(pvclust)

###Creating a phyloseq obj###

# Read in OTU table
otu_table_file = args[1]
otu_table_in <- read.csv(otu_table_file, sep = "\t", row.names = 1)
otu_table_in <- as.matrix(otu_table_in)
otu_table_in

# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy_file = args[2]
taxonomy <- read.csv(taxonomy_file, sep = "\t", row.names = 1, check.names = FALSE)
taxonomy <- as.matrix(taxonomy)
taxonomy

# Read in metadata
metadata_file = args[3]
metadata <- read.csv(metadata_file, sep = ",", row.names = 1)
metadata


# Import all as phyloseq objects
OTU <- otu_table(otu_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)

#Combine data
ps <- phyloseq(OTU, TAX, META)

# Creat tree
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

# Finally merge!
ps <- phyloseq(OTU, TAX, META, random_tree)
ps

# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)

# Same sample names
sample_names(OTU)
sample_names(META)

groupby = args[4]

cl = c("yellow","blue","red","cyan","green",
       "magenta","snow","mistyrose4","darkgreen","cyan4",
       "royalblue1","darkorange4","palegreen1","chocolate1","hotpink1",
       "chartreuse","navy","darkmagenta","springgreen","purple1",
       "lightskyblue","lightgoldenrod1","deeppink2","slateblue4","olivedrab1",
       "lavenderblush3","olivedrab","limegreen","brown3","turquoise",
       "mediumpurple1","seagreen3","deepskyblue2","lightsalmon2","gold3",
       "steelblue","orchid1","hotpink4","magenta3","springgreen4",
       "darkred","paleturquoise1","orangered","darkseagreen","mediumpurple3",
       "goldenrod1","aquamarine4","olivedrab3","hotpink3","aquamarine",
       "seagreen1","darkorange3","lightseagreen","purple3","khaki3",
       "skyblue3","salmon3","dodgerblue4","bisque","maroon1",
       "plum1","chartreuse3","springgreen3","darkseagreen1","lightpink1",
       "orchid3","dodgerblue3","brown1","rosybrown3","aquamarine3",
       "chartreuse4","sienna","deeppink4","darkolivegreen1","lightblue3",
       "blue3","green3","red3","violetred3","brown4",
       "forestgreen","dodgerblue","gold4","coral1","lightsteelblue2",
       "lightgoldenrod4","orangered3","orange","darkviolet","mediumorchid4",
       "slateblue1","plum","steelblue4","darkolivegreen","burlywood2",
       "mediumorchid1","cyan3","cadetblue","steelblue1","honeydew2",
       "tan1","deepskyblue3","darkolivegreen3","gold","royalblue3",
       "darkgoldenrod","violetred2","green4","lightskyblue4","midnightblue",
       "darkorchid1","maroon","darkseagreen4","lightpink4","cornflowerblue",
       "seagreen","purple4","darkorange2","cadetblue1","firebrick",
       "goldenrod","slateblue3","mediumorchid3")
###Verificar as ids das variaveis

rank_names(ps)

sample_variables(ps)


###Clustering analysis


taxon_matrix <- otu_table(ps)
taxon_matrix  <- as.matrix(taxon_matrix)
taxon_matrix

# Extract current column names (sample names)
sample_names <- colnames(taxon_matrix)

# Ensure metadata column names are characters
metadata$SampleName <- as.character(rownames(metadata))

# Only rename if groupby is NOT "False"
if (groupby != "False") {
  # Check if the column exists in metadata
  if (groupby %in% colnames(metadata)) {
    # Create a named vector for mapping SampleName -> Chosen Grouping Column
    sample_group_map <- setNames(metadata[[groupby]], metadata$SampleName)
    
    # Rename columns by appending the selected metadata column value
    new_colnames <- paste0(sample_names, "_", sample_group_map[sample_names])
    
    # Assign the new column names to taxon_matrix
    colnames(taxon_matrix) <- new_colnames
  } else {
    stop(paste("Error: Column", groupby, "not found in metadata"))
  }
}


result <- pvclust(as.matrix(taxon_matrix),
                  method.dist =  "correlation", 
                  method.hclust="ward.D2", nboot=1000)

svg("cluster_bootstrap1000.svg", width = 15, height = 8)

p = plot(result,
     main="Cluster analysis. Bootstrap = 1000"
) + pvrect(result, alpha=0.95, border = 2, lwd= 1, xpd = TRUE)


dev.off()
