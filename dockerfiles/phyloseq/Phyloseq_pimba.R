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
#library(RColorBrewer)
#library(colorRamps)
library(scales)
library(plyr)
library(Biostrings)
library(vegan)
#library(diveRsity)
library(dplyr)
library(grid)
library(reshape2)
library(gridExtra)
library(xfun)
library(pvclust)
library(dendextend) 

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

# Plot the phylum composition
ps_phylum <- ps %>%
  tax_glom(taxrank = "phylum") %>%                     # agglomerate at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(phylum)                                      # Sort data frame alphabetically by Phylum

# library(randomcoloR)
# set.seed(8675309)
# n <- length(levels(ps_phylum$phylum))
# palette <- distinctColorPalette(n)

# library(RColorBrewer)  ### Analysis (https://vaulot.github.io/tutorials/Phyloseq_tutorial.html)
# set.seed(8675309)
# n <- length(levels(ps_phylum$phylum))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

n <- length(unique(ps_phylum$phylum))
cl0=cl
while (length(cl) <= n){
  if(n-length(cl) >= 239){
    cl = c(cl,cl0)
  }  else{
    cl = c(cl,cl0[0:(n-length(cl)+1)])
  }
}

palette = cl

p = ggplot(ps_phylum, aes(x = Sample, y = Abundance, fill = phylum)) + 
  geom_bar(stat = "identity", position="fill") +
  #facet_grid(~local) +
  #facet_wrap(~metadata$, scales="free_x", nrow=1)+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1, size=12)) +
  ylab("Relative Abundance (Phylum)") +
  #ggtitle("phylum") + 
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 30, hjust = 1, size=11))+
  theme(legend.text=element_text(size=rel(1.2))) + ggtitle("Phylum barplot")

if(groupby!=FALSE){
  groupby=paste("~",groupby)
  p = p + facet_grid(groupby, scale="free")
}

svg("phylum_barplots.svg", width = 15, height = 8)
print(p)
dev.off()


# Plot the class composition
ps_class <- ps %>%
  tax_glom(taxrank = "class") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(class)                                      # Sort data frame alphabetically by Class

# set.seed(8675309)
# n <- length(levels(ps_class$class))
# palette <- distinctColorPalette(n)

# set.seed(8675309)
# n <- length(levels(ps_class$class))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# 
n <- length(unique(ps_class$class))

cl0=cl
while (length(cl) <= n){
  if(n-length(cl) >= 239){
    cl = c(cl,cl0)
  }  else{
    cl = c(cl,cl0[0:(n-length(cl)+1)])
  }
}
palette <- cl

p = ggplot(ps_class, aes(x = Sample, y = Abundance, fill = class)) + 
  geom_bar(stat = "identity", position="fill") +
  #facet_grid(~local) +
  #facet_wrap(~area, scales="free_x", nrow=1)+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, size=12, ncol=2)) +
  ylab("Relative Abundance (Class)") +
  #ggtitle("genus") + 
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 30, hjust = 1, size=11))+
  theme(legend.text=element_text(size=rel(1.2)))+
  theme (strip.text = element_text(size = 12)) + ggtitle("Class barplot")

if(groupby!=FALSE){
  #groupby=paste("~",groupby)
  p = p + facet_grid(groupby, scale="free")
}

svg("class_barplots.svg", width = 15, height = 8)
print(p)
dev.off()


# Plot the order composition
ps_order <- ps %>%
  tax_glom(taxrank = "order") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(order)                                      # Sort data frame alphabetically by Class

# set.seed(8675309)
# n <- length(levels(ps_order$order))
# palette <- distinctColorPalette(n)

# set.seed(8675309)
# n <- length(levels(ps_order$order))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

n <- length(unique(ps_order$order))

cl0=cl
while (length(cl) <= n){
  if(n-length(cl) >= 239){
    cl = c(cl,cl0)
  }  else{
    cl = c(cl,cl0[0:(n-length(cl)+1)])
  }
}
palette <- cl

p = ggplot(ps_order, aes(x = Sample, y = Abundance, fill = order)) + 
  geom_bar(stat = "identity", position="fill") +
  #facet_grid(~local) +
  #facet_wrap(~area, scales="free_x", nrow=1)+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, size=12, ncol=2)) +
  ylab("Relative Abundance (Order)") +
  #ggtitle("genus") + 
  scale_fill_manual(values = palette) +  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=11))+
  theme(legend.text=element_text(size=rel(1.2)))+
  theme (strip.text = element_text(size = 12)) + ggtitle("Order barplot")

if(groupby!=FALSE){
  #groupby=paste("~",groupby)
  p = p + facet_grid(groupby, scale="free")
}

svg("order_barplots.svg", width = 15, height = 8)
print(p)
dev.off()


# Plot the family composition
ps_family <- ps %>%
  tax_glom(taxrank = "family") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(family)                                      # Sort data frame alphabetically by Genus

# set.seed(8675309)
# n <- length(levels(ps_family$family))
# palette <- distinctColorPalette(n)

# set.seed(8675309)
# n <- length(levels(ps_family$family))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

n <- length(unique(ps_family$family))
cl0=cl
while (length(cl) <= n){
  if(n-length(cl) >= 239){
    cl = c(cl,cl0)
  }  else{
    cl = c(cl,cl0[0:(n-length(cl)+1)])
  }
}
palette <- cl


p = ggplot(ps_family, aes(x = Sample, y = Abundance, fill = family)) + 
  geom_bar(stat = "identity", position="fill") +
  #facet_grid(~local) +
  #facet_wrap(~area, scales="free_x", nrow=1)+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, size=12, ncol=3)) +
  ylab("Relative Abundance (Family)") +
  #ggtitle("genus") + 
  scale_fill_manual(values = palette) +  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=11))+
  theme(legend.text=element_text(size=rel(1.2)))+
  theme (strip.text = element_text(size = 12)) + ggtitle("Family barplot")

if(groupby!=FALSE){
  #groupby=paste("~",groupby)
  p = p + facet_grid(groupby, scale="free")
}

svg("family_barplots.svg", width = 15, height = 8)
print(p)
dev.off()

# Plot the genus composition

ps_genus <- ps %>%
  tax_glom(taxrank = "genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(genus)                                      # Sort data frame alphabetically by Genus

# set.seed(8675309)
# n <- length(levels(ps_genus$genus))
# palette <- distinctColorPalette(n)

# set.seed(8675309)
# n <- length(levels(ps_genus$genus))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# palette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cl0 = cl
n <- length(unique(ps_genus$genus))
while (length(cl) <= n){
  if(n-length(cl) >= 239){
    cl = c(cl,cl0)
  }  else{
    cl = c(cl,cl0[0:(n-length(cl)+1)])
  }
}
palette <- cl

p = ggplot(ps_genus, aes(x = Sample, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity", position="fill") +
  #facet_grid(~local) +
  #facet_wrap(~area, scales="free_x", nrow=1)+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, size=12, ncol=3)) +
  ylab("Relative Abundance (Genus)") +
  #ggtitle("genus") + 
  scale_fill_manual(values = palette) +  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=11))+
  theme(legend.text=element_text(size=rel(1.2)))+
  theme (strip.text = element_text(size = 12)) + ggtitle("Genus barplot")


if(groupby!=FALSE){
  #groupby=paste("~",groupby)
  p = p + facet_grid(groupby, scale="free")
}

svg("genus_barplots.svg", width = 15, height = 8)
print(p)
dev.off()

### Alpha diversity analysis ###############
p = plot_richness(ps, measures=c("Observed", "Chao1", "Shannon"))
svg("alpha_diversity_dotplot.svg",width = 15, height = 8)
print(p)
dev.off()

### Beta diversity analysis ###############
otu_dataFC <-t(otu_table(ps))

S <- specnumber(otu_dataFC)
(raremax <- min(rowSums(otu_dataFC)))
Srare <- rarefy(otu_dataFC, raremax)

lty <- c("solid", "dashed", "longdash", "dotdash")

pars <- expand.grid(col = palette, lty = lty, stringsAsFactors = FALSE)
head(pars)
svg("rarefaction_curve2.svg", width = 15, height = 8)
out <- with(pars[1:26, ],
            rarecurve(otu_dataFC, step = 50, sample = raremax, col = col,
                      lty = lty, label = TRUE)) 

#dev.copy(svg, "rarefaction_curve.svg", width = 15, height = 8)
dev.off()


###Clustering analysis

groupby = args[4]

hell.tip.labels <- as(get_variable(ps, groupby), "character")
# This is the actual hierarchical clustering call, specifying average-linkage clustering
d <- phyloseq::distance(ps, method="bray", type="samples")
hell.hclust     <- hclust(d, method="average")
svg("cluster_dendogram.svg", width = 15, height = 8)

if(groupby!=FALSE){
  order=hell.hclust$order
  hell.tip.labels <- as(get_variable(ps, groupby), "character")
  labels(hell.hclust) <- paste(hell.tip.labels[order], sample_names(META)[order], sep="-")
}

# dend <- as.dendrogram(hell.hclust)
# 
# dend %>%
#   set("leaves_pch", as.numeric(hell.tip.labels[order])) %>%
#   plot()

p = plot(hell.hclust)


dev.off()





taxon_matrix <- otu_table(ps)
taxon_matrix  <- as.matrix(taxon_matrix)
taxon_matrix

result <- pvclust(as.matrix(taxon_matrix),
                  method.dist =  "correlation", 
                  method.hclust="ward.D2", nboot=1000)

svg("cluster_bootstrap1000.svg", width = 15, height = 8)

if(groupby!=FALSE){
  order=result$hclust$order
  hell.tip.labels <- as(get_variable(ps, groupby), "character")
  labels(result$hclust) <-paste(hell.tip.labels[order], sample_names(META)[order], sep="-")
}






p = plot(result,
     main="Cluster analysis. Bootstrap = 1000"
) + pvrect(result, alpha=0.95, border = 2, lwd= 1, xpd = TRUE)


dev.off()

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps, method="unifrac", weighted=F)
ordination = ordinate(ps, method="PCoA", distance=wunifrac_dist)

svg("PCoA_unifrac.svg", width = 15, height = 8)

if(groupby!=FALSE){
  sample_data(ps)[,groupby] <- sapply(sample_data(ps)[,groupby], as.factor)
  plot_ordination(ps, ordination, color=groupby, label="SampleName") + 
    theme(aspect.ratio=1) +  scale_colour_manual(values=palette) +
    ggtitle("PCoA with unweighted UniFrac distance")
} else{
  plot_ordination(ps, ordination, label="SampleName") + 
    theme(aspect.ratio=1) +  scale_colour_manual(values=palette) +
    ggtitle("PCoA with unweighted UniFrac distance")
}

dev.off()


#Do multivariate analysis based on Bray-Curtis distance and NMDS ordination.
ps.ord <- ordinate(ps, "NMDS", "bray")

svg("NMDS_bray.svg", width = 15, height = 8)
if(groupby!=FALSE){
  sample_data(ps)[,groupby] <- sapply(sample_data(ps)[,groupby], as.factor)
  plot_ordination(ps, ps.ord, type="sample", color=groupby, 
                  title="OTUs", label="SampleName") + geom_point(size=2) +
    scale_colour_manual(values=palette) +  ggtitle("NMDS with Bray-Curtis distance")
}else{
  plot_ordination(ps, ps.ord, type="sample",
                  title="OTUs", label="SampleName") + geom_point(size=2) +
    scale_colour_manual(values=palette) +  ggtitle("NMDS with Bray-Curtis distance")
}

dev.off()


