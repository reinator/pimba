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

# Create color range (https://htmlcolorcodes.com/color-picker/) for multiple colors 

# col <- c(
#   "#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
#   "#A6761D","#666666","#A6CEE3","#1F78B4", "#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
#   "#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9",
#   "#FFF2AE","#F1E2CC","#CCCCCC","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#CCEBC5","#FFED6F",
#   "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
#   "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#B3DE69","#FCCDE5",
#   "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#87D3CB", "#599861", "#F3A392", "#72A04A", "#70CB1F", "#38CDBE","#D9D9D9","#BC80BD",
#   "#878CD3", "#2E3589", "#9C89BC", "#3498DB", "#A867C9", "#D362E8", "#F1C40F", "#CF51A1", "#B8358A", "#9A5A84", 
#   "#B899AD", "#B83548", "#B8999D", "#9A5A64", "#B86335", "#B8A499", "#9A705A", "#B8A535", "#9A905A", "#B8B399", "#ADB899", "#9DB899"
# )
# 
# cl <- colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# cl <- cl[grep('\\whitesmoke\\b', cl, invert = T)]
# cl <- cl[grep('\\black\\b', cl, invert = T)]
# cl <- cl[grep('\\white\\b', cl, invert = T)]
# cl <- cl[grep('\\aliceblue\\b', cl, invert = T)]
# cl <- cl[grep('1', cl, invert = T)]
# cl <- cl[grep('2', cl, invert = T)]
# cl <- cl[grep('3', cl, invert = T)]
# cl <- cl[grep('azure', cl, invert = T)]
# cl <- cl[grep('bisque', cl, invert = T)]
# cl <- cl[grep('\\cadetblue3\\b', cl, invert = T)]
# cl <- cl[grep('cornsilk', cl, invert = T)]
# cl <- cl[grep('darkseagreen', cl, invert = T)]
# cl <- cl[grep('forestgreen', cl, invert = T)]
# cl <- cl[grep('honeydew', cl, invert = T)]
# cl <- cl[grep('hotpink', cl, invert = T)]
# cl <- cl[grep('ivory', cl, invert = T)]
# cl <- cl[grep('indianred', cl, invert = T)]
# cl <- cl[grep('lavender', cl, invert = T)]
# cl <- cl[grep('lemonchiffon', cl, invert = T)]
# cl <- cl[grep('lightblue', cl, invert = T)]
# cl <- cl[grep('lightgreen', cl, invert = T)]
# cl <- cl[grep('lightyellow', cl, invert = T)]
# cl <- cl[grep('lightgoldenrod', cl, invert = T)]
# cl <- cl[grep('lightcyan', cl, invert = T)]
# cl <- cl[grep('lightsteelblue', cl, invert = T)]
# cl <- cl[grep('\\linen\\b', cl, invert = T)]
# cl <- cl[grep('maroon', cl, invert = T)]
# cl <- cl[grep('mediumaquamarin', cl, invert = T)]
# cl <- cl[grep('mediumorchid', cl, invert = T)]
# cl <- cl[grep('mediumpurple', cl, invert = T)]
# cl <- cl[grep('mediumslateblue', cl, invert = T)]
# cl <- cl[grep('\\slateblue3\\b', cl, invert = T)]
# cl <- cl[grep('\\mintcream\\b', cl, invert = T)]
# cl <- cl[grep('mistyrose', cl, invert = T)]
# cl <- cl[grep('\\oldlace\\b', cl, invert = T)]
# cl <- cl[grep('papaya', cl, invert = T)]
# cl <- cl[grep('paleturquoise', cl, invert = T)]
# cl <- cl[grep('peachpuff', cl, invert = T)]
# cl <- cl[grep('plum', cl, invert = T)]
# cl <- cl[grep('\\rosybrown\\b', cl, invert = T)]
# cl <- cl[grep('seashell', cl, invert = T)]
# cl <- cl[grep('sienna', cl, invert = T)]
# cl <- cl[grep('snow', cl, invert = T)]
# cl <- cl[grep('wheat', cl, invert = T)]

# cl <- colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# cl <- cl[grep('black', cl, invert = T)]
# cl <- cl[grep('white', cl, invert = T)]
# palette = c("yellow", "blue", "red", "cyan", "green",  "magenta")
# max_dist=1
# while(max_dist>0.08){
#   palette_min = c()
#   palette_min_dist=c()
#   cl = setdiff(cl, palette)
#   
#   for (i in 1:length(cl)){
#     min_dist= 1
#     col1 = t(col2rgb(cl[i]))
#     for(j in 1:length(palette)){
#       col2 =  t(col2rgb(palette[j]))
#       d=sqrt((col2[1]-col1[1])^2+(col2[2]-col1[2])^2+(col2[3]-col1[3])^2)
#       p=d/sqrt((255)^2+(255)^2+(255)^2)
#       if(p <= min_dist){
#         min_dist = p
#       }
#     }
#     palette_min = c(palette_min, cl[i])
#     palette_min_dist = c(palette_min_dist, min_dist)
#     
#   }
#   
#   max_dist = 0
#   max_i = 0
#   for(k in 1:length(palette_min)){
#     if(palette_min_dist[k] > max_dist){
#       max_dist = palette_min_dist[k]
#       new_color = palette_min[k]
#     }
#   }
#   palette = c(palette, new_color)
#   cl = setdiff(cl, new_color)
# }
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
p = plot(hell.hclust)


dev.off()





taxon_matrix <- otu_table(ps)
taxon_matrix  <- as.matrix(taxon_matrix)
taxon_matrix

result <- pvclust(as.matrix(taxon_matrix),
                  method.dist =  "correlation", 
                  method.hclust="ward.D2", nboot=1000)

svg("cluster_bootstrap1000.svg", width = 15, height = 8)
p = plot(result,
     main="Cluster analysis. Bootstrap = 1000"
) + pvrect(result, alpha=0.95, border = 2, lwd= 1, xpd = TRUE)


dev.off()

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps, method="unifrac", weighted=F)
ordination = ordinate(ps, method="PCoA", distance=wunifrac_dist)

svg("PCoA_unifrac.svg", width = 15, height = 8)
plot_ordination(ps, ordination, color=groupby, label="SampleName") + 
  theme(aspect.ratio=1) +  scale_colour_manual(values=palette) +
  ggtitle("PCoA with unweighted UniFrac distance")


dev.off()


#Do multivariate analysis based on Bray-Curtis distance and NMDS ordination.
ps.ord <- ordinate(ps, "NMDS", "bray")

svg("NMDS_bray.svg", width = 15, height = 8)
plot_ordination(ps, ps.ord, type="sample", color=groupby, 
                title="OTUs", label="SampleName") + geom_point(size=2) +
  scale_colour_manual(values=palette) +  ggtitle("NMDS with Bray-Curtis distance")


dev.off()


