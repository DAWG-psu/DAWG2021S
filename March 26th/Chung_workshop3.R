# Install Required packages
install.packages("devtools")
library(devtools)
install.packages("igraph")
library(igraph)
install.packages("qgraph")
library(qgraph)
install.packages("vegan")
library(vegan)
install.packages("MCL")
library(MCL)
install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
install.packages("intergraph")
library(intergraph)
install.packages("ggnet")
library(ggnet)


setwd("~/Desktop/DAWG/2021/")
#Read RDS. By using <-, we are restoring the RDS under the same name we used for our Phyloseq object
ps.bac <- readRDS("Phyloseq-Bacteria-SoilSteaming.rds")
####Now, let's start Filtering, Subsetting, and Performing Richnness, Ordination, and Multivariate analysis!
#First, show available taxa ranks in the dataset and create a table with number of features for each taxa.
rank_names(ps.bac)
table(tax_table(ps.bac)[, "Phylum"], exclude = NULL) #For Phylum

#1,549 features were annotated with a Phylum NA, these features might be artifacts and should be filtered. 
#Other Phyla have only one features and should be remove too.
table(tax_table(ps.bac)[, "Kingdom"], exclude = NULL) #For Domain/Kingdom

####Remove Phyla with low number of features and the Domain Archea. 
filterPhyla <- c("Dadabacteria", "Thermotogota", "WS2", "Zixibacteria", "NA", "Cyanobacteria") #Notice that we removed Cyanobacteria, as it might be chloroplasts (released during rhizodeposition)
filterDomain <- c("Archaea", "<NA>")

#Filtering
ps.bac.filt <- subset_taxa(ps.bac, !Phylum %in% filterPhyla)
table(tax_table(ps.bac.filt)[, "Phylum"], exclude = NULL) #to verify if it worked

ps.bac.filt <- subset_taxa(ps.bac.filt, !Kingdom %in% filterDomain)
table(tax_table(ps.bac.filt)[, "Kingdom"], exclude = NULL)

####Let's remove features with ambiguous Phylum annotation, such as "uncharacterized".
ps.bac.filt@tax_table #to view our tax_table

#Great, there are no ambigous Phylum annotations. If there were some, we could have used the code below
ps.bac.filt <- subset_taxa(ps.bac.filt, !is.na(Phylum) & !Phylum %in% c(" ", "uncharacterized"))

###Network analysis ###
ps.bac.filt



##data processing##
ps_reduced <- prune_taxa(taxa_sums(ps.bac.filt) > 1000, ps.bac.filt)
ps_reduced_f <- microbiomeutilities::format_to_besthit(ps_reduced)
otu.c <- t(otu_table(ps_reduced_f)@.Data) #extract the otu table from phyloseq object
tax.c <- as.data.frame(tax_table(ps_reduced_f)@.Data)#extract the taxonomy information
head(otu.c)
head(tax.c)


##spiec.easi##
set.seed(1244)
net.c <- spiec.easi(otu.c, method='mb', icov.select.params=list(rep.num=5)) # reps have to increases for real data
n.c <- symBeta(getOptBeta(net.c))
colnames(n.c) <- rownames(n.c) <- colnames(otu.c)
vsize <- log2(apply(otu.c, 2, mean)) # add log abundance as properties of vertex/nodes.

ig.mb <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
ig.mb # we can see all the attributes and weights

E(ig.mb)[weight > 0]$color<-"steelblue" #now color the edges based on their values positive is steelblue
E(ig.mb)[weight < 0]$color<-"orange"  #now color the edges based on their values
coords.fdr = layout_with_fr(ig.mb)
plot(ig.mb, layout=coords.fdr, vertex.size = 2, vertex.label.cex = 0.5)

mb.network <- asNetwork(ig.mb)
network::set.edge.attribute(mb.network, "color", ifelse(mb.network %e% "weight" > 0, "steelblue", "orange"))

phyla <- map_levels(colnames(otu.c), from = "best_hit", to = "Phylum", tax_table(ps_reduced_f))
mb.network %v% "Phylum" <- phyla
mb.network %v% "nodesize" <- vsize

mycolors <- scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",'red',"blue","green","purple"))

p1 <- ggnet2(mb.network, node.color = "Phylum", 
            label = TRUE, node.size = "nodesize", 
            label.size = 2, edge.color = "color") + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors

p1

stl.mb <- degree.distribution(ig.mb)
plot(0:(length(stl.mb)-1), stl.mb, ylim=c(0,.08), type='p', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# we will look at only taxa connect more than 10 others
p1 <- ggnet2(mb.network, node.color = "Phylum", 
            label = TRUE, 
            label.size = 3, edge.color = "color",
            size = "degree", size.min = 50) + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors
p1


betaMat=as.matrix(symBeta(getOptBeta(net.c)))
# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 
total=length(betaMat[betaMat!=0])/2 
print(list(positive, negative, total))

###### End of SPEIC-EASI tutorial ######

##### SparCC tutorial - from SPEIC-EASI package####
##Original SparCC was developed under python environment###

set.seed(1248)
net.c_sparcc <- sparcc(otu.c)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(net.c_sparcc$Cor) >= 0.7
diag(sparcc.graph) <- 0
library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
ig.sparcc <- adj2igraph(sparcc.graph)
library(igraph)
plot(ig.sparcc, vertex.size=vsize, vertex.label=NA, main="sparcc")

sparcc.network <- asNetwork(ig.sparcc)
network::set.edge.attribute(sparcc.network, "color", ifelse(sparcc.network %e% "weight" > 0, "steelblue", "orange"))

phyla <- map_levels(colnames(otu.c), from = "best_hit", to = "Phylum", tax_table(ps_reduced_f))
sparcc.network %v% "Phylum" <- phyla
sparcc.network %v% "nodesize" <- vsize

mycolors <- scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",'red',"blue","green","purple"))

p <- ggnet2(sparcc.network, node.color = "Phylum", 
            label = TRUE, node.size = "nodesize", 
            label.size = 2, edge.color = "color") + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors

p 

stl.mb <- degree.distribution(sparcc.network)
plot(0:(length(stl.mb)-1), stl.mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# we will look at only taxa connect more than 20 others
p <- ggnet2(sparcc.network, node.color = "Phylum", 
            label = TRUE, 
            label.size = 3, edge.color = "color",
            size = "degree", size.max = 30) + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolors
p

##End of sparCC tutorial##
## CoNet turorial using Cytoscape
