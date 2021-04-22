# How to integrate environmental/climate variables with sequencing data

library(corrplot)
library(microbiome)
library(phyloseq)
library(vegan)

# How many samples do you have? How many variables do you have?
# if number of variables > number of samples, RDA/CCA are not appropriate

ps.bac <- subset_taxa(ps, Kingdom == "Bacteria")

ps.rel  = transform_sample_counts(ps.bac, function(x) x / sum(x) )

# first remove samples that do not have the ancillary data
ps.chem <- subset_samples(ps.rel, Manganese_PPM != "NA")
# let's start with a PCoA ####

bac.bray = ordinate(ps.chem, "PCoA", "bray")
p <- plot_ordination(ps.chem, bac.bray, "samples", color="Depth", shape = "Time")

# assess whether any chemical variables are associated with the distribution along the x-axis
A1 <- p$data$Axis.1

chemistry <- data.frame(sample_data(ps.chem))
chemistry[1:8] <- NULL
chemistry$CEC.H_. <- NULL
chemistry$Limestone_. <- NULL
chemistry$Moisture_Analysis <- NULL
chemistry$A1 <- A1
chemistry$LimeReq_Tons.AF <- NULL
chemistry$GypReq_Calc_Tons.AF <- NULL
chemistry$CEC.Na_. <- NULL
M <- cor(chemistry, use = "pairwise.complete.obs",method = "spearman")

corrplot(cor(chemistry[, 1:27])[27, 1:27, drop=FALSE], cl.pos='n', method = "number", number.cex = 0.6, tl.col = "black")

# check the entire correlation matrix
P <- cor.mtest(chemistry, conf.level = 0.95, method = "spearman")
corrplot(M, type = "upper", p.mat = P$p, method = "number", tl.col = "black", insig = "blank", tl.cex = 0.5, number.cex = 0.5)

# rda with stepwise selection
filt.asv <- data.frame(otu_table(ps.chem))
# check axis lengths (if > 4, transform via hellinger or use CCA)
decorana(filt.asv)
chemistry$A1 <- NULL
chemistry.log  <- log(chemistry+1)
rda.1 <- rda(decostand(filt.asv, 'hellinger') ~ ., chemistry.log, scaling = 'sites')
envfit(rda.1 ~., chemistry)

# let's keep significant variables
chemistry.log[1:2] <- NULL
chemistry.log[3:7] <- NULL
chemistry.log[5:9] <- NULL
chemistry.log$CEC_meq.100g <- NULL
chemistry.log[7:10] <- NULL
chemistry.log$Moisture_. <- NULL

rda.2 <- rda(decostand(filt.asv, 'hellinger') ~ ., chemistry.log, scaling = 'sites')

safe_colorblind <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                     "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

metadata <- data.frame(sample_data(ps.chem))
pch2 <- c(15, 16)

plot(rda.2, type="n", display="sites",xlab = paste("PCA1"),
     ylab = paste("PCA2"), cex.lab = 1.5, cex.axis = 1.2)
points(rda.2, display="bp", lwd=2, col="black")
text(rda.2, display="bp", col="black", font=2, cex = 0.5)
points(rda.2, display="sites", col=safe_colorblind[metadata$Time],cex=3, font=1, pch = pch2[metadata$Depth])
legend("bottomright", legend = levels(metadata$Time), col = safe_colorblind[metadata$Time], pch = 16, cex = 1, title = "Time")
legend("bottomright", legend = levels(metadata$Depth), col = "black", pch = c(15,16), cex = 1, title = "Depth", inset = c(0.37, 0.00))

# Mantel Tests - Bacteria ####

x.dist<- vegdist(chemistry$Sulfate_meq.L, "euclidean")
y.dist <- vegdist(filt.asv, "bray")

mantel(x.dist, y.dist, method ="spearman", permutations = 999)

lm.all <- lm(y.dist~x.dist)
plot(y.dist~x.dist, ylab = "Bray-Curtis Dissimilarity (16S)", xlab = "Sulfate (meq)", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )

x.dist<- vegdist(chemistry, "euclidean")
y.dist <- vegdist(filt.asv, "bray")

mantel(x.dist, y.dist, method ="spearman", permutations = 999)

lm.all <- lm(y.dist~x.dist)
plot(y.dist~x.dist, ylab = "Bray-Curtis Dissimilarity (16S)", xlab = "Soil chemistrt dissimilarity", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )

# integrate into alpha-divesity estimates ####
ps.bac.chem <- subset_samples(ps.bac, Chemistry_data == "Y")
ps.bac.chem.1 = filter_taxa(ps.bac.chem, function(x) sum(x) > 0, TRUE)

ordered(sample_sums(ps.bac.chem.1))
set.seed(500)

ps.bac.rare <-rarefy_even_depth(ps.bac.chem.1, sample.size = 22000)

## Relative abundance calculation##

# rarefied to 5000 sequences per sample 
alpha.div.bac <-estimate_richness(ps.bac.rare, measures=c("Shannon", "Observed"))

even <- evenness(ps.bac.rare, 'pielou')

chemistry$Evenness <- even$pielou
chemistry$Shannon <- alpha.div.bac$Shannon
chemistry$Richness <- alpha.div.bac$Observed

M <- cor(chemistry, use = "pairwise.complete.obs",method = "spearman")
P <- cor.mtest(chemistry, conf.level = 0.95, method = "spearman")

corrplot(M, type = "upper", p.mat = P$p, method = "number", tl.col = "black", insig = "blank", tl.cex = 0.5, number.cex = 0.5)

corrplot(cor(chemistry[, 1:29])[27:29, 1:29, drop=FALSE], cl.pos='n', method = "number", number.cex = 0.6, tl.col = "black")

# We can do the same thing with our fungal dataset ####
ps.fun <- readRDS("ps_fun.rds")

ps.fun.rel  = transform_sample_counts(ps.fun, function(x) x / sum(x) )

# first remove samples that do not have the ancillary data 
ps.fun.chem <- subset_samples(ps.fun.rel, Chemistry_data == "Y")
ps.fun.chem.1 = filter_taxa(ps.fun.chem, function(x) sum(x) > 0, TRUE)

# Ordinations: PCoA ASV-Level ####

fun.bray = ordinate(ps.fun.chem.1, "PCoA", "bray")
plot_ordination(ps.fun.chem.1, fun.bray, "biplot", color="Depth", shape = "Time")
p <- plot_ordination(ps.fun.chem.1, fun.bray, "sites", color="Depth", shape = "Time")
p

# assess whether any chemical variables are associated with the distribution along the x-axis
A1 <- p$data$Axis.1

fun.chemistry <- data.frame(sample_data(ps.fun.chem))
fun.chemistry[1:8] <- NULL
fun.chemistry$CEC.H_. <- NULL
fun.chemistry$Limestone_. <- NULL
fun.chemistry$Moisture_Analysis <- NULL
fun.chemistry$A1 <- A1
fun.chemistry$LimeReq_Tons.AF <- NULL
fun.chemistry$GypReq_Calc_Tons.AF <- NULL
fun.chemistry$CEC.Na_. <- NULL
M <- cor(fun.chemistry, use = "pairwise.complete.obs",method = "spearman")

corrplot(cor(fun.chemistry[, 1:27])[27, 1:27, drop=FALSE], cl.pos='n', method = "number", number.cex = 0.6, tl.col = "black")

# check the entire correlation matrix
P <- cor.mtest(fun.chemistry, conf.level = 0.95, method = "spearman")
corrplot(M, type = "upper", p.mat = P$p, method = "number", tl.col = "black", insig = "blank", tl.cex = 0.5, number.cex = 0.5)

# rda with stepwise selection
filt.fun.asv <- data.frame(otu_table(ps.fun.chem.1))
# check axis lengths (if > 4, transform via hellinger or use CCA)
decorana(filt.fun.asv)
fun.chemistry$A1 <- NULL
fun.chemistry.log  <- log(fun.chemistry+1)
rda.1 <- rda(decostand(filt.fun.asv, 'hellinger') ~ ., fun.chemistry.log, scale = TRUE)
envfit(rda.1 ~., chemistry)


fun.chemistry.log[1:2] <- NULL
fun.chemistry.log[3:7] <- NULL
fun.chemistry.log[5:9] <- NULL
fun.chemistry.log[6:7] <- NULL
fun.chemistry.log$pH <- NULL
fun.chemistry.log$Saturation_. <- NULL

rda.2 <- rda(decostand(filt.fun.asv, 'hellinger') ~ ., fun.chemistry.log)


fun.metadata <- data.frame(sample_data(ps.fun.chem))

plot(rda.2, type="n", display="sites",xlab = paste("PCA1"),
     ylab = paste("PCA2"), cex.lab = 1.5, cex.axis = 1.2)
points(rda.2, display="bp", lwd=2, col="black")
text(rda.2, display="bp", col="black", font=2, cex = 0.5)
points(rda.2, display="sites", col=safe_colorblind[fun.metadata$Time],cex=3, font=1, pch = pch2[fun.metadata$Depth])
legend("bottomright", legend = levels(fun.metadata$Time), col = safe_colorblind[fun.metadata$Time], pch = 16, cex = 0.7, title = "Time")
legend("bottomright", legend = levels(fun.metadata$Depth), col = "black", pch = c(15,16), cex = 0.7, title = "Depth", inset = c(0, 0.2))


# mantel test

x.dist<- vegdist(fun.chemistry$Sulfate_meq.L, "euclidean")
y.dist <- vegdist(filt.fun.asv, "bray")

mantel(x.dist, y.dist, method ="spearman", permutations = 999)

lm.all <- lm(y.dist~x.dist)
plot(y.dist~x.dist, ylab = "Bray-Curtis Dissimilarity (ITS)", xlab = "Sulfate (meq)", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )


x.dist<- vegdist(fun.chemistry, "euclidean")
y.dist <- vegdist(filt.fun.asv, "bray")

mantel(x.dist, y.dist, method ="spearman", permutations = 999)

lm.all <- lm(x.dist~y.dist)
plot(y.dist~x.dist, ylab = "Bray-Curtis Dissimilarity (ITS)", xlab = "Soil chemistry distance", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )


# fungal ordinations at genus-level ####
ps.fun.gen <- tax_glom(ps.fun.chem, 'Genus')

fun.bray = ordinate(ps.fun.gen, "PCoA", "bray")
p <- plot_ordination(ps.fun.gen, fun.bray, "sites", color="Depth", shape = "Time")
p

# rda with stepwise selection
filt.fun.gen <- data.frame(otu_table(ps.fun.gen))
# check axis lengths (if > 4, transform via hellinger or use CCA)
decorana(filt.fun.gen)
fun.chemistry.log  <- log(fun.chemistry+1)
rda.1 <- rda(decostand(filt.fun.gen, 'hellinger') ~ ., fun.chemistry.log, scale = TRUE)
envfit(rda.1 ~., chemistry)


fun.chemistry.log[1:3] <- NULL
fun.chemistry.log$Exch.Ca_PPM <- NULL
fun.chemistry.log$Exch.Mg_PPM<- NULL
fun.chemistry.log$Exch.Na_PPM <- NULL
fun.chemistry.log[6:10] <- NULL
fun.chemistry.log[7:8] <- NULL
fun.chemistry.log$pH <- NULL
fun.chemistry.log[8:12]<- NULL

rda.2 <- rda(decostand(filt.fun.gen, 'hellinger') ~ ., fun.chemistry.log, scaling = 'sites')

plot(rda.2, type="n", display="sites",xlab = paste("PCA1"),
     ylab = paste("PCA2"), cex.lab = 1.5, cex.axis = 1.2)
points(rda.2, display="bp", lwd=2, col="black")
text(rda.2, display="bp", col="black", font=2, cex = 0.5)
points(rda.2, display="sites", col=safe_colorblind[fun.metadata$Time],cex=3, font=1, pch = pch2[fun.metadata$Depth])
text(rda.2, display = 'species', col = "#332288", cex = 0.5)
legend("bottomright", legend = levels(fun.metadata$Time), col = safe_colorblind[fun.metadata$Time], pch = 16, cex = 0.7, title = "Time")
legend("bottomright", legend = levels(fun.metadata$Depth), col = "black", pch = c(15,16), cex = 0.7, title = "Depth", inset = c(0, 0.2))


# exporting for funguild ####
# requires at least genus classification
ps.fun.no.gen <- subset_taxa(ps.fun, Genus != "NA")
fun.asv<- data.frame(otu_table(ps.fun.no.gen))
fun.taxa <- data.frame(tax_table(ps.fun.no.gen))
fun.taxa$Taxonomy <- unite(fun.taxa, Taxonomy, sep = ';')
fun.asv.t <- data.frame(t(fun.asv))

fun.taxa$Genus.2 <- gsub('g__', '', fun.taxa$Genus)
fun.taxa$Species.2 <- gsub('s__', '', fun.taxa$Species)

fun.taxa$Genus.sp <- paste(fun.taxa$Genus.2, fun.taxa$Species.2, sep = "_")
fun.taxa$Genus.sp <- paste0("|", fun.taxa$Genus.sp, "|")
fun.taxa$Genus.2 <- NULL
fun.taxa$Species.2 <- NULL
fun.taxa$taxonomy <- paste(fun.taxa$Genus.sp, fun.taxa$Taxonomy$Taxonomy, sep ="")

fun.asv.t$taxonomy <- fun.taxa$taxonomy
write.table(fun.asv.t, file = "Fungal-Steam-ASV-Table1.txt")


  