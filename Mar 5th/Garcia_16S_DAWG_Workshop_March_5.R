###First, create a working directory to save scripts and data.
getwd()
setwd("/Users/raymondo.garcia")

####Install and Load Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq","Biostrings","RDPutils","ggplot2","dplyr", "ggpubr", "vegan", "emmeans"))

#In case you recieved the following message "is not available for R version 4.0.2", use the code below
install.packages("remotes")
remotes::install_github("jfq3/RDPutils")
devtools::install_github("tidyverse/ggplot2")

#Load Packages
library(phyloseq)
library(Biostrings)
library(RDPutils)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(emmeans)

#Alternatively you can use p_load to load the packages installed, it can also try to install the packages if they are not install
#For the latter, include the: ..., install = TRUE
install.packages("pacman")
pacman::p_load(phyloseq, Biostrings, RDPutils, ggplot2, dplyr, ggpubr, vegan, emmeans)

####Import Metadata or Mapping File

#Read Metadata or Mapping file using read.table (). You could also use "Import Dataset"
Meta <- read.table("/Users/raymondo.garcia/Desktop/P. syringae/DAWG/Metadata.txt", header = TRUE, row.names = 1, stringsAsFactors = F, na.strings = c(".", "-"))

#Verify Metadata or Mapping file structure. We want to make sure that R assigned the variables correctly
str(Meta)

is.factor(Meta$Depth)
Meta$Depth <- as.factor(Meta$Depth) #changing to factor
is.numeric(Meta$Nitrate.N_PPM)
is.factor(Meta$Time)
Meta$Time <- as.factor(Meta$Time)
Meta$Replicate <- as.factor(Meta$Replicate)
Meta$CEC.Na_. <- as.numeric(Meta$CEC.Na_.) # nas produced because there were no values for certain entries
Meta$GypReq_Calc_Tons.AF <- as.numeric(Meta$GypReq_Calc_Tons.AF)

####Create a Phyloseq object. Here we merged the otu table, taxonomy table, and Metadata/Mapping File.
ps.bac <- merge_phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa.spp), sample_data(Meta))

#DADA2 outputs a OTU/ASV table in the sample (rows) by OTU/ASV (columns) format. Mnay packages might required a OTU/ASV (rows) by sample (columns) format, so to transpose use t(otu_table(nameoffile))
#Also, each column is a sequence, but we want to exchange the sequence for an 'ASV' label, to make things easier downstream
ps.bac@otu_table #to view your OTU table

####Duplicate the phyloseq object to preserve original
ps.bac1 <- ps.bac

####Change row names so that they are ASVs rather than individual nucleotide sequences
DNA <- Biostrings::DNAStringSet(taxa_names(ps.bac))
names(DNA) <- taxa_names(ps.bac)
ps.bac <- merge_phyloseq(ps.bac, DNA)
taxa_names(ps.bac) <- paste0("ASV", seq(ntaxa(ps.bac)))

ps.bac@otu_table #check to make sure it worked

####Order samples by the number of reads or counts or library size or depth of coverage
ordered(sample_sums(ps.bac))

#####Change the sample names so that they are shorter and are more descriptive
sample_names(ps.bac) <- ps.bac@sam_data$Name 

ps.bac@sam_data ##now check that it worked

####Save RDS and Read RDS
saveRDS(ps.bac, "Phyloseq-Bacteria-SoilSteaming.rds")

#Read RDS. By using <-, we are restoring the RDS under the same name we used for our Phyloseq object
ps.bac <- readRDS("~/Phyloseq-Bacteria-SoilSteaming.rds")

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

####Let's verify the prevalence of each taxon, in other words, the number of samples in which a taxon appears at least once.

#Prevalence of each feature and store as data.frame
prevdf <- apply(X = otu_table(ps.bac.filt),
               MARGIN = ifelse(taxa_are_rows(ps.bac.filt), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps.bac.filt),
                    tax_table(ps.bac.filt))

#Average prevalence and total read counts of the features in Phyla
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Plot the Prevalence vs Total Abundance or Counts
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps.bac.filt),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + #Include a guess for parameter for prevalence treshold, here is 5%
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#Define a prevalence treshold based on the plot, here we used 5% of total samples
prevalenceThreshold <- 0.05*nsamples(ps.bac.filt)

#Finally, let's filter based on the prevalence treshold established
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps.bac.filt.2 <- prune_taxa(keepTaxa, ps.bac.filt)

#Let's see how many taxa after filtering
length(get_taxa_unique(ps.bac.filt.2, taxonomic.rank = "Genus")) #587 Genera after filtering

####Agglomerate taxa at the Genus level and convert to long format for plotting.
ps.bac.filt.agglo <- tax_glom(ps.bac.filt.2, "Genus", NArm = TRUE) 

ps.bac.filt.agglo.long <- ps.bac.filt.agglo %>% #this is a function of the dplyr package not R base
  tax_glom(taxrank = "Phylum") %>%                     
  psmelt() %>%                                         
  arrange(Phylum) 

#####Create a stacked bar plot to observe community composition at the Phylum Level. For this we used counts, not relative abundance. 
phylum_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", 
                   "#5E738F","#D1A33D", "#8A7C64", "#599861","#999999","#E69F00","#56B4E9") 

Phylum_Stacked_Bar_Plot <- ggplot(ps.bac.filt.agglo.long, aes(x = SampleID, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position="stack") +
  theme_grey()+
  scale_fill_manual (values = phylum_colors) +
  theme() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Abundance") +
  xlab("Sample") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Phylum Composition")
print(Phylum_Stacked_Bar_Plot)

#Steaming seems to shfit the bacterial community, depleting the Phylum Proteobacteria and favoring Firmicutes
#Thus, let's subset for the Phylum Firmicutes
ps.bac.filt.agglo.firmicutes = subset_taxa(ps.bac.filt.agglo, Phylum=="Firmicutes")                 
ps.bac.filt.agglo.firmicutes@tax_table ###just to verify that it was subset appropiately 

ps.bac.filt.agglo.firmicutes.long <- ps.bac.filt.agglo.firmicutes %>%
tax_glom(taxrank = "Genus") %>%                     
  psmelt() %>%                                         
  arrange(Genus) 

#Create a stacked bar plot
Firmi_Stacked_Bar_Plot <- ggplot(ps.bac.filt.agglo.firmicutes.long, aes(x = SampleID, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position="stack") +
  theme_grey()+
  scale_color_brewer (palette = "Spectral") + #notice that here I change to scale_color_brewer instead of scale_fill_manual because there are more Families than colors in the vector Phylum colors I created earlier
  theme() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Abundance") +
  xlab("Sample") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle(" ")
print(Firmi_Stacked_Bar_Plot)

#the family Bacillaceae, is the most abundant post-treatment, which kind of make sense since this taxon holds a bunch of endospore forming bacteria, thus more resistant to abiotic stresses such as heat.

####Now, let's explore the alpha Diversity and have a better understanding of species richnness and evenness
#One option is to use the function plot_richness
plot_richness(ps.bac.filt.2) #not so pretty, right

#Let's choose some of the alpha diveristy indexes. I like ACE because it has a standard error bar; Chao1 estimate diversity from abundance data, highlights the importance of rare ASVs
plot_richness(ps.bac.filt.2, measures = c("Shannon", "Chao1", "Observed", "ACE"), x = "SampleID", color = "SampleID", title = "") #you could add shape = "", for example shape = "Depth"

#Alternatively we can output estimates of alpha diveristy/richness and creating a data frame.
Richness_DF <- estimate_richness(ps.bac.filt.2) 
Richness_DF$SampleID <- ps.bac.filt.2@sam_data$SampleID

##Create a box plot for each of the richness estimates.
A <- ggplot(Richness_DF, aes(x=SampleID, y=Fisher, fill=SampleID)) + geom_boxplot() + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5))
B <- ggplot(Richness_DF, aes(x=SampleID, y=Shannon, fill=SampleID)) + geom_boxplot() + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5))

install.packages("ggpubr") #arrange multiple ggplots over multiple pages
library(ggpubr)

Richness_Boxplot_Figure <- ggarrange(A, B, nrow = 1, ncol = 2, widths = 2, heights = 2, label.x = 0,
                                   labels = c("A","B"), common.legend = FALSE, legend = FALSE)
print(Richness_Boxplot_Figure)

#Here we verify the structure of the data with an histogram

#in order to decide if we need to use parametric or non-parametric analyses.
hist(Richness_DF$Fisher, main = "Fisher", xlab = "", breaks=10)
hist(Richness_DF$Shannon, main = "Shannon", xlab = "", breaks=10)

#if you want to observe all histograms simultaneously you can use
par(mfrow = c(4, 3)) #create a 4x3 plot environment
par(mfrow = c(1,1)) #to cancel it

#Perform a Shapiro-Wilk normality test.
shapiro.test(Richness_DF$Fisher) #The null hypothesis for this test is that the data are normally distributed, since results are significant we reject the null hypothesis.
shapiro.test(Richness_DF$Shannon) #normally distributed

#Analysis for non-normally distributed metrics
kruskal.test(Fisher ~ SampleID, data = Richness_DF) ##significant (p-value = 0.0011)

#Analysis for normally distributed metrics, here we will be using lm, but you could use other models such as lme(), glm()
linear.model <- lm(Shannon ~ SampleID, data = Richness_DF)

#Check for model assumptions (i.e. - Normality of residuals, equal variance, independency)
Richness_DF$Fitted <- fitted(linear.model)
Richness_DF$Residuals <- residuals(linear.model, type = "pearson")
qqnorm(Richness_DF$Residuals); qqline(Richness_DF$Residuals)

#Multiple comparison
compare.means = emmeans(linear.model, specs = "SampleID")
multcomp::cld(compare.means, Letters = c("abcdef"))

####Now, let's normalize data. Abundance value transformation to relative abundance or simple proportions. 
####Count proportions outperformed rarefied counts in most simulations due to better sensitivity, 
####but also suffered from a higher rate of false positives at larger values of effect size. 
ps.bac.filt.2.normalized = transform_sample_counts(ps.bac.filt.2, function (x) {x/ sum(x)}) ###counts are divided by the total library size. 
ps.bac.filt.2.normalized@otu_table #to make sure it worked or
View(ps.bac.filt.2.normalized)

###First, we fix sample_ID in sample data and then we ordinate at the samples level using PCoA and Bray-Curtis distances.  
sample_data(ps.bac.filt.2.normalized)$SampleID <- factor (sample_data(ps.bac.filt.2.normalized)$SampleID, levels = c("2cm.1D","2cm.1M","2cm.2M","2cm.5M","15cm.1D","15cm.1M","15cm.2M","15cm.5M","Negative_control"))

###Unconstrained Ordination and plot ordination.
PCoA <- ordinate (ps.bac.filt.2.normalized, method =  "PCoA", distance = "bray")

PCoA_Plot <- plot_ordination(physeq = ps.bac.filt.2.normalized, ordination = PCoA, color = "SampleID", title = ""
) + scale_fill_manual (values = phylum_colors) + geom_point(aes(color = SampleID), alpha = 0.7, size =4)
print(PCoA_Plot)

#Interestingly, samples are clustering based on time. So, let's plot time as color and depth as shape
PCoA_Plot_2 <- plot_ordination(physeq = ps.bac.filt.2.normalized, ordination = PCoA, color = "Time", shape = "Depth", title = ""
) + scale_fill_manual (values = phylum_colors) + geom_point(aes(color = Time), alpha = 0.7, size =4)
print(PCoA_Plot_2) #indeed, samples are grouping based on time after treatment

#Used this if you want to plot with ellipse of 0.95 confidence intervals. However, in this case we have to few points to calculate an ellipse
PCoA_Plot + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t") + theme_bw() 

set.seed(1) #to get repeatable p-values

#Calculate bray curtis distance matrix (i.e. - betadiversity)
Bray_Curtis <- phyloseq::distance(ps.bac.filt.2.normalized, method = "bray")

#Make a data frame from the sample data
Sample_DF <- data.frame(sample_data(ps.bac.filt.2.normalized))

install.packages("vegan")
library(vegan)

#Perform a multivariate analysis, specifically a PERMANOVA, using the function adonis from the vegan package. 
adonis(Bray_Curtis ~ SampleID, data = Sample_DF) #p-value = 0.001, so we reject the null hypothesis that our samples/treatments have similar bacterial community composition

#Beta-dispersion measures the compositional variation of the microbiome among a group of samples
Beta_Dispersion <- betadisper(Bray_Curtis, Sample_DF$SampleID)
permutest(Beta_Dispersion) #p-value = 0.831, meaning that we cannot reject the null hypothesis that our groups have the same dispersions. 
                           #Basically is telling us that our results from the PERMANOVA are real and not due to differences in disperions

####Constrained Ordination and plot ordination
#Subset samples to remove the negative control, it has missing values and the analysis doesn't work
ps.bac.filt.2.normalized.CAP <- subset_samples(ps.bac.filt.2.normalized, SampleID %in% c("2cm.1D","2cm.1M","2cm.2M","2cm.5M","15cm.1D","15cm.1M","15cm.2M","15cm.5M"))

#Subset taxa, let's use Phylum Firmicutes
ps.bac.filt.2.normalized.CAP.Firmi <- subset_taxa(ps.bac.filt.2.normalized.CAP, Phylum=="Firmicutes")

#Calculate the Bray-Curtis distances
Bray_Curtis_CAP <- phyloseq::distance(ps.bac.filt.2.normalized.CAP.Firmi, method = "bray")

#Ordinate
CAP <- ordinate(
  physeq = ps.bac.filt.2.normalized.CAP.Firmi, #you can use CCA or others constrained ordination analysis
  method = "CAP",
  distance = Bray_Curtis_CAP,
  formula = ~ Time + Depth) #our response vairiable is ASV

#To select the environmental factors or predictors that contribute to a greater extent to the variance observed on the bacterial community. 
#Here we will plot both regardless, but if you have multiple environmental covariates and you only want to model significant ones use this.
ordistep(CAP)

#Plot constrained ordination
CAP_Plot <- plot_ordination(
  physeq = ps.bac.filt.2.normalized.CAP.Firmi,
  ordination = CAP, color = "SampleID", title = "",
  axes = c(1,2)
) + scale_fill_manual(values = phylum_colors) #default is type = "samples", but if you include type = "split" you will create a split plot and you could include taxa (color = "Genus") in one side and samples (shape = "SampleID") in the other
                                              #if you include type = "biplot" it will plot taxa and samples in the same plot
#Add the environmental variables as arrows
Arrowenv <- vegan::scores(CAP, display = 'bp') 

#Add labels, make a data.frame
Arrowdf <- data.frame(labels = rownames(Arrowenv), Arrowenv)

#Define the arrow aesthetic mapping
Arrowmap <- aes(xend = CAP1,
                yend = CAP2,
                x = 0,
                y = 0,
                shape = NULL,
                color = NULL,
                label = labels)

#To make labels the same angle of the arrow. 
Arrowdf = transform(Arrowdf, radians = atan(CAP2/CAP1))
Arrowdf = transform(Arrowdf, angle = 360 * radians/(2 * pi))
#Quadrants II, III, IV
Arrowdf$quad234 <- apply(Arrowdf[, c("CAP1", "CAP2")], 1, function(x) {
  any(x < 0)
})
Arrowdf$quad4 <- apply(Arrowdf[, c("CAP1", "CAP2")], 1, function(x) {
  all(x < 0)
})
#If quadrant II, III, IV, add 180 to angle
if (any(Arrowdf$quad234)) {
  Arrowdf[Arrowdf$quad234, "angle"] <- Arrowdf[Arrowdf$quad234, "angle"] + 
    180
}
#If quadrant IV, add additional 180
if (any(Arrowdf$quad4)) {
  Arrowdf[Arrowdf$quad4, "angle"] <- Arrowdf[Arrowdf$quad4, "angle"] + 180
}
#For printing text, we want to flip if its greater than 90
if (any(Arrowdf$angle > 90)) {
  Arrowdf[Arrowdf$angle > 90, "angle"] <- Arrowdf[Arrowdf$angle > 90, "angle"] - 
    180
}

Labelmap <- aes(x = 1.3 * CAP1,
                y = 1.3 * CAP2,
                shape = NULL,
                color = NULL,
                label = labels,
                angle = angle)

Arrowhead = arrow(length = unit(0.01, "npc"))

#Make a new graphic
CAP_Plot + 
  geom_segment(
    mapping = Arrowmap,
    size = .5,
    data = Arrowdf,
    color = "black",
    arrow = Arrowhead
  ) + 
  geom_text(
    mapping = Labelmap,
    size = 2,
    data = Arrowdf,
    show.legend = FALSE
  )

#Lastly, let's do the same but this time plotting taxa. The result is messy because of the number of genera but I just wanted to provide you with an example
CAP_Plot_Taxa <- plot_ordination(
  physeq = ps.bac.filt.2.normalized.CAP.Firmi,
  ordination = CAP, type = "taxa", color = "Genus", title = "", #notice that we used type = "taxa"
  axes = c(1,2)
) + scale_fill_manual(values = phylum_colors)

#Add the environmental variables as arrows
Arrowenv <- vegan::scores(CAP, display = 'bp') 

#Add labels, make a data.frame
Arrowdf <- data.frame(labels = rownames(Arrowenv), Arrowenv)

#Define the arrow aesthetic mapping
Arrowmap <- aes(xend = CAP1,
                yend = CAP2,
                x = 0,
                y = 0,
                shape = NULL,
                color = NULL,
                label = labels)

#To make labels the same angle of the arrow. 
Arrowdf = transform(Arrowdf, radians = atan(CAP2/CAP1))
Arrowdf = transform(Arrowdf, angle = 360 * radians/(2 * pi))
#Quadrants II, III, IV
Arrowdf$quad234 <- apply(Arrowdf[, c("CAP1", "CAP2")], 1, function(x) {
  any(x < 0)
})
Arrowdf$quad4 <- apply(Arrowdf[, c("CAP1", "CAP2")], 1, function(x) {
  all(x < 0)
})
#If quadrant II, III, IV, add 180 to angle
if (any(Arrowdf$quad234)) {
  Arrowdf[Arrowdf$quad234, "angle"] <- Arrowdf[Arrowdf$quad234, "angle"] + 
    180
}
#If quadrant IV, add additional 180
if (any(Arrowdf$quad4)) {
  Arrowdf[Arrowdf$quad4, "angle"] <- Arrowdf[Arrowdf$quad4, "angle"] + 180
}
#For printing text, we want to flip if its greater than 90
if (any(Arrowdf$angle > 90)) {
  Arrowdf[Arrowdf$angle > 90, "angle"] <- Arrowdf[Arrowdf$angle > 90, "angle"] - 
    180
}

Labelmap <- aes(x = 1.3 * CAP1,
                y = 1.3 * CAP2,
                shape = NULL,
                color = NULL,
                label = labels,
                angle = angle)

Arrowhead = arrow(length = unit(0.01, "npc"))

#Make a new graphic
CAP_Plot_Taxa + 
  geom_segment(
    mapping = Arrowmap,
    size = .5,
    data = Arrowdf,
    color = "black",
    arrow = Arrowhead
  ) + 
  geom_text(
    mapping = Labelmap,
    size = 2,
    data = Arrowdf,
    show.legend = FALSE
  )
