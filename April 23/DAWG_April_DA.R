###DESeq2###

# Classify OTUs into Phylum
library('dplyr')
library('DESeq2')
library("ggplot2")

ps.bac_phylum <- tax_glom(ps, taxrank = "Phylum")
ps.bac_phylum <- subset_samples(ps.bac_phylum, fastq_prefix != "ASV25")

# Carry out binomial test for bacterial phyla using DESeq2 package.

Depth <- phyloseq_to_deseq2(ps.bac_phylum, ~ Depth)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(Depth), 1, gm_mean)
Depth = estimateSizeFactors(Depth, geoMeans = geoMeans)
Depth = DESeq(Depth, fitType="local")

res1 <- results(Depth, pAdjustMethod = "BH")
res1 <- res1[order(res1$padj, na.last=NA),]
alpha=0.05
sigtab1 = res1[(res1$padj < alpha),]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ps.bac_phylum)[rownames(sigtab1), ], "matrix"))
sigtab1 <- sigtab1[order(sigtab1$log2FoldChange),]

x = tapply(sigtab1$log2FoldChange, sigtab1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab1$Phylum = factor(as.character(sigtab1$Phylum), levels=names(x))
d3 <-ggplot(sigtab1, aes(y=reorder(Phylum,log2FoldChange), x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) + 
  scale_color_manual(values= c("Blue", "Red")) +
  theme_set(theme_bw()) + 
  labs(x ="log-2 Fold Change",y="Phylum")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=16),axis.text.y= element_text(size=16), legend.position="none") +
  theme(axis.title= element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d3

### ALDex2 - ANOVA-like differential expression tool 
###for high throughtput sequencing data

BiocManager::install("ALDEx2")
library(ALDEx2)

ASV_table <- otu_table(ps.bac_phylum)
aldex.all <- aldex(t(ASV_table), 
               conditions = sample_data(ps.bac_phylum)$Depth,
               mc.samples = 128, test= "t",
               effect =TRUE, denom ="all")
par(mfrow=c(1,2))
aldex.plot(aldex, type="MA", test="welch")
aldex.plot(aldex, type="MW", test="welch")

x <- aldex.clr(t(ASV_table), 
              sample_data(ps.bac_phylum)$Depth, mc.samples=128, denom="all", verbose=F)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
x.kw <- aldex.kw(x)

covariates <- data.frame("A" = Meta$Depth,
                         "B" = Meta$Time_binary)
mm <- model.matrix(~ A + B, covariates)
x <- aldex.clr(t(ASV_table), mm, mc.samples=128, denom="all")
glm.test <- aldex.glm(x, mm)

par(mfrow=c(1,2))
plot(glm.test[,"model.A2cm Pr(>|t|).BH"], aldex.all$we.eBH, log="xy",
     xlab="glm model A", ylab="Welch's t-test")
plot(glm.test[,"model.BOld Pr(>|t|).BH"], aldex.all$we.eBH, log="xy",
     xlab="glm model A", ylab="Welch's t-test")


