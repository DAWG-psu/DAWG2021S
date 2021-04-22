###First, create a working directory to save scripts and data.
getwd()
setwd("/Users/raymondo.garcia")

###Read Metadata or Mapping file using read.table (). You could also use "Import Dataset"
Meta <- read.table("/Users/raymondo.garcia/Desktop/P. syringae/DAWG/Metadata.txt", header = TRUE, row.names = 1, stringsAsFactors = F, na.strings = c(".", "-"))

####Install and Load Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("nlme","ggplot2","dplyr","ggpubr","emmeans"))

install.packages("pacman")
pacman::p_load(nlme, ggplot2, dplyr, ggpubr, emmeans)

###Verify Metadata or Mapping file structure. We want to make sure that R assigned the variables correctly
str(Meta)

is.factor(Meta$Depth)
Meta$Depth <- as.factor(Meta$Depth) #changing to factor
is.numeric(Meta$Nitrate.N_PPM)
is.factor(Meta$Time)
Meta$Time <- as.factor(Meta$Time)
Meta$Replicate <- as.factor(Meta$Replicate)
Meta$CEC.Na_. <- as.numeric(Meta$CEC.Na_.) # nas produced because there were no values for certain entries
Meta$GypReq_Calc_Tons.AF <- as.numeric(Meta$GypReq_Calc_Tons.AF)

###Create a Box Plot for pH and Nitrate values.
pH_BoxPlot <- ggplot(Meta, aes(x=SampleID, y=pH, fill=SampleID)) + geom_boxplot() + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5))
print(pH_BoxPlot)

Nitrate_BoxPlot <- ggplot(Meta, aes(x=SampleID, y=Nitrate.N_PPM, fill=SampleID)) + geom_boxplot() + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5))
print(Nitrate_BoxPlot)

###Here we verify the structure of the data with an histogram in order to know if we could use a parametric or non-parametric analyses or perhaps transform our data
hist(Meta$Nitrate.N_PPM, main = "Nitrate.N_PPM", xlab = "", breaks=5)

###Perform a Shapiro-Wilk normality test.
shapiro.test(Meta$Nitrate.N_PPM) #The null hypothesis for this test is that the data is normally distributed, since results are not significant we can't reject the null hypothesis.

###If your data was non-normally distributed, some options to analyze it are the non-parametric kruskal-wallis (for more than two independent variables) and wilcoxon rank sum test (for two independent variables)
kruskal.test(Nitrate.N_PPM ~ SampleID, data = Meta)
wilcox.test()

###Analysis for normally distributed metrics, here we will be using lm, but you could use more complex models, such as lme()
linear.model <- lm(Nitrate.N_PPM ~ SampleID, data = Meta)

ANOVA_1 <- aov(linear.model, data = Meta) ##pretty simple. Assuming there were significant differences, ANOVA would only tell you that there were signficant differences, but not between what independent variables.
summary(ANOVA_1)

summary(linear.model) ##more informative

###Check for model assumptions (i.e. - Normality of residuals, equal variance, independency, etc.)
Meta$Fitted <- fitted(linear.model)
Meta$Residuals <- residuals(linear.model, type = "pearson")
qqnorm(Meta$Residuals); qqline(Meta$Residuals)

###Multiple comparison
compare.means = emmeans(linear.model, specs = "SampleID")
multcomp::cld(compare.means, Letters = c("abcdef"))

###If your residuals were non-normally distributed then you could have used Generalized Linear Models, more flexible distributions, such as Negative binomial, overall fit microbiome data better.
glm()
