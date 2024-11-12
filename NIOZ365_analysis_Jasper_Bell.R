######           Phyloseq analysis of CASCABEL first run output                   #####

# Imports
library(dada2)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(readxl)
#library(writexl)
library(dplyr)
library(vegan)
library(DESeq2)
library(reshape2)
library(lme4)
library(forcats)
library(tidyverse)
library(RColorBrewer)
set.seed(5123)

#Loading datasets, change file locations to your own
samdf <- read.csv("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/NIOZ365_mapping_file.txt", sep="\t") #Matrix with samples as rows and metadata as columns. Mapping file from CASCABEL
seqtab.nochim <- read.csv("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/asv_table.txt", sep="\t") #Matrix with asvs (numbers) as rows and samples as columns. Numbers represent read counts. Taken from CASCABEL
taxa <- read.csv("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/representative_seq_set_tax_assignments.txt", sep = "\t") #matrix with asvs (numbers) as rows and taxonomic levels as columns. Taken from CASCABEL
diet <- read_excel("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/Diet_NIOZ365.xlsx") #Diet fraction (columns) per sample (rows). From Jan 
asv_proportions <- data.frame(read_excel("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/asv_proportions.xlsx")) #Read counts recalculated as fraction per sample. Based on seqtab.nochim, dataframe made in "Setting up PCA"
seqtab.long <- read_excel("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/seqtab_long.xlsx") #Long format of seqtab.nochim. Each row is a combination of a sample and an asv, has read count and phylum information. Made in "Exploring ASV table"
core_microbiome <- read_excel("~/Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/Core_microbiome.xlsx") #Only ASVs present in at least 50% of samples. ASVs as rows, samples as columns. Made in "Defining the core microbiome"

# Making sure row names are correct in all dataframes and that columns have the right data types
seqtab.nochim = seqtab.nochim %>% mutate_if(is.integer, as.numeric)
rownames(seqtab.nochim) = seqtab.nochim[,1]
seqtab.nochim <- seqtab.nochim[,-1]

rownames(taxa) = taxa[,1]
taxa <- taxa[,-1]

rownames(asv_proportions) <- rownames(taxa)

rownames(samdf) = samdf[,1]
samdf <- samdf[,-1]
samdf$Age <- as.factor(samdf$Age)

rownames(asv_proportions) <- rownames(seqtab.nochim)

rownames(core_microbiome) <- core_microbiome$ASV

#Adding a column called Bill_class where bill length <35 is short and >=35 is long
for (bird in 1: nrow(samdf)){
  if (!(is.na(samdf[bird, "Bill"]))) {
    if (samdf[bird, "Bill"] < 35 ) {
      samdf[bird, "Bill_class"] = "Short"
    }
    if (samdf[bird, "Bill"] >= 35 ) {
      samdf[bird, "Bill_class"] = "Long"
    }
  }
}

# Constructing a phyloseq object from dad2 outputs 
## I Don't understand why changing to a matrix makes this work, but it does
ps <- phyloseq(otu_table(as.matrix(seqtab.nochim), taxa_are_rows = TRUE),
               sample_data(samdf),
               tax_table(as.matrix(taxa)))

#Adding the proportion of seagrass, loripes and dosinia in the diet of each bird to the metadate of ps.rarefied. 
## This is so complicated because I only have both sample and ringnumber together in samdf, very smart of me to not add ringnumbers directly to the mapping file 
for (bird in 1:nrow(sample_data(ps))){
  ringnr <- rownames(sample_data(ps))[bird]
  ringnr <- samdf[ringnr,"Description"]
  if (ringnr %in% diet$Ring){
    loripes_perc <- diet[diet$Ring ==ringnr,][[2]]
    seagrass_perc <- diet[diet$Ring ==ringnr,][[4]]
    otherBivalves_perc <- diet[diet$Ring ==ringnr,][[6]] 
    sample_data(ps)[bird, "Seagrass"] <- seagrass_perc
    sample_data(ps)[bird, "Loripes"] <- loripes_perc
    sample_data(ps)[bird, "OtherBivalves"] <- otherBivalves_perc
  }
}

# Removing neg extraction controls and empty samples ("NIOZ365.0059.0068" & "NIOZ365.0195.0196")
## I also exclude "NIOZ365.0191.0008" as this is a weird sample that should not exist. Probably formed due to a sequencing error combined with unintentional inclusion of the barcode sequence in the mapping file
## This leaves me with a total of 83 samples with at least 1 read
neg_controls <- c("NIOZ365.0013.0022", "NIOZ365.0033.0042", "NIOZ365.0055.0064", "NIOZ365.0083.0092", "NIOZ365.0105.0114", "NIOZ365.0135.0144", "NIOZ365.0163.0172", "NIOZ365.0183.0192", "", "NIOZ365.0191.0008", "NIOZ365.0059.0068", "NIOZ365.0195.0196")
ps <- prune_samples(!(sample_names(ps) %in% neg_controls), ps)

summarize_phyloseq(ps)

otu_tab <-t(abundances(ps))

# Plotting Rarefaction plot, 
## hashtagged because takes some time to run on full dataset
#p <- vegan::rarecurve(otu_tab, step = 50, label = FALSE, sample = min(rowSums(otu_tab), Col ="red", cex = 0.5))

## Plotting the frequency distribution of number of reads in samples
ggplot(data.frame(rowSums(otu_tab)), aes(x=rowSums.otu_tab.)) + geom_histogram(bins = 35)

#Adding Alpha diversity measures to sample_data
erich = estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
sample_data(ps)$Observed<-erich$Observed
sample_data(ps)$Shannon<-erich$Shannon
sample_data(ps)$Simpson<-erich$Simpson
## Plotting the distribution of number of species per sample
ggplot(sample_data(ps), aes(x=Observed)) + geom_histogram(bins = 35)


#####                     Exploring asv table                                     #####
#Number of occurrences of each class in total dataset
ggplot(data=taxa, aes(y=fct_rev(fct_infreq(Phylum)), fill=Phylum)) +geom_bar() + 
  theme(legend.position="none", text=element_text(size=14), axis.text = element_text(size=12)) + xlab("Number of ASVs") + ylab("Phylum")# + ggtitle("Number of ASVs per phylum in dataset")
ggplot(data=taxa, aes(y=fct_rev(fct_infreq(Class)), fill=Phylum)) +geom_bar() #Seperating deeper than Class becomes really ugly

#Plotting average abundance per ASV
## Adding taxa isn't necessary, but its nice for the fill of the bars, though maybe the fill also distracts from the point of the plot
ASV.abundance.taxa <- cbind(asv_proportions, taxa)
ASV.abundance.taxa$Mean_abundance <- rowMeans(asv_proportions)
ASV.abundance.taxa$ASV <- rownames(ASV.abundance.taxa)
ggplot(data = ASV.abundance.taxa[1:83,], aes(y=Mean_abundance, x=reorder(ASV, - Mean_abundance), fill=Phylum)) + geom_col() + 
  theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5, size=12)) + xlab("ASV") + ylab("Mean relative abundance")

#Percent bar chart all samples, colored based on phylum
## Hashtagged because it takes a long time to run
## Resulting file of this block is already imported as "seqtab.long" from seqtab_long.xlsx
#seqtab.nochim2 <- seqtab.nochim[,!(names(seqtab.nochim) %in% neg_controls)] #Using only samples that are non-controls and have >0 reads
#seqtab.long <- data.frame(t(seqtab.nochim2))
#seqtab.long$Sample <- rownames(seqtab.long)
#seqtab.long <- pivot_longer(data=seqtab.long, cols=asv.1:asv.9835, names_to = "ASV", values_to="Abundance")
#for (row in 1:nrow(seqtab.long)){
#  asv.name <- seqtab.long[row, "ASV"]
#  asv.index <- match(asv.name, rownames(taxa))
#  seqtab.long[row, "Phylum"] <- taxa[asv.index, "Phylum"] 
#}

##Hashtagged because it takes some time to plot
#ggplot(data=seqtab.long, aes(x=Sample, y=Abundance, fill=Phylum)) + 
 # geom_bar(position="fill", stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
  #scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(seqtab.long$Phylum)), "Paired"))(length(unique(seqtab.long$Phylum))))
#####################################################################################
#####                             Alpha div analysis                         ########
#####################################################################################
#####     Plotting alpha diversity for different bird qualities all data          ##### 
plot_richness(ps, x="Age", measures=c("Observed", "Shannon", "Simpson"), color="Age") + geom_boxplot(alpha=0.6)+theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=1,vjust=1,size=12))
plot_richness(subset_samples(ps, Sex != ""), x="Sex", measures=c("Observed", "Shannon", "Simpson"), color="Sex") + geom_boxplot(alpha=0.6) + theme(axis.title=element_text(size = 14), legend.position="none", axis.text.x=element_text(angle=0, hjust=1,vjust=1,size=14), axis.text.y=element_text(size=14))
plot_richness(ps, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Expedition") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps, x="Bill_class", measures=c("Observed", "Shannon", "Simpson"), color="Bill_class") + geom_boxplot(alpha=0.6)+theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=0.5, vjust=1,size=12))

#Alpha diversity vs seagrass diet component
ggplot(data=as.data.frame(sample_data(ps)), aes(x=Seagrass, y=Simpson, color=Bill_class)) + geom_point(size=3) + 
  geom_smooth(method = "lm", formula= y~x) + 
  scale_color_discrete(name = "Bill length", labels = c("Long: > 35 mm", "Short: < 35 mm")) + 
  theme(axis.title=element_text(size = 14), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.title = element_text(size =13), legend.text = element_text(size=11)) +
  ylab("Simspon diversity score") + xlab("Zostera proportion in diet")

#####            Alpha diversity for full dataset                                 #####
#sample_data has to be extracted from the phyloseq object, otherwise lm() wont recognize it
data.alpha <- cbind(sample_data(ps))
#Making a linear model of the observed number of ASVs when bill length increases (richness vs Bill)
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness) #Based on Simpson measure
#Difference between Bill <35mm and Bill > 35mm
wilcox.test((subset(sample_data(data.alpha), Bill_class =="Long"))$Simpson, (subset(sample_data(data.alpha), Bill_class =="Short"))$Simpson, paired=FALSE)
#Influence of sex
wilcox.test((subset(sample_data(data.alpha), Sex =="M"))$Simpson, (subset(sample_data(data.alpha), Sex =="F"))$Simpson, paired=FALSE)
#Influence of age
wilcox.test((subset(sample_data(data.alpha), Age ==1))$Simpson, (subset(sample_data(data.alpha), Age ==3))$Simpson, paired=FALSE)
#Influence of expedition
kruskal.test(Simpson ~ as.factor(Expedition), data.alpha)
#Correlation with seagrass proportion
Simpson.Seagrass.lm <- lm(Simpson ~ Seagrass, data = data.alpha)
summary(Simpson.Seagrass.lm)

#####           Splitting dataset based on expedition                             #####  
# Making separate phyloseq objects for each expidition for ease of use downstream
ps.2019 <- subset_samples(ps, Expedition =="Mauritania_autumn_2019")
ps.2021 <- subset_samples(ps, Expedition =="Mauritania_autumn_2021")
ps.2022 <- subset_samples(ps, Expedition =="Mauritania_spring_2023") #For some reason dates from 2022 have been marked as 2023 in the expedition name. Does not impact analysis any further
# KEEP IN MIND: in the following 3 analyses I overwrite "data.alpha" multiple times, always make sure to update to the correct versions before use

# =========================  Autumn 2019 analysis                                 =================================================
plot_richness(ps.2019, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.2019, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))

#Bill length as an explanation for observed Simpson diversity
data.alpha <- cbind(sample_data(ps.2019))
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness)
#Difference between Bill <35mm and Bill > 35mm
wilcox.test((subset(sample_data(ps.2019), Bill_class =="Long"))$Simpson, (subset(sample_data(ps.2019), Bill_class =="Short"))$Simpson, paired=FALSE)
#Influence of sex --> Is significant, male and female means also caluclated
wilcox.test((subset(sample_data(ps.2019), Sex =="M"))$Simpson, (subset(sample_data(ps.2019), Sex =="F"))$Simpson, paired=FALSE)
mean(subset(sample_data(ps.2019), Sex =="M")$Simpson) #Though significant, difference in mean Simspon score is minimal (especialy compared to sd() )
mean(subset(sample_data(ps.2019), Sex =="F")$Simpson)
#Influence of age
wilcox.test((subset(sample_data(ps.2019), Age ==1))$Simpson, (subset(sample_data(ps.2019), Age ==3))$Simpson, paired=FALSE)
#Correlation with seagrass proportion
Simpson.Seagrass.lm <- lm(Simpson ~ Seagrass, data = data.alpha)
summary(Simpson.Seagrass.lm)

# =========================   Autumn 2021 analysis                                =================================================
plot_richness(ps.2021, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.2021, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))

#Bill length as an explanation for observed Simpson diversity
data.alpha <- cbind(sample_data(ps.2021))
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness)
#Difference between Bill <35mm and Bill > 35mm
wilcox.test((subset(sample_data(ps.2021), Bill_class =="Long"))$Simpson, (subset(sample_data(ps.2021), Bill_class =="Short"))$Simpson, paired=FALSE)
mean(subset(sample_data(ps.2021), Bill_class =="Short")$Simpson) - mean(subset(sample_data(ps.2021), Bill_class =="Long")$Simpson) #Again, difference is so small it is irrelevant
#Influence of sex
wilcox.test((subset(sample_data(ps.2021), Sex =="M"))$Simpson, (subset(sample_data(ps.2021), Sex =="F"))$Simpson, paired=FALSE)
#Influence of age
wilcox.test((subset(sample_data(ps.2021), Age ==1))$Simpson, (subset(sample_data(ps.2021), Age ==3))$Simpson, paired=FALSE)
#Correlation with seagrass proportion
Simpson.Seagrass.lm <- lm(Simpson ~ Seagrass, data = data.alpha)
summary(Simpson.Seagrass.lm)

# =========================   Spring 2022 analysis                                =================================================
plot_richness(ps.2022, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.2022, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))

#Bill length as an explanation for observed Simpson diversity
data.alpha <- cbind(sample_data(ps.2022))
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness)
#Difference between Bill <35mm and Bill > 35mm
wilcox.test((subset(sample_data(ps.2022), Bill_class =="Long"))$Simpson, (subset(sample_data(ps.2022), Bill_class =="Short"))$Simpson, paired=FALSE)
#Influence of sex
wilcox.test((subset(sample_data(ps.2022), Sex =="M"))$Simpson, (subset(sample_data(ps.2022), Sex =="F"))$Simpson, paired=FALSE)
#Influence of age --> does not work becasue there are nog juveniles in this subset
wilcox.test((subset(sample_data(ps.2022), Age ==1))$Simpson, (subset(sample_data(ps.2022), Age ==3))$Simpson, paired=FALSE)
#Correlation with seagrass proportion
Simpson.Seagrass.lm <- lm(Simpson ~ Seagrass, data = data.alpha)
summary(Simpson.Seagrass.lm)

#####             Mixed effects model, Batch as random effect                     #####

hist(sample_data(ps)$Simpson)
mixed.eff.data <- data.frame(sample_data(ps))
mixed.eff.data <- mixed.eff.data[mixed.eff.data$Sex != "",]
mixed.eff.model <- lmer(Simpson ~ Bill + Age + Sex + Expedition +(1|Extraction_batch), data = data.frame(sample_data(ps)))
summary(mixed.eff.model)
# Extraction batch accounts for 9.7% of the variance explained by random factors (batch and residual). Since this 

#####################################################################################
######                           Beta div analysis                              #####
#####################################################################################
#####                   Bray-curtis distance and ordination plots                 #####
#Determine Bray Curtis distance and ordination
## Distances between samples are calculated, and then ordered (ordinated) on who is closest to who
ps.rel <- microbiome::transform(subset_samples(ps, !(is.na(Seagrass))), "compositional")  #I need the extra check on having a seagrass proportion since one sample doesn't, and that breaks the PERMANOVA
ord.nmds.bray <- ordinate(ps.rel, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel, ord.nmds.bray, color="Sex", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel, ord.nmds.bray, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel, ord.nmds.bray, color="Expedition", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel, ord.nmds.bray, color="Bill_class") +stat_ellipse() + geom_point(size=3) + 
  theme(axis.title=element_text(size = 14), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), legend.title = element_text(size =13), legend.text = element_text(size=11)) +
  labs(color = "Bill length") + scale_color_discrete(name = "Bill length", labels = c("Long: > 35 mm", "Short: < 35 mm"))

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel)) # create a data frame with sample data 
ps.rel.bray <- phyloseq::distance(ps.rel, method = "bray") #Calculate distance between each pair of samples
adonis2(ps.rel.bray ~ Seagrass+Bill_class+Sex+Age+Expedition, data=subset(sampledf, (!is.na(sampledf$Seagrass) )), by = "margin", permutations=10000)  #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices

#Plotting the distribution of beta diversities
hist(x = ps.rel.bray, xlim = range(0,1), ylim = range(0, 1000), xlab ="Bray-Curtis dissimilarity score", ylab = "Frequency", main = "") 


#####                   Beta diversity per expedition                             #####
# ==================                  Autumn 2019                                 =========================
ps.rel.2019 <- microbiome::transform(ps.2019, "compositional")
ord.nmds.bray.2019 <- ordinate(ps.rel.2019, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.2019, ord.nmds.bray.2019, color="Sex", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2019, ord.nmds.bray.2019, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2019, ord.nmds.bray.2019, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel.2019)) # create a data frame with sample data 
ps.rel.bray.2019 <- phyloseq::distance(ps.rel.2019, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.2019 ~ Bill, data = sampledf) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices

# ==================                  Autumn 2021                                 ================================
ps.rel.2021 <- microbiome::transform(ps.2021, "compositional")
ord.nmds.bray.2021 <- ordinate(ps.rel.2021, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.2021, ord.nmds.bray.2021, color="Sex", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2021, ord.nmds.bray.2021, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2021, ord.nmds.bray.2021, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel.2021)) # create a data frame with sample data 
ps.rel.bray.2021 <- phyloseq::distance(ps.rel.2021, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.2021 ~ Bill, data = sampledf) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices

# ==================                  Spring 2022                                 ================================
ps.rel.2022 <- microbiome::transform(ps.2022, "compositional")
ord.nmds.bray.2022 <- ordinate(ps.rel.2022, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.2022, ord.nmds.bray.2022, color="Sex", title = "Bray NMDS") +stat_ellipse() 
plot_ordination(ps.rel.2022, ord.nmds.bray.2022, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2022, ord.nmds.bray.2022, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel.2022)) # create a data frame with sample data 
ps.rel.bray.2022 <- phyloseq::distance(ps.rel.2022, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.2022 ~ Bill, data = sampledf) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices    


#####################################################################################
#####                                   PCA                                     #####
#####################################################################################
#####                           setting up PCA                                    #####
#Making the otu table into proportions for comparison, as samples are not rarefied
## Hashtagged because it takes a long time to run. Generated file is loaded at the start of the script as asv_proportions
##In hindsight, I could have easily done this with "microbiome::transform" on "compositional", would have been way quicker. Oh well

#asv_proportions <- data.frame(otu_table(ps))
#for (j in 1:ncol(otu_table(ps))) {
  #sample_sum = sum(otu_table(ps)[,j])
  #for (i in 1:nrow(otu_table(ps))){
    #asv_proportions[i,j] = otu_table(ps)[i,j]/sample_sum
  #}
#}

#now that I have the abundances of each ASV per sample as a fraction (sum = 1), I need to limit the number of ASVs I use to 83 (The number of samples I have)
#adding row means as a column, then order dataframe rows on this column, remove column again, transpose so that now I take the 83 ASVs with highest mean abundance
ASV_means <- rowMeans(asv_proportions)
asv_proportions <- asv_proportions[order(-ASV_means),]
asv_proportions_top83 <- t(asv_proportions[1:83,]) 

#Running actual PCA
pca.model <- princomp(asv_proportions_top83)

#Adding PCA score of each sample on each PC to the sample_data (metadata) such that I can plot the principle components vs bird metadata
if (!(colnames(pca.model$scores)[1] %in% colnames(sample_data(ps)))){ #If statement makes sure I cannot cbind twice, which would duplicate columns and column names
  sample_data(ps) <- cbind(sample_data(ps), pca.model$scores)  
}

#####                     General stats of PCA model                              #####
summary(pca.model)

tail(sort(pca.model$loadings[,1])) #Most dominant ASVs in PC1 are ASV 1, 3, 5, 2, 4, 10
head(sort(pca.model$loadings[,1]))

tail(sort(pca.model$loadings[,2])) #Most dominant ASVs in PC2 are ASV 7, 8, 14, 24, 22, 19 (loaded more negatively than positively)
head(sort(pca.model$loadings[,2]))


#####           Defining ploting function that makes PC loadings plot             #####
PC_loadings_plot <- function(PC){
  asv.loadings = cbind(taxa[colnames(asv_proportions_top83),], "Loads"=pca.model$loadings[,PC], "ASV"=colnames(asv_proportions_top83)) #Making dataframe with columns of all phylogenetic levels, the loading on the PC given and the ASV number
  ggplot(data=asv.loadings, aes(x=reorder(ASV, -Loads), fill=Genus)) + geom_col(aes(y=Loads)) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5,size=12), text=element_text(size=14)) +
    xlab("ASV") + 
    ylab("Load") +
    #ggtitle(paste("Load per ASV in", PC)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(asv.loadings$Genus)), "Paired"))(length(unique(asv.loadings$Genus))))
}

#####                              Scree plot                                     #####
#making a scree plot of explained variance per Principle component
prop.of.variance = pca.model$sdev^2/sum(pca.model$sdev^2) #The variance explained is hard to get from the summary object, but can easily be calculated based on the st. dev. like this
prop.of.variance = data.frame(prop.of.variance)
prop.of.variance$PC <- rownames(prop.of.variance)
#Adding the cumulative proportions of variance explained to the dataframe 
prop.of.variance["Comp.1", "Cumulative.proportion"] <- prop.of.variance['Comp.1', "prop.of.variance"] #Adding first cumulative proportion by hand such that next loop is easy
for (row in 2:nrow(prop.of.variance)){ #Looping over each row, adding resp. proportion of row to cumulative proportion in the previous 
  prop.of.variance[row, "Cumulative.proportion"] = prop.of.variance[row-1, "Cumulative.proportion"] + prop.of.variance[row, "prop.of.variance"]
}
#The final scree plot
##Variance explained is almost 0 by PC50, so I only plot the first 50 PCs for a nicer figure
ggplot(data = data.frame(prop.of.variance[1:50,]), aes(x=reorder(PC, -prop.of.variance))) + 
  geom_col(aes(y=prop.of.variance)) +
  geom_line(aes(y=Cumulative.proportion, group=1)) + 
  scale_y_continuous(
    name = "Variance Explained", 
    sec.axis= sec_axis (trans= ~.*1 ,name = "Cumulative variance explained")) +
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5,size=12), text=element_text(size=16)) +
  xlab("Principal Component")

#####                     Other PCA related plots                                 ######
#Plotting PC1 and PC2, labelled with zostera fraction in diet
plot(Comp.2~Comp.1, data.frame(sample_data(ps)), cex=0)+ text(Comp.2~Comp.1, data.frame(sample_data(ps)), label=format(round(Seagrass, 2), nsmall = 2))

#Plotting fraction of Seagrass in diet vs PC1
##Jans version of the plot
model0=lm(Seagrass~Comp.1, data.frame(sample_data(ps)))
summary(model0)
plot(Seagrass~Comp.1, data.frame(sample_data(ps)), pch=21, bg="yellow", cex=3) +
  abline(model0) +
  text(Seagrass~Comp.1, data.frame(sample_data(ps)), label=Extraction_batch)

#GGplot versions of the plot
## Plots give warning of removed row, this is because NIOZ365.0125.0134 (Ringnr Z099517) has no input for the amount of seagrass
## Should probably fix this with Jan, as it is not in the diet file that Jan sent me.
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) + geom_point(aes(color=as.factor(Extraction_batch)), size=3) + geom_smooth(method = "loess", formula = y ~ x) 
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) + geom_point(aes(color=Age), size=3) + geom_smooth(method = "lm", formula = y ~ x) 
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) + geom_point(aes(color=Sex), size=3)+ geom_smooth(method = "lm", formula = y ~ x)
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) + geom_point(aes(color=as.factor(Expedition)), size=3)+ geom_smooth(method = "lm", formula = y ~ x)
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) + geom_point(aes(color=Bill_class), size=3)+ geom_smooth(method = "lm", formula = y ~ x)

#####                           PCs deep dive                                     #####
#=============                        PC1                                         #####
# Loading per ASv on PC1 plot
PC_loadings_plot("Comp.1")

#I think PC1 correlates best with Diet (either seagrass or Otherbivalves)
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)        #Seagrass
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Loripes)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)         #Loripes
ggplot(data =sample_data(ps), aes(x=Comp.1, y=OtherBivalves)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)   #OtherBivalves
ggplot(data =sample_data(ps), aes(x=Comp.1, y=Bill)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)            #Bill
ggplot(data =sample_data(ps), aes(x=Sex, y=Comp.1)) + geom_boxplot()                                                               #Sex
ggplot(data =sample_data(ps), aes(x=Age, y=Comp.1)) + geom_boxplot()                                                               #Age
ggplot(data =sample_data(ps), aes(x=as.factor(Expedition), y=Comp.1)) + geom_boxplot()

#linear model having PC1 score explain Seagrass fraction per sample 
model <- lm(Seagrass ~ Comp.1, data.frame(sample_data(ps)))
summary(model)
#=============                        PC2                                         #####
# plotting load per ASv on PC2 with function defined above
PC_loadings_plot("Comp.2")

#I think PC2 also correlates best with Diet (either seagrass or Otherbivalves)
ggplot(data =sample_data(ps), aes(x=Comp.2, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)        #Seagrass
ggplot(data =sample_data(ps), aes(x=Comp.2, y=Loripes)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)         #Loripes
ggplot(data =sample_data(ps), aes(x=Comp.2, y=OtherBivalves)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)   #OtherBivalves
ggplot(data =sample_data(ps), aes(x=Comp.2, y=Bill)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)            #Bill
ggplot(data =sample_data(ps), aes(x=Sex, y=Comp.2)) + geom_boxplot()                                                               #Sex
ggplot(data =sample_data(ps), aes(x=Age, y=Comp.2)) + geom_boxplot()                                                               #Age
ggplot(data =sample_data(ps), aes(x=as.factor(Expedition), y=Comp.2)) + geom_boxplot()                                             #Expedition

#linear model having PC2 score explain Seagrass fraction per sample 
model <- lm(Seagrass ~ Comp.2, data.frame(sample_data(ps)))
summary(model)

#=============                        PC3                                         #####
# plotting load per ASv on PC2 with function defined above
PC_loadings_plot("Comp.3")

#I think PC3 also still correlates best with Diet (either seagrass or Otherbivalves)
ggplot(data =sample_data(ps), aes(x=Comp.3, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)        #Seagrass
ggplot(data =sample_data(ps), aes(x=Comp.3, y=Loripes)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)         #Loripes
ggplot(data =sample_data(ps), aes(x=Comp.3, y=OtherBivalves)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)   #OtherBivalves
ggplot(data =sample_data(ps), aes(x=Comp.3, y=Bill)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)            #Bill
ggplot(data =sample_data(ps), aes(x=Sex, y=Comp.3)) + geom_boxplot()                                                               #Sex
ggplot(data =sample_data(ps), aes(x=Age, y=Comp.3)) + geom_boxplot()                                                               #Age
ggplot(data =sample_data(ps), aes(x=as.factor(Expedition), y=Comp.3)) + geom_boxplot()                                             #Expedition

#linear model having PC3 score explain Seagrass fraction per sample 
model <- lm(Seagrass ~ Comp.3, data.frame(sample_data(ps)))
summary(model)

#=============                        PC4                                         #####
# plotting load per ASv on PC2 with function defined above
PC_loadings_plot("Comp.4")

#I think PC4 moight correlate best with Sex though it is very vague and unclear
ggplot(data =sample_data(ps), aes(x=Comp.4, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)        #Seagrass
ggplot(data =sample_data(ps), aes(x=Comp.4, y=Loripes)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)         #Loripes
ggplot(data =sample_data(ps), aes(x=Comp.4, y=OtherBivalves)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)   #OtherBivalves
ggplot(data =sample_data(ps), aes(x=Comp.4, y=Bill)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)            #Bill
ggplot(data =sample_data(ps), aes(x=Sex, y=Comp.4)) + geom_boxplot()                                                               #Sex
ggplot(data =sample_data(ps), aes(x=Age, y=Comp.4)) + geom_boxplot()                                                               #Age
ggplot(data =sample_data(ps), aes(x=as.factor(Expedition), y=Comp.4)) + geom_boxplot()                                             #Expediti

ggplot(data =sample_data(ps), aes(x=Comp.4, y=Seagrass, color=Sex)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)        


##### ASV specific search
#First, I want to determine which PC associates best with amount of Zostera in diet
ggplot(data=sample_data(ps), aes(x=Comp.1, y=Comp.2)) + geom_point(aes(color=Seagrass), size=3) +scale_color_gradient(low="red", high="green")

ggplot(data =sample_data(ps), aes(x=Comp.2, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)
ggplot(data =sample_data(ps), aes(x=Comp.3, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)

ggplot(data =sample_data(ps)[sample_data(ps)$Comp.1 < 0.2], aes(x=Comp.1, y=Seagrass)) + geom_point(size=3) + geom_smooth(method = "lm", formula = y ~ x)
model1=lm(Seagrass~Comp.1, data.frame(sample_data(ps)[sample_data(ps)$Comp.1 < 0.2]))
summary(model1)

ggplot(data =sample_data(ps), aes(x=Comp.1, y=Seagrass)) +geom_point()

#####          Searching for a PC significant with Seagrass fraction in diet      #####
#Setting up a dataframe to store p-values of linear models based for each PCs correlation with the amount of seagrass in the diet
p_scores_seagrass <- data.frame(matrix(ncol=2, nrow=length(colnames(sample_data(ps))[26:ncol(sample_data(ps))])))
colnames(p_scores_seagrass) <- c("Principle.component", "p_value")
rownames(p_scores_seagrass) <- colnames(sample_data(ps))[26:ncol(sample_data(ps))]

#Populating P-value dataframe
for (PC in colnames(sample_data(ps))[26:ncol(sample_data(ps))]) {
  model = lm(paste("Seagrass ~",PC[[1]]), data.frame(sample_data(ps)))   #For each PC, a linear model is built on its correlation with the fraction of seagrass in the diet
  p_scores_seagrass[PC,"Principle.component"] <- PC                               #Adding the name of the PC as a column as wel as rowname for later plottnig
  p_scores_seagrass[PC,"p_value"] <- summary(model)$coefficients[8]               # Adding the P-value of the lm of the PC with seagrass to the p-value dataframe
}

p_scores_seagrass$p_value_corrected <- p.adjust(p_scores_seagrass$p_value, method="BH")   #Performing  Benjamin and hochberg correction on calculatd p-values
min(p_scores_seagrass$p_value_corrected)
##Sadly, nothing is even close to significant after correction

#####          Searching for a PC significant with expedition     !UNFINISHED!    #####
#Setting up a dataframe to store p-values of linear models based for each PCs correlation with the amount of seagrass in the diet
p_scores_expedition <- data.frame(matrix(ncol=2, nrow=length(colnames(sample_data(ps))[26:ncol(sample_data(ps))])))
colnames(p_scores_expedition) <- c("Principle.component", "p_value")
rownames(p_scores_expedition) <- colnames(sample_data(ps))[26:ncol(sample_data(ps))]


#Expedition werkt natuurlijk niet met een LM, want dat is categoriaal
expedition_dataframe <- data.frame(sample_data(ps))
expedition_dataframe$Expedition <- as.factor(expedition_dataframe$Expedition)
for (PC in colnames(sample_data(ps))[26:ncol(sample_data(ps))]) {
  model = kruskal.test(paste("Expedition ~",PC[[1]]), expedition_dataframe)   #For each PC, a linear model is built on its correlation with the fraction of seagrass in the diet
  p_scores_expedition[PC,"Principle.component"] <- PC                               #Adding the name of the PC as a column as wel as rowname for later plottnig
  p_scores_expedition[PC,"p_value"] <- model[3][[1]]               # Adding the P-value of the lm of the PC with seagrass to the p-value dataframe
}
p_scores_expedition$p_value_corrected <- p.adjust(p_scores_expedition$p_value, method="BH")   #Performing  Benjamin and hochberg correction on calculatd p-values

##Sadly, nothing is even close to significant after correction
min(p_scores_expedition$p_value_corrected)
test <- kruskal.test(as.factor(Expedition)~Comp.1, data.frame(sample_data(ps)))
test[3][[1]]
             
#####################################################################################
#####                             Core microbiome                               #####
#####################################################################################
#####                       Defining the core microbiome                          #####

# Core microbiome already imported at the start of the script. Codeblock below was used to make this subset

#I use the definition for core microbiome as "found in >50% of sampled cases", so >0 relative abundance in >50% of birds
## I can calculate this using the relative abundance table "asv_proportions". Per row, I look check how many values in that row are >0, and divide the number of times this is true by the length of the row
## I also already have a data frame that has both the relative abundances and the taxonomic values of all asv s (ASV.abundance.taxa), So i can subset rows from there to a new data frame
## !! This code only works if both asv_proportions and ASV.abundance.taxa have their rows in the same order !!  
#core_microbiome <- data.frame()
#for(row in 1:nrow(asv_proportions)){
 # if ( (sum(asv_proportions[row,] >0) / length(asv_proportions[row,])) > 0.5 ) {
  #  core_microbiome <- rbind(core_microbiome, ASV.abundance.taxa[row,])
  #}
#}

#Percentage of ASVs that make it into the core microbiome
nrow(core_microbiome) / nrow(asv_proportions) *100
#Number of ASVs that make it into the core microbiome
nrow(core_microbiome)
sum(colnames(asv_proportions_top83) %in% rownames(core_microbiome)) #60 of the 83 mean most abundant ASVs are in the core microbiome 

#####                         Plotting core microbiome                            #####
#Number of occurrences of each class in total dataset
ggplot(data=core_microbiome, aes(y=fct_rev(fct_infreq(Phylum)), fill=Phylum)) +geom_bar() + 
  theme(legend.position="none", text=element_text(size=14), axis.text = element_text(size=12)) + xlab("Number of ASVs") + ylab("Phylum") + ggtitle("Number of ASVs per phylum in dataset")
ggplot(data=core_microbiome, aes(y=fct_rev(fct_infreq(Class)), fill=Phylum)) +geom_bar()

ggplot(data=core_microbiome, aes(y=fct_rev(fct_infreq(Genus)), fill=Phylum)) +geom_bar()+ 
  theme(legend.position="none", text=element_text(size=14), axis.text = element_text(size=12)) +
  xlab("Number of ASVs in core microbiome") + ylab("Genus")

ggplot(data=core_microbiome, aes(y=fct_rev(fct_infreq(Genus)), fill=Class)) +geom_bar()

sort(table(core_microbiome$Phylum), decreasing = TRUE)
sort(table(core_microbiome$Class), decreasing = TRUE)
sort(table(core_microbiome$Genus), decreasing = TRUE)
#####################################################################################
#####                       Miscelanious and unfinished                         #####
#####################################################################################
#####                       Bill length vs diet fractions                         #####
#Plotting bill length as a function of seagrass diet, loripes and otherBivalve fraction
ggplot(sample_data(ps), aes(Seagrass, Bill)) + geom_point(aes(color=Expedition), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps), aes(Loripes, Bill)) + geom_point(aes(color=Expedition), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps), aes(OtherBivalves, Bill)) + geom_point(aes(color=Expedition), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)

ggplot(sample_data(ps), aes(Loripes, Bill)) + geom_point(aes(color=Age), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps), aes(Loripes, Bill)) + geom_point(aes(color=Sex), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)


#Plotting bill length as a function of seagrass fraction per expedition
ggplot(sample_data(ps.2019), aes(Seagrass, Bill)) + geom_point(aes(color=Sex), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps.2021), aes(Seagrass, Bill)) + geom_point(aes(color=Sex), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps.2022), aes(Seagrass, Bill)) + geom_point(aes(color=Sex), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)

#Plotting alpha div as function of seagrassdiet, loripes and otherBivalve fraction
ggplot(sample_data(ps), aes(Seagrass, Simpson)) + geom_point(aes(color=Expedition), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps), aes(Loripes, Simpson)) + geom_point(aes(color=Expedition), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)
ggplot(sample_data(ps), aes(OtherBivalves, Simpson)) + geom_point(aes(color=Expedition), size=3) + theme(text=element_text(size=17))+ geom_smooth(method = "loess", formula = y~x)

##### Stacked bar chart of diet composition per expedition per sample             #####

for (bird in 1:nrow(diet)){
  ringnr <- diet[bird, "Ring"][[1]]
  if (ringnr %in% sample_data(ps)$Description){
    expedition <- sample_data(ps)[sample_data(ps)$Description == ringnr,][["Expedition"]]
    diet[diet$Ring == ringnr, "Expedition"] <- expedition
  }
}


diet2 <- melt(diet, id.vars = c("Ring", "Expedition"), measure.vars=c("Loripes", "Seagrass(+Abra)", "OtherBivalves"))
diet2 <- na.omit(diet2)

ggplot(diet2, aes(fill=variable, y=value, x=Ring)) + geom_bar(position = "fill", stat = "identity") + facet_wrap(~ Expedition, scales ="free_x") + theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1))

#####################################################################################
#####################################################################################
#####################################################################################
########                      Rarefied alpha div analysis                    ######## 
#####################################################################################
#####                     Rarefying and splitting dataset                         #####
#Rarefying dataset such that all samples have the same number of reads as the lowest one
ps.rarefied = rarefy_even_depth(ps)
summarize_phyloseq(ps.rarefied)

#Adding Alpha diversity measures to sample_data
erich = estimate_richness(ps.rarefied, measures = c("Observed", "Shannon", "Simpson"))
sample_data(ps.rarefied)$Observed<-erich$Observed
sample_data(ps.rarefied)$Shannon<-erich$Shannon
sample_data(ps.rarefied)$Simpson<-erich$Simpson

###### Plotting alpha diversity for different bird qualities all data rarefied    #####
plot_richness(ps.rarefied, x="Age", measures=c("Observed", "Shannon", "Simpson"), color="Age") + geom_boxplot(alpha=0.6)+theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=1,vjust=1,size=12))
plot_richness(subset_samples(ps.rarefied, Sex != ""), x="Sex", measures=c("Observed", "Shannon", "Simpson"), color="Sex") + geom_boxplot(alpha=0.6) + theme(axis.title=element_text(size = 14), legend.position="none", axis.text.x=element_text(angle=0, hjust=1,vjust=1,size=14), axis.text.y=element_text(size=14))
plot_richness(ps.rarefied, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Expedition") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.rarefied, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.rarefied, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.rarefied, x="Age", measures=c("Observed", "Shannon", "Simpson"), color="Age") + geom_boxplot(alpha=0.6)+theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=1,vjust=1,size=12))
plot_richness(ps.rarefied, x="Bill_class", measures=c("Observed", "Shannon", "Simpson"), color="Bill_class") + geom_boxplot(alpha=0.6)+theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=0.5, vjust=1,size=12))


#####            Alpha Diversity comparisons for full dataset rarefied            #####
data.alpha <- cbind(sample_data(ps.rarefied))
#Making a linear model of the observed number of ASVs when bill length increases (richness vs Bill)
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness) #Based on Simpson measure
#Influence of sex
wilcox.test((subset(sample_data(data.alpha), Sex =="M"))$Simpson, (subset(sample_data(data.alpha), Sex =="F"))$Simpson, paired=FALSE)
#Influence of age
wilcox.test((subset(sample_data(data.alpha), Age ==1))$Simpson, (subset(sample_data(data.alpha), Age ==3))$Simpson, paired=FALSE)
#Influence of expedition
levels(as.factor(sample_data(data.alpha)$Expedition))
kruskal.test(Simpson ~ as.factor(Expedition), data.alpha)


######                    Splitting data based on expedition rarefied             #####

# Making separate phyloseq objects for each expidition for ease of use downstream
ps.rarefied.2019 <- subset_samples(ps.rarefied, Expedition =="Mauritania_autumn_2019")
ps.rarefied.2021 <- subset_samples(ps.rarefied, Expedition =="Mauritania_autumn_2021")
ps.rarefied.2022 <- subset_samples(ps.rarefied, Expedition =="Mauritania_spring_2023") #For some reason dates from 2022 have been marked as 2023 in the expidition name. Does not impact analysis any further
# KEEP IN MIND: in the following 3 analyses I overwrite "data.alpha" multiple times, always make sure to update to the correct versions before use

# =========================   Autumn 2019 analysis rarefied                       =================================================
plot_richness(ps.rarefied.2019, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.rarefied.2019, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))

#Bill length as an explanation for observed Simpson diversity
data.alpha <- cbind(sample_data(ps.rarefied.2019))
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness)
#Influence of sex --> Is significant, male and female means also caluclated
wilcox.test((subset(sample_data(ps.rarefied.2019), Sex =="M"))$Simpson, (subset(sample_data(ps.rarefied.2019), Sex =="F"))$Simpson, paired=FALSE)
mean(subset(sample_data(ps.rarefied.2019), Sex =="M")$Simpson)
mean(subset(sample_data(ps.rarefied.2019), Sex =="F")$Simpson)
#Influence of age
wilcox.test((subset(sample_data(ps.rarefied.2019), Age ==1))$Simpson, (subset(sample_data(ps.rarefied.2019), Age ==3))$Simpson, paired=FALSE)


# =========================   Autumn 2021 analysis rarefied                       =================================================
plot_richness(ps.rarefied.2021, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.rarefied.2021, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))

#Bill length as an explanation for observed Simpson diversity
data.alpha <- cbind(sample_data(ps.rarefied.2021))
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness)
#Influence of sex --> Is significant, male and female means also caluclated
wilcox.test((subset(sample_data(ps.rarefied.2021), Sex =="M"))$Simpson, (subset(sample_data(ps.rarefied.2021), Sex =="F"))$Simpson, paired=FALSE)
#Influence of age
wilcox.test((subset(sample_data(ps.rarefied.2021), Age ==1))$Simpson, (subset(sample_data(ps.rarefied.2021), Age ==3))$Simpson, paired=FALSE)


# =========================   Spring 2022 analysis rarefied                       =================================================
plot_richness(ps.rarefied.2022, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Sex") +geom_point(size=2) + theme(text=element_text(size=17))
plot_richness(ps.rarefied.2022, x="Bill", measures=c("Observed", "Shannon", "Simpson"), color="Age") +geom_point(size=2) + theme(text=element_text(size=17))


#Bill length as an explanation for observed Simpson diversity
data.alpha <- cbind(sample_data(ps.rarefied.2022))
model.richness <- (lm(Simpson ~ Bill, data.alpha))
summary(model.richness)
#Influence of sex --> Is significant, male and female means also caluclated
wilcox.test((subset(sample_data(ps.rarefied.2022), Sex =="M"))$Simpson, (subset(sample_data(ps.rarefied.2022), Sex =="F"))$Simpson, paired=FALSE)
#Influence of age
wilcox.test((subset(sample_data(ps.rarefied.2022), Age ==1))$Simpson, (subset(sample_data(ps.rarefied.2022), Age ==3))$Simpson, paired=FALSE)






#####################################################################################
######                     Rarefied  Beta div analysis                          #####
#####################################################################################
#####               Rarefied Bray-curtis distance and ordination plots            #####
#Determine Bray Curtis distance and ordination
## Distances between samples are calculated, and then ordered (ordinated) on who is closest to who
ps.rel.rarefied <- microbiome::transform(subset_samples(ps.rarefied, !(is.na(Seagrass))), "compositional")
ord.nmds.bray <- ordinate(ps.rel, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.rarefied, ord.nmds.bray, color="Sex", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.rarefied, ord.nmds.bray, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.rarefied, ord.nmds.bray, color="Expedition", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.rarefied, ord.nmds.bray, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf.rarefied <- data.frame(sample_data(ps.rel.rarefied)) # create a data frame with sample data 
ps.rel.bray.rarefied <- phyloseq::distance(ps.rel.rarefied, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.rarefied ~ Bill, data = sampledf.rarefied) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices

#Plotting the distribution of beta diversities
hist(x = ps.rel.bray.rarefied, xlim = range(0,1), xlab ="Bray-Curtis dissimilarity score", ylab = "Frequency", main = "Frequency distirbution of Bray-Curtis dissimilarities") 
bray_distances.rarefied <- matrix(ps.rel.bray.rarefied)
ggplot(data=data.frame(ps.rel.bray.rarefied), aes(x=ps.rel.bray.rarefied)) + geom_bar()

######                    BETA DIVERSITY PER EXPEDITION                           #####
# ==================            Autumn 2019 rarefied                              =========================
ps.rel.2019.rarefied <- microbiome::transform(ps.rarefied.2019, "compositional")
ord.nmds.bray.2019.rarefied <- ordinate(ps.rel.2019.rarefied, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.2019.rarefied, ord.nmds.bray.2019.rarefied, color="Sex", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2019.rarefied, ord.nmds.bray.2019.rarefied, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2019.rarefied, ord.nmds.bray.2019.rarefied, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel.2019.rarefied)) # create a data frame with sample data 
ps.rel.bray.2019.rarefied <- phyloseq::distance(ps.rel.2019.rarefied, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.2019.rarefied ~ Bill, data = sampledf) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices

# ==================            Autumn 2021 rarefied                              ================================
ps.rel.2021.rarefied <- microbiome::transform(ps.rarefied.2021, "compositional")
ord.nmds.bray.2021.rarefied <- ordinate(ps.rel.2021.rarefied, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.2021.rarefied, ord.nmds.bray.2021.rarefied, color="Sex", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2021.rarefied, ord.nmds.bray.2021.rarefied, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2021.rarefied, ord.nmds.bray.2021.rarefied, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel.2021.rarefied)) # create a data frame with sample data 
ps.rel.bray.2021.rarefied <- phyloseq::distance(ps.rel.2021.rarefied, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.2021.rarefied ~ Bill, data = sampledf) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices


# ==================            Spring 2022 rarefied                              ================================
ps.rel.2022.rarefied <- microbiome::transform(ps.rarefied.2022, "compositional")
ord.nmds.bray.2022.rarefied <- ordinate(ps.rel.2022, method="NMDS", distance="bray")

#Plotting ordination results
plot_ordination(ps.rel.2022.rarefied, ord.nmds.bray.2022.rarefied, color="Sex", title = "Bray NMDS") +stat_ellipse() 
plot_ordination(ps.rel.2022.rarefied, ord.nmds.bray.2022.rarefied, color="Age", title = "Bray NMDS") +stat_ellipse()
plot_ordination(ps.rel.2022.rarefied, ord.nmds.bray.2022.rarefied, color="Bill_class", title = "Bray NMDS") +stat_ellipse()

# Calculating effect size of Bill on community similarity
sampledf <- data.frame(sample_data(ps.rel.2022.rarefied)) # create a data frame with sample data 
ps.rel.bray.2022.rarefied <- phyloseq::distance(ps.rel.2022.rarefied, method = "bray") # generate a distance matrix using Bray distances
adonis2(ps.rel.bray.2022.rarefied ~ Bill, data = sampledf) #ADONIS is another name for PERMANOVA --> analysis of variance using distance matrices    