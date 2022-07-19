library(ggplot2)
library(phyloseq)
library(vegan)
library("ALDEx2")
library(dplyr)
library(ggpubr)
library(otuSummary)


##########TEST
taxa_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/taxonomy_complete.csv", sep=';', row.names=1)
abundace_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/abundances.csv", row.names=1)
metadata = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/metadata_tidy.csv", sep=';', row.names=1)

#some extra columns were inserted at the end of the file and need to be deleted

taxa_table <- taxa_table[1:(length(taxa_table)-9)]

#Remove mitochondria and chloroplasts

#Phyloseq object
ps <- phyloseq(otu_table(abundace_table, taxa_are_rows = FALSE),
               tax_table(as.matrix(taxa_table)), sample_data(metadata))

#Unforunatly we lost 4 samples during the pre-processing so we need to also remove them from the metadata
lost_samples <- c('N97', 'N116', 'N122', 'N53')
metadata <- subset(metadata, !(row.names(metadata) %in% lost_samples))

#need to remove “junk” sequences from the taxonomy table such as Chloroplasts and Mitochondria

is.chloroplast <- taxa_table[,"Order"] %in% "Chloroplast"
seqtable.nochloro <- abundace_table[,!is.chloroplast]
dim(seqtable.nochloro) 
taxa_table.nochloro <- taxa_table[!is.chloroplast,]
dim(taxa_table.nochloro)
#31 sequence containing chloroplasts removed

is.mitochondria <- taxa_table.nochloro[,"Family"] %in% "Mitochondria"
seqtable.nomito <- seqtable.nochloro[,!is.mitochondria]
taxonomy.nomito <- taxa_table.nochloro[!is.mitochondria,]
dim(seqtable.nomito) #121 sequences containing mitochodnria removed
dim(seqtable.nochloro)

seqtable_nochim <- seqtable.nomito #abundance table
taxa_table <- taxonomy.nomito #taxonomy table
metadata <- metadata #metadata table

write.csv(seqtable_nochim, file="/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/asv-no-chloro-mito.csv")
write.csv(taxa_table, file="/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/taxa-no-chloro-mito.csv")
write.csv(metadata, file="/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/metadata-final.csv")

#Rarefraction

ps <- phyloseq(otu_table(seqtable_nochim, taxa_are_rows = FALSE),
               tax_table(as.matrix(taxa_table)), sample_data(metadata))

#Creating a rarefraction curve with a vegan package function "rarecurve". 
#This is ONLY to visualise the distribution of ASVs amongst the samples - NEVER RARIFY!
rarecurve(seqtable_nochim, step=100, lwd=2, ylab="ASVs", label=TRUE)
abline(v=(min(rowSums(abundace_table))))
max(sample_sums(ps)) #43619
min(sample_sums(ps)) #5553

rarefied <-subset_samples(ps,sample_sums(ps)>9999) #834 OTUs lost 
rarefied_data <-rarefy_even_depth(rarefied) #rarify without replacement 
sample_sums(rarefied_data) # I can see that all samples now have been rarified to 10258  

#save rarefied dataset
write.csv(otu_table(rarefied_data),"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/rarefied-asv.csv", row.names = TRUE)
write.csv(tax_table(rarefied_data),"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/rarefied-taxa.csv", row.names = TRUE)
write.csv(sample_data(rarefied_data),"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/rarefied-meta.csv", row.names = TRUE)

meta <- data.frame(sample_data(rarefied_data)) #Accessing my sample information from the ps object containing rarefied data
otus <- data.frame(otu_table(rarefied_data, taxa_are_rows = TRUE)) #Acessing an otu table from the ps object containing rarefied data
taxa <- data.frame(tax_table(rarefied_data)) #Acessing an otu table from the ps object containing rarefied data

#Read it
otus <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/rarefied-asv-numbered-correct.csv", sep=";")
otus <- otus[,-1]
rownames(otus) <- otus$X.1
otus <- otus[,-1]

meta <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/rarefied-meta-numbered-correct.csv", sep=";")
meta <- meta[,-1]
rownames(meta) <- meta$X.1
meta <- meta[,-1]

taxa <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/rarefied-taxa.csv", row.names = 1)

#De-replicate at species level

#ASVs with species
species_vec <- c()
for(i in 1:length(taxa$Species)){
  if(!is.na(taxa$Species[i])){
    species_vec <- append(species_vec, rownames(taxa)[i])
  }
}

otus_taxa <- otus[,colnames(otus) %in% species_vec]
taxa_taxa <- taxa[rownames(taxa) %in% species_vec,]

ps <- phyloseq(otu_table(otus_taxa, taxa_are_rows = FALSE),
               tax_table(as.matrix(taxa_taxa)))

ps <- tax_glom(ps, taxrank=rank_names(ps)[7], NArm=FALSE)

aggregated_otus <- as.data.frame(otu_table(ps))
aggregated_taxa <- as.data.frame(tax_table(ps))

#replace otus with species with aggregated ones
otus <- otus[,!colnames(otus) %in% colnames(otus_taxa)]
otus <- cbind(otus, aggregated_otus)

#replace taxa with species with aggregated ones
taxa <- taxa[!rownames(taxa) %in% rownames(taxa_taxa),]
taxa <- rbind(taxa, aggregated_taxa)

#Replacing sampling source with pollution gradient:
# Drinking Water Resservoir  	1
# Nature resservoir	2
# Regular pond		3
# Polluted pond		4
# WWTP		5

for(i in 1:length(rownames(meta))){
  if(meta$Type[i] == "Drinking water Reservoir"){
    meta$Type[i] <- 1
  }else if(meta$Type[i] == "Nature reserve"){
    meta$Type[i] <- 2
  }else if(meta$Type[i] == "Regular pond"){
    meta$Type[i] <- 3
  }else if(meta$Type[i] == "Polluted Pond"){
    meta$Type[i] <- 4
  }else if(meta$Type[i] == "WWTP "){
    meta$Type[i] <- 5
  }
}

#1
length(which(meta$Type == "1")) #28
#2
length(which(meta$Type == "2")) #54
#3
length(which(meta$Type == "3")) #27
#4
length(which(meta$Type == "4")) #35
#5
length(which(meta$Type == "5")) #9


#Save the data
write.csv(otus,"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/final-asvs-for-ordination.csv", row.names = TRUE)
write.csv(taxa,"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/final-taxa-for-ordination.csv", row.names = TRUE)
write.csv(meta,"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/final-meta-for-ordination.csv", row.names = TRUE)

#Chao1 Richness
sum(otus[1,]) #quick check if it worked, the sum gives back a number of 10258 so its good to go
rarecurve((otus), step=100, lwd=2, ylab="ASVs", label=TRUE) #visualising the rarefraction curve after rarefying 
otus <- t(otus)

# Getting singletons and doubletons for Chao richness
meta$Richness <- colSums(otus != 0) #Calculating species richness by summing all non-zero values of my data
singletons <- colSums(otus == 1) #Summing all values that appeared only once 
doubletons <- colSums(otus == 2) # Summing all values that apperad twice
rares = (singletons/(2*doubletons))
rares[doubletons==0] <- 0 
meta$Chao1 = meta$Richness + rares

#Turning counts into proportions
sums <- apply(otus, 2, sum) # get the column sums (sample sums)
norm_counts <- otus # generate new data frame for normalized data
for (i in 1:ncol(otus)){ # divide per sample by sample total
  norm_counts [,i] <- otus[,i]/sums[i]
}

#Shannon Diversity
shannon_div <- function(vector){
  vector <- vector*log(vector)
  # the vector has NA values if the species proportions are 0
  vectorsum <- sum(na.omit(vector)) 
  return((-1)*vectorsum)
}

# Calculate diversities
meta$Shannon = apply(norm_counts, 2, shannon_div)
meta$Pielou = meta$Shannon / log(meta$Richness)
meta$Simpson = colSums(norm_counts * norm_counts)
meta$InvSimpson = 1 / meta$Simpson

colours <- c("#c96d44",
             "#777acd")

long_data <- reshape2::melt(meta, id.vars="experiment")
long_data = long_data[long_data$variable %in% c('Chao1', 'Shannon', 'Pielou', 'InvSimpson', 'Richness'),] #selecting indices
long_data$value <- as.numeric(long_data$value)

library(ggpubr)
#Plotting indices
ggboxplot(long_data, x="experiment", y = "value", fill="experiment")+
  scale_fill_manual(values = colours)+
  stat_compare_means()+
  theme(legend.position="none",strip.text = element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(name="", labels = scales::comma)+
  facet_wrap(. ~ variable, scales = "free", ncol = 5)+
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Daphnia

meta_daphnia <- meta[meta$experiment == 'Daphnia',] #Filter out only Daphnia data

#  h0= There is no difference between different diversity indices for sampling locations
print(kruskal.test(Chao1~Type, data=meta_daphnia)) # p-value = 0.515 H0 accepted
print(kruskal.test(Shannon~Type, data=meta_daphnia)) # p-value = 0.8162 H0 accepted
print(kruskal.test(Pielou~Type, data=meta_daphnia)) #p-value = 0.1664 H0 accepted
print(kruskal.test(InvSimpson~Type, data=meta_daphnia)) #p-value = 0.4883 H0 accepted
print(kruskal.test(Richness~Type, data=meta_daphnia)) #p-value = 0.4938 H0 accepted
# The diversity indices between the samples are not significantly different, H0 accepted

level_order <- c('1', '2', '3', '4') #this vector might be useful for other plots/analyses
colours <- c("#c96d44",
             "#777acd",
             "#7aa456",
             "#c65999")

#cleaning up metadata
meta_daphnia <- meta_daphnia[,-c(5:8)]
meta_daphnia <- meta_daphnia[,-3]

long_data <- reshape2::melt(meta_daphnia, id.vars="Type")
long_data = long_data[long_data$variable %in% c('Chao1', 'Shannon', 'Pielou', 'InvSimpson', 'Richness'),] #selecting indices
long_data$value <- as.numeric(long_data$value)

#Plotting indices
ggboxplot(long_data, x="Type", y = "value", fill="Type")+
  scale_fill_manual(values = colours)+
  stat_compare_means()+
  theme(legend.position="none",strip.text = element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(name="", labels = scales::comma)+
  facet_wrap(. ~ variable, scales = "free", ncol = 5)+
  theme_pubclean()


##Bacterioplankton
meta_bac <- meta[meta$experiment == 'Bacterioplankton',] #Filter out only Daphnia data

#  h0= There is no difference between different diversity indices for sampling locations
print(kruskal.test(Chao1~Type, data=meta_bac)) # p-value = 0.01226 H0 rejected
print(kruskal.test(Shannon~Type, data=meta_bac)) # p-value = 0.000831 H0 accepted
print(kruskal.test(Pielou~Type, data=meta_bac)) #p-value = 0.07394 H0 accepted
print(kruskal.test(InvSimpson~Type, data=meta_bac)) #p-value = 0.01051 H0 accepted
print(kruskal.test(Richness~Type, data=meta_bac)) #p-value = 0.01148 H0 accepted
# The diversity indices between the samples are not significantly different, H0 accepted

level_order <- c('1', '2', '3', '4', '5') #this vector might be useful for other plots/analyses
colours <- c("#71a557",
             "#ab62c0",
             "#c47d3c",
             "#648ace",
             "#ca566f")

#cleaning up metadata
meta_bac <- meta_bac[,-c(6:8)]
meta_bac <- meta_bac[,-3]

long_data <- reshape2::melt(meta_bac, id.vars="Type")
long_data = long_data[long_data$variable %in% c('Chao1', 'Shannon', 'Pielou', 'InvSimpson', 'Richness'),] #selecting indices
long_data$value <- as.numeric(long_data$value)

#Plotting indices for Type
ggboxplot(long_data, x="Type", y = "value", fill="Type")+
  scale_fill_manual(values = colours)+
  stat_compare_means()+
  theme(legend.position="none",strip.text = element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(name="", labels = scales::comma)+
  facet_wrap(. ~ variable, scales = "free", ncol = 5)+
  theme_pubclean()

long_data <- reshape2::melt(meta_bac, id.vars="Layer")
long_data = long_data[long_data$variable %in% c('Chao1', 'Shannon', 'Pielou', 'InvSimpson', 'Richness'),] #selecting indices
long_data$value <- as.numeric(long_data$value)

colours <- c("#71a557",
             "#ab62c0",
             "#c47d3c",
             "#648ace",
             "#ca566f")

#Empty --> Mid
for(i in 1:length(long_data$Layer)){
  if(long_data$Layer[i] ==""){
    long_data$Layer[i] <- "Mid"
  }
}


#Plotting indices for Layer
ggboxplot(long_data, x="Layer", y = "value", fill="Layer")+
  scale_fill_manual(values = colours)+
  stat_compare_means()+
  theme(legend.position="none",strip.text = element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(name="", labels = scales::comma)+
  facet_wrap(. ~ variable, scales = "free", ncol = 5)+
  theme_pubclean()


#modifications to the tables
meta$Joint <- paste(meta$experiment, meta$Type) #making a join column for experiment and Type
taxa$combined <- paste(taxa$Family, taxa$Genus, taxa$Species, sep=" ")

#####################ps#############################################################

#incorporating all the elements into ps elements:
ps <- phyloseq(otu_table(otus, taxa_are_rows = FALSE),
               tax_table(as.matrix(taxa)), sample_data(meta))

###############Filtering out the zeros###############################################
filt_ps <- prevFilter(ps, 0.1) #Filter taxa that are not in at least 10% of samples
otu <- as.data.frame(otu_table(filt_ps))

####################Ordination########################################################

#Bray distance matrix
bray <- vegdist(otus, method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, distance='bray', na.action='na.omit')

#Turn Eigenvalues into %
explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

#PERMANOVA test with adonis

####################Experiment######################################

permanova = adonis2(bray ~ meta$experiment)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ meta$experiment)
# Df SumOfSqs      R2      F Pr(>F)    
# meta$experiment   1    7.348 0.11917 20.429  0.001 ***
#   Residual        151   54.310 0.88083                  
# Total           152   61.657 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#0.11917 = 12% of the variation in distances is explained by the grouping being tested = experiment

##########Dispersion
dispr <- vegan::betadisper(bray, meta$experiment) #test of the equality of variances.
permutest <- anova(dispr) #pvalue = 0.03457 ** samples are not homogenous

#####PCoA plot with experiment == Daphnia vs Bacterioplankton

#Clusters for the elipses
H_CLustering=hclust(bray)

#Colour for the points
color_origin <- c("Daphnia"="#96568b","Bacterioplankton"="#99a666")

#Setting treatment variable
treatment <- meta$experiment

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 1, bg = color_origin[treatment], pch= 21, lwd = 0.5, col= "black")
legend("bottomright", inset = c(-0.4, 0.01), legend=c("Daphnia","Bacterioplankton"), col=c("#96568b","#99a666"), pch=19, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

#Pairwise -- Joint

ps <- phyloseq(otu_table(t(abundace_table), taxa_are_rows = TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(metadata))


#Bray distance matrix
bray <- vegdist(abundace_table, method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, distance='bray', na.action='na.omit')

#Turn Eigenvalues into %
explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

permanova = adonis2(bray ~ metadata$Joint) #0.001 ***

##########Dispersion
dispr <- vegan::betadisper(bray, metadata$Joint) #test of the equality of variances.
permutest <- anova(dispr) 

#Clusters for the elipses
H_CLustering=hclust(bray)

#assign colours to each treatment
colours <- c("#c95292",
             "#7ecf5b",
             "#8449c2",
             "#cdb353",
             "#4d304a",
             "#8ec9b4",
             "#c25942",
             "#9695c1",
             "#58643d")

groups <- unique(metadata$Joint)
names(colours) = groups

H_CLustering=hclust(bray)

#Setting treatment variable
treatment <- unique(meta_bac$Joint)

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
#ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
#ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 3, bg = colours, pch= 21, lwd = 0.5, col= "black")
legend("topright", inset = c(-0.4, 0.01), legend= groups, col= colours, pch=19, cex=0.5, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

#################################Pairwise Location########## 01/07/2022
library(devtools)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

new <- cbind(abundace_table,metadata$Joint)
new_meta <- t(subset(metadata, select = Joint))

#Pairwise adonis
pair.mod <- pairwise.adonis2(abundace_table ~ Joint, data = new)
result <- dplyr::bind_rows(test[2:37], .id= "id")

write.csv(result,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Taxa-figures/extras/pairwise-adonis-full", row.names = TRUE)

#cleaing the table
clean_result <- result
clean_result <- clean_result[-c(2:5),]

clean_result <- clean_result[,grepl('Joint', colnames(clean_result))]
clean_result <- t(clean_result)

new_df <- matrix(ncol = 9, nrow = 9)
colnames(new_df) <- unique(metadata$Joint)
rownames(new_df) <- unique(metadata$Joint)


#terrible double nested for loop but does the job --> matrix with p-values
for(i in 1:length(rownames(clean_result))){
  for(j in 1:length(rownames(new_df))){
    for(k in 1:length(colnames(new_df))){
      foo <- unlist(strsplit(clean_result[i,1], split = "_vs_"))
      if(all(foo[1] == colnames(new_df)[j] & foo[2] == rownames(new_df)[k]) | all(foo[2] == colnames(new_df)[j] & foo[1] == rownames(new_df)[k])){
        new_df[j,k] <- clean_result[i,2]
      }}}}


write.csv(new_df,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Taxa-figures/extras/pairwise-adonis-p-vals.csv", row.names = TRUE)

####Biplots

#colours
treatment <- c("Daphnia"="#96568b","Bacterioplankton"="#99a666")
otu_taxa <- abundace_table
taxa_table$combined <- paste(taxa_table$Family, taxa_table$Genus, taxa_table$Species, sep=" ")

#replace ASVs with taxa-names
colnames(otu_taxa) <- taxa_table$combined[match(colnames(otu_taxa), rownames(taxa_table))]
colnames(otu_taxa) <- gsub("NA", "Unknown", colnames(otu_taxa))
#Delete taxa that only has unknowns (it is not identified)
colnames(otu_taxa) <- colnames(otu_taxa)[!grepl("Unknown Unknown Unknown", colnames(otu_taxa))] 
#Delete all unknowns
colnames(otu_taxa) <- gsub("Unknown", " ", colnames(otu_taxa))

#Preparing metadata
#Loading full metadata with added parameters ie. DO, temp etc. for the treatment ennvfit

#["pH","Redox","Temperature...C.","LDO..mg.L.","total.Chla...µg.l..", "Conductivity"]

meta = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/TRIAL1.csv", row.names=1, sep=";", dec=".")
meta_full <- as.data.frame(sapply(meta, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- rownames(meta)

meta_full=assignMetadataTypes(meta_full,categoric=c("experiment","Name", "Layer", "Joint")) #Set as categorical

#Change blanks to NA
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

groups=as.vector(metadata$experiment)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=20) #Replace NAs with group's means


#Plot
seqPCoA(t(otu_taxa), metadata=meta_full,groups=groups, pAdjMethod="fdr", arrowFactor= 0.000105, metadataFactor= 0.000000000000005, hideGroupLegend= TRUE, drawEllipse = TRUE, ellipseConf=0.95, colors= treatment[groups], xlim= c(-0.40, 0.45), ylim=c(-0.40, 0.40), taxonColor= "#2E2E2E")
legend("topright", legend=c("Daphnia","Bacterioplankton"), col= c("#96568b", "#99a666"),pch = 19, cex=0.5, pt.cex=1.5)
abline(h=-10:10, v=-10:10, col = "gray", lty = "dotted", lwd = 0.7)
points(abundace_table,cex= 1, bg = color_origin[treatment], pch = 21, lwd = 0.5, col= "black")
text(x= -0.3 , y= 0.2 ,paste("PERMANOVA R2: 0.0359, p-value: 0.001 ", pval, sep = ""), cex = 1)

#################Bacterioplankton###############################################

meta_bac <- metadata[metadata$experiment == 'Bacterioplankton',] #Filter out only Bacterioplankton data
otu_bac <- abundace_table %>% 
  filter(rownames(abundace_table) %in% rownames(meta_bac)) #Filter out only Bacterioplankton data

ps <- phyloseq(otu_table(t(otu_bac), taxa_are_rows =TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(meta_bac))

filt_ps_bac <- prevFilter(ps, 0.1) #Filter taxa that are not in at least 0.1 of samples
otu_bac <- as.data.frame(otu_table(filt_ps_bac))

#Bray distance matrix
bray <- vegdist(t(otu_bac), method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

permanova <- adonis2(bray ~ meta_bac$Type)# 0.049 * significant but just

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ meta_bac$Type)
# Df SumOfSqs      R2      F Pr(>F)  
# meta_bac$Type  1   0.5727 0.01824 1.6724  0.049 *
#   Residual      90  30.8174 0.98176                
# Total         91  31.3900 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 

##########Dispersion
dispr <- vegan::betadisper(bray, meta_bac$Type) #test of the equality of variances.
permutest <- anova(dispr) #pvalue = 0.6502 ** samples are homogenous

#####PCoA plot with experiment == Daphnia vs Bacterioplankton

#Clusters for the elipses
H_CLustering=hclust(bray)

#Colour for the points
color_origin <- c("1"="#a94754","2"="#8eca66", "3"="#8a51b1","4"="#b98a46","5"="#728a8f")

#Setting treatment variable
treatment <- meta_bac$Type

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 1, bg = color_origin[treatment], pch= 21, lwd = 0.5, col= "black")
legend("bottomright", legend=c("1","2","3","4","5"), col=c("#a94754","#8eca66","#8a51b1","#b98a46","#728a8f"), pch=19, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

###########Triplot###################
#Preparing metadata
#Loading full metadata with added parameters ie. DO, temp etc. for the treatment ennvfit

#["pH","Redox","Temperature...C.","LDO..mg.L.","total.Chla...µg.l..", "Conductivity"]

meta = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/TRIAL1.csv", row.names=1, sep=";", dec=".")
meta_full <- as.data.frame(sapply(meta, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- rownames(meta)

meta_full=assignMetadataTypes(meta_full,categoric=c("experiment","Name", "Layer", "Joint")) #Set as categorical

#Change blanks to NA
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

#Keep only bacterioplankton

meta_full <- meta_full[meta_full$experiment == 'Bacterioplankton',] #Filter out only Bacterioplankton data
otu_taxa_bac <- otu_taxa[rownames(otu_taxa) %in% rownames(meta_full),]

groups=as.vector(meta_full$Type)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=20) #Replace NAs with group's means

treatment <- c("1"="#8a51b2",
               "2"="#8ec765",
               "3"="#a84655",
               "4"="#71888d",
               "5"="#bf8847")

#Plot
seqPCoA(t(otu_taxa_bac), metadata=meta_full,groups=groups, pAdjMethod="fdr", arrowFactor= 0.000105, metadataFactor= 0.000000000000005, hideGroupLegend= TRUE, drawEllipse = TRUE, ellipseConf=0.95, colors= treatment[groups], xlim= c(-0.40, 0.45), ylim=c(-0.40, 0.40), taxonColor= "#2E2E2E")
legend("bottomright", legend=c("1","2","3","4","5"), col= c("#8a51b2","#8ec765","#a84655","#71888d","#bf8847"),pch = 19, cex=0.5, pt.cex=1.5)
abline(h=-10:10, v=-10:10, col = "gray", lty = "dotted", lwd = 0.7)
points(abundace_table,cex= 1, bg = color_origin[treatment], pch = 21, lwd = 0.5, col= "black")
text(x= -0.3 , y= 0.2 ,paste("PERMANOVA R2: 0.0551, p-value: 0.0809 ", sep = ""), cex = 1)

# [1] "0 significant numeric metadata found, in order of significance:"
# [1] "0 significant categoric metadata found, in order of significance:"
# 'adonis' will be deprecated: use 'adonis2' instead
# [1] "Adonis to test for significant difference in group compositions"
# [1] "Adonis R2: 0.0551, p-value: 0.0889"
# [1] "Cluster quality index silhouette"
# [1] -0.03787475
# [1] "Among the top 10 covarying taxa, 8 are significant."
# [1] "T34    "
# [1] "Moraxellaceae Acinetobacter  "
# [1] "Flavobacteriaceae Flavobacterium terrigena"
# [1] "Mycoplasmataceae    "
# [1] "Comamonadaceae Limnohabitans  "
# [1] "Sphingomonadaceae Sphingomonas  "
# [1] "Mycoplasmataceae Candidatus Bacilloplasma  "
# [1] "Microbacteriaceae Candidatus Limnoluna"
# 

#################Name######01/07/2022###Bacterioplankton

meta_bac <- metadata[metadata$experiment == 'Bacterioplankton',] #Filter out only Bacterioplankton data
otu_bac <- abundace_table %>% 
  filter(rownames(abundace_table) %in% rownames(meta_bac)) #Filter out only Bacterioplankton data

ps <- phyloseq(otu_table(t(otu_bac), taxa_are_rows =TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(meta_bac))

#Bray distance matrix
bray <- vegdist(otu_bac, method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue


#assign colours to each treatment
colours <- c("#cb50bf",
             "#84d14f",
             "#7042c9",
             "#ccb749",
             "#552d6c",
             "#79d2a1",
             "#cd4a74",
             "#536c39",
             "#727bcc",
             "#d25b34",
             "#77b0c3",
             "#7d3f34",
             "#c595bb",
             "#3b3b43",
             "#cdaf8c")

groups <- unique(meta_bac$Name)
names(colours) = groups

H_CLustering=hclust(bray)

#Setting treatment variable
treatment <- unique(meta_bac$Name)

#Statistical test
permanova <- adonis2(bray ~ meta_bac$Name)
dispr <- vegan::betadisper(bray, meta_bac$Name) #test of the equality of variances.
permutest <- anova(dispr) 

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
#ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
#ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 3, bg = colours, pch= 21, lwd = 0.5, col= "black")
legend("topright", inset = c(-0.4, 0.01), legend= groups, col= colours, pch=19, cex=0.5, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)


#################Daphnia###############################################

meta_daphnia <- metadata[metadata$experiment == 'Daphnia',] #Filter out only Daphnia data
otu_daphnia <- abundace_table %>% 
  filter(rownames(abundace_table) %in% rownames(meta_daphnia)) #Filter out only Daphnia data

ps <- phyloseq(otu_table(t(otu_daphnia), taxa_are_rows =TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(meta_daphnia))

filt_ps_dap <- prevFilter(ps, 0.1) #Filter taxa that are not in at least 0.1 of samples
otu_dap <- as.data.frame(otu_table(filt_ps_dap))

#Bray distance matrix
bray <- vegdist(t(otu_dap), method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

permanova <- adonis2(bray ~ meta_daphnia$Type)# 0.049 * significant but just

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ meta_daphnia$Type)
# Df SumOfSqs      R2      F Pr(>F)  
# meta_daphnia$Type  1   0.5103 0.02822 1.7134  0.069 .
# Residual          59  17.5717 0.97178                
# Total             60  18.0820 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########Dispersion
dispr <- vegan::betadisper(bray, meta_daphnia$Type) #test of the equality of variances.
permutest <- anova(dispr) #pvalue = 0.6502 ** samples are homogenous

#####PCoA plot with experiment == Daphnia vs Bacterioplankton

#Clusters for the elipses
H_CLustering=hclust(bray)

#Colour for the points
color_origin <- c("1"="#a94754","2"="#8eca66", "3"="#8a51b1","4"="#b98a46")

#Setting treatment variable
treatment <- meta_daphnia$Type

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 2.3, bg = color_origin[treatment], pch= 21, lwd = 0.5, col= "black")
legend("bottomright", legend=c("1","2","3","4"), col=c("#a94754","#8eca66","#8a51b1","#b98a46"), pch=19, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

###########Triplot###################
#Preparing metadata
#Loading full metadata with added parameters ie. DO, temp etc. for the treatment ennvfit

#["pH","Redox","Temperature...C.","LDO..mg.L.","total.Chla...µg.l..", "Conductivity"]

meta = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/TRIAL1.csv", row.names=1, sep=";", dec=".")
meta_full <- as.data.frame(sapply(meta, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- rownames(meta)

meta_full=assignMetadataTypes(meta_full,categoric=c("experiment","Name", "Layer", "Joint")) #Set as categorical

#Change blanks to NA
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

#Keep only bacterioplankton

meta_full <- meta_full[meta_full$experiment == 'Daphnia',] #Filter out only Bacterioplankton data
otu_taxa_daphnia <- otu_taxa[rownames(otu_taxa) %in% rownames(meta_full),]

groups=as.vector(meta_full$Type)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=20) #Replace NAs with group's means

treatment <- c("1"="#8a51b2",
               "2"="#8ec765",
               "3"="#a84655",
               "4"="#71888d")

#Plot
seqPCoA(t(otu_taxa_daphnia), metadata=meta_full,groups=groups, pAdjMethod="fdr", arrowFactor= 0.000105, metadataFactor= 0.000000000000005, hideGroupLegend= TRUE, drawEllipse = TRUE, ellipseConf=0.95, colors= treatment[groups], xlim= c(-0.40, 0.45), ylim=c(-0.40, 0.40), taxonColor= "#2E2E2E")
legend("bottomright", legend=c("1","2","3","4"), col= c("#8a51b2","#8ec765","#a84655","#71888d"),pch = 19, cex=0.5, pt.cex=1.5)
abline(h=-10:10, v=-10:10, col = "gray", lty = "dotted", lwd = 0.7)
points(abundace_table,cex= 1, bg = color_origin[treatment], pch = 21, lwd = 0.5, col= "black")
text(x= -0.3 , y= 0.2 ,paste("PERMANOVA R2: 0.0619, p-value: 0.1299 ", sep = ""), cex = 1)

#################Name######01/07/2022###Daphnia
meta_daphnia <- metadata[metadata$experiment == 'Daphnia',] #Filter out only Daphnia data
otu_daphnia <- abundace_table %>% 
  filter(rownames(abundace_table) %in% rownames(meta_daphnia)) #Filter out only Daphnia data

ps <- phyloseq(otu_table(t(otu_daphnia), taxa_are_rows =TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(meta_daphnia))

filt_ps_dap <- prevFilter(ps, 0.1) #Filter taxa that are not in at least 0.1 of samples
otu_dap <- as.data.frame(otu_table(ps))

#Bray distance matrix
bray <- vegdist(t(otu_dap), method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

#assign colours to each treatment
colours <- c("#8accb0",
             "#7e4cc7",
             "#8ed055",
             "#c954a1",
             "#cba553",
             "#46345e",
             "#d35345",
             "#9e9dc5",
             "#4c623c",
             "#824740")

groups <- unique(meta_daphnia$Name)
names(colours) = groups

H_CLustering=hclust(bray)

#Setting treatment variable
treatment <- unique(meta_daphnia$Name)

#Statistical test
permanova <- adonis2(bray ~ meta_daphnia$Name)# 0.049 * significant but just
dispr <- vegan::betadisper(bray, meta_daphnia$Name) #test of the equality of variances.
permutest <- anova(dispr) #pvalue = 0.6502 ** samples are homogenous

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
#ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
#ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 3, bg = colours, pch= 21, lwd = 0.5, col= "black")
legend("topright", inset = c(-0.4, 0.01), legend= groups, col= colours, pch=19, cex=0.5, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

#####################ASSOCIATIONS#####################
#Data preparation

#Reloading the data:
#rarefied:
taxa_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/rarefied-taxa.csv", row.names=1)
abundace_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/rarefied-asv.csv", row.names=1)
metadata = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/rarefied-meta.csv", row.names=1)

#metadata :: full 
meta = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/TRIAL1.csv", row.names=1, sep=";", dec=".")
meta_full <- as.data.frame(sapply(meta, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- rownames(metadata)

meta_full=assignMetadataTypes(meta_full,categoric=c("experiment","Name", "Layer", "Joint")) #Set as categorical

#Change blanks to NA based on the group average == Type
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

groups=as.vector(metadata$experiment)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=15) #Replace NAs with group's means

#modifications to the tables
taxa_table$combined <- paste(taxa_table$Family, taxa_table$Genus, taxa_table$Species, sep=" ")

ps <- phyloseq(otu_table(t(abundace_table), taxa_are_rows = TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(meta_full))

#Filter out taxa that are not in less than 10% of all samples
filt_ps_daphnia <- prevFilter(ps, 0.1) 

taxa_table <- as.data.frame(tax_table(filt_ps_daphnia))
otus <- as.data.frame(t(otu_table(filt_ps_daphnia)))
metadata <- as.data.frame(sample_data(filt_ps_daphnia))

#Cleaning taxa table

#Replacing ASVs with taxa names
taxa_table$combined <- gsub("NA", "Unknown", taxa_table$combined)
#Taxa with 3x Unknowns will be replace with their ASVs

for(i in 1:length(taxa_table$combined)){
  if(taxa_table$combined[i] == "Unknown Unknown Unknown"){
    taxa_table$combined[i] <- paste(taxa_table$Kingdom[i],taxa_table$Phylum[i],taxa_table$Class[i], sep=" ")
  }
}

#Delete all unknowns
taxa_table$combined <- gsub("Unknown", " ", taxa_table$combined)
colnames(otus) <- taxa_table$combined[match(colnames(otus), rownames(taxa_table))]

#Modfying metadata
metadata <- meta_full[,-c(5:11)]
metadata <- metadata[,-c(2,4)]

#experiment
metadata$Daphnia <- metadata[,1]
metadata$Bacterioplankton <- metadata[,1]

#Daphnia column
for(i in 1:length(metadata$Daphnia)){
  if(metadata$Daphnia[i] == "Daphnia"){
    metadata$Daphnia[i] <- "1"
  }else{
    metadata$Daphnia[i] <- "0"
  }
}

#Bacterioplankton column
for(i in 1:length(metadata$Bacterioplankton)){
  if(metadata$Bacterioplankton[i] == "Bacterioplankton"){
    metadata$Bacterioplankton[i] <- "1"
  }else{
    metadata$Bacterioplankton[i] <- "0"
  }
}

#Type
#make empty columnds first
metadata$DrinkingWaterReservoir <- 0
metadata$NatureReserve <- 0
metadata$RegularPond <- 0
metadata$PollutedPond <- 0
metadata$WTTP <- 0

for(i in 1:length(metadata$Type)){
  if(metadata$Type[i] == 1){
    metadata$DrinkingWaterReservoir[i] <- 1
  }else if(metadata$Type[i] == 2){
    metadata$NatureReserve[i] <- 1
  }else if(metadata$Type[i] == 3){
    metadata$RegularPond[i] <- 1
  }else if(metadata$Type[i] == 4){
    metadata$PollutedPond[i] <- 1
  }else if(metadata$Type[i] == 5){
    metadata$WTTP[i]  <- 1
  }
}

#Taxa names cannot repeat:
colnames(otus) <- make.unique(colnames(otus), sep="_")

write.csv(meta_full,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/network/network-meta.csv", row.names = TRUE)
write.csv(taxa_table,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/network/netowrk-taxa.csv", row.names = TRUE)
write.csv(otus,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/network/network-otus.csv", row.names = TRUE)

#separately: 
##Daphnia
##Bacterioplankton 
###Metadata: Pollution gradient


##############################Pairwise location 

taxa_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/rarefied-taxa.csv", row.names=1)
abundace_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/rarefied-asv.csv", row.names=1)
metadata = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Cleaned_tables/rarefied-meta.csv", row.names=1)


#modifications to the tables
metadata$Joint <- paste(metadata$experiment, metadata$Name) #making a join column for experiment and Type
taxa_table$combined <- paste(taxa_table$Family, taxa_table$Genus, taxa_table$Species, sep=" ")

#####################ps#############################################################
#incorporating all the elements into ps elements:
ps <- phyloseq(otu_table(t(abundace_table), taxa_are_rows = TRUE),
               tax_table(as.matrix(taxa_table)), sample_data(metadata))

#Bray distance matrix
bray <- vegdist(abundace_table, method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, distance='bray', na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

#Significance
permanova = adonis2(bray ~ metadata$Joint)
dispr <- vegan::betadisper(bray, metadata$Joint) #test of the equality of variances.
permutest <- anova(dispr) #Not significant

#Clusters for the elipses
H_CLustering=hclust(bray)

#Stat test
pval <- permanova$`Pr(>F)`[1]

#Dispersion
disp_test <- permutest$`Pr(>F)`[1]

#assign colours to each treatment
colours <- c("#de8296",
             "#d93fb7",
             "#e0b2c6",
             "#de3864",
             "#9a6982",
             "#e34993",
             "#944759",
             "#dd83bc",
             "#bb476b",
             "#a94382",
             "#686d51",
             "#9be83f",
             "#636627",
             "#59d955",
             "#a7a881",
             "#dbdd41",
             "#57814b",
             "#d2ea8f",
             "#489037",
             "#d8deb6",
             "#81b837",
             "#97c386",
             "#88892d",
             "#81db80",
             "#bcb656")

groups1 <- unique(metadata$Joint)
names(colours) = groups1


#x and y labels
xlabel= paste("PC1:", explainedvar1, "%")
ylabel= paste("PC2:", explainedvar2, "%")

par(mar = c(4.1,3, 4.1, 2.1))
plot(PCoA,main="PCoA: Bray-Curtis Dissimilarity", type="n", xlab="", ylab="",ylim=c(-2, 2), xlim=c(-2, 2), cex.axis=1, tck = -0.01, mgp = c(2, 0.5, 0),
     xaxp  = c(-4, 4, 2), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=1.5, cex.lab=1)
title(xlab=xlabel, line=1.5, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
ordiellipse(PCoA, groups = treatment,  kind = "sd",conf= 0.95, draw ="polygon", col = "azure", lwd = 0.5, lty=2)
par(lty=2)
ordicluster(PCoA, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
points(PCoA,cex= 1.5, bg = colours, pch= 21, lwd = 0.5, col= "black")
legend("bottomright", legend= metadata$Joint, col= colours, pch= symbols, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)


######Pairwise Name
library(pairwiseAdonis)

new <- cbind(abundace_table,metadata$Joint)
new_meta <- t(subset(metadata, select = Joint))

#Pairwise adonis
pair.mod <- pairwise.adonis2(new ~ Joint, data = new)
result <- dplyr::bind_rows(test[2:37], .id= "id")

write.csv(result,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Taxa-figures/extras/pairwise-adonis-full", row.names = TRUE)

#cleaing the table
clean_result <- result
clean_result <- clean_result[-c(2:5),]

clean_result <- clean_result[,grepl('Joint', colnames(clean_result))]
clean_result <- t(clean_result)

new_df <- matrix(ncol = 9, nrow = 9)
colnames(new_df) <- unique(metadata$Joint)
rownames(new_df) <- unique(metadata$Joint)


#terrible double nested for loop but does the job --> matrix with p-values
for(i in 1:length(rownames(clean_result))){
  for(j in 1:length(rownames(new_df))){
    for(k in 1:length(colnames(new_df))){
      foo <- unlist(strsplit(clean_result[i,1], split = "_vs_"))
      if(all(foo[1] == colnames(new_df)[j] & foo[2] == rownames(new_df)[k]) | all(foo[2] == colnames(new_df)[j] & foo[1] == rownames(new_df)[k])){
        new_df[j,k] <- clean_result[i,2]
      }}}}


###########ALDEX####################################


