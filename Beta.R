library(ggplot2)
library(phyloseq)
library(vegan)
library("ALDEx2")
library(dplyr)
library(ggpubr)
library(otuSummary)

#####################BETA DIVERSITIES AND ORDINATION#############################################################
#LOAD THE DATA
otus <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/final-asvs-for-ordination.csv", row.names = 1)
taxa <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/final-taxa-for-ordination.csv", row.names = 1)
meta <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/final-meta-for-ordination.csv", row.names = 1)

####modifications to the tables
meta$Joint <- paste(meta$experiment, meta$Type) #making a join column for experiment and Type
taxa$combined <- paste(taxa$Family, taxa$Genus, taxa$Species, sep=" ")

#Replacing ASVs with taxa names
taxa$combined <- gsub("NA", "Unknown", taxa$combined)
#Taxa with 3x Unknowns will be replace with their ASVs

for(i in 1:length(taxa$combined)){
  if(taxa$combined[i] == "Unknown Unknown Unknown"){
    taxa$combined[i] <- paste(taxa$Kingdom[i],taxa$Phylum[i],taxa$Class[i], sep=" ")
  }
}

#Delete all unknowns
taxa$combined <- gsub("Unknown", " ", taxa$combined)

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
points(PCoA,cex= 1, bg = color_origin[treatment], pch= 21, lwd = 0.5, col= "black")
legend("bottomright", inset = c(-0.4, 0.01), legend=c("Daphnia","Bacterioplankton"), col=c("#96568b","#99a666"), pch=19, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

#################Bacterioplankton###############################################

meta_bac <- meta[meta$experiment == 'Bacterioplankton',] #Filter out only Bacterioplankton data
otu_bac <- otus %>% 
  filter(rownames(otus) %in% rownames(meta_bac)) #Filter out only Bacterioplankton data

#Bray distance matrix
bray <- vegdist(otu_bac, method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

permanova <- adonis2(bray ~ meta_bac$Type)# 0.001 ***

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ meta_bac$Type)
# Df SumOfSqs      R2      F Pr(>F)    
# meta_bac$Type  1    2.548 0.07707 7.2655  0.001 ***
#   Residual      87   30.515 0.92293                  
# Total         88   33.063 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########Dispersion
dispr <- vegan::betadisper(bray, meta_bac$Type) #test of the equality of variances.
permutest <- anova(dispr) #pvalue = 0.1027892 ** samples are homogenous

#####PCoA plot 

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
points(PCoA,cex= 1.5, bg = color_origin[treatment], pch= 21, lwd = 0.5, col= "black")
legend("bottomright",inset = c(-0.4, 0.01), legend=c("1","2","3","4","5"), col=c("#a94754","#8eca66","#8a51b1","#b98a46","#728a8f"), pch=19, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

###########Triplot###################
#Preparing metadata
#Loading full metadata with added parameters ie. DO, temp etc. for the treatment ennvfit

#["pH","Redox","Temperature...C.","LDO..mg.L.","total.Chla...µg.l..", "Conductivity"]

meta_full = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/TRIAL1.csv", row.names=1, sep=";")
all(rownames(meta_full)==rownames(meta))
meta_full <- as.data.frame(sapply(meta_full, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- rownames(meta)

meta_full=assignMetadataTypes(meta_full,categoric=c("Name")) #Set as categorical

#Change blanks to NA
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

#one hot encoding
library(mltools)
library(data.table)

meta_full <- one_hot(as.data.table(meta_full), cols = "Name")
rownames(meta_full) <- rownames(meta)

#Keep only bacterioplankton
meta_full <- meta_full[meta_full$Daphnia == 0,] #Filter out only Bacterioplankton data
otu_taxa_bac <- otus[rownames(otus) %in% rownames(meta_full),]

groups=as.vector(meta_full$Type)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=8) #Replace NAs with group's means

treatment <- c("1"="#8a51b2",
               "2"="#8ec765",
               "3"="#a84655",
               "4"="#71888d",
               "5"="#bf8847")

#Rename ASVs to taxa
colnames(otu_bac) <- taxa$combined[match(colnames(otu_bac), rownames(taxa))]

#Plot
seqPCoA(t(otu_bac), metadata=meta_full,groups=groups, pAdjMethod="fdr", arrowFactor= 0.0003, metadataFactor= 1.1, hideGroupLegend= TRUE, drawEllipse = TRUE, ellipseConf=0.95, colors= treatment[groups], xlim= c(-0.40, 0.45), ylim=c(-0.40, 0.40), taxonColor= "#2E2E2E", topTaxa= 10)
legend("topleft", legend=c("1","2","3","4","5"), col= c("#8a51b2","#8ec765","#a84655","#71888d","#bf8847"),pch = 19, cex=0.5, pt.cex=1.5)
abline(h=-10:10, v=-10:10, col = "gray", lty = "dotted", lwd = 0.7)
text(x= -0.3 , y= 0.2 ,paste("PERMANOVA R2: 0.2233, p-value: 0.001", sep = ""), cex = 1)

# 
# [1] "17 significant numeric metadata found, in order of significance:"
# [1] "Name_Evangelie Boom"
# [1] "Name_Harelbeke"
# [1] "Name_Meer van Rotselaar"
# [1] "Name_Olsene"
# [1] "Name_St. Donatus Park"
# [1] "Type"
# [1] "pH"
# [1] "Oxygen...."
# [1] "Temperature...C."
# [1] "LDO..mg.L."
# [1] "Name_MVR"
# [1] "Redox"
# [1] "Name_Leuven"
# [1] "Conductivity"
# [1] "Name_Kluizen"
# [1] "Name_LRF"
# [1] "Name_Kortrijk Blauwwe Poort Park"
# [1] "0 significant categoric metadata found, in order of significance:"
# 'adonis' will be deprecated: use 'adonis2' instead
# [1] "Adonis to test for significant difference in group compositions"
# [1] "Adonis R2: 0.2233, p-value: 0.001"
# [1] "Cluster quality index silhouette"
# [1] 0.08285213
# [1] "Among the top 10 covarying taxa, 6 are significant."
# [1] "Spirosomaceae Pseudarcicella  "
# [1] "NS11-12 marine group    "
# [1] "Armatimonadaceae Armatimonas  "
# [1] "Microbacteriaceae Candidatus Limnoluna"
# [1] "Microbacteriaceae Aurantimicrobium minutum"
# [1] "Burkholderiaceae Polynucleobacter cosmopolitanus"

######Pairwise############

######Pairwise Name
library(pairwiseAdonis)

new <- cbind(otu_bac,meta_bac$Type)
new_meta <- t(subset(metadata, select = Joint))

#Pairwise adonis
pair.mod <- pairwise.adonis2(new ~ Type, data = new)
result <- dplyr::bind_rows(pair.mod[2:37], .id= "id")

write.csv(result,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Taxa-figures/extras/pairwise-adonis-full", row.names = TRUE)

#cleaing the table
clean_result <- result
clean_result <- clean_result[,-c(2:5)]

for(i in 1:length(clean_result$`Pr(>F)`)){
  if(is.na(clean_result$`Pr(>F)`)[i]){
    clean_result <- clean_result[-i,]
  }
}

new_df <- matrix(ncol = 5, nrow = 5)
colnames(new_df) <- unique(meta_bac$Type)
rownames(new_df) <- unique(meta_bac$Type)


#terrible double nested for loop but does the job --> matrix with p-values
for(i in 1:length(rownames(clean_result))){
  for(j in 1:length(rownames(new_df))){
    for(k in 1:length(colnames(new_df))){
      foo <- unlist(strsplit(clean_result[i,1], split = "_vs_"))
      if(all(foo[1] == colnames(new_df)[j] & foo[2] == rownames(new_df)[k]) | all(foo[2] == colnames(new_df)[j] & foo[1] == rownames(new_df)[k])){
        new_df[j,k] <- clean_result[i,2]
      }}}}


write.csv(new_df,"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/pairwise-adonis-p-vals-BACTERIOPLANKTON.csv", row.names = TRUE)

#################Daphnia###############################################

meta_dap <- meta[meta$experiment == 'Daphnia',] #Filter out only Bacterioplankton data
otu_dap <- otus %>% 
  filter(rownames(otus) %in% rownames(meta_dap)) #Filter out only Bacterioplankton data

#Bray distance matrix
bray <- vegdist(otu_dap, method="bray") #creating dissimilarity matrix with vegdist using bray curtis
#PCoA ordination
PCoA <- capscale(bray~1, na.action='na.omit')

#Turn Eigenvalues into %

explainedvar1 <- round(eigenvals(PCoA)[[1]] / sum(eigenvals(PCoA)), 2) * 100 #21% of diversity explained by the first eigenvalue
explainedvar2 <- round(eigenvals(PCoA)[[2]] / sum(eigenvals(PCoA)), 2) * 100 #10% explained by the second eigenvalue

permanova <- adonis2(bray ~ meta_dap$Type)# 0.001 ***

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ meta_dap$Type)
# Df SumOfSqs      R2      F Pr(>F)    
# meta_dap$Type  1   2.4083 0.11335 7.9261  0.001 ***
#   Residual      62  18.8380 0.88665                  
# Total         63  21.2463 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 


##########Dispersion
dispr <- vegan::betadisper(bray, meta_dap$Type) #test of the equality of variances.
permutest <- anova(dispr) #pvalue = 0.004522 ** samples are not homogenous

#Clusters for the elipses
H_CLustering=hclust(bray)

#Colour for the points
color_origin <- c("1"="#a94754","2"="#8eca66", "3"="#8a51b1","4"="#b98a46")

#Setting treatment variable
treatment <- meta_dap$Type

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
points(PCoA,cex= 1.5, bg = color_origin[treatment], pch= 21, lwd = 0.5, col= "black")
legend("bottomright",inset = c(-0.4, 0.01), legend=c("1","2","3","4","5"), col=c("#a94754","#8eca66","#8a51b1","#b98a46","#728a8f"), pch=19, cex=0.8, pt.cex=1.5 )
text(x= -3 , y= 1.5 ,expression(paste("Permanova")), cex = 1)
text(x= -3 , y= 1.4 ,paste("p-value = ", pval, sep = ""), cex = 1)
text(x= -3 , y= 1.2 ,expression(paste("Betadisper")), cex = 1)
text(x= -3 , y= 1.1 ,paste("p-value = ", disp_test, sep = ""), cex = 1)

###########Triplot###################
#Preparing metadata
#Loading full metadata with added parameters ie. DO, temp etc. for the treatment ennvfit

#["pH","Redox","Temperature...C.","LDO..mg.L.","total.Chla...µg.l..", "Conductivity"]

meta_full = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/TRIAL1.csv", row.names=1, sep=";")
all(rownames(meta_full)==rownames(meta))
meta_full <- as.data.frame(lapply(meta_full, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- rownames(meta)

meta_full=assignMetadataTypes(meta_full,categoric=c("Name")) #Set as categorical

#Change blanks to NA
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

#one hot encoding
library(mltools)
library(data.table)

meta_full <- one_hot(as.data.table(meta_full), cols = "Name")
rownames(meta_full) == rownames(meta)
rownames(meta_full) <- rownames(meta)

#Keep only Daphnia
meta_full <- meta_full[meta_full$Daphnia == 1,] #Filter out only Bacterioplankton data
meta_full <- as.data.frame(meta_full)
rownames(meta_full) <- rownames(meta)

groups=as.vector(meta_full$Type)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=8) #Replace NAs with group's means

treatment <- c("1"="#8a51b2",
               "2"="#8ec765",
               "3"="#a84655",
               "4"="#71888d")

#Rename ASVs to taxa
colnames(otu_dap) <- taxa$combined[match(colnames(otu_dap), rownames(taxa))]

#Plot
seqPCoA(t(otu_dap), metadata=meta_full,groups=groups, pAdjMethod="fdr", arrowFactor= 0.00003, metadataFactor= 1.1, hideGroupLegend= TRUE, drawEllipse = TRUE, ellipseConf=0.95, colors= treatment[groups], xlim= c(-0.40, 0.45), ylim=c(-0.40, 0.40), taxonColor= "#2E2E2E", topTaxa= 10)
legend("topleft", legend=c("1","2","3","4"), col= c("#8a51b2","#8ec765","#a84655","#71888d"),pch = 19, cex=0.5, pt.cex=1.5)
abline(h=-10:10, v=-10:10, col = "gray", lty = "dotted", lwd = 0.7)
text(x= -0.3 , y= 0.2 ,paste("PERMANOVA R2: 0.2233, p-value: 0.001", sep = ""), cex = 1)

# [1] "12 significant numeric metadata found, in order of significance:"
# [1] "Name_LRF"
# [1] "Name_Reserve next to Kulak"
# [1] "Type"
# [1] "Redox"
# [1] "Oxygen...."
# [1] "LDO..mg.L."
# [1] "Name_Kortrijk Blauwwe Poort Park"
# [1] "Temperature...C."
# [1] "Name_De Bourgoyen"
# [1] "Conductivity"
# [1] "Name_St. Donatus Park"
# [1] "pH"
# [1] "0 significant categoric metadata found, in order of significance:"
# 'adonis' will be deprecated: use 'adonis2' instead
# [1] "Adonis to test for significant difference in group compositions"
# [1] "Adonis R2: 0.2137, p-value: 0.001"
# [1] "Cluster quality index silhouette"
# [1] 0.1658904
# [1] "Among the top 10 covarying taxa, 4 are significant."
# [1] "T34    "
# [1] "Mycoplasmataceae    "
# [1] "Flavobacteriaceae Flavobacterium terrigena"
# [1] "Candidatus Hepatincola  


######Pairwise############

######Pairwise Name
library(pairwiseAdonis)

new <- cbind(otu_dap,meta_dap$Type)

#Pairwise adonis
pair.mod <- pairwise.adonis2(new ~ Type, data = new)
result <- dplyr::bind_rows(pair.mod[2:37], .id= "id")

write.csv(result,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Taxa-figures/extras/pairwise-adonis-full", row.names = TRUE)

#cleaing the table
clean_result <- result
clean_result <- clean_result[,-c(2:5)]

#two iterations needed
for(i in 1:length(clean_result$`Pr(>F)`)){
  if(is.na(clean_result$`Pr(>F)`)[i]){
    clean_result <- clean_result[-i,]
  }
}

new_df <- matrix(ncol = 4, nrow = 4)
colnames(new_df) <- unique(meta_dap$Type)
rownames(new_df) <- unique(meta_dap$Type)


#terrible double nested for loop but does the job --> matrix with p-values
for(i in 1:length(rownames(clean_result))){
  for(j in 1:length(rownames(new_df))){
    for(k in 1:length(colnames(new_df))){
      foo <- unlist(strsplit(clean_result[i,1], split = "_vs_"))
      if(all(foo[1] == colnames(new_df)[j] & foo[2] == rownames(new_df)[k]) | all(foo[2] == colnames(new_df)[j] & foo[1] == rownames(new_df)[k])){
        new_df[j,k] <- clean_result[i,2]
      }}}}


write.csv(new_df,"/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/pairwise-adonis-p-vals-DAPHNIA.csv", row.names = TRUE)









