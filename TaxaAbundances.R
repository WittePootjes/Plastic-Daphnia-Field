library(ggplot2)
library(phyloseq)
library(vegan)
library("ALDEx2")
library(dplyr)
library(ggpubr)
library(otuSummary)

otus <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/asv-no-chloro-mito.csv", row.names=1)
taxa <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/taxa-no-chloro-mito.csv", row.names=1)
meta <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/metadata-final.csv", row.names=1)


######Taxonomic graphs
options(scipen = 999)
# Normalize counts by the sample sums
sums <- apply(otus, 2, sum)

for (i in 1:ncol(otus)){
  otus[,i] <- otus[,i]/sums[i]
}

#Phyloseq object
ps <- phyloseq(otu_table(otus, taxa_are_rows = FALSE),
               tax_table(as.matrix(taxa)), sample_data(meta))


#Filtering out taxa
prevFilter = function(phy, prev){
  #' @title Prevalence filter
  #'
  #' @description Filters taxa present in fewer samples than the specified fraction. 
  #' @details Filtered taxa are retained in the Bin taxon to preserve sample sums. 
  #'
  #' @param phy Phyloseq object
  #' @param prev Minimum fraction of samples
  #'
  prev0 = apply(X = otu_table(phy),
                MARGIN = ifelse(taxa_are_rows(phy), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
  prevalenceThreshold = prev * nsamples(phy)
  nonrares = prune_taxa((prev0 > prevalenceThreshold), phy)
  rares = prune_taxa((prev0 < prevalenceThreshold), phy)
  rares = merge_taxa(rares, taxa_names(rares))
  otus = data.frame(otu_table(nonrares))
  otus = cbind(otus, data.frame(otu_table(rares)))
  tax = data.frame(tax_table(nonrares), stringsAsFactors = FALSE)
  tax = rbind(tax, rep("Other", 7))
  otus <- t(otus)
  rownames(otus) <- c(taxa_names(nonrares), 'Bin')
  rownames(tax) <- c(taxa_names(nonrares), 'Bin')
  newphy = phyloseq(otu_table(otus, taxa_are_rows = TRUE), sample_data(phy), tax_table(as.matrix(tax)))
  return(newphy)
}

#Bacterioplankton
metadata_bac <- meta[which(meta$experiment == 'Bacterioplankton'),names(meta)]
ps_bac = prune_samples(rownames(metadata_bac), ps)

prev60 <-prevFilter(ps_bac, 0.6) #Filtering at 60% prevalence 

tax <-data.frame(tax_table(prev60))
otus <-data.frame(otu_table(prev60))
otus <- head(otus, - 1)  
tax <- head(tax, - 1) 

long_data <-cbind(tax, otus)
long_data <- reshape2::melt(long_data)

#Adding treatments

#Adding Type
long_data$Type <- as.character(long_data$variable)
for (i in 1:length(long_data$Type)){
  sample_name <- long_data$variable[[i]]
  long_data$Type[[i]] <- as.character(metadata_bac$Type[rownames(metadata_bac) == sample_name])
}


#Adding Layer
meta_layer <- metadata_bac[!(is.na(metadata_bac$Layer) | metadata_bac$Layer==""), ]
ps_bac_layer = prune_samples(rownames(meta_layer), ps)

prev60 <-prevFilter(ps_bac_layer, 0.6) #Filtering at 60% prevalence 
tax <-data.frame(tax_table(prev60))
otus <-data.frame(otu_table(prev60))
otus <- head(otus, - 1)  
tax <- head(tax, - 1) 

long_data <-cbind(tax, otus)
long_data <- reshape2::melt(long_data)

long_data$Layer <- as.character(long_data$variable)
for (i in 1:length(long_data$Layer)){
  sample_name <- long_data$variable[[i]]
  long_data$Layer[[i]] <- as.character(meta_layer$Layer[rownames(meta_layer) == sample_name])
}

colourCount = length(rownames(tax))
mycolors <- c("#a5495f",
              "#50c150",
              "#c650ba",
              "#81b748",
              "#7b63d4",
              "#b5bf3a",
              "#5f6cb8",
              "#e68320",
              "#659ad7",
              "#dbaf3f",
              "#c590da",
              "#498b2a",
              "#da3e7f",
              "#56d1a2",
              "#d53f41",
              "#3cbed0",
              "#c05727",
              "#4cac93",
              "#994e89",
              "#65b573",
              "#e17fa6",
              "#3c8446",
              "#b05148",
              "#297250",
              "#e98877",
              "#52641b",
              "#cd8a37",
              "#7f8d48",
              "#94632e",
              "#aab56e",
              "#827223",
              "#d3a36c",
              "#9d972e",
              "#77703b")

colourCount = length(rownames(tax))
mycolors <- colorRampPalette(brewer.pal(6,"Set3"))(6)

#############Creating a plot with filtered out genus at over 60% frequency
ggplot(data=long_data, aes(x=variable, y=value, fill=Phylum)) + # base ggplot2 data
  geom_bar(position='fill', stat='identity') + # defines stacked bar plot
  labs(x='', y='Relative abundance for Phylum at over 60% frequency') + 
  theme_minimal() + # adds titles and minimal look
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position="bottom") + 
  facet_wrap(vars(Type), ncol=2, scales='free')+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  scale_fill_manual(values = mycolors)

#############Creating a plot with filtered out Phylum at over 60% frequency
ggplot(data=long_data, aes(x=variable, y=value, fill=Genus)) + # base ggplot2 data
  geom_bar(position='fill', stat='identity') + # defines stacked bar plot
  labs(x='', y='Relative abundance for Genus at over 60% frequency') + 
  theme_pubclean() + # adds titles and minimal look
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position="bottom") + 
  facet_wrap(vars(Layer), ncol=2, scales='free')+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  scale_fill_manual(values = mycolors)

##########Only Daphnia 
metadata_dap <- meta[which(meta$experiment == 'Daphnia'),names(meta)]
ps_dap = prune_samples(rownames(metadata_dap), ps)

prev60 <-prevFilter(ps_dap, 0.4) #Filtering at 60% prevalence 
tax <-data.frame(tax_table(prev60))
otus <-data.frame(otu_table(prev60))
otus <- head(otus, - 1)  
tax <- head(tax, - 1) 

long_data <-cbind(tax, otus)
long_data <- reshape2::melt(long_data)
long_data
length(long_data)

#Adding treatments

#Adding Type
long_data$Type <- as.character(long_data$variable)
for (i in 1:length(long_data$Type)){
  sample_name <- long_data$variable[[i]]
  long_data$Type[[i]] <- as.character(metadata_dap$Type[rownames(metadata_dap) == sample_name])
}

colourCount = length(rownames(tax))
mycolors <- c("#a5495f",
               "#50c150",
               "#c650ba",
               "#81b748",
               "#7b63d4",
               "#b5bf3a",
               "#5f6cb8",
               "#e68320",
               "#659ad7",
               "#dbaf3f",
               "#c590da",
               "#498b2a",
               "#da3e7f",
               "#56d1a2",
               "#d53f41",
               "#3cbed0",
               "#c05727",
               "#4cac93",
               "#994e89",
               "#65b573",
               "#e17fa6",
               "#3c8446",
               "#b05148",
               "#297250",
               "#e98877",
               "#52641b",
               "#cd8a37",
               "#7f8d48",
               "#94632e",
               "#aab56e",
               "#827223",
               "#d3a36c",
               "#9d972e",
               "#77703b")

mycolors <- colorRampPalette(brewer.pal(6,"Set3"))(6)
#############Creating a plot with filtered out genus at over 40% frequency

ggplot(data=long_data, aes(x=variable, y=value, fill=Phylum)) + # base ggplot2 data
  geom_bar(position='fill', stat='identity') + # defines stacked bar plot
  labs(x='', y='Relative abundance for Phylum at over 40% frequency') + 
  theme_minimal() + # adds titles and minimal look
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position="bottom") + 
  facet_wrap(vars(Type), ncol=2, scales='free')+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  scale_fill_manual(values = mycolors)

#####Plastic

otus <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/asv-no-chloro-mito.csv", row.names=1)
taxa <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/taxa-no-chloro-mito.csv", row.names=1)
meta <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/metadata-final.csv", row.names=1)

#Phyloseq object
ps <- phyloseq(otu_table(otus, taxa_are_rows = FALSE),
               tax_table(as.matrix(taxa)), sample_data(meta))

plast_bac <- c("Flavobacterium","Streptomyces","Pseudomonas", "Alcaligenes","Ralstonia","Acidovorax","Comamonas", "Bdellovibrio","Amycolatopsis","Ideonella", "Saccharomonospora","Thermobifida","Cryptococcus","Rhodococcus", "Schlegelella","Geobacillus","Caldimonas","Marinobacter","Diaphorobacter","Cupriavidus","Azotobacter","Aspergillus","Burkholderia","Paenibacillus","Rhizopus","Pseudozyma","Fusarium","Undibacterium","Paracoccus","Roseateles","Kocuria","Thermobifida","Thermomonospora","Oleispira")

#Filter out Daphnia
Burk_daphnia <- subset_samples(ps, experiment == "Daphnia")
#Filter out Bakterioplankton
Burk_bac <- subset_samples(ps, experiment == "Bacterioplankton")

#Filter out plastic genera
Burk = subset_taxa(Burk_daphnia, Genus %in% plast_bac)
Burk = subset_taxa(Burk_bac, Genus %in% plast_bac)

tax <-data.frame(tax_table(Burk))
otus <-data.frame(otu_table(Burk))
otus <- t(otus)

long_data <-cbind(tax, otus)
long_data <- long_data[,-c(1:5)]
long_data <- long_data[,-2]
long_data <- reshape2::melt(long_data, id.vars="Genus")

#Adding Type
long_data$Type <- as.character(long_data$variable)
for (i in 1:length(long_data$Type)){
  sample_name <- long_data$variable[[i]]
  long_data$Type[[i]] <- as.character(meta$Type[rownames(meta) == sample_name])
}

colourCount = length(plast_bac)
mycolors <- c("#a5495f",
               "#50c150",
               "#c650ba",
               "#81b748",
               "#7b63d4",
               "#b5bf3a",
               "#5f6cb8",
               "#e68320",
               "#659ad7",
               "#dbaf3f",
               "#c590da",
               "#498b2a",
               "#da3e7f",
               "#56d1a2",
               "#d53f41",
               "#3cbed0",
               "#c05727",
               "#4cac93",
               "#994e89",
               "#65b573",
               "#e17fa6",
               "#3c8446",
               "#b05148",
               "#297250",
               "#e98877",
               "#52641b",
               "#cd8a37",
               "#7f8d48",
               "#94632e",
               "#aab56e",
               "#827223",
               "#d3a36c",
               "#9d972e",
               "#77703b")

ggplot(data=long_data, aes(x=Type, y=value, fill=Genus)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x='Source', y='Abundance') + 
  theme_minimal()+
  facet_wrap(~ Genus, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=0.1), legend.position="bottom")+
  coord_flip()+
  theme_minimal()



