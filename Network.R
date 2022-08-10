library(phyloseq)
library(vegan)
library(dplyr)


#####################ASSOCIATIONS#####################
#Data preparation

#Reloading the data:
#rarefied:
taxa_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/rarefied-taxa.csv", row.names=1)
abundace_table = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/rarefied-asv-numbered-correct.csv", sep=";",row.names=1)
rownames(abundace_table) <- abundace_table[,1]
abundace_table <- abundace_table[,-1]
meta_full = read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Results/TRIAL1.csv", row.names=1, sep=";")
meta_names <- rownames(meta_full)

#metadata :: full 
meta_full <- as.data.frame(sapply(meta_full, gsub, pattern = ",", replacement= ".")) #change comma separator to a dot
rownames(meta_full) <- meta_names

meta_full=assignMetadataTypes(meta_full,categoric=c("Name", "Layer")) #Set as categorical

#Change blanks to NA based on the group average == Type
meta_full[] <- lapply(meta_full, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

groups=as.vector(meta_full$Type)
meta_full <- setNAToGroupMean(meta_full, groups=groups, na.threshold=15) #Replace NAs with group's means

#modifications to the tables
taxa_table$combined <- paste(taxa_table$Family, taxa_table$Genus, taxa_table$Species, sep=" ")

ps <- phyloseq(otu_table(abundace_table, taxa_are_rows = FALSE),
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
meta_full$DrinkingWaterReservoir <- 0
meta_full$NatureReserve <- 0
meta_full$RegularPond <- 0
meta_full$PollutedPond <- 0
meta_full$WTTP <- 0

for(i in 1:length(meta_full$Type)){
  if(meta_full$Type[i] == 1){
    meta_full$DrinkingWaterReservoir[i] <- 1
  }else if(meta_full$Type[i] == 2){
    meta_full$NatureReserve[i] <- 1
  }else if(meta_full$Type[i] == 3){
    meta_full$RegularPond[i] <- 1
  }else if(meta_full$Type[i] == 4){
    meta_full$PollutedPond[i] <- 1
  }else if(meta_full$Type[i] == 5){
    meta_full$WTTP[i]  <- 1
  }
}

#Taxa names cannot repeat:
colnames(otus) <- make.unique(colnames(otus), sep="_")

#one hot-encoding

#one hot encoding for Name
library(mltools)
library(data.table)

meta_full <- one_hot(as.data.table(meta_full), cols = "Name")
rownames(meta_full) <- meta_names



write.csv(meta_full,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Results/Network/network-meta.csv", row.names = TRUE)
write.csv(otus,"/Users/u0145079/Desktop/Daphnia_MiSeq/Results/Results/Network/network-otus.csv", row.names = TRUE)
