library(ggplot2)
library(phyloseq)
library(vegan)
library("ALDEx2")
library(dplyr)
library(ggpubr)
library(otuSummary)

otus <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/asv-no-chloro-mito.csv", row.names=1)
taxa <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/taxa-no-chloro-mito.csv", row.names=1)
meta <- read.csv("/Users/u0145079/Desktop/Daphnia_MiSeq/Preprocessing check/metadata-final.csv", row.names=1)


#############For ALDEX2 we need our original dataset (raw_counts)
abund_table <- t(otus)

ps <- phyloseq(otu_table(abund_table, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxa)), sample_data(meta))

#Define prevalace filter again
prevFilter =function(phy, prev){
  prev0 = apply(X = otu_table(phy),
                MARGIN = ifelse(taxa_are_rows(phy), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
  prevalenceThreshold = prev * nsamples(phy)
  nonrares = prune_taxa((prev0 > prevalenceThreshold), phy)
  rares = prune_taxa((prev0 < prevalenceThreshold), phy)
  rares = merge_taxa(rares, taxa_names(rares))
  otus = data.frame(otu_table(nonrares))
  otus = rbind(otus, data.frame(otu_table(rares)))
  tax = data.frame(tax_table(nonrares), stringsAsFactors = FALSE)
  tax = rbind(tax, rep("Other", 7))
  rownames(otus) <- c(taxa_names(nonrares), 'Bin')
  rownames(tax) <- c(taxa_names(nonrares), 'Bin')
  newphy = phyloseq(otu_table(otus, taxa_are_rows = TRUE), sample_data(phy), tax_table(as.matrix(tax)))
  return(newphy)
}

#Only keep those ASVs that are in more than 10%
filt_ps <- prevFilter(ps, 0.1) #Only 287 ASVs left

#Reassign the variables
abund_table <- as.data.frame(otu_table(filt_ps))
taxa_table <- as.data.frame(tax_table(filt_ps))

########Experiment
#Define groups for ALDEX2
groups <- as.factor(meta$experiment)
#CLR tranformation
ald.clr <- aldex.clr(abund_table, groups, mc.samples=128, verbose=TRUE)
#T-test (two groups bakterioplanton vs Daphnia)
ald.t <- aldex.ttest(ald.clr, groups, paired.test=FALSE)
#ald.kw <- aldex.kw(ald.clr)
#size effect calculation
ald.effect <- aldex.effect(ald.clr, groups, verbose=FALSE, CI=TRUE, include.sample.summary=FALSE)

ald.all <- data.frame(ald.t,ald.effect)
sig_by_both <- which(ald.all$we.ep < 0.05 & ald.all$wi.ep < 0.05) #42
sig_by_both_fdr <- which(ald.all$we.eBH < 0.05 & ald.all$wi.eBH < 0.05) #21

ald.all.significant <- ald.all[as.array(sig_by_both_fdr),]

#Aldex graphs
par(mfrow=c(1,2))
aldex.plot(ald.all, type="MA", test="welch")
aldex.plot(ald.all, type="MW", test="welch")

#Replace ASVs with Genus name
taxa_table$combined <- paste(taxa_table$Family,taxa_table$Genus,taxa_table$Species, sep=" " )
ald.all.significant["taxa_names"] <- taxa_table$combined[match(rownames(ald.all.significant), rownames(taxa_table))]
ald.all.significant$taxa_names <- gsub("NA", "Unknown", ald.all.significant$taxa_names)
#Delete taxa that only has unknowns (it is not identified)
ald.all.significant <- ald.all.significant[!grepl("Unknown Unknown Unknown", ald.all.significant$taxa_names),] 
#Delete all unknowns
ald.all.significant$taxa_names <- gsub("Unknown", " ", ald.all.significant$taxa_names)


#Significant taxa histogram
p<-ggplot(data=ald.all.significant, aes(x=taxa_names, y=effect)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  xlab("Median effect size")+ # for the x axis label
  ylab("Taxa")+
  theme(axis.text = element_text(size = 7))# for the y axis label
p

# Horizontal bar plot
p + coord_flip()

######## Low pollution vs High Pollution - Daphnia
#Define groups for ALDEX2
meta$Joint <- paste(meta$experiment, meta$Type)

groups_of_interest <- c("Daphnia Nature reserve", "Daphnia Polluted Pond")

#Retain only groups of interest
meta_Daphnia <- meta %>%
  filter(Joint %in% groups_of_interest)

#Retain only groups of interest
abund_table <- t(abund_table)
abun_Daphnia <- abund_table %>%
  subset(rownames(abund_table) %in% rownames(meta_Daphnia))

#match rows between otus and meta tables
abun_Daphnia <- abun_Daphnia[match(rownames(meta_Daphnia), rownames(abun_Daphnia)), ]

#create condition vector
groups <- as.factor(meta_Daphnia$Joint)
#CLR tranformation
ald.clr <- aldex.clr(t(abun_Daphnia), groups, mc.samples=128, verbose=TRUE)
#T-test (two groups bakterioplanton vs Daphnia)
ald.t <- aldex.ttest(ald.clr, groups, paired.test=FALSE)
#ald.kw <- aldex.kw(ald.clr)
#size effect calculation
ald.effect <- aldex.effect(ald.clr, groups, verbose=FALSE, CI=TRUE, include.sample.summary=FALSE)

ald.all <- data.frame(ald.t,ald.effect)
sig_by_both <- which(ald.all$we.ep < 0.05 & ald.all$wi.ep < 0.05) #4 (before p-value adjustment)
sig_by_both_fdr <- which(ald.all$we.eBH < 0.05 & ald.all$wi.eBH < 0.05) #0 (after adjustment)

ald.all.significant <- ald.all[as.array(sig_by_both_fdr),]

#Aldex graphs
par(mfrow=c(1,2))
aldex.plot(ald.all, type="MA", test="welch")
aldex.plot(ald.all, type="MW", test="welch")

#Replace ASVs with Genus name
taxa_table$combined <- paste(taxa_table$Family,taxa_table$Genus,taxa_table$Species, sep=" " )
ald.all.significant["taxa_names"] <- taxa_table$combined[match(rownames(ald.all.significant), rownames(taxa_table))]
ald.all.significant$taxa_names <- gsub("NA", "Unknown", ald.all.significant$taxa_names)
#Delete taxa that only has unknowns (it is not identified)
ald.all.significant <- ald.all.significant[!grepl("Unknown Unknown Unknown", ald.all.significant$taxa_names),] 
#Delete all unknowns
ald.all.significant$taxa_names <- gsub("Unknown", " ", ald.all.significant$taxa_names)


#Significant taxa histogram
p<-ggplot(data=ald.all.significant, aes(x=taxa_names, y=effect)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  xlab("Median effect size")+ # for the x axis label
  ylab("Taxa")+
  theme(axis.text = element_text(size = 7))# for the y axis label
p

# Horizontal bar plot
p + coord_flip()

######## Low pollution vs High Pollution - Bakterioplankton

#Define groups for ALDEX2
groups_of_interest <- c("Bacterioplankton Nature reserve", "Bacterioplankton WWTP ")

#Retain only groups of interest
meta_Bact <- meta %>%
  filter(Joint %in% groups_of_interest)

#Retain only groups of interest
abun_Bact <- abund_table %>%
  subset(rownames(abund_table) %in% rownames(meta_Bact))

#create condition vector
groups <- as.factor(meta_Bact$Joint)
#CLR tranformation
ald.clr <- aldex.clr(t(abun_Bact), groups, mc.samples=128, verbose=TRUE)
#T-test (two groups bakterioplanton vs Daphnia)
ald.t <- aldex.ttest(ald.clr, groups, paired.test=FALSE)
#ald.kw <- aldex.kw(ald.clr)
#size effect calculation
ald.effect <- aldex.effect(ald.clr, groups, verbose=FALSE, CI=TRUE, include.sample.summary=FALSE)

ald.all <- data.frame(ald.t,ald.effect)
sig_by_both <- which(ald.all$we.ep < 0.05 & ald.all$wi.ep < 0.05) #0 (before p-value adjustment)
sig_by_both_fdr <- which(ald.all$we.eBH < 0.05 & ald.all$wi.eBH < 0.05) #0 (after adjustment)

ald.all.significant <- ald.all[as.array(sig_by_both_fdr),]

#Aldex graphs
par(mfrow=c(1,2))
aldex.plot(ald.all, type="MA", test="welch")
aldex.plot(ald.all, type="MW", test="welch")












