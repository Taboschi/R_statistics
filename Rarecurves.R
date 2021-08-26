#### Rarecurves ####

library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2") # das Package gibt es fuer meine R version nicht 

setwd("C:/Users/Tabby/Desktop/Master/R Daten vorläufig/")

EMB <- readRDS("EMB238_balanced_neu.rds")
EMB

samples_to_remove<-c("M1", "M2", "B1", "B2") # remove the mock samples

x=subset_samples(EMB, !(Name %in% samples_to_remove))

EMB238_nomockctrl=prune_taxa(taxa_sums(x)>0, x) #ASVs withouth input deleted

# remove Chloroplasts, Eukaryota, Mitochondria -> all not of interest
EMB238_nochloroplast <- subset_taxa(EMB238_nomockctrl, !(Order %in% c("Chloroplast"))) 
EMB238_noeu <- subset_taxa(EMB238_nochloroplast, !(Kingdom %in% c("Eukaryota")))
EMB238_nomito <- subset_taxa(EMB238_noeu, !(Family %in% c("Mitochondria")))
EMB238_clean <- subset_taxa(EMB238_nomito, !is.na(Kingdom))

?subset_samples
# for subsetting it is important to connect the arguments with & and use == instead of =
y<- subset_samples(EMB238_clean, St_EMB238 == "8_5"& Core=="3")
y

#rarecurve((otu_table(y)), step=50, cex=0.5) # einfacher ohne t kein transponieren 

otu_table(y)
tax_table(y)

rarecurve(t(sam_data(y)), step=50, cex=0.5) #  the t() transpones the data
# with dim() you can watch if rows and zeilen are correct or if you need to transpone
view(y)

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(y), "matrix")
# transpose if necessary
if(taxa_are_rows(y)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

rarecurve(OTUdf)

# Print the metadata using the phyloseq function
sample_data(EMB)

rarecurve(t(otu_table(EMB)), step=50, cex=0.5) 
?rarecurve


#### Stationen und Kerne Example ####

# 10_4
St_10_4_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "10_4"& Core=="1")
rarecurve((otu_table(St_10_4_C_1)), step=50, cex=0.5)

St_10_4_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "10_4"& Core=="2")
rarecurve((otu_table(St_10_4_C_2)), step=50, cex=0.5)

St_10_4_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "10_4"& Core=="3")
rarecurve((otu_table(St_10_4_C_3)), step=50, cex=0.5)



