#### Rarecurves ####

library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2") # das Package gibt es fuer meine R version nicht 

setwd("C:/Users/Tabby/Desktop/Master/R Daten vorläufig/")

EMB <- readRDS("EMB238_balanced_neu.rds")
EMB

samples_to_remove<-c("M1", "M2", "B1", "B2")

x=subset_samples(EMB, !(Name %in% samples_to_remove))

EMB238_nomockctrl=prune_taxa(taxa_sums(x)>0, x)

EMB238_nochloroplast <- subset_taxa(EMB238_nomockctrl, !(Order %in% c("Chloroplast")))
EMB238_noeu <- subset_taxa(EMB238_nochloroplast, !(Kingdom %in% c("Eukaryota")))
EMB238_nomito <- subset_taxa(EMB238_noeu, !(Family %in% c("Mitochondria")))
EMB238_clean <- subset_taxa(EMB238_nomito, !is.na(Kingdom))

?subset_samples
y<- subset_samples(EMB238_clean, St_EMB238 == "8_5"& Core=="3")
y

#rarecurve((otu_table(y)), step=50, cex=0.5) # einfacher ohne t kein transponieren 

otu_table(y)
tax_table(y)

rarecurve(t(sam_data(y)), step=50, cex=0.5)
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


#### Stationen und Kerne ####

# 10_4
St_10_4_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "10_4"& Core=="1")
rarecurve((otu_table(St_10_4_C_1)), step=50, cex=0.5)

St_10_4_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "10_4"& Core=="2")
rarecurve((otu_table(St_10_4_C_2)), step=50, cex=0.5)

St_10_4_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "10_4"& Core=="3")
rarecurve((otu_table(St_10_4_C_3)), step=50, cex=0.5)

# 13_6
St_13_6_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "13_6"& Core=="1")
rarecurve((otu_table(St_13_6_C_1)), step=50, cex=0.5)

St_13_6_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "13_6"& Core=="2")
rarecurve((otu_table(St_13_6_C_2)), step=50, cex=0.5)

St_13_6_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "13_6"& Core=="3")
rarecurve((otu_table(St_13_6_C_3)), step=50, cex=0.5)

# 15_5
St_15_5_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "15_5"& Core=="1")
rarecurve((otu_table(St_15_5_C_1)), step=50, cex=0.5)

St_15_5_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "15_5"& Core=="2")
rarecurve((otu_table(St_15_5_C_2)), step=50, cex=0.5)

St_15_5_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "15_5"& Core=="3")
rarecurve((otu_table(St_15_5_C_3)), step=50, cex=0.5)

# 17_6
St_17_6_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "17_6"& Core=="1")
rarecurve((otu_table(St_17_6_C_1)), step=50, cex=0.5)
# empty rows removed
St_17_6_C_1_15min_RNA_Biomics <- subset_samples(St_17_6_C_1, Beating!="6 min"&
                                               Nucleic_Acid!="DNA"&
                                                Beating!="30 sec")

St_17_6_C_1_6min_DNA_Quick <- subset_samples(St_17_6_C_1, Beating!="15 min"&
                                                Beating!="30 sec")

St_17_6_C_1
View(St_17_6_C_1)
rarecurve((otu_table(St_17_6_C_1_6min_RNA_Quick)), step=50, cex=0.5)

St_17_6_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "17_6"& Core=="2")

St_17_6_C_2_15min_RNA_Biomics <- subset_samples(St_17_6_C_2, Beating!="6 min"&
                                                 Nucleic_Acid!="DNA"&
                                                 Beating!="30 sec")
rarecurve((otu_table(St_17_6_C_2_15min_RNA_Biomics)), step=50, cex=0.5)

St_17_6_C_2_6min_DNA_Quick <- subset_samples(St_17_6_C_2, Beating!="15 min"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_17_6_C_2_6min_DNA_Quick)), step=50, cex=0.5)

St_17_6_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "17_6"& Core=="3")

St_17_6_C_3_15min_RNA_Quick <- subset_samples(St_17_6_C_3, Beating!="6 min"&
                                                  Nucleic_Acid!="DNA"&
                                                  Beating!="30 sec")
rarecurve((otu_table(St_17_6_C_3_15min_RNA_Quick)), step=50, cex=0.5)

St_17_6_C_3_6min_DNA_Quick <- subset_samples(St_17_6_C_3, Beating!="15 min"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_17_6_C_3_6min_DNA_Quick)), step=50, cex=0.5)

St_17_6_C_1_15min_DNA_Biomics <- subset_samples(St_17_6_C_1, Beating!="6 min"&
                                                  Nucleic_Acid!="RNA"&
                                                  Beating!="30 sec")
rarecurve((otu_table(St_17_6_C_1_15min_DNA_Biomics)), step=50, cex=0.5)

# 18_6
St_18_6_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "18_6"& Core=="1")
rarecurve((otu_table(St_18_6_C_1)), step=50, cex=0.5)
# auch mit beadbeating getrennt

St_18_6_C_1_30sec_RNA_Quick <- subset_samples(St_18_6_C_1, Beating!="6 min"&
                                                Nucleic_Acid!="DNA"&
                                                Beating!="15 min")
rarecurve((otu_table(St_18_6_C_1_30sec_RNA_Quick)), step=50, cex=0.5)

St_18_6_C_1_15min_RNA_Biomics <- subset_samples(St_18_6_C_1, Beating!="6 min"&
                                                Nucleic_Acid!="DNA"&
                                                Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_1_15min_RNA_Biomics)), step=50, cex=0.5)

St_18_6_C_1_15min_DNA_Biomics <- subset_samples(St_18_6_C_1, Beating!="6 min"&
                                                  Nucleic_Acid!="RNA"&
                                                  Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_1_15min_DNA_Biomics)), step=50, cex=0.5)

St_18_6_C_1_6min_DNA_Quick <- subset_samples(St_18_6_C_1, Beating!="15 min"&
                                                  Nucleic_Acid!="RNA"&
                                                  Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_1_6min_DNA_Quick)), step=50, cex=0.5)

St_18_6_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "18_6"& Core=="2")

St_18_6_C_2_15min_RNA_Quick <- subset_samples(St_18_6_C_2, Beating!="6 min"&
                                                  Nucleic_Acid!="DNA"&
                                                  Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_2_15min_RNA_Quick)), step=50, cex=0.5)

St_18_6_C_2_6min_DNA_Quick <- subset_samples(St_18_6_C_2, Beating!="15 min"&
                                               Nucleic_Acid!="RNA"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_2_6min_DNA_Quick)), step=50, cex=0.5)

St_18_6_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "18_6"& Core=="3")

St_18_6_C_3_15min_RNA_Quick <- subset_samples(St_18_6_C_3, Beating!="6 min"&
                                                Nucleic_Acid!="DNA"&
                                                Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_3_15min_RNA_Quick)), step=50, cex=0.5)

St_18_6_C_3_6min_DNA_Quick <- subset_samples(St_18_6_C_3, Beating!="15 min"&
                                               Nucleic_Acid!="RNA"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_18_6_C_3_6min_DNA_Quick)), step=50, cex=0.5)

# 5_5
St_5_5_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "5_5"& Core=="1")
rarecurve((otu_table(St_5_5_C_1)), step=50, cex=0.5)

St_5_5_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "5_5"& Core=="2")
rarecurve((otu_table(St_5_5_C_2)), step=50, cex=0.5)

St_5_5_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "5_5"& Core=="3")
rarecurve((otu_table(St_5_5_C_3)), step=50, cex=0.5)

# 8_5
St_8_5_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "8_5"& Core=="1")
rarecurve((otu_table(St_8_5_C_1)), step=50, cex=0.5)

St_8_5_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "8_5"& Core=="2")
rarecurve((otu_table(St_8_5_C_2)), step=50, cex=0.5)

St_8_5_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "8_5"& Core=="3")
rarecurve((otu_table(St_8_5_C_3)), step=50, cex=0.5)

# 2_4
St_2_4_C_1 <- subset_samples(EMB238_clean, St_EMB238 == "2_4"& Core=="1")

St_2_4_C_1_15min_RNA_Quick <- subset_samples(St_2_4_C_1, Beating!="2 min"&
                                               Kit!="unknown"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_2_4_C_1_15min_RNA_Quick)), step=50, cex=0.5)

St_2_4_C_1_30sec_RNA_Quick <- subset_samples(St_2_4_C_1, Beating!="2 min"&
                                               Kit!="unknown"&
                                               Beating!="15 min")
rarecurve((otu_table(St_2_4_C_1_30sec_RNA_Quick)), step=50, cex=0.5)

St_2_4_C_1_2min_RNA_Quick <- subset_samples(St_2_4_C_1, Beating!="30 sec"&
                                               Kit!="unknown"&
                                               Beating!="15 min")
rarecurve((otu_table(St_2_4_C_1_2min_RNA_Quick)), step=50, cex=0.5)

St_2_4_C_1_unknown <- subset_samples(St_2_4_C_1, Beating!="30 sec"&
                                       Beating!="2 min"&
                                              Beating!="15 min")
rarecurve((otu_table(St_2_4_C_1_unknown)), step=50, cex=0.5)

St_2_4_C_2 <- subset_samples(EMB238_clean, St_EMB238 == "2_4"& Core=="2")

St_2_4_C_2_15min_RNA_Quick <- subset_samples(St_2_4_C_2, Beating!="2 min"&
                                               Kit!="unknown"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_2_4_C_2_15min_RNA_Quick)), step=50, cex=0.5)

St_2_4_C_2_30sec_RNA_Quick <- subset_samples(St_2_4_C_2, Beating!="2 min"&
                                               Kit!="unknown"&
                                               Beating!="15 min")
rarecurve((otu_table(St_2_4_C_2_30sec_RNA_Quick)), step=50, cex=0.5)

St_2_4_C_2_2min_RNA_Quick <- subset_samples(St_2_4_C_2, Beating!="30 sec"&
                                              Kit!="unknown"&
                                              Beating!="15 min")
rarecurve((otu_table(St_2_4_C_2_2min_RNA_Quick)), step=50, cex=0.5)

St_2_4_C_3 <- subset_samples(EMB238_clean, St_EMB238 == "2_4"& Core=="3")

St_2_4_C_3_15min_RNA_Quick <- subset_samples(St_2_4_C_3, Beating!="2 min"&
                                               Kit!="unknown"&
                                               Beating!="30 sec")
rarecurve((otu_table(St_2_4_C_3_15min_RNA_Quick)), step=50, cex=0.5)



