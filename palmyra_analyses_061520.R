#library(BiocManager)
#BiocManager::install("phyloseq")
#BiocManager::install("decontam")
#library(devtools)
#devtools::install_github("jbisanz/qiime2R")

library(phyloseq)
library(ggplot2)
library(decontam)
library(qiime2R)
library(vegan)
library(psadd)
library(ape)
library(cluster)
library(tidyverse)
library(ranacapa)
library(picante)
library(nVennR)


##qiime2R is needed to import qiime2 artifacts directly into R. all of these should be in the qiime2 output, just need to specify file paths

palmyra_physeq<-qza_to_phyloseq(features="~/Documents/GitHub/palmyra_edna/raw-data/table.qza",tree="~/Documents/GitHub/palmyra_edna/raw-data/rooted-tree.qza",taxonomy="~/Documents/GitHub/palmyra_edna/raw-data/16S-taxonomy-silva138-V1V3.qza",metadata="~/Documents/GitHub/palmyra_edna/raw-data/metadata_qubit_sampinfo.tsv")

#remove 17B, 17A, 19B for sure
#also removing 17B115 and 17B515 for now bc lack qubit
#potentially remove in the future: 19A, 19A1, 19A5, others with no site info!

palmyra_prune = subset_samples(palmyra_physeq, sample_names(palmyra_physeq) != "17B")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "17A")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19B")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "17B115")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "17B515")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19A")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19A1")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19A1b")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19A5")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19A5b")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19Ab")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "19B")
palmyra_prune = subset_samples(palmyra_prune, sample_names(palmyra_prune) != "PAL19")

## number of reads per sample
sample_sums(palmyra_prune)

## after pruning, remove taxa with no occurrences
palmyra_prune = subset_samples_no_zero(palmyra_prune)

#identify contaminants using frequency method

contamdf.freq <- isContaminant(palmyra_prune, method="frequency", conc="quant_reading")
head(contamdf.freq)
head(which(contamdf.freq$contaminant))
which(contamdf.freq$contaminant)

#compare prevalence of contaminants to non-contaminants

mean(contamdf.freq$prev[which(contamdf.freq$contaminant == FALSE)])
mean(contamdf.freq$prev[which(contamdf.freq$contaminant)])

ks.test(contamdf.freq$prev[which(contamdf.freq$contaminant == FALSE)],contamdf.freq$prev[which(contamdf.freq$contaminant)])

## remove contaminants!

palmyra_noncontam <- prune_taxa(!contamdf.freq$contaminant, palmyra_prune)
palmyra_noncontam

##############
###revising###
##############

####bray####

################# Ordinations ###################
#Ordination
PCOA <- ordinate(
  physeq = palmyra_noncontam, 
  method = "PCoA", 
  distance = "bray"
)

#PERMANOVA to test if sites differ
rare<-rarefy_even_depth(physeq = palmyra_noncontam)
# Calculate Bray distance matrix
palmyra_bray <- phyloseq::distance(rare, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(rare))

sampledf$PC1=PCOA$vectors[,1]
sampledf$PC2=PCOA$vectors[,2]

plot(sampledf$PC1,sampledf$PC2,col=as.integer(factor(sampledf$Sample_type)))

legend(x=-0.1,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c(1,2,3,4),pch=1)

#by substrate
substrate<-adonis(palmyra_bray ~ Sample_type, data = sampledf)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Sample_type  3    4.9039 1.63463  6.7018 0.48912  0.001 ***
#Residuals   21    5.1221 0.24391         0.51088           
#Total       24   10.0260                 1.00000 

#By site, just water samples
water<-subset_samples(rare, Sample_type=="Water")
water_bray <- phyloseq::distance(water, method = "bray")
water.df<- data.frame(sample_data(water))
adonis(water_bray ~ Site, data = water.df)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Site       1   0.27698 0.27698  1.0295 0.17075  0.389
#Residuals  5   1.34519 0.26904         0.82925       
#Total      6   1.62217                 1.00000   

#Just sessile samples
sessile<-subset_samples(rare, Sample_type=="Sessile")
sessile_bray <- phyloseq::distance(sessile, method = "bray")
sessile.df<- data.frame(sample_data(sessile))
adonis(sessile_bray ~ Site, data = sessile.df)
#Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
#Site       1   0.38478 0.38478  5.2982 0.56981 0.01111 *
#Residuals  4   0.29050 0.07262         0.43019          
#Total      5   0.67528                 1.00000    

#Sediment 1
sed1<-subset_samples(rare, Sample_type=="Sediment (100 - 500 um)")
sed1_bray <- phyloseq::distance(sed1, method = "bray")
sed1.df<- data.frame(sample_data(sed1))
adonis(sed1_bray ~ Site, data = sed1.df)
#Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
#Site       1   0.61410 0.61410  4.6845 0.53941 0.07778 .
#Residuals  4   0.52437 0.13109         0.46059          
#Total      5   1.13847                 1.00000    

#Sediment 2
sed2<-subset_samples(rare, Sample_type=="Sediment (500 um - 2 mm)")
sed2_bray <- phyloseq::distance(sed2, method = "bray")
sed2.df<- data.frame(sample_data(sed2))
adonis(sed2_bray ~ Site, data = sed2.df)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Site       1   0.78643 0.78643  3.4681 0.46439 0.1556
#Residuals  4   0.90704 0.22676         0.53561       
#Total      5   1.69347                 1.00000     

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3861 taxa and 25 samples ]
# sample_data() Sample Data:       [ 25 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 3861 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3861 tips and 3856 internal nodes ]


###Jaccard

PCOA <- ordinate(
  physeq = palmyra_noncontam, 
  method = "PCoA", 
  distance = "jaccard"
)


# Calculate Jaccard distance matrix
palmyra_jaccard <- phyloseq::distance(rare, method = "jaccard")

sampledf$PC1=PCOA$vectors[,1]
sampledf$PC2=PCOA$vectors[,2]

plot(sampledf$PC1,sampledf$PC2,col=as.integer(factor(sampledf$Sample_type)))

legend(x=-0.1,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c(1,2,3,4),pch=1)

#By substrate
substrate<-adonis(palmyra_jaccard ~ Sample_type, data = sampledf)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Sample_type  3    4.1392 1.37972    4.43 0.38758  0.001 ***
#Residuals   21    6.5404 0.31145         0.61242           
#Total       24   10.6795                 1.00000 

#Just water samples
water<-subset_samples(rare, Sample_type=="Water")
water_jaccard <- phyloseq::distance(water, method = "jaccard")
water.df<- data.frame(sample_data(water))
adonis(water_jaccard ~ Site, data = water.df)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Site       1    0.3572 0.35720  1.0318 0.17106  0.381
#Residuals  5    1.7309 0.34618         0.82894       
#Total      6    2.0881                 1.00000   

#Just sediment samples
sessile<-subset_samples(rare, Sample_type=="Sessile")
sessile_jaccard <- phyloseq::distance(sessile, method = "jaccard")
sessile.df<- data.frame(sample_data(sessile))
adonis(sessile_jaccard ~ Site, data = sessile.df)
#Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
#Site       1   0.56953 0.56953  4.3721 0.52222 0.01111 *
#Residuals  4   0.52106 0.13026         0.47778          
#Total      5   1.09059                 1.00000    

sed1<-subset_samples(rare, Sample_type=="Sediment (100 - 500 um)")
sed1_jaccard <- phyloseq::distance(sed1, method = "jaccard")
sed1.df<- data.frame(sample_data(sed1))
adonis(sed1_jaccard ~ Site, data = sed1.df)
#Df SumsOfSqs MeanSqs F.Model     R2  Pr(>F)  
#Site       1   0.73608 0.73608  3.7654 0.4849 0.06667 .
#Residuals  4   0.78194 0.19548         0.5151          
#Total      5   1.51801                 1.0000     

sed2<-subset_samples(rare, Sample_type=="Sediment (500 um - 2 mm)")
sed2_jaccard <- phyloseq::distance(sed2, method = "jaccard")
sed2.df<- data.frame(sample_data(sed2))
adonis(sed2_jaccard ~ Site, data = sed2.df)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Site       1   0.87101 0.87101  3.5349 0.46914 0.1556
#Residuals  4   0.98560 0.24640         0.53086       
#Total      5   1.85661                 1.00000     

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3861 taxa and 25 samples ]
# sample_data() Sample Data:       [ 25 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 3861 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3861 tips and 3856 internal nodes ]




## create dataset merging replicate samples

palmyra_sampmerge<-merge_samples(palmyra_noncontam,as.data.frame(sample_data(palmyra_noncontam))$Sample)

df <- as.data.frame(sample_data(palmyra_sampmerge)) # Put sample_data into a ggplot-friendly data.frame
df$SampleType=sub('..','',rownames(df))
sample_data(palmyra_sampmerge)<-df

## create dataset merging samples

palmyra_typemerge<-merge_samples(palmyra_noncontam,as.data.frame(sample_data(palmyra_noncontam))$Sample_type)

## count # of phyla

phys=unique(data.frame(tax_table(palmyra_sampmerge)[,2])$Phylum)
length(phys)
#70 "phyla" - 3 = 67 phyla
length(which(data.frame(tax_table(palmyra_sampmerge)[,2])$Phylum == "uncultured"))
#14 ASVs identified as "uncultured" at phylum level
length(which(is.na(data.frame(tax_table(palmyra_sampmerge)[,2])$Phylum)==TRUE))
#764 ASVs identified as unknown
length(which(data.frame(tax_table(palmyra_sampmerge)[,2])$Phylum == "Incertae_Sedis"))
#24 ASVs identified as uncertain placement

ords=unique(data.frame(tax_table(palmyra_sampmerge)[,4])$Order)
length(ords)
#234-3 (uncultured/NA/incertae)
length(which(data.frame(tax_table(palmyra_sampmerge)[,4])$Order == "uncultured"))
#38
length(which(is.na(data.frame(tax_table(palmyra_sampmerge)[,4])$Order)==TRUE))
#1146
length(which(data.frame(tax_table(palmyra_sampmerge)[,4])$Order == "Incertae_Sedis"))
#46

species=unique(data.frame(tax_table(palmyra_sampmerge)[,7])$Species)
length(species)
#411
length(which(is.na(data.frame(tax_table(palmyra_sampmerge)[,7])$Species)==TRUE))
#2149


### Colors for bar-plotting
library(sp)
colours<- bpy.colors(n = 10, cutoff.tails = 0.2, alpha = 0.75)
c25 <- c(
  "#E31A1C", # red
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)




#### Plotting 10 most common orders(absolute frequency) at sample level
palmyra_sampmerge_order = tax_glom(palmyra_sampmerge, "Order")

relative_ord = transform_sample_counts(palmyra_sampmerge_order, function(OTU) OTU / sum(OTU))

orderfreqsum=data.frame(tax_table(relative_ord))
orderfreqsum$freqs = taxa_sums(relative_ord)
orderfreqsort=orderfreqsum[order(orderfreqsum$freqs,decreasing=TRUE),]

head(orderfreqsort,10)

write.csv(orderfreqsort,file="~/Documents/GitHub/palmyra_edna/results/palmyra_ordertable_samples.csv")

Top10ords_relative = names(sort(taxa_sums(relative_ord), TRUE)[1:10])
comparetop10ords_relative = prune_taxa(Top10ords_relative, relative_ord)

plot_bar(physeq = comparetop10ords_relative, fill = "Order") + scale_fill_manual(values=c25)

ggsave("~/Documents/GitHub/palmyra_edna/figures/barplot_orders.pdf",device="pdf")

# Plotting by relative frequency

plotto = transform_sample_counts(comparetop10ords_relative, function(OTU) OTU / sum(OTU))

plot_bar(physeq = plotto, fill = "Order") + scale_fill_manual(values=c25)

ggsave("~/Documents/GitHub/palmyra_edna/figures/barplot_orders_relative.pdf",device="pdf")

#### Plotting 10 most common orders(absolute frequency) at type level
palmyra_typemerge_order = tax_glom(palmyra_typemerge, "Order")

relative_ord = transform_sample_counts(palmyra_typemerge_order, function(OTU) OTU / sum(OTU))

orderfreqsum=data.frame(tax_table(relative_ord))
orderfreqsum$freqs = taxa_sums(relative_ord)
orderfreqsort=orderfreqsum[order(orderfreqsum$freqs,decreasing=TRUE),]

head(orderfreqsort,10)

write.csv(orderfreqsort,file="~/Documents/GitHub/palmyra_edna/results/palmyra_ordertable_types.csv")

Top10ords_relative = names(sort(taxa_sums(relative_ord), TRUE)[1:10])
comparetop10ords_relative = prune_taxa(Top10ords_relative, relative_ord)

plot_bar(physeq = comparetop10ords_relative, fill = "Order") + scale_fill_manual(values=c25)

ggsave("~/Documents/GitHub/palmyra_edna/figures/barplot_orders_type.pdf",device="pdf")

# Plotting by relative frequency

plotto = transform_sample_counts(comparetop10ords_relative, function(OTU) OTU / sum(OTU))

plot_bar(physeq = plotto, fill = "Order") + scale_fill_manual(values=c25) + xlab("Substrate")

ggsave("~/Documents/GitHub/palmyra_edna/figures/barplot_orders_type_relative.pdf",device="pdf")



#### Plotting 10 most common orders(absolute frequency) at type level
palmyra_typemerge_family = tax_glom(palmyra_typemerge, "Family")

relative_fam = transform_sample_counts(palmyra_typemerge_family, function(OTU) OTU / sum(OTU))

famfreqsum=data.frame(tax_table(relative_fam))
famfreqsum$freqs = taxa_sums(relative_fam)
famfreqsort=famfreqsum[order(famfreqsum$freqs,decreasing=TRUE),]

head(famfreqsort,10)

write.csv(famfreqsort,file="~/Documents/GitHub/palmyra_edna/results/palmyra_famtable_types.csv")

Top10fams_relative = names(sort(taxa_sums(relative_fam), TRUE)[1:10])
comparetop10fams_relative = prune_taxa(Top10fams_relative, relative_fam)

plot_bar(physeq = comparetop10fams_relative, fill = "Family") + scale_fill_manual(values=c25)

ggsave("~/Documents/GitHub/palmyra_edna/figures/barplot_orders_type.pdf",device="pdf")

# Plotting by relative frequency

plotto = transform_sample_counts(comparetop10ords_relative, function(OTU) OTU / sum(OTU))

plot_bar(physeq = plotto, fill = "Order") + scale_fill_manual(values=c25)

ggsave("~/Documents/GitHub/palmyra_edna/figures/barplot_orders_type_relative.pdf",device="pdf")



#### Rarefaction curves (ASV-level) by sample
options(scipen=10000)
p<-ggrare(palmyra_sampmerge, step = 1000, color = "SampleType", label = "Site", se = FALSE) + ylab("ASVs")
plot(p)

#### Rarefaction curves (order-level) by sample
p <- ggrare(palmyra_sampmerge_order, step = 1000, color = "SampleType", label = "Site", se = FALSE) + ylab("Orders")
plot(p)


### Faith's Phylogenetic Diversity and observed species richness (OTUs)

ASV <- as(otu_table(palmyra_typemerge), "matrix")
tree=phy_tree(palmyra_typemerge)
pd_type=pd(ASV,tree)
write.csv(pd_type,file="~/Documents/GitHub/palmyra_edna/results/palmyra_pd_asv_types.csv")

ASV <- as(otu_table(palmyra_sampmerge), "matrix")
tree=phy_tree(palmyra_sampmerge)
pd_samp=pd(ASV,tree)
write.csv(pd_type,file="~/Documents/GitHub/palmyra_edna/results/palmyra_pd_asv_types.csv")


### Other diversity metrics for different taxonomic/classification schemes

div_samp_phy=estimate_richness(palmyra_sampmerge_phylum,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_samp_phy,file="~/Documents/GitHub/palmyra_edna/results/div_samp_phy.csv")
div_samp_ord=estimate_richness(palmyra_sampmerge_order,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_samp_ord,file="~/Documents/GitHub/palmyra_edna/results/tables/div_samp_ord.csv")
div_samp_asv=estimate_richness(palmyra_sampmerge,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_samp_asv,file="~/Documents/GitHub/palmyra_edna/results/div_samp_asv.csv")

div_type_phy=estimate_richness(palmyra_typemerge_phylum,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_type_phy,file="~/Documents/GitHub/palmyra_edna/results/div_type_phy.csv")
div_type_ord=estimate_richness(palmyra_typemerge_order,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_type_ord,file="~/Documents/GitHub/palmyra_edna/results/div_type_ord.csv")
div_type_asv=estimate_richness(palmyra_typemerge,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_type_asv,file="~/Documents/GitHub/palmyra_edna/results/div_type_asv.csv")


### do Venn on OTUs + orders
### find out what orders proportion of orders found in different combos 

otudf=data.frame(otu_table(palmyra_typemerge_order))

sed1tax=colnames(otudf[which(otudf[1,]>0)])

sed2tax=colnames(otudf[which(otudf[2,]>0)])

sesstax=colnames(otudf[which(otudf[3,]>0)])

watertax=colnames(otudf[which(otudf[4,]>0)])

ordvenn<-plotVenn(list(Sediment1=sed1tax,Water=watertax,Sessile=sesstax,Sediment2=sed2tax),outFile="~/Documents/GitHub/palmyra_edna/figures/ordervenn.svg",setColors=c("purple","blue","green","gold"), borderWidth=3, opacity=0.2)




# Venn diagram for asvs


asvdf=data.frame(otu_table(palmyra_typemerge))

sed1tax=colnames(asvdf[which(asvdf[1,]>0)])

sed2tax=colnames(asvdf[which(asvdf[2,]>0)])

sesstax=colnames(asvdf[which(asvdf[3,]>0)])

watertax=colnames(asvdf[which(asvdf[4,]>0)])

asvvenn<-plotVenn(list("Sediment (100-500 um)"=sed1tax,Water=watertax,Sessile=sesstax,"Sediment (500 um - 2 mm)"=sed2tax),outFile="~/Desktop/asvvenn.svg",setColors=c("purple","blue","green","gold"), borderWidth=3, opacity=0.2)

asvvenn<-plotVenn(list("Sediment (100-500um)"=sed1tax,"Sediment (500um-2mm)"=sed2tax,Sessile=sesstax,Water=watertax),outFile="~/Desktop/asvvenn.svg",setColors=c("purple","gold","green","blue"), borderWidth=2, opacity=0.2)

ggsave("~/Documents/GitHub/palmyra_edna/figures/asv_venn.pdf",device="pdf")




### Plotting PCA
### ASV level

rare_asv<-rarefy_even_depth(physeq = palmyra_sampmerge)
# Calculate distance matrix
palmyra_bray_asv <- phyloseq::distance(rare_asv, method = "bray")
# make a data frame from the sample_data
sampledf_asv <- data.frame(sample_data(rare_asv))

PCOA_asv <- ordinate(
  physeq = rare_asv, 
  method = "PCoA", 
  distance = "bray"
)


sampledf_asv$PC1=PCOA_asv$vectors[,1]
sampledf_asv$PC2=PCOA_asv$vectors[,2]

plot(sampledf_asv$PC1,sampledf_asv$PC2,col=c("purple","gold","green","blue","blue","blue","blue","blue","purple","gold","green","purple","gold","green","blue","blue"),xlab="PC1",ylab="PC2",pch=19,main="Bray-Curtis")

legend(x=0.1,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c("purple","gold","green","blue"),pch=1,box.col="white")



PCOA_asv <- ordinate(
  physeq = rare_asv, 
  method = "PCoA", 
  distance = "jaccard"
)


sampledf_asv$PC1=PCOA_asv$vectors[,1]
sampledf_asv$PC2=PCOA_asv$vectors[,2]

plot(sampledf_asv$PC1,sampledf_asv$PC2,col=c("purple","gold","green","blue","blue","blue","blue","blue","purple","gold","green","purple","gold","green","blue","blue"),xlab="PC1",ylab="PC2",pch=19,main="Jaccard")

legend(x=0.1,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c("purple","gold","green","blue"),pch=1,box.col="white")
