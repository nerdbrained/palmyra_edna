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

## create dataset merging sample types

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

## glom OTUs (using pairwise distance of 0.03 as threshold)
palmyra_h03_tipglom=tip_glom(palmyra_sampmerge,h=0.03)
## 2612 OTUs
## OTUs 
palmyra_h03_tipglom_type=tip_glom(palmyra_typemerge,h=0.03)



dd = as.dist(cophenetic.phylo(phy_tree(palmyra_sampmerge)))
psclust = cutree(as.hclust(agnes(dd)), h = 0.03)
taxGlomKey <- data.frame(key = psclust) %>% rownames_to_column
taxGlomKey$spid=data.frame((tax_table(palmyra_sampmerge)))[,7]
nGrp = max(taxGlomKey$key) 
h03_grps = vector("list", nGrp)
for (loc_key in 1:nGrp){
      h03_grps[[loc_key]] = taxGlomKey %>% filter(key == loc_key) %>% pull(spid)
}

h03_nona = lapply(h03_grps, function(x) x[!is.na(x)])
h03_id = h03_nona[lapply(h03_nona,length)>0]
length(h03_id)
## 1164 OTUs with some sort of ID. Some identified as same sp split!


### Colors for bar-plottin'
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

# Plotting by relative frequency

plotto = transform_sample_counts(comparetop10ords_relative, function(OTU) OTU / sum(OTU))

plot_bar(physeq = plotto, fill = "Order") + scale_fill_manual(values=c25)

#### Rarefaction curves (ASV-level) by sample
options(scipen=10000)
p<-ggrare(palmyra_sampmerge, step = 1000, color = "SampleType", label = "Site", se = FALSE) + ylab("ASVs")
plot(p)
#### Rarefaction curves (OTU-level) by sample
p <- ggrare(palmyra_h03_tipglom, step = 1000, color = "SampleType", label = "Site", se = FALSE) + ylab("OTUs")
plot(p)
#### Rarefaction curves (order-level) by sample
p <- ggrare(palmyra_sampmerge_order, step = 1000, color = "SampleType", label = "Site", se = FALSE) + ylab("Orders")
plot(p)


### Faith's Phylogenetic Diversity and observed species richness (OTUs)
 
OTU <- as(otu_table(palmyra_h03_tipglom), "matrix")
tree=phy_tree(palmyra_h03_tipglom)

pd_sample=pd(OTU,tree)

write.csv(pd_sample,file="~/Documents/GitHub/palmyra_edna/results/palmyra_pd_samples.csv")


OTU <- as(otu_table(palmyra_h03_tipglom_type), "matrix")
tree=phy_tree(palmyra_h03_tipglom_type)

pd_type=pd(OTU,tree)

write.csv(pd_type,file="~/Documents/GitHub/palmyra_edna/results/palmyra_pd_types.csv")


### Other diversity metrics for different taxonomic/classification schemes

div_samp_phy=estimate_richness(palmyra_sampmerge_phylum,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_samp_phy,file="~/Documents/GitHub/palmyra_edna/results/div_samp_phy.csv")
div_samp_ord=estimate_richness(palmyra_sampmerge_order,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_samp_ord,file="~/Documents/GitHub/palmyra_edna/results/tables/div_samp_ord.csv")
div_samp_otu=estimate_richness(palmyra_h03_tipglom,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_samp_otu,file="~/Documents/GitHub/palmyra_edna/results/div_samp_otu.csv")


div_type_phy=estimate_richness(palmyra_typemerge_phylum,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_type_phy,file="~/Documents/GitHub/palmyra_edna/results/div_type_phy.csv")
div_type_ord=estimate_richness(palmyra_typemerge_order,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_type_ord,file="~/Documents/GitHub/palmyra_edna/results/div_type_ord.csv")
div_type_otu=estimate_richness(palmyra_h03_tipglom_type,measures=c("Observed","Chao","Shannon","ACE"))
write.csv(div_type_otu,file="~/Documents/GitHub/palmyra_edna/results/div_type_otu.csv")


### do Venn on OTUs + orders?
### find out what orders proportion of orders found in different combos 

otudf=data.frame(otu_table(palmyra_typemerge_order))

sed1tax=colnames(otudf[which(otudf[1,]>0)])

sed2tax=colnames(otudf[which(otudf[2,]>0)])

sesstax=colnames(otudf[which(otudf[3,]>0)])

watertax=colnames(otudf[which(otudf[4,]>0)])

ordvenn<-plotVenn(list(Sediment1=sed1tax,Water=watertax,Sessile=sesstax,Sediment2=sed2tax),outFile="~/Desktop/ordervenn.svg",setColors=c("purple","blue","green","gold"), borderWidth=3, opacity=0.2)



otudf=data.frame(otu_table(palmyra_h03_tipglom_type))

sed1tax=colnames(otudf[which(otudf[1,]>0)])

sed2tax=colnames(otudf[which(otudf[2,]>0)])

sesstax=colnames(otudf[which(otudf[3,]>0)])

watertax=colnames(otudf[which(otudf[4,]>0)])

otuvenn<-plotVenn(list(Sediment1=sed1tax,Water=watertax,Sessile=sesstax,Sediment2=sed2tax),outFile="~/Desktop/otuvenn.svg",setColors=c("purple","blue","green","gold"), borderWidth=3, opacity=0.2)





### COMMUNITAYS - ordinate at OTU level
rare_OTU<-rarefy_even_depth(physeq = palmyra_h03_tipglom)
# Calculate distance matrix
palmyra_bray_OTU <- phyloseq::distance(rare_OTU, method = "bray")
# make a data frame from the sample_data
sampledf_OTU <- data.frame(sample_data(rare_OTU))

PCOA_OTU <- ordinate(
  physeq = rare_OTU, 
  method = "PCoA", 
  distance = "bray"
)


sampledf_OTU$PC1=PCOA_OTU$vectors[,1]
sampledf_OTU$PC2=PCOA_OTU$vectors[,2]

plot(sampledf_OTU$PC1,sampledf_OTU$PC2,col=as.integer(factor(sampledf_OTU$SampleType)),xlab="PC1",ylab="PC2",pch=19,main="OTU, Bray-Curtis")

legend(x=-0.2,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c(1,2,3,4),pch=19)



PCOA_OTU <- ordinate(
  physeq = rare_OTU, 
  method = "PCoA", 
  distance = "jaccard"
)


sampledf_OTU$PC1=PCOA_OTU$vectors[,1]
sampledf_OTU$PC2=PCOA_OTU$vectors[,2]

plot(sampledf_OTU$PC1,sampledf_OTU$PC2,col=as.integer(factor(sampledf_OTU$SampleType)),xlab="PC1",ylab="PC2",pch=19,main="OTU, Jaccard")

legend(x=-0.2,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c(1,2,3,4),pch=19)


#### Order level

rare_ord<-rarefy_even_depth(physeq = palmyra_sampmerge_order)
# Calculate distance matrix
palmyra_bray_ord <- phyloseq::distance(rare_ord, method = "bray")
# make a data frame from the sample_data
sampledf_ord <- data.frame(sample_data(rare_ord))

PCOA_ord <- ordinate(
  physeq = rare_ord, 
  method = "PCoA", 
  distance = "bray"
)


sampledf_ord$PC1=PCOA_ord$vectors[,1]
sampledf_ord$PC2=PCOA_ord$vectors[,2]

plot(sampledf_ord$PC1,sampledf_ord$PC2,col=as.integer(factor(sampledf_ord$SampleType)),xlab="PC1",ylab="PC2",pch=19,main="Order, Bray-Curtis")

legend(x=-0.2,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c(1,2,3,4),pch=1)

PCOA_ord <- ordinate(
  physeq = rare_ord, 
  method = "PCoA", 
  distance = "jaccard"
)


sampledf_ord$PC1=PCOA_ord$vectors[,1]
sampledf_ord$PC2=PCOA_ord$vectors[,2]

plot(sampledf_ord$PC1,sampledf_ord$PC2,col=as.integer(factor(sampledf_ord$SampleType)),xlab="PC1",ylab="PC2",pch=19,main="Order, Jaccard")

legend(x=-0.2,y=0.4,legend=c("Sediment (500 um - 2 mm)","Sediment (100 um - 500 um)","Sessile","Water"), col=c(1,2,3,4),pch=19)



#### PERMANOVA

palmyra_jaccard_ord <- phyloseq::distance(rare_ord, method = "jaccard")

palmyra_jaccard_OTU <- phyloseq::distance(rare_OTU, method = "jaccard")


#TYPE - OTU
perm_type_otu<-adonis(palmyra_bray_OTU ~ SampleType, data = sampledf_OTU)

#TYPE - Order
perm_type_order<-adonis(palmyra_bray_ord ~ SampleType, data = sampledf_ord)

#TYPE - OTU
perm_type_jotu<-adonis(palmyra_jaccard_OTU ~ SampleType, data = sampledf_OTU)

#TYPE - Order
perm_type_jorder<-adonis(palmyra_jaccard_ord ~ SampleType, data = sampledf_ord)


