#!/bin/env Rscript

# Calculate population genetic stats

library(PopGenome)
genome.metro_CDS <- readData("./CDS",include.unknown=T) # ./CDS is a directory of fasta alignments of the CDS of each gene

# pop1 = Mn, pop2 = Ms
#SET POPULATION. INdividual names correspond to sequence names in the input fasta files
genome.metro_CDS <- set.populations(genome.metro_CDS, list(c("Mn_01","Mn_02","Mn_03","Mn_04","Mn_05","Mn_06","Mn_07","Mn_08","Mn_09","Mn_10"),c("Ms_01","Ms_02","Ms_03","Ms_04","Ms_05","Ms_06","Ms_07","Ms_08","Ms_09","Ms_10")))
genome.metro_CDS <- set.filter(genome.metro_CDS ,missing.freqs=TRUE,miss.lower.bound=0,miss.upper.bound=0.5)

#GET FST
genome.metro_CDS <- F_ST.stats(genome.metro_CDS , mode = "nucleotide", detail = TRUE)
F_ST.stats.metro_CDS <- genome.metro_CDS@nucleotide.F_ST
colnames(F_ST.stats.metro_CDS) <- "FST"
F_ST.stats.metro_CDS <- as.data.frame(F_ST.stats.metro_CDS)
F_ST.stats.metro_CDS[which(F_ST.stats.metro_CDS$FST < 0),] <- 0

#GET DXY
genome.metro_CDS <- diversity.stats.between(genome.metro_CDS)
dxy.metro_CDS <- as.data.frame(genome.metro_CDS@nuc.diversity.between)
colnames(dxy.metro_CDS) <- "DXY"
dxy.metro_CDS$DXY_persite <- dxy.metro_CDS$DXY / genome.metro_CDS@n.sites

# GET neutrality stats  (for Tajimas D and Watterson's theta)
genome.metro_CDS <- neutrality.stats(genome.metro_CDS)
dtaj.metro_CDS <- as.data.frame(genome.metro_CDS@Tajima.D)
watterson.metro_CDS <- as.data.frame(genome.metro_CDS@theta_Watterson)
colnames(dtaj.metro_CDS) <- c("dtaj_Mn", "dtaj_Ms")
colnames(watterson.metro_CDS) <- c("Wattersons_Theta_Mn","Wattersons_Theta_Ms")

#GET PI
genome.metro_CDS <- diversity.stats(genome.metro_CDS)
Pi.metro_CDS <- as.data.frame(genome.metro_CDS@nuc.diversity.within)
colnames(Pi.metro_CDS) <- c("Pi.Mn","Pi.Ms")
Pi.metro_CDS$Pi_persite_Mn <- Pi.metro_CDS$Pi.Mn / n.sites.metro_CDS[[1]]
Pi.metro_CDS$Pi_persite_Ms <- Pi.metro_CDS$Pi.Ms / n.sites.metro_CDS[[1]]

#GET FIXED, SHARED,MONOMORPHIC SITES, AND OTHER SUMMARY DATA
genome.metro_CDS <- calc.fixed.shared(genome.metro_CDS)
fixed.shared.metro_CDS <- cbind(genome.metro_CDS@n.fixed.sites,genome.metro_CDS@n.shared.sites,genome.metro_CDS@n.monomorphic.sites)
colnames(fixed.shared.metro_CDS) <- c("n.fixed.sites","n.shared.sites","n.monomorphic.sites","n.segregating.sites_Mn","n.segregating.sites_Ms")
genedata.metro_CDS <- get.sum.data(genome.metro_CDS)

# Outlier detection
library(MINOTAUR)

input <- read.csv("CDS_minotaur_input.csv", sep=",", header=T) # "CDS_minotaur_input.csv" is a csv file with various population genetic statistics calculated above for each gene. FST is in column 8, DXY is in column 9, PiMn is in column 14, PiMs is in column 15

### compute rank-based P values
FST_DXY_Pvals <- stat_to_pvalue(dfv=input,column.nums = c(8,9), two.tailed = c(FALSE,FALSE),right.tailed = c(TRUE,TRUE))
FST_PiMn_Pvals <- stat_to_pvalue(dfv=input,column.nums = c(8,14), two.tailed = c(FALSE,FALSE),right.tailed = c(TRUE,FALSE))
FST_PiMs_Pvals <- stat_to_pvalue(dfv=input,column.nums = c(8,15), two.tailed = c(FALSE,FALSE),right.tailed = c(TRUE,FALSE))

### compute Mahalanobis distance
Md_FST_DXY  <- Mahalanobis(FST_DXY_Pvals,M=c(1,1))
Md_FST_PiMn  <- Mahalanobis(FST_PiMn_Pvals,M=c(1,1))
Md_FST_PiMs  <- Mahalanobis(FST_PiMs_Pvals,M=c(1,1))

### Get 95th and 99th quantiles for each dataset
Q95_FST_DXY <- quantile(Md_FST_DXY , 0.95)
Q95_FST_PiMn <- quantile(Md_FST_PiMn , 0.95)
Q95_FST_PiMs <- quantile(Md_FST_PiMs , 0.95)

Q99_FST_DXY <- quantile(Md_FST_DXY , 0.99)
Q99_FST_PiMn <- quantile(Md_FST_PiMn , 0.99)
Q99_FST_PiMs <- quantile(Md_FST_PiMs , 0.99)


### Supplementary plots

# FST vs DXY dotplot
FST_DXY_95 <- subset(input, Md_FST_DXY > Q95_FST_DXY)
FST_DXY_99 <- subset(input, Md_FST_DXY > Q99_FST_DXY)

plot(input$FST,input$DXY,col=rgb(0,0,0,alpha=10,maxColorValue = 100),pch=21,bg=rgb(0,0,0,alpha=10,maxColorValue = 100),ylab="dXY",xlab="FST")
points(FST_DXY_95$FST,FST_DXY_95$DXY,col=rgb(100,55,0,alpha=50,maxColorValue = 100),pch=21,bg=rgb(100,55,0,alpha=50,maxColorValue = 100))
points(FST_DXY_99$FST,FST_DXY_99$DXY,col=rgb(100,0,0,alpha=50,maxColorValue = 100),pch=21,bg=rgb(100,0,0,alpha=50,maxColorValue = 100))
legend("topleft", legend = c("non-outlier","outlier","strong outlier"),col=c(rgb(0,0,0,alpha=10,maxColorValue = 100),rgb(100,55,0,alpha=50,maxColorValue = 100),rgb(100,0,0,alpha=50,maxColorValue = 100)),pt.bg=c(rgb(0,0,0,alpha=10,maxColorValue = 100),rgb(100,55,0,alpha=50,maxColorValue = 100),rgb(100,0,0,alpha=50,maxColorValue = 100)),pch=21,cex=1.5,pt.cex=1)


# Pi Density plot
par(mfrow=c(1,1))
den_PiMn <- density(input$Pi.Mn)
den_PiMs <- density(input$Pi.Ms)

plot(den_PiMn,type="l",xlab="Pi",main="",xlim=c(0,0.002),cex.lab=2,cex.axis=1.5,ylim=c(0,3100),col="lightcoral",lwd=3)
lines(den_PiMs,col="cornflowerblue",lwd=3)
legend("topright",lty=1, lwd=3,col=c("lightcoral","cornflowerblue"),legend=c("M. nervulosa", "M. sclerocarpa"),text.font=3,cex=2)

### GO analysis, (example for FST,DXY outliers, enrichment for all other sets of outliers was conducted in the same way)

library("topGO")

# ME_TAIR_geneid2go.txt has GO terms for all genes in TopGO format
geneID2GO <- readMappings(file = "ME_TAIR_geneid2go.txt") 
# Minotaur_results.csv has outlier status coded as 0 for not outlier, 1 for outlier and 2 for strong outlier
outlier_results= read.csv("Minotaur_results.csv")

# functions to select genes
selector95 <- function(theScore) {
+ return (theScore > 0)}
selector99 <- function(theScore) {
+ return (theScore > 1)}


FST_DXY_data = outlier_results$Md_FST_DXY_Outlier_Status
names(FST_DXY_data) = outlier_results$Transcript.ID
FST_DXY_data = na.omit(FST_DXY_data)

FST_DXY_GOdata_95 = new("topGOdata",description = "",ontology = "BP",allGenes = FST_DXY_data,geneSel = selector95,annot = annFUN.gene2GO,nodeSize = 10,gene2GO = geneID2GO)
FST_DXY_resultFisher.weig_95<- runTest(FST_DXY_GOdata_95, algorithm = "weight", statistic = "fisher")
FST_DXY_allRes_95 = GenTable(FST_DXY_GOdata_95, weight_fisher_P = FST_DXY_resultFisher.weig_95,topNodes = 50, numChar = 1000)
FST_DXY_allRes_sub_95 <- subset(FST_DXY_allRes_95, Significant > 1 & weight_fisher_P  < 0.05)
FST_DXY_allRes_sub_95$test <- "FST_DXY"
FST_DXY_allRes_sub_95$sig.level <- 95

FST_DXY_GOdata_99 = new("topGOdata",description = "",ontology = "BP",allGenes = FST_DXY_data,geneSel = selector99,annot = annFUN.gene2GO,nodeSize = 10,gene2GO = geneID2GO)
FST_DXY_resultFisher.weig_99<- runTest(FST_DXY_GOdata_99, algorithm = "weight", statistic = "fisher")
FST_DXY_allRes_99 = GenTable(FST_DXY_GOdata_99, weight_fisher_P = FST_DXY_resultFisher.weig_99,topNodes = 50, numChar = 1000)
FST_DXY_allRes_sub_99 <- subset(FST_DXY_allRes_99, Significant > 1 & weight_fisher_P  < 0.05)
FST_DXY_allRes_sub_99$test <- "FST_DXY"
FST_DXY_allRes_sub_99$sig.level <- 99
