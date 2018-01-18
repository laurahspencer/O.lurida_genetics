

# Have we sampled enough loci? 
NFH.2016.acccurve <- genotype_curve(NFH.2016, sample=1000, quiet=T) 
NFW.2017.acccurve <- genotype_curve(NFW.2017, sample=1000, quiet=T)
NF.acccurve <- genotype_curve(NF, sample=1000, quiet=T)

# Summary stats on # alleles
NFH.2016.table <- locus_table(NFH.2016) #mean # alleles = 18.67
NFW.2017.table <- locus_table(NFW.2017) #mean # alleles = 19.33
NF.table <- locus_table(NF)

info_table(NF, type="missing", plot=TRUE) #see how the missing data is distributed over the 2 populations
mlg.table(NF) #genotype eveness. Result is N=199; MLG=199
NF.pop <- poppr(NF) #summary stats on each population
(NF.pop$N / (NF.pop$N - 1)) * NF.pop$lambda #corrected simpson's index (N/(N-1)) #all different genotypes

#---- Abbreviation	Statistic ---#
# Pop: Population name.
# N: Number of individuals observed.
# MLG:	Number of multilocus genotypes (MLG) observed.
# eMLG:	The number of expected MLG at the smallest sample size ≥ 10 based on rarefaction
# SE:	Standard error based on eMLG.
# H:	Shannon-Wiener Index of MLG diversity (Shannon, 2001).
# G:	Stoddart and Taylor’s Index of MLG diversity (Stoddart & Taylor, 1988).
# lambda:	Simpson’s Index (Simpson, 1949). 0 = no genotypes are differet; 1 = all genotypes are different
# E.5:	Evenness, E5E5 (Pielou, 1975; Ludwig & Reynolds, 1988; Grünwald et al., 2003).
# Hexp:	Nei’s unbiased gene diversity (Nei, 1978).
# Ia:	The index of association, IAIA (Brown, Feldman & Nevo, 1980; Smith et al., 1993).
# rbarD:	The standardized index of association, r¯dr¯d [@].

# Are our populations in Hardy-Weinberg equilibrium? 
# Hardy‐Weinberg Assumptions
# • infinite population
# • discrete generations
# • random mating
# • no selection
# • no migration in or out of population
# • no mutation
# • equal initial genotype frequencies in the two sexes
# Equilibrium is reached after one generation of mating under the
# Hardy‐Weinberg assumptions...Genotype frequencies remain the same from generation to generation.


library("pegas")
NF.HW <- seppop(NF) %>% lapply(hw.test, B=1000) #all P-values >0.05; do not reject the null that these populations are under HWE. 
NF.HW.P <- sapply(test, "[", i=TRUE, j=3) #pvalues of HW chi-squared test for all loci, both pops

# Are populations in linkage disequilibrium? The null hypothesis tested is that alleles observed at different loci are not linked if populations are sexual while alleles recombine freely into new genotypes during the process of sexual reproduction
# IA =VO/VE -1 
# ... where V0 is the observed variance of K and VE is the expected variance of K, where K is the number of loci at which two individuals differ.

library("magrittr")
NF.ia.H <- ia(popsub(NF, "NFH-2016"), sample=999)
NF.ia.W <- ia(popsub(NF, "NFW-2017"), sample=999)
#we find significant support for the hypothesis that alleles are linked across loci with P<0.001

#which loci are linked? run pairwise assessment: 
NF.W2017.pair <- pair.ia(popsub(NF, "NFW-2017"))
NF.H2016.pair <- pair.ia(popsub(NF, "NFH-2016"))
pair.range <- range(c(NF.W2017.pair, NF.H2016.pair), na.rm=TRUE)
plot(NF.W2017.pair, limits=pair.range)
plot(NF.H2016.pair, limits=pair.range)
NF.W2017.pair
NF.H2016.pair

# it looks like loci 13, 15 & 19 are possibly linked

NF.freq <- rraf(NF, by_pop=TRUE)
NF.freq.t <- t(NF.freq)
plot(NF.freq.t)

#---- Generate statistics using a different tutorial / R package ------# 
# http://popgen.nescent.org/startMicrosatellite.html
library("adegenet")
library("pegas")
library("hierfstat")
NF.summary <- summary(NF)
NFW.2017.summary <- summary(popsub(NF, "NFW-2017"))
NFW.2017.summary
plot(NFW.2017.summary$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus, NFW 2017")
plot(NFW.2017.summary$Hobs, NFW.2017.summary$Hexp, xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected ~ Observed Heterozygosity per locus, NFW 2017")
bartlett.test(list(NFW.2017.summary$Hexp, NFW.2017.summary$Hobs)) #indicates no difference between mean observed and expected heterozygosity 

NFH.2016.summary <- summary(popsub(NF, "NFH-2016"))
NFH.2016.summary
plot(NFH.2016.summary$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus, NFH 2016")
plot(NFH.2016.summary$Hobs, NFH.2016.summary$Hexp, xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected ~ Observed Heterozygosity per locus, NFH 2016")
bartlett.test(list(NFH.2016.summary$Hexp, NFH.2016.summary$Hobs)) #indicates no difference between mean observed and expected heterozygosity 

# ------Calculate and plot allelic frequencies-----#

#Used the following script: http://www.molecularecologist.com/wp-content/uploads/2012/03/Allelefrequency_calculations2.txt

NF.csv <- read.csv(file="Data/Oly2016NFH+2017NFW_Merged.csv", stringsAsFactors = F)
names(NF.csv) <- NF.csv[2,]
NF.csv <- NF.csv[-1:-2,]
NFW.df <- subset(NF.csv, NF.csv$Population == "NFW-2017")
NFH.df <- subset(NF.csv, NF.csv$Population == "NFH-2016")
NFW.df.1 <- NFW.df[,-1:-2]
NFH.df.1 <- NFH.df[,-1:-2]

allelic.freq <- function(df) {
  L=ncol(df)   #how many columns are there?
  locus_positions=(2*(unique(round((1:(L-2))/2)))+1)   #find the starting column number for each locus
  lnames=colnames(df)                          #locus names, from the header
  OUT=list()        #create a null dataset to append allele freqs to
  
  for (x in locus_positions) {                       #begin for loop, to calculate frequencies for each locus
    alleles=c(df[,x],df[,x+1])        #For example, combine columns 1 and 2 for locus 1 (two columns because they are diploid)
    alleles2=as.data.frame(table(alleles))             #count each allele at locus x
    missing=alleles2[which(alleles2[,1]==0),2]          #count missing data at locus x, entered as '0' in this dataset (not used further for simplicity)
    alleles3=alleles2[which(!alleles2[,1]==0),]          #remove missing data (otherwise 0 would be counted in total number of alleles)
    alleles4=cbind(alleles3,alleles3[,2]/sum(alleles3[,2])) #calculate frequencies
    output=cbind(x,lnames[x],alleles4)                        #combine x, locusname, and frequencies
    OUT[[x]] <- output
  }  
  OUT.1 <- do.call(rbind, OUT)
  colnames(OUT.1) <- c("Number","Locus","allele","count","frequency") #add column headers
  return(OUT.1)
}

NFH.afreq <- allelic.freq(NFW.df.1)[,-1]
NFW.afreq <- allelic.freq(NFH.df.1)[,-1]
write.table(NFH.afreq,file="Analyses/NFH-2016-Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(NFW.afreq,file="Analyses/NFW-2017-Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)

Freq.Plots <- function(wanted_locus, frequency_table, population) {
  Locus=frequency_table[which(frequency_table[,1]==wanted_locus),]
  plot(as.numeric(as.character(Locus[,2])),as.numeric(as.character(Locus[,4])),xlab="Allele",ylab="Frequency",main=paste(population, "_", "Locus_",Locus[1,1]),pch=21,bg="blue",cex=1.5)
  plot(1:length(Locus[,2]),sort(as.numeric(as.character(Locus[,4])),decreasing=TRUE),xlab="Allele (orderd by frequency)",ylab="Frequency",main=paste(population, "_", "Locus_",Locus[1,1], sep=""),pch=21,bg="blue",cex=1.5)
}

# Generate plots for NFH-2016 population
loci <- unique(NFH.afreq$Locus)
for (i in 1:length(loci)) {
  path <- file.path(paste("Analyses/", "FreqPlots_NFH-2016_", loci[i], ".pdf", sep = ""))
  pdf(file=path)
  Freq.Plots(loci[i], NFH.afreq, "NFH-2016")
  dev.off()
}

# Generate plots for NFW-2017 population
loci <- unique(NFW.afreq$Locus)
for (i in 1:length(loci)) {
  path <- file.path(paste("Analyses/", "FreqPlots_NFW-2017_", loci[i], ".pdf", sep = ""))
  pdf(file=path)
  Freq.Plots(loci[i], NFW.afreq, "NFW-2017")
  dev.off()
}



