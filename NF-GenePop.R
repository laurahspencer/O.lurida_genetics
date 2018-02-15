# Oly genetic stats with GenePop program

install.packages("genepop")
library(genepop)

setwd("~/Documents/Roberts Lab/O.lurida_genetics")

# ===== NF 2016/2017 =======#

# Get summary info 
basic_info(inputFile="Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-Info.txt", verbose=T)

# Hardy Weinberg test - null is that alleles are in HW equilibrium 
test_HW(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", which="Proba", outputFile = "Analyses/NF-HWE.txt", enumeration = FALSE, dememorization = 10000, batches = 100, iterations = 1000, verbose = interactive())

# Loci disequilibrium test - null is that alleles are in equilibrium 
test_LD(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-LD.txt", dememorization = 10000, batches = 100, iterations = 1000, verbose = TRUE)

# Check for null alleles
nulls(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-null.txt", nullAlleleMethod = "B96", CIcoverage = 0.95, verbose = TRUE)

# Run exact test for genotypic differential; null hypothesis is there is not difference
test_diff(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-Diff.txt", genic=FALSE, pairs=TRUE, dememorization = 10000, batches = 100, iterations = 1000, verbose = TRUE)

# Calculate FST, the proportion of the total genetic variance contained in a subpopulation (the S subscript) relative to the total genetic variance (the T subscript). Values can range from 0 to 1. High FST implies a considerable degree of differentiation among populations. 
Fst(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-Fst.txt", sizes=F, pairs=TRUE, dataType="Diploid", verbose = TRUE)

#Evaluates Fis and gene diversities. FIS (inbreeding coefficient) is the proportion of the variance in the subpopulation contained in an individual. High FIS implies a considerable degree of inbreeding. 
genedivFis(inputFile="Data/Oly2016NFH+2017NFW_Merged.txt", sizes=FALSE, outputFile = "Analyses/NF-Fis.txt", dataType = "Diploid", verbose=interactive())

#Power Analysis
# How? # 

# ------Calculate and plot allelic frequencies-----#

#Used the following script: http://www.molecularecologist.com/wp-content/uploads/2012/03/Allelefrequency_calculations2.txt

# Read in .csv of hatchery & wild NF data to generate allele frequencies 
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


# ======== 2014 SS & PGB data =========#

basic_info(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", outputFile = "Analyses/2014-Info.txt", verbose=T)
#

test_HW(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", which="Proba", outputFile = "Analyses/2014-HWE.txt", enumeration = FALSE, dememorization = 10000, batches = 100, iterations = 1000, verbose = interactive())
# All loci in HWE for all pops 

test_LD(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", outputFile = "Analyses/2014-LD.txt", dememorization = 10000, batches = 100, iterations = 1000, verbose = TRUE)

nulls(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", outputFile = "Analyses/2014-null.txt", nullAlleleMethod = "B96", CIcoverage = 0.95, verbose = TRUE)

test_diff(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", outputFile = "Analyses/2014-Diff.txt", genic=FALSE, pairs=TRUE, dememorization = 10000, batches = 100, iterations = 1000, verbose = TRUE)

Fst(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", outputFile = "Analyses/2014-Fst.txt", sizes=F, pairs=TRUE, dataType="Diploid", verbose = TRUE)

genedivFis(inputFile="Data/oly_msats_2015/2015-SS-PG-GenePop.txt", sizes=FALSE, outputFile = "Analyses/2014-Fis.txt", dataType = "Diploid", verbose=interactive())

# Power Analysis? 

# Read in .csv of hatchery & wild NF data to generate allele frequencies 
SS.PG.2014.csv <- read.csv(file="Data/oly_msats_2015/oly_msats_tandem_edited_LHS.csv", stringsAsFactors = F)
# names(SS.GB.2014.csv) <- SS.GB.2014.csv[1,]
PG.W.2014 <- subset(SS.PG.2014.csv, SS.GB.2014.csv$Population == "PG-W")
PG.H.2014 <- subset(SS.PG.2014.csv, SS.GB.2014.csv$Population == "PG-H")
SS.W.2014 <- subset(SS.PG.2014.csv, SS.GB.2014.csv$Population == "SS-W")
SS.H.2014 <- subset(SS.PG.2014.csv, SS.GB.2014.csv$Population == "SS-H")
PG.W.2014.1 <- PG.W.2014[,-1:-2]
PG.H.2014.1 <- PG.H.2014[,-1:-2]
SS.W.2014.1 <- SS.W.2014[,-1:-2]
SS.H.2014.1 <- SS.H.2014[,-1:-2]
PG.W.2014.afreq <- allelic.freq(PG.W.2014.1)[,-1]
PG.H.2014.afreq <- allelic.freq(PG.H.2014.1)[,-1]
SS.W.2014.afreq <- allelic.freq(SS.W.2014.1)[,-1]
SS.H.2014.afreq <- allelic.freq(SS.H.2014.1)[,-1]
write.table(PG.W.2014.afreq,file="Analyses/PGW-2014-Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(PG.H.2014.afreq,file="Analyses/PGH-2014-Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(SS.W.2014.afreq,file="Analyses/SSW-2014-Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(SS.H.2014.afreq,file="Analyses/SSH-2014-Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)

# Allele frequency dataframes 
NFH.afreq$Population <- "NFH.2016"
NFW.afreq$Population <- "NFW-2017"
PG.W.2014.afreq$Population <- "PGW.2014"
PG.H.2014.afreq$Population <- "PGH.2014"
SS.W.2014.afreq$Population <- "SSW.2014"
SS.H.2014.afreq$Population <- "SSH.2014"

#Merge allele frequency dataframes from each population (2014, 2016) into 1 
allele.frequencies <- rbind(NFH.afreq, NFW.afreq, PG.W.2014.afreq, PG.H.2014.afreq, SS.W.2014.afreq, SS.H.2014.afreq)
allele.frequencies$allele <- as.numeric(as.character(unlist(allele.frequencies$allele)), na.action=na.omit)

# Generate bar plots
library(plotly)
Olur10 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur10"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur10 Allele Frequency",
         yaxis = list(title = 'Frequency'),
         xaxis = list(title = 'Allele Size (bp)'),
         legend = list(x=.95, y=.95))

Olur11 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur11"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur11 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

Olur12 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur12"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur12 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

Olur13 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur13"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur13 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

Olur15 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur15"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur15 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

Olur17 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur17"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur17 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

Olur18 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur18"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur18 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

Olur19 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Locus == "Olur19"), x = ~allele, y = ~count, type="bar", color=~Population, hovertext=~count) %>%  #generate plotly plot
  layout(title="Olur19 Allele Frequency",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

allele.frequencies$type <- allele.frequencies$Population
allele.frequencies$type <- gsub("NFH.2016", "F1", allele.frequencies$type)
allele.frequencies$type <- gsub("NFW-2017", "Wild", allele.frequencies$type)
allele.frequencies$type <- gsub("PGH.2014", "F1", allele.frequencies$type)
allele.frequencies$type <- gsub("PGW.2014", "Wild", allele.frequencies$type)
allele.frequencies$type <- gsub("SSH.2014", "F1", allele.frequencies$type)
allele.frequencies$type <- gsub("SSW.2014", "Wild", allele.frequencies$type)

OlurAll.NF.2016 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Population == "NFH.2016" | allele.frequencies$Population == "NFW-2017"), x = ~allele, y = ~count, type="bar", color=~type, hovertext=~count) %>%  #generate plotly plot
  layout(title="Fidalgo Bay Allele Frequency, All Loci (2016)",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

OlurAll.SS.2014 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Population == "SSH.2014" | allele.frequencies$Population == "SSW.2014"), x = ~allele, y = ~count, type="bar", color=~type, hovertext=~count) %>%  #generate plotly plot
  layout(title="South Sound Allele Frequency, All Loci (2014)",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))

OlurAll.PG.2014 <- plot_ly(data = subset(allele.frequencies, allele.frequencies$Population == "PGH.2014" | allele.frequencies$Population == "PGW.2014"), x = ~allele, y = ~count, type="bar", color=~type, hovertext=~count) %>%  #generate plotly plot
  layout(title="Port Gamble Allele Frequency, All Loci (2014)",
         xaxis = list(title = 'Allele Size (bp)'),
         yaxis = list(title = 'Frequency'),
         legend = list(x=.95, y=.95))
