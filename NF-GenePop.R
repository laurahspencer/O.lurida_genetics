# Oly genetic stats with GenePop program

install.packages("genepop")
library(genepop)

setwd("~/Documents/Roberts Lab/O.lurida_genetics")

# ===== NF 2016/2017 =======#
basic_info(inputFile="Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-Info.txt", verbose=T)

test_HW(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", which="Proba", outputFile = "Analyses/NF-HWE.txt", enumeration = FALSE, dememorization = 10000, batches = 100, iterations = 1000, verbose = interactive())

test_LD(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-LD.txt", dememorization = 10000, batches = 100, iterations = 1000, verbose = TRUE)

nulls(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-null.txt", nullAlleleMethod = "B96", CIcoverage = 0.95, verbose = TRUE)

test_diff(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-Diff.txt", genic=FALSE, pairs=TRUE, dememorization = 10000, batches = 100, iterations = 1000, verbose = TRUE)

Fst(inputFile = "Data/Oly2016NFH+2017NFW_Merged.txt", outputFile = "Analyses/NF-Fst.txt", sizes=F, pairs=TRUE, dataType="Diploid", verbose = TRUE)

genedivFis(inputFile="Data/Oly2016NFH+2017NFW_Merged.txt", sizes=FALSE, outputFile = "Analyses/NF-Fis.txt", dataType = "Diploid", verbose=interactive())


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

