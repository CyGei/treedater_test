annot <- read.csv("data/annot.csv")
dna <-
  adegenet::fasta2DNAbin(file = "data/MEGAall_sequences_aligned_trimmed.fas")
rownames(dna) <- sub("\\s.*", "", rownames(dna))
dna <- dna[rownames(dna) %in% annot$fasta_id, ]



dna

#"""Building the ML tree with phangorn:"""
tre_ini <-  ape::nj(ape::dist.dna(dna, model = "TN93"))
phydat_dna <- phangorn::as.phyDat(dna)
fit_ini <- phangorn::pml(tre_ini, phydat_dna, k = 4)

# optimize tree topology (optNni = TRUE), base frequencies (optBf = TRUE),
# the rates of all possible subtitutions (optQ = TRUE), and use a gamma distribution
# to model variation in the substitution rates across sites (optGamma = TRUE):
fit <- phangorn::optim.pml(
  fit_ini,
  optNni = TRUE,
  optBf = TRUE,
  optQ = TRUE,
  optGamma = TRUE
)

anova(fit_ini, fit)
AIC(fit_ini)
AIC(fit)

# """Rooting the tree with the Wuhan strain"""
ml_tree <-
  ape::root(fit$tree,
            outgroup = c("WUHAN", "ALPHA"),
            resolve.root = TRUE)
ml_tree <- ape::ladderize(ml_tree)

plot(ml_tree)

#"""# Tree Dater"""
library(treedater)
#make sure we match the tree tip label with the correct fasta ID
annot <- annot[match(ml_tree$tip.label, annot$fasta_id), ]

#sample times
library( lubridate )
sts <- as.Date(annot$sampling_date) |> decimal_date()
names(sts) <- annot$fasta_id
sts

#datetree
dtr <- treedater::dater(ml_tree, sts, s = 29375)

dtr
dtr$mean.rate 


plot(dtr)

#Red and black respectively indicate sample and internal nodes.
treedater::rootToTipRegressionPlot( dtr )

ot<- treedater::outlierTips( dtr , alpha = 0.20)
#remove 2 outliers not including outgroup: 
ml_tree2 <- drop.tip( ml_tree, c('OQ866550.1', 'OQ866547.1'))
#redo
(dtr2 <- dater( ml_tree2, sts, s = 29375))
# Note: Initial guess of substitution rate not provided. Will attempt to guess starting conditions. Provide initial guesses of the rate using *omega0* parameter. 
# Note: *dater* called with non binary tree. Will proceed after resolving polytomies.
# Note: Minimum temporal branch length  (*minblen*) set to 0.00602739726027415. Increase *minblen* in the event of convergence failures. 
# Tree is rooted. Not estimating root position.
# Initial guesses of substitution rate: 0.000285044944986525,0.000360833609735679,0.000577784875073646,0.000668354765635498,0.000870524805160766,0.000975875921535318 
# 
# Phylogenetic tree with 52 tips and 51 internal nodes.
# 
# Tip labels:
#   WUHAN, ALPHA, OQ800894.1, OQ800895.1, OQ800896.1, OQ800897.1, ...
# 
# Rooted; includes branch lengths.
# 
#  Time of common ancestor 
# 2014.00245114339 
# 
#  Time to common ancestor (before most recent sample) 
# 9.13179543195247 
# 
#  Weighted mean substitution rate (adjusted by branch lengths) 
# 0.000127042983380417 
# 
#  Unadjusted mean substitution rate 
# 0.000127042983380417 
# 
#  Clock model  
# strict 
# 
#  Coefficient of variation of rates 
# 0 

treedater::relaxedClockTest( dtr2, sts, s = 29375, ncpu = 4) 
#favours uncorrelated rc

#redo
(dtr3 = dater( ml_tree2, sts, s = 29375, clock = 'uncorrelated', omega0=0.001, numStartConditions = 0))
dtr3
# 
# NOTE: The estimated coefficient of variation of clock rates is high (>1). Sometimes this indicates a poor fit and/or a problem with the data.
# 		
# 
# The following steps may help to fix or alleviate common problems: 
# * Check that the vector of sample times is correctly named and that the units are correct. 
# * If passing a rooted tree, make sure that the root position was chosen correctly, or estimate the root position by passing an unrooted tree (e.g. pass ape::unroot(tree))
# * The root position may be poorly estimated. Try increasing the _searchRoot_ parameter in order to test more lineages as potential root positions. 
# * The model may be fitted by a relaxed or strict molecular clock. Try changing the _clock_ parameter 
# * A poor fit may be due to a small number of lineages with unusual / outlying branch lengths which can occur due to sequencing error or poor alignment. Try the *outlierTips* command to identify and remove these lineages. 
# * Check that there is adequate variance in sample times in order to estimate a molecular clock by doing a root-to-tip regression. Try the *rootToTipRegressionPlot* command. If the clock rate can not be reliably estimated, you can fix the value to a range using the _meanRateLimits_ option which would estimate a time tree given the previous estimate of clock rates. 
# 
# Phylogenetic tree with 52 tips and 51 internal nodes.
# 
# Tip labels:
#   WUHAN, ALPHA, OQ800894.1, OQ800895.1, OQ800896.1, OQ800897.1, ...
# 
# Rooted; includes branch lengths.
# 
#  Time of common ancestor 
# 2019.93447001563 
# 
#  Time to common ancestor (before most recent sample) 
# 3.19977655971434 
# 
#  Weighted mean substitution rate (adjusted by branch lengths) 
# 0.000366419608010347 
# 
#  Unadjusted mean substitution rate 
# 0.00268129512911975 
# 
#  Clock model  
# uncorrelated 
# 
#  Coefficient of variation of rates 
# 2.87891336416268 

plot(dtr3)

pb3 <- parboot( dtr3, ncpu = 4, nreps = 200)
pb3
#                            pseudo ML        2.5 %       97.5 %
# Time of common ancestor 2.019934e+03 2.017026e+03 2019.9879452
# Mean substitution rate  2.681295e-03 7.100237e-04    0.0101255
# 
#  For more detailed output, $trees provides a list of each fit to each simulation 
