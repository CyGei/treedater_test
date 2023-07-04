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
sts <- as.Date(annot$sampling_date) - min(as.Date(annot$sampling_date))
names(sts) <- annot$fasta_id
sts

#datetree
dtr <- treedater::dater(ml_tree, sts, s = 29375)

dtr
# rate * 365 days
dtr$mean.rate * 365

#number of subs per year over the whole genome:
dtr$mean.rate * 365 * 29375

plot(dtr)

#Red and black respectively indicate sample and internal nodes.
treedater::rootToTipRegressionPlot( dtr )

treedater::outlierTips( dtr , alpha = 0.20)
treedater::relaxedClockTest( dtr, sts, s = 29375, ncpu = 1) #error

