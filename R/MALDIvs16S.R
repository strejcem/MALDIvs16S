## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----MSdata--------------------------------------------------------------
library(here)
library(tidyverse)
library(cowplot)
# load convenience R functions for MS data manipulations
source(here("MALDIbacteria.R"))

## import MS and Bruker Biotyper data
# Ms spectra from MzMl format

suppressMessages(
  s <- importMzMl(path = here("mzML"))
)

# assign name to spectra from file names
fullNames <- metaParam(s, "file")
fullNames <- str_replace(fullNames, pattern = "^.*\\\\(.*)\\.mzML$", "\\1")
s <- chngMetaParam(s, "fullName", fullNames)

# sort spectra by their full names
s <- s[order(metaParam(s, "fullName"))]

# read MS spectra metadata with BioTyper identification data
metadata <- read_tsv(here("metadata.tsv"))

# incorporate metadata into 'MassSpectrum' objetcs
identical(metadata$fullName, metaParam(s, "fullName"))

for (param in names(metadata)) {
  s <- chngMetaParam(s, param, pull(metadata, param))
}

## ----screen, warning=5---------------------------------------------------
# process MS data with convenience function 'classicMaldi()' that wraps the whole
MS <- classicMALDI(s, range = c(2000, 20000))

# to print mass dedrogram of all mass spectra (Supplementary Figure 1) you can
# use PDFplot() function.
PDFplot(MS$hc, file = "dendrogram.pdf", cutoffs = c(1-0.79, 1-0.92))

# plot distribution of average cosine similarity of each spectrum to its technical replicates
qplot(MS$ACSreps, geom = "histogram", bins = 100)

## ----outliers------------------------------------------------------------
s <- s[-(which(MS$ACSreps < 0.85))]

## ----taxonomy16S---------------------------------------------------------
EzT <- read_tsv(here("EzTaxon.tsv"))

species <- plyr::mapvalues(metaParam(s, "culture"),
                           from = EzT$culture,
                           to = EzT$species
                           )

## ----sda, warning=F------------------------------------------------------
library(sda)
MS <- classicMALDI(s, range = c(2000, 20000))
ldar <- sda.ranking(Xtrain=MS$fm, L = as.factor(species),
                    ranking.score = "entropy", fdr = TRUE, diagonal = F)

# number of predictive signals
sum(ldar[,"lfdr"] < 0.2)

# construct data frame from sda data for plotting purposes
ldadf <- data.frame(ldar[,])
ldadf$peaks <- as.numeric(rownames(ldadf))
ldadf$cutoff <- ifelse(ldadf$lfdr < 0.2, "predictive", "non-predictive")

# plot predictive vs non-predictive peaks
ggplot(ldadf, aes(peaks)) + 
  geom_histogram(aes(fill = cutoff), alpha = 1, binwidth = 250,
                 color = "black", position = "identity") +
  geom_vline(xintercept = c(4000, 10000), linetype = 5) +
  labs(x=expression(paste(italic("m/z"), "[Da]")),
       y = "count", fill = element_blank()) +
  annotate("label", x = c(4000, 10000), y = 50, label = c(4000, 10000)) +
  theme(legend.position = c(0.76, 0.9))

## ----sdaValidation-------------------------------------------------------

library(crossval)

#define prediction function
predfun <- function(Xtrain, Ytrain, Xtest, Ytest,
                      numVars, diagonal=F) {
    # estimate ranking and determine the best numVars variables
    ra <- sda.ranking(Xtrain, Ytrain,
                      verbose=FALSE, diagonal=diagonal, fdr=FALSE)
    selVars <- ra[,"idx"][1:numVars]
    # fit and predict
    sda.out <- sda(Xtrain[, selVars, drop=FALSE], Ytrain,
                   diagonal=diagonal, verbose=FALSE)
    ynew <- predict(sda.out, Xtest[, selVars, drop=FALSE],
                    verbose=FALSE)$class
    # compute accuracy
    acc <- mean(Ytest == ynew)
    return(acc)
  }

# extract 'culture' metada attribute from all sprectra to a vector
cultures <- metaParam(s, "culture")
cv <- crossval(predfun, X=MS$fm, Y=factor(cultures), diagonal=F, verbose=F,
               K=10, B=10,  # 10x10 cross-validation
               numVars=150)  # number of top predicting peaks to test

## Acuracy
cv$stat

## ----reproducibility-----------------------------------------------------
# Interval for Ms ranges to be calculated
unit <- 1000
intervals <- seq(2000, 20000, unit)
iteration <- seq_along(intervals)
iteration <- iteration[-length(iteration)]

MSranges <- NULL
cumulInt <- 0
fm <- MS$fm
techReps <- metaParam(s, "replicate")
cultures <- metaParam(s, "culture") # biological replicates
fmRange <- as.numeric(colnames(fm))

for(iter in iteration) {
  bot <- intervals[iter]
  upper <- intervals[iter+1]

  predictives <- sum(ldadf$peaks >= bot & ldadf$peaks < upper &
                       ldadf$cutoff == "predictive")
  
  fmSubset <- fm[, fmRange >= bot & fmRange < upper]
  numPeaks <- ncol(fmSubset)
  Int <- sum(fmSubset)
  
  avgCosRange <- avgCos(fmSubset, techReps)
  avgCosRange <- mean(avgCosRange$meanACS, na.rm = TRUE)
  
  avgCosRangeS <- avgCos(fmSubset, cultures)
  avgCosRangeS <- mean(avgCosRangeS$meanACS, na.rm = TRUE)

  MSranges <- as.data.frame(rbind(MSranges,
                                  c(bot = bot,up = upper, numPeaks = numPeaks,
                                    Int = Int, predPeaks = predictives,
                                    ACSR = avgCosRange, ACSS = avgCosRangeS)))
}
MSranges

## ----intervalPlot--------------------------------------------------------
intenzity <- data.frame(mz = rep(MSranges$bot, each = 2),
                        Int = rep(MSranges$Int, each = 2),
                        ACSR = rep(MSranges$ACSR, each = 2),
                        ACSS = rep(MSranges$ACSS, each = 2))
intenzity$mz <- intenzity$mz + c(unit/2-unit/4, unit/2+unit/4)

# generate the plot
ggplot(MSranges, aes(bot+unit/2)) +
  geom_area(data = intenzity, aes(x = mz, y = Int/max(Int)),
            alpha = 0.4, fill = "green") +
  geom_col(aes(y = numPeaks/max(numPeaks)), width = unit/2, fill = 1,
           color = "black", alpha = 0.3) +
  geom_col(aes(y = numPeaks/max(numPeaks)), width = unit/2, fill = NA,
           color = "black", alpha = 1, size = 0.05) +
  geom_col(aes(y = predPeaks/max(numPeaks)), width = unit/2, fill = "blue",
           color = "black", alpha = 0.8, size = 0.1) +
  geom_text(aes(y = numPeaks/max(numPeaks), label = numPeaks), size = 3,
            position = position_nudge(y = 0.02)) +
  geom_text(aes(y = predPeaks/max(numPeaks), 
                label = scales::percent(predPeaks/max(numPeaks))), size = 3,
            position = position_nudge(y = 0.07), col = "black", angle = 90) +
  geom_line(data = intenzity, aes(x = mz, y = ACSS), col = "red", size = 1) +
  labs(x = expression(italic("m/z")), y = "") +
  scale_x_continuous(name = expression(italic("m/z")),
                     breaks = seq(2000, 20000, 2000),
                     minor_breaks = seq(2000, 20000, 1000)) +
  scale_y_continuous(name = "", labels = scales::percent)

## ----data16S, message=F--------------------------------------------------
# read 16S data
library(DECIPHER)
dna <- readDNAStringSet(here("16S.fas"))

simil16S <- matrix(100, nrow = length(dna), ncol = length(dna), 
                  dimnames = list(names(dna), names(dna)))

# To reduce computational time
validM <- lower.tri(simil16S, diag = F)

# can be also parallelized
for (x in seq_along(dna)) {
  for (y in seq_along(dna)) {
    if(!validM[x,y]) next
    simil16S[x,y] <- pid(pairwiseAlignment(dna[x], dna[y], type = "local"),
                         type = "PID2")
  }
}

simil16S <- simil16S / 100
simil16S[upper.tri(simil16S)] <- simil16S[lower.tri(simil16S)]

# UPGMA dendrogram
plot(hclust(as.dist(1-simil16S), method = "ave"), hang=-1)

## ----16Spair-------------------------------------------------------------
# to discard the duplicated entries later
simil16S[upper.tri(simil16S, diag = F)] <- -1

# transform full matrix into 2-column matrix
df16S <- reshape2::melt(simil16S, varnames = c("row", "col"), as.is = TRUE)

# now discard duplicated entries
df16S <- df16S[df16S$value >= 0, ]

df16S$pair <- paste(df16S$row, df16S$col, sep = "-")

## ----MSpair, warning=5---------------------------------------------------
# the default MS range is 4000-10000, here it is just to stress it out 
MS <- classicMALDI(s, range = c(4000, 10000), labels = "fullName")
similMS <- 1-as.matrix(MS$d, dimnames = list(labels(MS$d), labels(MS$d)))
similMS[upper.tri(similMS, diag = F)] <- -1
dfMS <- reshape2::melt(similMS, varnames = c("row", "col"), as.is = TRUE)
dfMS <- dfMS[dfMS$value >= 0, ]

pairUnique <- paste(dfMS$row, dfMS$col, sep = "-")
dfMS$pair <- str_replace_all(pairUnique, "_[0-9]_[0-9]", "")

## ----pairData------------------------------------------------------------
pairData <- merge(dfMS, df16S, by = "pair")
names(pairData) <- c("name", "nameAunique", "nameBunique", "similMS",
                          "nameA", "nameB", "simil16S") 

## ----taxonomyPairData----------------------------------------------------
for (taxLevel in names(EzT)[-1]) {
  pairData[paste0(taxLevel,"A")] <-
    plyr::mapvalues(pairData$nameA, EzT$culture, EzT[[taxLevel]])

  pairData[paste0(taxLevel,"B")] <-
    plyr::mapvalues(pairData$nameB, EzT$culture, EzT[[taxLevel]])
}

# identify the lowest shared taxonomy level between the pair
pairData$taxRelation <- case_when(
    pairData$speciesA == pairData$speciesB ~ "Species",
    pairData$genusA == pairData$genusB ~ "Genus",
    pairData$familyA == pairData$familyB ~ "Family",
    pairData$classA == pairData$classB ~ "Class",
    pairData$orderA == pairData$orderB ~ "Order",
    pairData$phylumA == pairData$phylumB ~ "Phylum",
    TRUE ~ "Domain")


# add alpha chanel value for ploting purposes
pairData$alpha <- plyr::mapvalues(pairData$taxRelation,
                c("Species", "Genus", "Family", "Class", "Order", "Phylum",
                  "Domain"),
                c(0.5, 0.5, 0.3, 0.3, 0.2, 0.2, 0.1))



## ----similarityRelations-------------------------------------------------
# group data by culture
dfGrouped <- group_by(pairData, name) %>%
  summarise(avgMS = mean(similMS), avg16S = mean(simil16S), sd = sd(similMS),
            relation = unique(taxRelation), alpha = as.numeric(unique(alpha)))
dfGrouped$relation <- factor(dfGrouped$relation,
         levels = c("Species", "Genus", "Family", "Class", "Order", "Phylum",
                  "Domain"))

# generate plot
ggplot(dfGrouped, aes(avgMS, avg16S*100, col = relation, alpha = alpha)) +
  geom_errorbarh(aes(xmax= avgMS + sd, xmin = avgMS - sd), alpha = 0.2) +
  scale_colour_manual(name = element_blank(),
                      values = c("red", "blue", "green", "purple", "orange",
                                 "pink", "grey")) +
  geom_point() +
  scale_alpha(guide = 'none') +
  geom_hline(yintercept = 98.65, linetype = 5) +
  geom_vline(xintercept = c(0.79, 0.92), linetype = c(2,1)) +
  annotate("label", x=c(0.13,0.79, 0.92), y=c(98.65, 78, 78),
           label = c(98.65, 0.79, 0.92))  +
  labs(x="Mass spectra cosine similarity",
       y="16S rRNA gene sequence similarity [%]") +
  theme(legend.position = c(0.5, 0.4),
        legend.background = element_rect(
          fill = "white", linetype = 1, color = 1))


## ----F1score, message=F--------------------------------------------------
library(BiocParallel)
set.seed(123)
B <- 1
K <- 2
beta <- 1 # this defines the F-beta score

# parallelized F1 score calculation for the closest type strain (phylotype)
# classification
closestTypeF1 <- bplapply(1:B, function(B, K, beta, pairData) {
  
  source('MALDIbacteria.R')
  resultsB <- NULL
  fold <- sample(1:K, size = nrow(pairData), replace = T)
  for (h in seq(0, 0.99, 0.01)) {
    resultsK <- lapply(1:K, function(K, h) {
      Xval <- pairData[fold != K, ]
      # base thruth is: the closest type strain of culture A (speciesA) is the same
      # as the closest type strain of culture B (speciesB)
      F1train <- FscPrecRecAcc(Xval$similMS >= h,
                               Xval$speciesA == Xval$speciesB, beta = beta)
      names(F1train) <- paste0("train.", names(F1train))
      
      Xval <- pairData[fold == K, ]
      F1test <- FscPrecRecAcc(Xval$similMS >= h, Xval$speciesA == Xval$speciesB,
                              beta = beta)
      names(F1test) <- paste0("test.", names(F1test))
      
      output <- cbind(F1train[3], F1test[c(1,2,4)])
      return(output)
    }, h=h)
    
    resultsK <- do.call(rbind, resultsK)
    resultsK <- apply(resultsK, 2, mean)
    resultsK <- as.data.frame(t(c(h=h, resultsK)))
    
    resultsB <- rbind(resultsB, resultsK)
  }
  maxF <- which.max(resultsB$test.Fscore)
  return(cbind(B=B, resultsB))
}, K=K, beta=beta, pairData=pairData)

closestTypeF1 <- do.call(rbind, closestTypeF1)
closestTypeF1 <- closestTypeF1 %>%
  select(-B) %>%
  group_by(h) %>%
  summarise_all(mean)
  

# similarity based classification
similarity9865F1 <- bplapply(1:B, function(B, K, beta, pairData) {
  source('MALDIbacteria.R')
  resultsB <- NULL
  fold <- sample(1:K, size = nrow(pairData), replace = T)
  for (h in seq(0, 0.99, 0.01)) {
    resultsK <- lapply(1:K, function(K, h) {
      Xval <- pairData[fold != K, ]
      # base thruth is: 16S rRNA gene similarity >= 98.65%
      F1train <- FscPrecRecAcc(Xval$similMS >= h, Xval$simil16S >= 0.9865,
                               beta = beta)
      names(F1train) <- paste0("train.", names(F1train))
      
      Xval <- pairData[fold == K, ]
      F1test <- FscPrecRecAcc(Xval$similMS >= h, Xval$simil16S >= 0.9865,
                              beta = beta)
      names(F1test) <- paste0("test.", names(F1test))
      
      output <- cbind(F1train[3], F1test[c(1,2,4)])
      return(output)
    }, h=h)
    resultsK <- do.call(rbind, resultsK)
    resultsK <- apply(resultsK, 2, mean)
    resultsK <- as.data.frame(t(c(h=h, resultsK)))
    
    resultsB <- rbind(resultsB, resultsK)
  }
  maxF <- which.max(resultsB$test.Fscore)
  return(cbind(B=B, resultsB))
}, K=K, beta=beta, pairData=pairData)

similarity9865F1 <- do.call(rbind, similarity9865F1)
similarity9865F1 <- similarity9865F1 %>%
  select(-B) %>%
  group_by(h) %>%
  summarise_all(mean)

## ----F1------------------------------------------------------------------
dfPrecRec <- merge(closestTypeF1, similarity9865F1, by="h")

names(dfPrecRec) <- c("h", "F1.type", "prec.type","rec.type", "acc.type",
                      "F1.seq","prec.seq", "rec.seq", "acc.seq")

dfPRF <- reshape2::melt(dfPrecRec, id.vars = "h")
dfPRF <- dfPRF %>% mutate(
  Reference =
    case_when(
      grepl("type", variable) ~ "Closest Type Strain",
      grepl("seq", variable) ~ "Similarity based"
    ),
  color =
    case_when(
      grepl("F1", variable) ~ "F1 score",
      grepl("prec", variable) ~ "Precision",
      grepl("rec", variable) ~ "Recall",
      grepl("acc", variable) ~ "Accuracy"
    ),
  lineType =
    case_when(
      grepl("seq", variable) ~ "Similarity based",
      grepl("type", variable) ~ "Closest Type Strain"
    )
  )

#not symetric dataset, skewed to many POSITIVEs, accuracy is not a good performance measure
dfPRF <- dfPRF %>% filter(color != "Accuracy")


# max F1
maxType <- dfPRF %>%
  filter(variable == "F1.type") %>%
  top_n(1, value) %>%
  select(h) %>%
  summarise(h = mean(h)) %>%
  as.numeric()

maxSimil <- dfPRF %>%
  filter(variable == "F1.seq") %>%
  top_n(1, value) %>%
  select(h) %>%
  summarise(h = mean(h)) %>%
  as.numeric()


# generate plot using ggplot2
ggplot(dfPRF, aes(h, value)) +
  geom_line(aes(color=color, linetype=lineType)) +
  geom_vline(xintercept = c(maxSimil, maxType), linetype = c(2,1)) +
  annotate("label", x = c(maxSimil, maxType), y = 0.1, label = c(maxSimil, maxType)) +
  labs(x="Cosine similarity threshold",
       y="value", color = element_blank(), linetype = element_blank()) +
  theme(legend.position = c(0, 0.55),
        legend.box.background = element_rect(
          fill = "white", size = 0.5, linetype = 1, color = 1)
        ) +
  scale_color_manual(values=c(4,3,2))
  

## ----versions------------------------------------------------------------
sessionInfo()

