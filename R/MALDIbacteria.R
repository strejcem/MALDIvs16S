# convenience tools for bacterial MS MALDI-TOF identification and classification

require("MALDIquant")
require("MALDIquantForeign")
require("XML")
require("coop")

message("v0.986 (2018-21-Mar)")

###############################################################################

parseBrukerSpectra <- function(path = ".", dimension = 10000,
                               repeatedMeasurements = FALSE) {

  ## 1) finds Bruker report files (.htm) in the specified path
  ## 2) import spectra from fid and adds taxonomy/score/replicates

  dirs <- list.dirs(path = path, recursive = F)
  message("Found ", length(dirs), " directories:") 

  ## should be faster
  parsedSpectra <- vector("list", dimension)
  ## Starting position in prepared list parsedSpectra
  pos_replicate2 <- 0 # for following replicate numbering
  pos2 <- 0 #position of imported object to the pre-created output data.frame
  okDirs <- 0 #processing statistics

  for(i in seq_along(dirs)) {
    message("Proccesing: ", i, "/", length(dirs))
    BT_report <- list.files(pattern = "\\.ht(m|ml)$", path = dirs[i], recursive = T)

    ## Check if report file exists
    if(length(BT_report) == 0) {
      message(dirs[i], " does not include Bruker report file (.htm)")
	      next()
    }

    ## Check for multiple Bruker report files
    if(length(BT_report) > 1) {
      message(dirs[i], "includes multiple Bruker report files;
                  Using: ", BT_report[1])
      BT_report <- BT_report[1]
    }

    ## Parsing HTML
    reportPath <- paste(dirs[i], BT_report, sep = "/")
                        
  	data <- readHTMLTable(reportPath, as.data.frame = F, header = T,
  	                      which = 2)[c("AnalyteID", "AnalyteName",
                          "ScoreValue", "Organism(best match)")]

    ## In some cases, HTMLparsing incorporated "\r\n +" into some fields (encoding)
  	data <- lapply(data, gsub, pattern = "\r\n +", replacement = "")


    ## Check if report file is in compatible format
  	if(is.null(data$"AnalyteName")) {
  	  message(dirs[i], "/", BT_report,
  	               " has incompatible format. Save the report as HTML!")
  	  next()
  	}

  	## processing parsed HTML & excluding empty spectra
  	## Variables with _R are from a report.
  	## Variables with _S are derived from spectra files.
    noPeaks <- data$"ScoreValue" == "< 0"
    data <- lapply(data, "[", !noPeaks)
    scoreValues_R <- as.numeric(data$ScoreValue)
    taxonomy_R <- data$`Organism(best match)`
      taxonomy_R[taxonomy_R == "not reliable identification"] <- "NA NA"
      taxSpecies_R <- strsplit(taxonomy_R, split = " ")
  	  taxGenus_R <- sapply(taxSpecies_R, "[[", 1)
  	  taxSpecies_R <- sapply(taxSpecies_R, "[[", 2)
    sampleNames_R <- data$AnalyteID
    
    #MALDIquant import changes special characters to "_"
    sampleNames_R <- gsub(pattern = "[[:punct:]]| ",
                          replace = "_", sampleNames_R, perl = T)
  	spots_R <- data$AnalyteName
  	spots_R <- sub(spots_R, pattern = " .*$", replacement = "")
  	
  	sampleSpots_R <- paste(sampleNames_R, spots_R, sep = ".")

  	s <- suppressWarnings(importBrukerFlex(path = dirs[i], verbose = F))
  	spots_S <- sapply(s, function(x)metaData(x)$spot)
  	sampleNames_S <- sapply(s, function(x)metaData(x)$sampleName)
  	sampleSpots_S <- sapply(s, function(x)metaData(x)$fullName)
  	replicates_S <- as.numeric(factor(sampleNames_S)) + pos2
  	
  	## Check if report and spectral data match
  	
  	if(!all(unique(sort(sampleSpots_S[spots_S %in% spots_R])) ==
  	      unique(sort(sampleSpots_R[spots_R %in% spots_S]))) ||
  	   length(sampleSpots_S[spots_S %in% spots_R]) == 0) {
  	  message(BT_report, " - Inconsistency between Report File and",
  	  "spectral data. Skipping directory!")
  	  next()
  	}
  	
  	## check for repeated measurements of the same spot. Bruker report shows
  	## only the best score. If "multiMeasurement = T" the best score will be
  	## assigned to multiple spectra.
  	if(length(spots_S) != length(unique(spots_S))) {
  	  
  	  if(!repeatedMeasurements) {
  	    message(dirs[i],"includes repeated measurements of the same spot(s).",
  	                "Skipping directory!",
  	                "If you want to assign the same (best) score to multiple",
  	                "spectra, use 'repeatedMeasurements = TRUE'")
  	    next()
  	  }
  	  repeated_S <- sampleSpots_S[duplicated(sampleSpots_S)]
  	  message(dirs[i], "contains repeated measurements of",
  	              paste(repeated_S, collapse = ", "),
  	              "samples! Scores will not reflect reality")
  	}

  	pos_replicate <- as.numeric(factor(sampleNames_S)) + pos_replicate2

	  for(y in seq_along(s)) {
	    metaData(s[[y]])$taxonomy <- taxonomy_R[which(spots_R == metaData(s[[y]])$spot)]
	    metaData(s[[y]])$taxGenus <- taxGenus_R[which(spots_R == metaData(s[[y]])$spot)]
	    metaData(s[[y]])$taxSpecies <- taxSpecies_R[which(spots_R == metaData(s[[y]])$spot)]
	    metaData(s[[y]])$scoreValue <- scoreValues_R[which(spots_R == metaData(s[[y]])$spot)]
      metaData(s[[y]])$replicate <- pos_replicate[y]
      metaData(s[[y]])$treatment <- sub(".*/", "", dirs[i], fixed = F)
	  }

  	## putting data in the final list
  	pos1 <- pos2 + 1
  	pos2 <- pos1 + length(s) - 1
  	parsedSpectra[pos1:pos2] <- s
  	
  	pos_replicate2 <- max(pos_replicate)
  	
  	## count succesfully parsed directories
  	okDirs <- okDirs + 1
  }

  ## Remove empty entries
  parsedSpectra[sapply(parsedSpectra, is.null)] <- NULL
  message("******************************************")
  message("done")
  message("******************************************")
  
  message("Succesfully proccessed directories: ", okDirs, "/",
               length(dirs)," (", length(parsedSpectra), " spectra)" )
  if(okDirs == 0) {
    message("!!! No data were proccessed !!!")
  }
  
  return(parsedSpectra)
}

###############################################################################

classicMALDI <- function(x, dist.method = "cos", range = c(4000, 10000),
                         tol = 0.002, SNR = 3, smooth.hw = 20, peak.hw = 20,
                         calib.method = "TIC", transMethod = "sqrt",
                         detPeak.method = "SuperSmoother", clust.method = "average",
                         labels = "sampleName", base.method = "SNIP", iter = 50,
                         end = c("cluster", "distance", "feature", "binpeaks",
                                 "peaks"), averageReps = NULL, alignment = F,
                         filterPeaks = NULL) {

  
  ## classic cluster analysis of spectra by calculating all-to-all correlation
  ## coefficients or cosines

  ## spectra pre-processing
  end <- match.arg(end, c("cluster", "distance", "feature", "binpeaks", "peaks"))
  
  spectra <- x
  if(transMethod != "none") {
    spectra <- transformIntensity(x, method=transMethod)
  }
  spectra <- MALDIquant::trim(spectra, range=range)
  spectra <- smoothIntensity(spectra, method="SavitzkyGolay",
                             halfWindowSize=smooth.hw)
  if(base.method == "SNIP") {
    spectra <- removeBaseline(spectra, base.method, iterations=iter)
  } else if (base.method == "TopHat") {
    spectra <- removeBaseline(spectra, base.method, halfWindowSize=iter)
  } else if(base.method == "ConvexHull") {
    spectra <- removeBaseline(spectra, base.method)
  } else if(base.method == "median") {
    spectra <- removeBaseline(spectra, base.method, halfWindowSize=iter)
  } else stop("Incorrect baseline corrction method!")
 
  ## spectral aligment
  spectra <- calibrateIntensity(spectra, method=calib.method)
  if(alignment) spectra <- alignSpectra(spectra)
  
  ## average
  if(!is.null(averageReps)) {
    labs <-
      sapply(lapply(spectra, function(x)metaData(x)[averageReps]), paste,
             collapse = " ")
    labs <- as.factor(labs)
    spectra <- averageMassSpectra(spectra, labels = labs, method = "sum")
  }
  
  ## run peak detection
  peaks <- detectPeaks(spectra, method=detPeak.method, halfWindowSize=peak.hw,
                       SNR=SNR)

  ## stop if only peak calling is desired
  if (end == "peaks") return(peaks)
  
  ## bin peaks
  peaks <- binPeaks(peaks, tolerance=tol, method = "strict")
  if (!is.null(filterPeaks)) {
    if (filterPeaks < 1 & filterPeaks > 0) {
    peaks <- filterPeaks(peaks, minFrequency = filterPeaks, 
                         labels = metaParam(x, "replicate"),
                         mergeWhitelists = T)
    } else stop("'filterPeaks' value is not in the 0-1 range!")
  } 

  
  ## stop if only binned peaks are desired
  if (end == "binpeaks") return(peaks)
  
  featureMatrix <- intensityMatrix(peaks, spectra)
  #featureMatrix[is.na(featureMatrix)] <-0

  rownames(featureMatrix) <- 
    sapply(lapply(spectra, function(x)metaData(x)[labels]), paste,
           collapse = " ")
  
  ## stop if only feature matrix is desired
  if (end == "feature") return(list("fm"=featureMatrix))
  
  ## faster implementation coop::pcor
  if(dist.method == "cor") {
    similarity <- coop::pcor(t(featureMatrix))
    d <- as.dist(sqrt(1-similarity))
  } else if(dist.method == "cos") {
    similarity <- coop::cosine(t(featureMatrix))
    d <- as.dist(1-similarity)
  } else if(dist.method == "euc") {
    d <- dist(featureMatrix)
  }
  
  #"quality" of each measurement, mean cosine of a point with all other replicates
  reps <- sapply(spectra, function(x) metaData(x)$replicate)
  if(is.null(averageReps) & class(reps) != "list") {
    qual <- sapply(seq_along(reps), function(idx) {
      if(sum(reps == reps[idx]) > 1){
        repIdx <- which(reps == reps[idx])
        otherIdx <- repIdx[!repIdx %in% idx]
        mean(similarity[idx, otherIdx])
      } else 0
    })
    
  } else qual <- NULL
    # number of peaks
    numPeaks <- sapply(peaks, function(x) sum(intensity(x)>0))
    
  ## stop if only distance matrix is desired
  if (end == "distance") return(list("d"=d, "fm"=featureMatrix, "ACSreps"=qual,
                                     "numPeaks" = numPeaks))
  
  hc <- hclust(d, method=clust.method)

  return(list("hc"=hc, "d"=d, "fm"=featureMatrix, "ACSreps"=qual,
              "numPeaks" = numPeaks))
}

###############################################################################

## PDF dendrogram

PDFplot <- function(hc, cutoffs = 0.17, color = NULL, file = NULL,
                    scaleToPrint = T, main = NULL) {
  
  if (class(hc) != "hclust") stop("Object is not of 'hclust' class")
  if (is.null(file)) file <-
      paste(format(Sys.time(), "%y%m%d-%H%M%S"), "pdf", sep = ".")
  
  ps <- ifelse(length(hc$labels) >= 50, length(hc$labels)/50, 1) # plot scale
  if(!scaleToPrint) ps <- 1
  
  pdfWidth <- if(11.68*ps > 200) 200 else 11.68*ps
  
  if(is.null(main)) main <- Sys.Date()
  
  #open pdf device
  pdf(file = file, width = pdfWidth, height = 7, pointsize = 8)

  plot(hc, hang = -1, main = main)
  clusText <- NULL
  textPosY <- NULL
  if(!is.null(cutoffs)) {
    for (idx in seq_along(cutoffs)) {
      colIdx <- ifelse(is.null(color), idx+1, color[idx]) 
      cluster <- cutree(hc, h = cutoffs[idx])[hc$order]
      clusNum <- sum(table(cluster) > 0, na.rm = TRUE)
      rect.hclust(hc, h = cutoffs[idx], border = colIdx)
      clusText <- c(clusText, paste("similarity:", 1-cutoffs[idx], " clusters:", clusNum))
      textPosY <- c(textPosY, max(hc$height) - 0.05*idx)
    }
    
    text(x = length(hc$labels), y = textPosY, labels = clusText, adj = 1,
         col = seq_along(cutoffs)+1)
  }
  # close pdf
  dev.off()
}

###############################################################################


## remove peaks

removePeaks <- function(x, remove) {
  
  if(class(x) != "MassPeaks") stop(
    "Object is not of 'MassPeaks{MALDIquant}' class")
  
  if(length(remove) == 0) warning("Zero peaks removed!")
  if(length(remove) > 0 ){
    
    #stop("Zero peaks to remove!") 
  
    x@snr <- x@snr[-remove]
    x@mass <- x@mass[-remove]
    x@intensity <- x@intensity[-remove]
  }
  x
}

###############################################################################

## get peaks

getPeaks <- function(x, peaks) {
  if(class(x) != "MassPeaks") stop(
    "Object is not of 'MassPeaks{MALDIquant}' class")
  
  if(length(remove) == 0) warning("Zero peaks selected!")
  if(length(remove) > 0 ){
    x@snr <- x@snr[peaks]
    x@mass <- x@mass[peaks]
    x@intensity <- x@intensity[peaks]
  }
  x
}

###############################################################################
## Selects top n peaks

topNPeaks <- function(x, n = 100) {
  
  if(class(x) != "MassPeaks") stop(
    "Object is not of 'MassPeaks{MALDIquant}' class")
  
  nPeaks <- ifelse(n > length(x@mass), length(x@mass), n)
  x[order(x@intensity, decreasing = T)][1:nPeaks]
}



###############################################################################
## Change meta data parametr 

chngMetaParam <- function(x, param, values) {

    f <- function(x, param, values) {
    metaData(x)[[param]] <- values
    x
  }
  
  if(class(x) == "list") {
    if(!inherits(x[[1]], c("MassPeaks", "MassSpectrum"))) {
      stop(
        "'x' must be an object (list) of class 'MassPeaks' or 'MassSpectrum'")
    }
    
    if(length(x) == length(values)) {
      mapply(f, x, values, MoreArgs = list(param=param))
    } else if(length(values) == 1 | is.null(values)) {
      sapply(x, f, param, values)
    } else stop("Length of 'values' must be the same as 'x' or 1")
  } else {
    if(!inherits(x, c("MassPeaks", "MassSpectrum"))) {
      stop("'x' must be an object (list) of class 'MassPeaks' or 'MassSpectrum'")
    }
    f(x, param, values)
  }
}

###############################################################################
## print meta data parametr 

metaParam <- function(x, param) {
    sapply(x, function(x) metaData(x)[[param]])
  }


###############################################################################
## averageCosine inside groups

avgCos <- function(fm, groups) {
  stopifnot(nrow(fm) == length(groups))
  if(!any(fm > 0)) return(NA)
  groupsU <- unique(groups)
  similarity <- coop::cosine(t(fm))

  output <- lapply(seq_along(groupsU), function(idx) {
    if(sum(groups == groupsU[idx]) > 1){
      groupIdx <- which(groups == groupsU[idx])
      
      groupDist <- as.dist(similarity[groupIdx, groupIdx])
      data.frame(sample = as.character(groupsU[idx]),
                 meanACS = mean(groupDist, na.rm = TRUE), 
                 sdACS = sd(groupDist, na.rm = TRUE))
    } else data.frame(sample = as.character(groupsU), meanACS = NA, sdACS = NA)
  })

  return(do.call(rbind, output))
}




###############################################################################
## averageEucDist inside groups

avgEucDist <- function(fm, groups) {
  stopifnot(nrow(fm) == length(groups))
  groupsU <- unique(groups)
  distance <- as.matrix(dist(fm))
  output <- sapply(seq_along(groupsU), function(idx) {
    if(sum(groups == groupsU[idx]) > 1){
      groupIdx <- which(groups == groupsU[idx])
      
      mean(as.dist(distance[groupIdx, groupIdx]),na.rm = T)
    } else NA
  })
  names(output) <- as.character(groupsU)
  return(output)
}

###############################################################################
## averagedist between groups

avgDistBetween <- function(d, groups) {
  d <- as.matrix(d)
  diag(d) <- NA
  stopifnot(dim(d) == length(groups))
  groupsU <- unique(groups)
  output <- matrix(data = 0, length(groupsU), length(groupsU), dimnames = list(groupsU, groupsU))
  validM <- lower.tri(output, diag = T)
    for (idx in seq_along(groupsU)) {
      for (idy in seq_along(groupsU)) {
        if(!validM[idx,idy]) next
        group1 <- which(groups == groupsU[idx])
        group2 <- which(groups == groupsU[idy])
        
        output[idx,idy] <- mean(d[group1, group2], na.rm=T)
      }
      output[upper.tri(output)] <- output[lower.tri(output)]
  }
  return(output)
}

###############################################################################
## Calculates Precision, Recall and F1 measure bsased on Kim et al. 2014
###############################################################################

FscPrecRecAcc <- function(pred, ref, beta = 1) {
  stopifnot(length(pred)==length(ref))
  TP <- as.numeric(sum(pred & ref))
  TN <- as.numeric(sum(!pred & !ref))
  FP <- as.numeric(sum(pred & !ref))
  FN <- as.numeric(sum(!pred & ref))
  prec <- TP/(TP+FP)
  rec <- TP/(TP+FN)
  acc <- (TP+TN)/(TP+TN+FP+FN)
  Fscore <- (1+beta^2)*prec*rec/(beta^2*prec+rec)
  MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  FPR <- FP/(FP+TN)
  FDR <- FP/(FP+TP)

  return(data.frame(prec=prec, rec=rec, Fscore=Fscore,
                    Acc=acc,MCC=MCC, FPR=FPR, FDR=FDR))
}