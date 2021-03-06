source("vars.R")

getSubsets <- function(raw, intervals_ctrl) {
  subsets <- list()
  for(s in names(intervals_ctrl)) {
    df <- raw[, intervals_ctrl[[s]]]
    row.names(df) <- df[, 1]
    df <- df[,-1]
    for(c in 1:ncol(df)) {
      # convert 0s to smallest non-zero floating-point number (to be able to get log2)
      df[df[, c]==0, c] <- .Machine$double.xmin
    }
    subsets[[s]] <- log2(df)
  }
  return(subsets)
}

subsets <- getSubsets(raw, intervals_ctrl)

######################## ORDINARY-T METHOD:

getSigGenes <- function(subsets, pval) {
  pval <- as.numeric(pval)
  sigs <- list()
  for(s in names(subsets)) {
    df <- subsets[[s]]
    sigGenes <- c()
    invalidRows <- c()
    for(r in 1:nrow(df)) {
      controlRange <- c(1:9)
      caseRange <- c(10:ncol(df))
      if(length(levels(factor(as.numeric(df[r, ])))) == 1) {
        #### eliminate the constant rows to overcome "t-test data are essentially constant" error
        invalidRows <- c(invalidRows, row.names(df)[r])
      }
      else if(t.test(df[r, controlRange], df[r, caseRange])$p.val <= pval) {
        sigGenes <- c(sigGenes, row.names(df)[r])
      }
    }
    cat("There are", length(invalidRows), "invalid rows:", invalidRows, "\n")
    sigs[[s]] <- sigGenes
  }
  return(sigs)
}

sigs <- getSigGenes(subsets, P_VAL)

getSigDFs <- function(subsets, sigs) {
  sigDFs <- list()
  for(s in names(subsets)) {
    df <- subsets[[s]]
    sigDFs[[s]] <- subset(df, row.names(df) %in% sigs[[s]])
  }
  return(sigDFs)
}

sigDFs <- getSigDFs(subsets, sigs)

###################################### FC cutoff:

getDEGs <- function(sigDFs, FC) {
  degs <- list()
  for(s in names(sigDFs)) {
    df <- sigDFs[[s]]
    controlRange <- c(1:9)
    caseRange <- c(10:ncol(df))
    tmpDegs <- c()
    for(r in 1:nrow(df)) {
      controlMean <- apply(df[r, controlRange], 1, mean)
      caseMean <- apply(df[r, caseRange], 1, mean)
      if(abs(controlMean-caseMean) >= FC) {
        tmpDegs <- c(tmpDegs, row.names(df)[r])
      }
    }
    degs[[s]] <- subset(df, row.names(df) %in% tmpDegs)
  }
  return(degs)
}

degs <- getDEGs(sigDFs, FC)

##################################### WRITE TO FILE:

writeDEGs <- function(degs, FC, P_VAL) {
  for(s in names(degs)) {
    df <- data.frame("Symbol"=rownames(degs[[s]]))
    fname <- paste("DEGs/", s, "_fc", Stringify(FC),
                   "_p", Stringify(P_VAL), ".csv", sep="")

    write.table(df, fname, row.names = FALSE)
    cat("Writing", fname, "\n")
  }
  cat("ALL DONE!\n")
}

writeDEGs(degs, FC, P_VAL)

######################## TREAT OR MODERATED-T METHOD:
# library(limma)
#
# designs <- list()
# designs[["mets"]] <- cbind(Intercept=1,Group=c(rep(0,9),rep(1,(length(subsets[["mets"]])-9))))
# designs[["t2d"]]  <- cbind(Intercept=1,Group=c(rep(0,9),rep(1,(length(subsets[["t2d"]])-9))))
# designs[["cad"]]  <- cbind(Intercept=1,Group=c(rep(0,9),rep(1,(length(subsets[["cad"]])-9))))
#
# fit.mets <- lmFit(subsets[["mets"]], designs[["mets"]])
# fit.t2d  <- lmFit(subsets[["t2d"]], designs[["t2d"]])
# fit.cad  <- lmFit(subsets[["cad"]], designs[["cad"]])

######################## TREAT METHOD:
# tfit.mets <- treat(fit.mets)
# tt.mets <- topTreat(tfit.mets, coef=2, number = nrow(df.mets))
# deg.mets <- tt.mets[which(abs(tt.mets$logFC) >= FC & tt.mets$adj.P.Val <= P_VAL),]
#
# tfit.cad <- treat(fit.cad)
# tt.cad <- topTreat(tfit.cad, coef=2, number = nrow(df.cad))
# deg.cad <- tt.cad[which(abs(tt.cad$logFC) >= FC & tt.cad$adj.P.Val <= P_VAL),]
#
# tfit.t2d <- treat(fit.t2d)
# tt.t2d <- topTreat(tfit.t2d, coef=2, number = nrow(df.t2d))
# deg.t2d <- tt.t2d[which(abs(tt.t2d$logFC) >= FC & tt.t2d$adj.P.Val <= P_VAL),]

######################## MODERATED-T METHOD:
# bfit.mets <- eBayes(fit.mets)
# tb.mets <- topTable(bfit.mets, coef=2, number = nrow(subsets[["mets"]]))
# deg.mets <- tb.mets[which(abs(tb.mets$logFC) >= FC & tb.mets$adj.P.Val <= P_VAL),]
#
# bfit.t2d <- eBayes(fit.t2d)
# tb.t2d <- topTable(bfit.t2d, coef=2, number = nrow(subsets[["cad"]]))
# deg.t2d <- tb.t2d[which(abs(tb.t2d$logFC) >= FC & tb.t2d$adj.P.Val <= P_VAL),]
#
# bfit.cad <- eBayes(fit.cad)
# tb.cad <- topTable(bfit.cad, coef=2, number = nrow(subsets[["t2d"]]))
# deg.cad <- tb.cad[which(abs(tb.cad$logFC) >= FC & tb.cad$adj.P.Val <= P_VAL),]
