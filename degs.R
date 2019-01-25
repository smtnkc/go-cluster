source("vars.R")

########################################################

createEntrezRaw <- function(s_raw) {
  e_raw <- s_raw
  map_raw <- as.data.frame(fread("map_raw_unique.csv", header = TRUE, sep = ','))
  e_raw$Entrez <- map_raw$ENTREZID
  e_raw <- na.omit(e_raw)
  row.names(e_raw) <- NULL
  e_raw[,1] <- e_raw[,37]
  colnames(e_raw)[1] <- "Entrez"
  e_raw <- e_raw[,1:36]
  e_raw[,1] <- as.character(e_raw[,1])
  write.csv(e_raw, "e_RAW.csv", row.names = FALSE)
  return(e_raw)
}

readRaws <- function(idTypes) {
  raws <- list()  
  for(i in idTypes) {
    fname <- paste(i, "_RAW.csv", sep="")
    
    if(!file.exists(fname)) {
      if(i == "s") {
        cat("ERROR:", fname, "does not exist!\n")
        return(NULL)
      } else {
        cat("Creating", fname, "...\n")
        raws[[i]] <- createEntrezRaw(raws[["s"]])
        cat("DONE!\n")
      }
    } else {
      cat("Reading", fname, "...\n")
      raws[[i]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
    }
  }
  return(raws)
}

raws <- readRaws(idTypes)

########################################################

getSubsets <- function(raws, intervals_ctrl) {
  subsets <- list()
  for(i in names(raws)) {
    for(s in names(intervals_ctrl)) {
      tmp <- raws[[i]][, intervals_ctrl[[s]]]
      row.names(tmp) <- tmp[,1]
      tmp <- tmp[,-1]
      for(c in 1:ncol(tmp)) {
        # convert 0s to smallest non-zero floating-point number (to be able to get log2)
        tmp[tmp[, c]==0, c] <- .Machine$double.xmin
      }
      subsets[[i]][[s]] <- log2(tmp)
    }
  }
  return(subsets)
}

subsets <- getSubsets(raws, intervals_ctrl)

######################## ORDINARY-T METHOD:

getSigGenes <- function(subsets, pval) {
  pval <- as.numeric(pval)
  sigs <- list()
  for(i in names(subsets)) {
    for(s in names(subsets[[i]])) {
      df <- subsets[[i]][[s]]
      SigGenes <- c()
      invalidRows <- c()
      for(r in 1:nrow(df)) {
        controlRange <- c(1:9)
        caseRange <- c(10:ncol(df))
        if(length(levels(factor(as.numeric(df[r, ])))) == 1) {
          #### eliminate the constant rows to overcome "t-test data are essentially constant" error
          invalidRows <- c(invalidRows, row.names(df)[r])
        }
        else if(t.test(df[r, controlRange], df[r, caseRange])$p.val <= pval) {
          SigGenes <- c(SigGenes, row.names(df)[r])
        }
      }
      cat("There are", length(invalidRows), "invalid rows:", invalidRows, "\n")
      sigs[[i]][[s]] <- SigGenes
    }
  }
  return(sigs)
}

sigs <- getSigGenes(subsets, P_VAL)

getSigDFs <- function(subsets, sigs) {
  sigDFs <- list()
  for(i in names(subsets)) {
    for(s in names(subsets[[i]])) {
      df <- subsets[[i]][[s]]
      sigDFs[[i]][[s]] <- subset(df, row.names(df) %in% sigs[[i]][[s]])  
    }
  }
  return(sigDFs)
}

sigDFs <- getSigDFs(subsets, sigs)

###################################### FC cutoff:

getDEGs <- function(sigDFs, FC) {
  degs <- list()
  for(i in names(sigDFs)) {
    for(s in names(sigDFs[[i]])) {
      df <- sigDFs[[i]][[s]]
      controlRange <- c(1:9)
      caseRange <- c(10:ncol(df))
      tmp_degs <- c()
      for(r in 1:nrow(df)) {
        controlMean <- apply(df[r, controlRange], 1, mean)
        caseMean <- apply(df[r, caseRange], 1, mean)
        if(abs(controlMean-caseMean) >= FC) {
          tmp_degs <- c(tmp_degs, row.names(df)[r])
        }
      }
      degs[[i]][[s]] <- subset(df, row.names(df) %in% tmp_degs)
    }
  }
  return(degs)
}

degs <- getDEGs(sigDFs, FC)

##################################### WRITE TO FILE:

writeDEGs <- function(degs, FC, P_VAL) {
  for(i in names(degs)) {
    for(s in names(degs[[i]])) {
      if(i == "s") id <- "Symbol" else id <- "Entrez"
      df <- data.frame(id=rownames(degs[[i]][[s]]))
      fname <- paste("DEGs/", i, "_", s, "_fc", Stringify(FC),
                     "_p", Stringify(P_VAL), ".csv", sep="")

      write.table(df, fname, row.names = FALSE)
      cat("Writing", fname, "\n")
    }
  }
  cat("ALL DONE!\n")
}

writeDEGs(degs, FC, P_VAL)




######################## TREAT OR MODERATED-T METHOD:
# design.mets <- cbind(Intercept=1,Group=c(rep(0,9),rep(1,(length(df.mets)-9))))
# design.cad  <- cbind(Intercept=1,Group=c(rep(0,9),rep(1,(length(df.cad)-9))))
# design.t2d  <- cbind(Intercept=1,Group=c(rep(0,9),rep(1,(length(df.t2d)-9))))
# 
# fit.mets <- lmFit(df.mets, design.mets)
# fit.cad  <- lmFit(df.cad, design.cad)
# fit.t2d  <- lmFit(df.t2d, design.t2d)

######################## TREAT METHOD:
# tfit.mets <- treat(fit.mets)
# tt.mets <- topTreat(tfit.mets, coef=2, number = nrow(df.mets))
# deg.mets <- tt.mets[which(abs(tt.mets$logFC) >= 1 & tt.mets$adj.P.Val <= 0.05),]
# 
# tfit.cad <- treat(fit.cad)
# tt.cad <- topTreat(tfit.cad, coef=2, number = nrow(df.cad))
# deg.cad <- tt.cad[which(abs(tt.cad$logFC) >= 1 & tt.cad$adj.P.Val <= 0.05),]
# 
# tfit.t2d <- treat(fit.t2d)
# tt.t2d <- topTreat(tfit.t2d, coef=2, number = nrow(df.t2d))
# deg.t2d <- tt.t2d[which(abs(tt.t2d$logFC) >= 1 & tt.t2d$adj.P.Val <= 0.05),]

######################## MODERATED-T METHOD:
# bfit.mets <- eBayes(fit.mets)
# tb.mets <- topTable(bfit.mets, coef=2, number = nrow(df.mets))
# deg.mets <- tb.mets[which(abs(tb.mets$logFC) >= 1 & tb.mets$adj.P.Val <= 0.05),]
# 
# bfit.cad <- eBayes(fit.cad)
# tb.cad <- topTable(bfit.cad, coef=2, number = nrow(df.cad))
# deg.cad <- tb.cad[which(abs(tb.cad$logFC) >= 1 & tb.cad$adj.P.Val <= 0.05),]
# 
# bfit.t2d <- eBayes(fit.t2d)
# tb.t2d <- topTable(bfit.t2d, coef=2, number = nrow(df.t2d))
# deg.t2d <- tb.t2d[which(abs(tb.t2d$logFC) >= 1 & tb.t2d$adj.P.Val <= 0.05),]
