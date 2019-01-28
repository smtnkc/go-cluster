source("vars.R")

readSlinks <- function(topologies) {
  slinks <- list()
  for(t in topologies) {
    fname <- paste("LINKS/slinks_", t, ".csv", sep="")
    cat("Reading", fname, "...\n")
    slinks[[t]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
  }
  return(slinks)
}

slinks <- readSlinks(topologies)

################################################

getNodes <- function(slinks) {
  nodes <- list()
  for(t in names(slinks)) {
    nodes[[t]] <- union(slinks[[t]][, 1], slinks[[t]][, 2])
  }
  return(nodes)
}

nodes <- getNodes(slinks)

################################################

readDegs <- function(FC, pval, subjects) {
  degs <- list()
  for(s in subjects) {
    fname <- paste("DEGs/", s, "_fc", Stringify(FC), 
                   "_p", Stringify(pval), ".csv", sep = "")
    cat("Reading", fname, "...\n")
    degs[[s]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))[,1]
  }
  return(degs)
}

degs <- readDegs(FC, P_VAL, subjects)

########## get raw vals for degs:

getDegVals <- function(degs, raw, intervals) {
  degVals <- list()
  for(s in names(degs)) {
    df <- raw[raw[, 1] %in% degs[[s]], intervals[[s]]]
    row.names(df) <- df[, 1] # set rownames as first column
    df <- df[, -1] # remove first column
    degVals[[s]] <- df
  }
  return(degVals)
}

degVals <- getDegVals(degs, raw, intervals)

####################################### Filter links

FilterLinks <- function(dfLinks, filterIds, cutoff) {
  cat("For cutoff =", cutoff,"\n")
  # Omit links of unmapped genes
  dfTmp <- dfLinks
  nAll <- nrow(dfTmp)
  cat("Total links =", nAll,"\n")
  
  dfTmp <- dfTmp[dfTmp[, 1] %in% filterIds, ]
  dfTmp <- dfTmp[dfTmp[, 2] %in% filterIds, ]
  nMapped <- nrow(dfTmp)
  cat("Mapped links =", nMapped, "   ", round(nMapped*100/nAll, 1), "%\n")
  
  # Omit insignificant links
  if(!missing(cutoff) && !is.null(cutoff))
    dfTmp <- dfTmp[dfTmp[, 3] >= cutoff, ]
  nSig <- nrow(dfTmp)
  cat("Significant links =", nSig, "   ", round(nSig*100/nAll,1), "%\n\n")
  
  row.names(dfTmp) <- NULL
  return(dfTmp)
}

getMsLinks <- function(slinks, raw, cutoff) {
  msLinks <- list()
  for(t in names(slinks)) {
    msLinks[[t]] <- FilterLinks(slinks[[t]], raw[, 1], cutoff[[t]])
  }
  return(msLinks)
}

msLinks <- getMsLinks(slinks, raw, cutoffs[[TOP_N]])

msNodes <- getNodes(msLinks)

################# MSDEGS

getMsDegVals <- function(msNodes, degVals) {
  msDegs <- list()
  for(t in names(msNodes)) {
    for(s in names(degVals)) {
      DV <- degVals[[s]]
      msDegs[[t]][[s]] <- DV[row.names(DV) %in% msNodes[[t]], ] 
    }
  }
  return(msDegs)
}

msDegVals <- getMsDegVals(msNodes, degVals)

######### get msdeglinks:

getMsDegLinks <- function(msLinks, msDegVals) {
  msDegLinks <- list()
  for(t in names(msDegVals)) {
    for(s in names(msDegVals[[t]])) {
      L <- msLinks[[t]]
      DV <- msDegVals[[t]][[s]]
      df <- L[L[,1] %in% row.names(DV), 1:2]
      df <- df[df[,2] %in% row.names(DV), ]
      msDegLinks[[t]][[s]] <- df
    }
  }
  return(msDegLinks)
}

msDegLinks <- getMsDegLinks(msLinks, msDegVals)

######### get msdegnodes:

getMsDegNodes <- function(msDegLinks) {
  msDegNodes <- list()
  for(t in names(msDegLinks)) {
    for(s in names(msDegLinks[[t]])) {
      df <- msDegLinks[[t]][[s]]
      genes <- unique(as.character(rbind(df[, 1], df[, 2])))
      msDegNodes[[t]][[s]] <- genes
    }
  }
  return(msDegNodes)
}

msDegNodes <- getMsDegNodes(msDegLinks)

######### write links:

WriteLinks <- function(msDegLinks, p, top) {
  for(t in names(msDegLinks)) {
    for(s in names(msDegLinks[[t]])) {
      outname = paste("RES/", t, "_", s, "_P",
                       Stringify(p), "_T", top, ".csv", sep = "")
      write.csv(msDegLinks[[t]][[s]], outname, row.names = FALSE)
    }
  }
}

WriteLinks(msDegLinks, p=P_VAL, top=TOP_N)
