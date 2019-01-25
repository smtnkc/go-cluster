source("vars.R")

##################################################

readRaws <- function(idTypes) {
  raws <- list()  
  for(i in idTypes) {
    fname <- paste(i, "_RAW.csv", sep="")
    raws[[i]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
  }
  return(raws)
}

raws <- readRaws(idTypes)

readLinks <- function(idTypes, topologies) {
  links <- list()
  for(i in idTypes) {
    for(t in topologies) {
      fname <- paste("LINKS/", i, "links_", t, ".csv", sep="")
      cat("Reading", fname, "...\n")
      links[[i]][[t]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
      
      # convert integer entrezID columns to character
      if(i == "e") {
        links[[i]][[t]][, 1] <- as.character(links[[i]][[t]][, 1])
        links[[i]][[t]][, 2] <- as.character(links[[i]][[t]][, 2])
      }
    }
  }
  return(links)
}

links <- readLinks(idTypes, topologies)

getNodes <- function(links) {
  nodes <- list()
  for(i in names(links)) {
    for(t in names(links[[i]])) {
      nodes[[i]][[t]] <- union(links[[i]][[t]][,1], links[[i]][[t]][,2])
    }
  }
  return(nodes)
}

nodes <- getNodes(links)

################################################

readDegs <- function(idTypes, FC, pval, subjects) {
  degs <- list()
  
  for(i in idTypes) {
    for(s in subjects) {
      fname <- paste("DEGs/", i, "_", s, "_fc", Stringify(FC), 
                     "_p", Stringify(pval), ".csv", sep="")
      cat("Reading", fname, "...\n")
      degs[[i]][[s]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))[,1]
    }
  }
  return(degs)
}

degs <- readDegs(idTypes, FC, P_VAL, subjects)

########## get raw vals for degs:

getDegVals <- function(degs, raws, intervals) {
  degVals <- list()
  for(i in names(degs)) {
    for(s in names(degs[[i]])) {
      tmp <- raws[[i]][raws[[i]][, 1] %in% degs[[i]][[s]], intervals[[s]]]
      row.names(tmp) <- tmp[, 1] # set rownames as first column
      tmp <- tmp[, -1] # remove first column
      degVals[[i]][[s]] <- tmp
    }
  }
  return(degVals)
}

degVals <- getDegVals(degs, raws, intervals)

####################################### Filter links

FilterLinks <- function(df.links, filterIds, cutoff) {
  cat("For cutoff =", cutoff,"\n")
  # Omit links of unmapped genes
  df.result <- df.links
  all <- nrow(df.result)
  cat("Total links =", all,"\n")
  
  df.result <- df.result[df.result[, 1] %in% filterIds, ]
  df.result <- df.result[df.result[, 2] %in% filterIds, ]
  mapped <- nrow(df.result)
  cat("Mapped links =", mapped, "   ", round(mapped*100/all, 1), "%\n")
  
  # Omit insignificant links
  if(!missing(cutoff) && !is.null(cutoff))
    df.result <- df.result[df.result[,3] >= cutoff, ]
  sig <- nrow(df.result)
  cat("Significant links =", sig, "   ", round(sig*100/all,1), "%\n\n")
  
  row.names(df.result) <- NULL
  return(df.result)
}

getMsLinks <- function(links, raws, cutoff) {
  ms_links <- list()
  for(i in names(links)) {
    for(t in names(links[[i]])) {
      ms_links[[i]][[t]] <- FilterLinks(links[[i]][[t]], raws[[i]][,1], cutoff[[i]][[t]])
    }
  }
  return(ms_links)
}

ms_links <- getMsLinks(links, raws, cutoffs[[TOP_N]])

ms_nodes <- getNodes(ms_links)

################# MSDEGS

getMsDegs <- function(idTypes, topologies, subjects, ms_nodes) {
  ms_degs <- list()
  for(i in idTypes) {
    for(t in topologies) {
      for(s in subjects) {
        D <- degVals[[i]][[s]]
        ms_degs[[i]][[t]][[s]] <- D[row.names(D) %in% ms_nodes[[i]][[t]], ] 
      }
    }
  }
  return(ms_degs)
}

ms_degs <- getMsDegs(idTypes, topologies, subjects, ms_nodes)

######### get msdeglinks:

GetMSDegLinks <- function(ms_links, ms_degs) {
  ms_deglinks <- list()
  for(i in names(ms_degs)) {
    for(t in names(ms_degs[[i]])) {
      for(s in names(ms_degs[[i]][[t]])) {
        L <- ms_links[[i]][[t]]
        D <- ms_degs[[i]][[t]][[s]]
        tmp <- L[L[,1] %in% row.names(D), 1:2]
        tmp <- tmp[tmp[,2] %in% row.names(D), ]
        ms_deglinks[[i]][[t]][[s]] <- tmp
      }
    }
  }
  return(ms_deglinks)
}

ms_deglinks <- GetMSDegLinks(ms_links, ms_degs)

######### get msdegnodes:

GetMSDegNodes <- function(ms_deglinks) {
  ms_degnodes <- list()
  for(i in names(ms_deglinks)) {
    for(t in names(ms_deglinks[[i]])) {
      for(s in names(ms_deglinks[[i]][[t]])) {
        df <- ms_deglinks[[i]][[t]][[s]]
        genes <- unique(as.character(rbind(df[, 1], df[, 2])))
        ms_degnodes[[i]][[t]][[s]] <- genes
      }
    }
  }
  return(ms_degnodes)
}

ms_degnodes <- GetMSDegNodes(ms_deglinks)

######### write links:

WriteLinks <- function(ms_deglinks, p, top) {
  for(i in names(ms_deglinks)) {
    for(t in names(ms_deglinks[[i]])) {
      for(s in names(ms_deglinks[[i]][[t]])) {
        out_file = paste("RES/", i, "_", t, "_", s, "_P", Stringify(p), "_T", top,".csv", sep="")
        write.csv(ms_deglinks[[i]][[t]][[s]], out_file, row.names = FALSE)
      }
    }
  }
}

WriteLinks(ms_deglinks, p=P_VAL, top=TOP_N)
