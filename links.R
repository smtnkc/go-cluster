source("vars.R")

##### GETMAPS: Protein to Symbol for STRING, Entrez to SYMBOL for INET

getMaps <- function(topologies) {
  maps <- list()
  for(t in topologies) {
    fname <- paste("map_", t, ".csv", sep="")
    df <- as.data.frame(fread(fname, header = TRUE, sep = ','))
    if(t == "inet") {
      df <- transform(df, ENTREZID = as.character(ENTREZID))
    }
    for(r in 1:nrow(df)) {
      maps[[t]][df[r, 1]] <- df[r, 2]
    }
  }
  maps <- lapply(maps, function(x) x[!is.na(x)]) # remove NA values
  return(maps)
}

maps <- getMaps(topologies)

##### Note that elinks_inet and plinks_string are downloaded files.

removeReverseDuplicates <- function(plinksDup) {
  plinks <- plinksDup
  start_time <- proc.time()
  #plinks <- plinks[!duplicated(apply(plinks, 1, function(x) paste(sort(x), collapse='_'))), ]
  plinks <- plinks[plinks[,1] < plinks[,2], ]
  elapsed_time <- proc.time() - start_time
  cat("UNIQUED IN", elapsed_time[[3]], "SECONDS\n")
  fname <- "LINKS/plinks_string.csv"
  cat("WRITING", fname, "...\n")
  write.csv(plinks, fname, row.names = FALSE)
  return(plinks)
}

readPlinks <- function() {
  fname <- "LINKS/plinks_string.csv"
  if(file.exists(fname)) {
    cat("READING", fname, "...\n")
    plinks <- as.data.frame(fread(fname, header = TRUE, sep = ','))
  } else {
    cat("CREATING", fname, "...\n")
    plinksDup <- as.data.frame(fread("LINKS/plinks_string_dup.csv", header = TRUE, sep = ' '))
    plinks <- removeReverseDuplicates(plinksDup)
  }
  return(plinks)
}

createLinks <- function(dfLinks, map) {
  portion <- 100000
  i <- 1
  n <- nrow(dfLinks)
  dfNewLinks <- dfLinks
  dfNewLinks[, c("symbol1", "symbol2")] <- NA

  while(i <= n) {
    j <- i + portion - 1
    if(j > n) {
      j <- n
    }
    cat("range", i, "to", j, "...\n")

    dfNewLinks[i:j, "symbol1"] <- map[as.character(dfNewLinks[i:j, 1])]
    dfNewLinks[i:j, "symbol2"] <- map[as.character(dfNewLinks[i:j, 2])]

    i <- i+portion
  }

  dfNewLinks <- dfNewLinks[, c("symbol1", "symbol2", "combined_score")]
  nOld <- nrow(dfNewLinks)
  dfNewLinks <- dfNewLinks[!is.na(dfNewLinks[, "symbol1"]), ] # remove NA rows
  dfNewLinks <- dfNewLinks[!is.na(dfNewLinks[, "symbol2"]), ] # remove NA rows
  nNew <- nrow(dfNewLinks)
  
  cat("*********", nOld-nNew, "NULL rows are removed!\n")
  return(dfNewLinks)
}

getSlinks <- function(maps) {
  slinks <- list()
  for(m in names(maps)) {
    fname = paste("LINKS/slinks_", m, ".csv", sep="")
    if(file.exists(fname)) {
      cat(paste("READING FROM -> ", fname, sep=""), "\n")
      slinks[[m]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
    } else {
      cat(paste("CREATING -> ", fname, sep = ""), "\n")
      if(m == "string") {
        xlinks <- readPlinks()
      } else {
        xlinks <- as.data.frame(fread("LINKS/elinks_inet.csv", header = TRUE, sep = ','))
      }
      slinks[[m]] <- createLinks(xlinks, maps[[m]])
      cat("WRITING", fname, "\n")
      write.csv(slinks[[m]], fname, row.names = FALSE)
    }
  }
  return(slinks)
}

slinks <- getSlinks(maps)
