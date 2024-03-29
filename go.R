source("vars.R")

readMsDegLinks <- function(topologies, subjects, p, top) {
  msDegLinks <- list()
  for(t in topologies) {
    for(s in subjects) {
      fname <- paste("RES/", t, "_", s, "_P", Stringify(p),
                     "_T", Stringify(top), ".csv", sep ="")
      msDegLinks[[t]][[s]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
    }
  }
  return(msDegLinks)
}

msDegLinks <- readMsDegLinks(topologies, subjects, P_VAL, TOP_N)

getOntologies <- function(ontTypes) {
  ontologies <- list()
  for(o in ontTypes) {
    ontologies[[o]] <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont = o)
  }
  return(ontologies)
}

ontologies <- getOntologies(ontTypes)

##############################################

createGoSim <- function(msDegLinks, measures, ontTypes) {
  gosim <- list()
  cat("topology,subject,measure,seconds\n")
  for(t in names(msDegLinks)) {
    for(s in names(msDegLinks[[t]])) {
      for(m in measures) {
        for(o in ontTypes) {
          start_time <- proc.time()
          dfGo <- msDegLinks[[t]][[s]]
          dfGo[, o] <- NA
          for(r in 1:nrow(dfGo)) {
            dfGo[r, o] <- geneSim(dfGo[r,1], dfGo[r,2],
                                  semData = ontologies[[o]],
                                  measure = m, combine = "BMA")[[1]]
          }
          gosim[[t]][[s]][[m]][[o]] <- dfGo
          elapsed_time <- proc.time() - start_time
          cat(paste(t, s, m, o, round(elapsed_time[[3]],3), sep=","), "\n")
        }
      }
    }
  }
  return(gosim)
}

addCombinedSimilarityScores <- function(gosimObj) {
  # Takes mean for information content-based methods (resnik, lin, rel, jiang)
  gosimComb <- gosimObj
  for(t in names(gosimComb)) {
    for(s in names(gosimComb[[t]])) {
      subj <- gosimObj[[t]][[s]]
      gosimComb[[t]][[s]][["Comb"]] <- subj[[1]] # just take a copy

      for(o in names(gosimComb[[t]][[s]][["Comb"]])) {
        gosimComb[[t]][[s]][["Comb"]][[o]][, o] <- NA # delete copied weights
        total <- 0
        for(m in names(subj)) {
          df <- gosimObj[[t]][[s]][[m]]
          total <- total + df[[o]][, 3]
        }
        gosimComb[[t]][[s]][["Comb"]][[o]][, 3] <- total / 4
      }
    }
  }
  return(gosimComb)
}

writeGosim <- function(gosimObj, spici, naming) {
  if(naming == "SYMBOL") {
    fdir <- "RES/GOSIM/bySYMBOL/"
  } else {
    fdir <- "RES/GOSIM/byPROBEID/"
  }
  if(spici) {
    fdir <- paste(fdir, "forSPICi/", sep="")
  }
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          df <- gosimObj[[t]][[s]][[m]][[o]]
          if(spici) {
            df <- na.omit(df) # Since NAs are invalid for SPICi
            row.names(df) <- NULL # Reset row ids
            df[df == 0] <- 0.00000001 # Since zero-weights are invalid for SPICi
            hasColNames <- FALSE # Since headers are invalid for SPICi
          } else {
            hasColNames <- TRUE
          }
          fname <- paste(fdir, t, "_", s, "_", m, "_", o, ".tsv", sep="")
          print(fname)
          if(!file.exists(fdir)) dir.create(fdir) # create directory
          write.table(df, fname, row.names = FALSE, col.names = hasColNames, sep="\t")
        }
      }
    }
  }
}

##############################################

generateTestLinks <- function(msDegLinks, size) {
  testLinks <- list()
  for(t in names(msDegLinks)) {
    for(s in names(msDegLinks[[t]])) {
      nAll <- nrow(msDegLinks[[t]][[s]])
      if(size > nAll) {
        cat("ERROR: Choose a sample size less than", nAll, "!")
        return(NULL)
      } else {
        sampleRows <- sample(1:nAll, size, replace = FALSE)
      }
      testLinks[[t]][[s]] <- msDegLinks[[t]][[s]][sampleRows,]
    }
  }
  return(testLinks)
}

# testLinks <- generateTestLinks(msDegLinks, size = 10)
# gosimTest <- createGoSim(testLinks, measures, ontTypes)
# gosimTestWithComb <- addCombinedSimilarityScores(gosimTest)

##############################################

# gosim <- createGoSim(msDegLinks, measures, ontTypes)
# gosimWithComb <- addCombinedSimilarityScores(gosim)

gosimWithComb <- readGosim(topologies, subjects, measures, ontTypes,
                           includeComb = TRUE, naming = "SYMBOL")

##########################################################

convertSymbols2Probes <- function(gosimObj) {
  gosimPr <- list()
  hguSym <- hgu133aSYMBOL
  mappedPr <- mappedkeys(hguSym)
  dfSymPr <- as.data.frame(hguSym[mappedPr])

  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          cat(paste("mapping ", t,"_",s,"_",m,"_",o," ...\n", sep=""))
          df <- gosimObj[[t]][[s]][[m]][[o]]
          ont <- colnames(df)[3]
          nodes <- sort(union(df[, 1], df[, 2]))
          map <- dfSymPr[dfSymPr$symbol %in% nodes, ]
          dfPr <- data.frame(probe1 = character(),
                             probe2 = character(),
                             ont = character(),
                             stringsAsFactors = FALSE)
          for(i in 1:nrow(df)) {
            fromPr <- map[map$symbol == df[i, "symbol1"], "probe_id"]
            toPr <- map[map$symbol == df[i, "symbol2"], "probe_id"]
            for(pr1 in fromPr) {
              for(pr2 in toPr) {
                dfPr <- rbind(dfPr, data.frame(probe1=pr1, probe2=pr2, ont=df[i,3],
                                               stringsAsFactors = FALSE))
              }
            }
          }
          gosimPr[[t]][[s]][[m]][[o]] <- dfPr
        }
      }
    }
  }
  return(gosimPr)
}

gosimPrWithComb <- convertSymbols2Probes(gosimWithComb)

writeGosim(gosimPrWithComb, spici = TRUE, naming = "PROBEID")
