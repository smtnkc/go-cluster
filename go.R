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
          cat(paste(t, ",", s, ",", m, ",", o, ",", round(elapsed_time[[3]],3), sep=""), "\n") 
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

writeGosim <- function(gosimObj, spici) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          df <- gosimObj[[t]][[s]][[m]][[o]]
          if(spici) {
            df <- na.omit(df) # Since NAs are invalid for SPICi
            row.names(df) <- NULL
            df[df == 0] <- 0.00000001 # Since zero-weights are invalid for SPICi
            hasColNames <- FALSE # Since headers are invalid for SPICi
            fdir <- "RES/SPICi/"
          } else {
            hasColNames <- TRUE
            fdir <- "RES/GOSIM/"
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
# writeGosim(gosimWithComb, spici = TRUE) # When exporting for SPIci, remove NAs, zero-weighths, and colnames

gosimWithComb <- readGosim(topologies, subjects, measures, ontTypes, includeComb = TRUE)
