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

#######################

createGoSimHelper <- function(subject, measure, ontologies) {
  dfGo <- subject
  for(o in names(ontologies)) {
    dfGo[, o] <- NA
    for(r in 1:nrow(dfGo)) {
      dfGo[r, o] <- geneSim(dfGo[r,1], dfGo[r,2], 
                             semData = ontologies[[o]],
                             measure = measure, combine = "BMA")[[1]]
    }
  }
  return(dfGo)
}

createGoSim <- function(msDegLinks, measures, ontologies) {
  gosim <- list()
  cat("topology,subject,measure,seconds\n")
  for(t in names(msDegLinks)) {
    for(s in names(msDegLinks[[t]])) {
      for(m in measures) {
        start_time <- proc.time()
        gosim[[t]][[s]][[m]] <- createGoSimHelper(msDegLinks[[t]][[s]], m, ontologies)
        elapsed_time <- proc.time() - start_time
        cat(paste(t, ",", s, ",", m, ",", round(elapsed_time[[3]],3), sep=""), "\n")
      }
    }
  }
  return(gosim)
}

addCombinedSimilarityScores <- function(gosim) {
  # Takes mean for information content-based methods (resnik, lin, rel, jiang)
  gosimComb <- gosim
  for(t in names(gosimComb)) {
    for(s in names(gosimComb[[t]])) {
      subj <- gosim[[t]][[s]]
      gosimComb[[t]][[s]][["Comb"]] <- subj[[1]] # just take a copy
      gosimComb[[t]][[s]][["Comb"]][, c(3:5)] <- NA # delete copied scores

      for(o in ontTypes) {
        total <- 0
        for(m in names(subj)[!names(subj) == "Wang"]) {
          df <- gosim[[t]][[s]][[m]]
          total <- total + df[, o]
        }
        gosimComb[[t]][[s]][["Comb"]][, o] <- total / 4
      }
    }
  }
  return(gosimComb)
}

writeGosim <- function(gosimObj, ontTypes, removeNA, hasColNames) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in ontTypes) {
          df <- gosimObj[[t]][[s]][[m]][, c("symbol1", "symbol2", o)]
          if(removeNA) {
            df <- na.omit(df)
            row.names(df) <- NULL
          }
          fdir <- paste("RES/GOSIM/", m, "/", sep = "")
          fname <- paste(fdir, t, "_", s, "_", o, ".tsv", sep="")
          print(fname)
          if(!file.exists(fdir)) dir.create(fdir) # create directory
          write.table(df, fname, row.names = FALSE, col.names = hasColNames, sep="\t")
        }
      }
    }
  }
}

readGosim <- function(topologies, subjects, measures, ontTypes, includeComb) {
  gosim <- list()
  if(includeComb) measures <- c(measures, "Comb")

  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/GOSIM/", m , "/", t, "_", s, "_", o, ".tsv", sep="")
          if(o == "MF") {
            gosim[[t]][[s]][[m]] <- as.data.frame(fread(fname, header = TRUE, sep = '\t'))
          } else {
            gosim[[t]][[s]][[m]][,o] <- as.data.frame(fread(fname, header = TRUE, sep = '\t'))[,3]
          }
        }
      }
    }
  }
  return(gosim)
}

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
# gosim <- createGoSim(testLinks, measures, ontologies)
# gosim <- createGoSim(msDegLinks, measures, ontologies)
# gosimWithComb <- addCombinedSimilarityScores(gosim)
# writeGosim(gosimWithComb, ontTypes, removeNA = FALSE, hasColNames = TRUE)

gosimWithComb <- readGosim(topologies, subjects, measures, ontTypes, includeComb = TRUE)
