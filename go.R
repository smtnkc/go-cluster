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

getDfGo <- function(subject, measure, ontologies) {
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

getGoSim <- function(msDegLinks, measures, ontologies) {
  gosim <- list()
  cat("topology,subject,measure,seconds\n")
  for(t in names(msDegLinks)) {
    for(s in names(msDegLinks[[t]])) {
      for(m in measures) {
        start_time <- proc.time()
        gosim[[t]][[s]][[m]] <- getDfGo(msDegLinks[[t]][[s]], m, ontologies)
        elapsed_time <- proc.time() - start_time
        cat(paste(t, ",", s, ",", m, ",", round(elapsed_time[[3]],3), sep=""), "\n")
      }
    }
  }
  return(gosim)
}

getTestLinks <- function(msDegLinks, size) {
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

#testLinks <- getTestLinks(msDegLinks, size = 10)
#gosim <- getGoSim(testLinks, measures, ontologies)

gosim <- getGoSim(msDegLinks, measures, ontologies)

writeGosim <- function(gosim, removeNA) {
  for(t in names(gosim)) {
    for(s in names(gosim[[t]])) {
      for(m in names(gosim[[t]][[s]])) {
        df <- gosim[[t]][[s]][[m]]
        fname <- paste("RES/GOSIM/ALL/", t, "_", s, "_", m, ".csv", sep="")
        print(fname)
        if(removeNA) {
          df <- na.omit(df)
          row.names(df) <- NULL
        }
        write.csv(df, fname, row.names = FALSE)
      }
    }
  }
}

writeGosim(gosim, removeNA = FALSE)

writeGosimByAnnotation <- function(gosim, removeNA, hasColNames) {
  for(t in names(gosim)) {
    for(s in names(gosim[[t]])) {
      for(m in names(gosim[[t]][[s]])) {
        for(a in c(3:5)) {
          df <- gosim[[t]][[s]][[m]]
          fname <- paste("RES/GOSIM/", colnames(df)[a], "/", t, "_", s, "_", m, ".csv", sep="")
          print(fname)
          df <- df[, c(1,2,a)]
          if(removeNA) {
            df <- na.omit(df)
          }
          write.table(df, fname, row.names = FALSE, col.names = hasColNames, sep="\t")
        }
      }
    }
  }
}

writeGosimByAnnotation(gosim, removeNA = TRUE, hasColNames = FALSE)
