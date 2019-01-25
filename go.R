source("vars.R")

##################################################

readMsDegLinks <- function(idTypes, topologies, subjects, P_VAL, TOP_N) {
  ms_deglinks <- list()
  for(i in idTypes) {
    for(t in topologies) {
      for(s in subjects) {
        fname <- paste("RES/", i, "_", t, "_", s, "_P", Stringify(P_VAL),
                       "_T", Stringify(TOP_N), ".csv", sep ="")
        ms_deglinks[[i]][[t]][[s]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
      }
    }
  }
  return(ms_deglinks)
}

ms_deglinks <- readMsDegLinks(idTypes, topologies, subjects, P_VAL, TOP_N)

getOntologies <- function(ontTypes) {
  ontologies <- list()
  ont <- list()
  for(o in ontTypes) {
    ont[[o]] <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont = o)
  }
  ontologies[["s"]] <- ont
  ont <- list()
  for(o in ontTypes) {
    ont[[o]] <- godata('org.Hs.eg.db', keytype = "ENTREZID", ont = o)
  }
  ontologies[["e"]] <- ont
  return(ontologies)
}

ontologies <- getOntologies(ontTypes)

#######################

getGoDf <- function(subject, measure, ontology) {
  df.go <- subject
  for(o in names(ontology)) {
    df.go[, o] <- NA
    for(r in 1:nrow(df.go)) {
      df.go[r, o] <- geneSim(df.go[r,1], df.go[r,2], 
                             semData = ontology[[o]],
                             measure = measure, combine = "BMA")[[1]]
    }
  }
  return(df.go)
}

getGoSim <- function(ms_deglinks, measures, ontologies) {
  gosim <- list()
  for(i in names(ms_deglinks)) {
    for(t in names(ms_deglinks[[i]])) {
      for(s in names(ms_deglinks[[i]][[t]])) {
        for(m in measures) {
          cat(paste(i, "_", t, "_", s, "_", m, sep=""))
          start_time <- proc.time()
          gosim[[i]][[t]][[s]][[m]] <- getGoDf(ms_deglinks[[i]][[t]][[s]], m, ontologies[[i]])
          elapsed_time <- proc.time() - start_time
          cat("\t-----", elapsed_time[[3]], "seconds\n")
        }
      }
    }
  }
  return(gosim)
}

getTestLinks <- function(ms_deglinks) {
  testLinks <- list()
  for(i in names(ms_deglinks)) {
    for(t in names(ms_deglinks[[i]])) {
      for(s in names(ms_deglinks[[i]][[t]])) {
        testLinks[[i]][[t]][[s]] <- ms_deglinks[[i]][[t]][[s]][1:10,]
      }
    }
  }
  return(testLinks)
}

testLinks <- getSampleLinks(ms_deglinks)

#gosim <- getGoSim(testLinks, measures, ontologies)
gosim <- getGoSim(ms_deglinks, measures, ontologies)






