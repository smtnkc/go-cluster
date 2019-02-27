source("vars.R")
library("clValid")

NAMING = "PROBEID"

readClusters <- function(type, topologies, subjects, measures, ontTypes,
                         includeComb, naming) {
  gosimMCL <- list()
  if(includeComb) measures <- c(measures, "Comb")
  
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/", type, "/CLUSTERS/by", naming, "/", t ,
                         "_", s, "_", m, "_", o, ".csv", sep="")
          gosimMCL[[t]][[s]][[m]][[o]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
        }
      }
    }
  }
  return(gosimMCL)
}

gosimMCL <- readClusters("MCL", topologies, subjects, measures, ontTypes,
                         includeComb = FALSE, naming = NAMING)

gosimSpici <- readClusters("SPICi", topologies, subjects, measures, ontTypes,
                           includeComb = FALSE, naming = NAMING)

################################################################

getBHIScores <- function(gosimClusters) {
  BHIScores <- list()
  for(t in names(gosimClusters)) {
    for(s in names(gosimClusters[[t]])) {
      for(m in names(gosimClusters[[t]][[s]])) {
        BHIScores[[t]][[s]][[m]] <- list()
        for(o in names(gosimClusters[[t]][[s]][[m]])) {
          df <- gosimClusters[[t]][[s]][[m]][[o]]
          clusters <- df$cluster
          names(clusters) <- df$node
          score <- BHI(clusters, annotation="hgu133a.db", names=names(clusters), category="all")
          print(score)
          BHIScores[[t]][[s]][[m]][[o]] <- score
        }
      }
    }
  }
  return(BHIScores)
}

# BHIScoresMCL <- getBHIScores(gosimMCL)
# BHIScoresSpici <- getBHIScores(gosimSpici)

writeBHIScores <- function(BHIScores, type) {
  fname  <- paste("RES/", type, "/SCORES.csv", sep = "")
  df <- data.frame(topology=character(),
                   subject=character(),
                   measure=character(),
                   ontology=character(),
                   score=double(),
                   stringsAsFactors = FALSE)
  
  for(t in names(BHIScores)) {
    for(s in names(BHIScores[[t]])) {
      for(m in names(BHIScores[[t]][[s]])) {
        for(o in names(BHIScores[[t]][[s]][[m]])) {
          df <- rbind(df, data.frame(topology=t,
                                     subject=s,
                                     measure=m,
                                     ontology=o,
                                     BHI=round(BHIScores[[t]][[s]][[m]][[o]],3),
                                     stringsAsFactors = FALSE))

        }
      }
    }
  }
  write.table(df, fname, row.names = FALSE, col.names = TRUE, sep=",")
}

# writeBHIScores(BHIScoresMCL, "MCL")
# writeBHIScores(BHIScoresSpici, "SPICi")

#################################################################

readBHIScores <- function(type, topologies, subjects, measures, ontTypes,
                          includeComb, naming) {
  BHIScores <- list()
  fname  <- paste("RES/", type, "/SCORES.csv", sep = "")
  df <- as.data.frame(fread(fname, header = TRUE, sep = ','))
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        BHIScores[[t]][[s]][[m]] <- list()
        for(o in ontTypes) {
          BHIScores[[t]][[s]][[m]][[o]] <- 
            df[df$topology==t & df$subject==s & df$measure==m & df$ontology==o, "BHI"]
        }
      }
    }
  }
  return(BHIScores)
}

BHIScoresMCL <- readBHIScores("MCL", topologies, subjects, measures, ontTypes,
                              includeComb, naming = NAMING)

BHIScoresSpici <- readBHIScores("SPICi", topologies, subjects, measures, ontTypes,
                              includeComb, naming = NAMING)
