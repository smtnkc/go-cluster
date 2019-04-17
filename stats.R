source("vars.R")
NAMING <- "SYMBOL"
gosim <- readGosim(topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)
removeNAs <- function(gosimObj) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          gosimObj[[t]][[s]][[m]][[o]] <- na.omit(gosimObj[[t]][[s]][[m]][[o]])
        }
      }
    }
  }
  return(gosimObj)
}
gosim <- removeNAs(gosim)

getGosimCounts <- function(gosim) {
  gosimCounts <- list()

  for(t in topologies) {
    for(s in subjects) {
      n <- 0  
      e <- 0
      i <- 0
      for(m in measures) {
        for(o in ontTypes) {
          n = n + length(union(gosim[[t]][[s]][[m]][[o]][,1], gosim[[t]][[s]][[m]][[o]][,2]))
          e = e + nrow(gosim[[t]][[s]][[m]][[o]])
          i <- i+1
        }
      }
      n <- n/i
      e <- e/i
      gosimCounts[[t]][[s]] <- paste(round(n,0), "/", round(e,0), sep="")
    }
  }
  return(gosimCounts)
}
gosimCounts <- getGosimCounts(gosim)