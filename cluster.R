source("vars.R")

readGosim <- function() {
  gosim <- list()
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        fname <- paste("RES/GOSIM/", t, "_", s, "_", m, ".csv", sep="")
        gosim[[t]][[s]][[m]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
      }
    }
  }
  return(gosim)
}

gosim <- readGosim()

# http://compbio.cs.princeton.edu/spici/
# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7300425
