source("vars.R")
NAMING <- "PROBEID"
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

##############################################

printGosimCounts <- function(gosim) {
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
      print(paste(t,"_",s,"      ",round(n,0), "/", round(e,0), sep=""))
    }
  }
}
printGosimCounts(gosim)

##############################################

readClusters <- function(type, topologies, subjects, measures, ontTypes, includeComb, naming) {
  gosimX <- list()
  if(includeComb) measures <- c(measures, "Comb")
  
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/", type, "/CLUSTERS/by", naming, "/", t ,
                         "_", s, "_", m, "_", o, ".csv", sep="")
          gosimX[[t]][[s]][[m]][[o]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
        }
      }
    }
  }
  return(gosimX)
}
gosimMCL <- readClusters("MCL", topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)
gosimSpici <- readClusters("SPICi", topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)
gosimLinkcomm <- readClusters("LINKCOMM", topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)

printClusterStats <- function(clusters) {
  for(c in names(clusters)) {
    for(t in names(clusters[[c]])) {
      for(s in names(clusters[[c]][[t]])) {
        nn<-0
        nc<-0
        for(m in names(clusters[[c]][[t]][[s]])) {
          for(o in names(clusters[[c]][[t]][[s]][[m]])) {
            df <- clusters[[c]][[t]][[s]][[m]][[o]]
            nn <- nn + length(unique(df[df$cluster != 9999,]$node))
            nc <- nc + length(unique(df[df$cluster != 9999,]$cluster))
          }
        }
        print(paste(c,"_",t,"_",s,"      ",round(nn/15,2),"/",round(nc/15,2),"=",round(nn/nc,2), sep=""))
      }
    }
  }
}
printClusterStats(list("mcl" = gosimMCL, "spici" = gosimSpici, "linkcomm" = gosimLinkcomm))

#########################################################################################

bhi.MCL <- as.data.frame(fread("RES/MCL/SCORES.csv", header = TRUE, sep = ','))
bhi.SPICi <- as.data.frame(fread("RES/SPICi/SCORES.csv", header = TRUE, sep = ','))
bhi.Linkcomm <- as.data.frame(fread("RES/LINKCOMM/SCORES.csv", header = TRUE, sep = ','))

bhi.MCL$clustering <- c(rep("MCL", nrow(bhi.MCL)))
bhi.SPICi$clustering <- c(rep("SPICi", nrow(bhi.SPICi)))
bhi.Linkcomm$clustering <- c(rep("Linkcomm", nrow(bhi.Linkcomm)))
bhi.ALL <- rbind(bhi.MCL, bhi.SPICi, bhi.Linkcomm)
library(dplyr)
printBHIStats <- function(bhi.ALL) {
  bhi.ALL %>% group_by(topology, subject, clustering) %>% summarise(BHI_mean=(mean(BHI)))
}
printBHIStats(bhi.ALL)


