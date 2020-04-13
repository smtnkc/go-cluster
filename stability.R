source("vars.R")

NAMING <- "PROBEID"

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
getGosimNodes <- function(gosimObj) {
  gosimNodes <- list()
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          df <- gosimObj[[t]][[s]][[m]][[o]]
          gosimNodes[[t]][[s]][[m]][[o]] <- sort(union(df[,1], df[,2]))
        }
      }
    }
  }
  return(gosimNodes)
}
getLC <- function(df) {
  lc <- getLinkCommunities(df, plot = FALSE, verbose = FALSE)
  t_end <- Sys.time()
  dfLC <- lc[["nodeclusters"]]
  dfLC <- transform(dfLC, node = as.character(node))
  dfLC <- transform(dfLC, cluster = as.integer(cluster))
  return(dfLC)
}
addNonClusteredNodes <- function(gosimNodes, gosimLinkcomm) {
  gosimClustersAll <- list()
  for(t in names(gosimLinkcomm)) {
    for(s in names(gosimLinkcomm[[t]])) {
      for(m in names(gosimLinkcomm[[t]][[s]])) {
        for(o in names(gosimLinkcomm[[t]][[s]][[m]])) {
          clusteredNodesDf <- gosimLinkcomm[[t]][[s]][[m]][[o]]
          nonClusteredNodesDf <- data.frame(node = setdiff(gosimNodes[[t]][[s]][[m]][[o]],
                                                           gosimLinkcomm[[t]][[s]][[m]][[o]]$node),
                                            cluster = 9999, stringsAsFactors = FALSE)
          allNodesDf <- rbind(clusteredNodesDf, nonClusteredNodesDf)
          gosimClustersAll[[t]][[s]][[m]][[o]] <- allNodesDf
        }
      }
    }
  }
  return(gosimClustersAll)
}
getSymbols <- function() {
  dfSymPr <- as.data.frame(hgu133aSYMBOL[mappedkeys(hgu133aSYMBOL)])
  gene_symbols_parsed <- c()
  for(r in dfSymPr$symbol) {
    for(s in strsplit(r, "-")) {
      gene_symbols_parsed <- c(gene_symbols_parsed, s)
    }
  }
  gene_symbols_parsed <- unique(gene_symbols_parsed)
  return(gene_symbols_parsed)
}
getPerturbedGraphs <- function(df, perturbationRatio) {
  pGraphs <- list()
  allNodes <- union(df[,1], df[,2])
  numOfNodesToBeRemoved <- floor(length(allNodes)*perturbationRatio)
  nodesNotRemovedYet <- allNodes
  for(i in 1:(1/perturbationRatio)) {
    nodesToBeRemoved <- sample(nodesNotRemovedYet, numOfNodesToBeRemoved, replace = FALSE)
    nodesNotRemovedYet <- setdiff(nodesNotRemovedYet, nodesToBeRemoved) # to make a different random selection each time
    coreNodes <- setdiff(allNodes, nodesToBeRemoved)
    pGraphs[[i]] <- df[df[,1] %in% coreNodes & df[,2] %in% coreNodes, ]
  }
  return(pGraphs)
}
getPerturbedClusters <- function(perturbedGraphs, alg) {
  cl <- list()
  for(i in 1:length(perturbedGraphs)) {
    cat("Clustering the perturbed graph #", i, "\n")
    if(alg == "linkcomm")
      cl[[i]] <- getLC(perturbedGraphs[[i]])
    else if(alg == "mcl")
      cl[[i]] <- getMCL(perturbedGraphs[[i]])
  }
  return(cl)
}
getAvgBestMatchRatio <- function(originalClusters, listOfPerturbedClusters) {
  avgBMRlist <- c()
  for(i in 1:length(listOfPerturbedClusters)) {
    pCl <- listOfPerturbedClusters[[i]]
    totalBMR <- 0
    for(clusterId in unique(pCl$cluster)) {
      cl <- pCl[pCl$cluster == clusterId,] # a single cluster
      BMR <- getBestMatchRatio(originalClusters, cl)
      totalBMR <- totalBMR + BMR
    }
    avgBMR <- totalBMR / length(unique(pCl$cluster))
    cat("Avg BMR for perturbed graph #", i, "is", avgBMR, "\n")
    avgBMRlist <- c(avgBMRlist, avgBMR)
  }
  BMRmean <- round(mean(avgBMRlist),4)
  BMRsd <- sd(avgBMRlist)
  return(c(BMRmean, BMRsd))
}
getBestMatchRatio <- function(multiClusterDf, singleClusterDf) {
  bestMatchRatio <- 0
  for(candidateClusterId in unique(multiClusterDf$cluster)) {
    nodes <- multiClusterDf[multiClusterDf$cluster == candidateClusterId, 1]
    overlap <- Reduce(intersect, list(v1=nodes, v2=singleClusterDf$node))
    matchRatio <- length(overlap) / max(length(nodes), length(singleClusterDf$node))
    if(matchRatio > bestMatchRatio)
      bestMatchRatio <- matchRatio
  }
  return(bestMatchRatio)
}

# FOR MCL
getAdjMatrix <- function(df) {
  nodes <- union(df[, 1], df[, 2])
  N <- length(nodes)
  adj <- matrix(0, N*N, nrow=N, ncol=N, dimnames = list(nodes, nodes))
  for(r in 1:nrow(df)) {
    adj[df[r,1], df[r,2]] <- df[r, 3]
    adj[df[r,2], df[r,1]] <- df[r, 3]
  }
  return(adj)
}
arrangeClusterIDs <- function(dfMCL) {
  # This function arranges all the cluster IDs to be consecutive numbers
  id <- 1
  currClusterId <- dfMCL[1, 2]
  prevClusterId <- currClusterId
  
  for(i in 1:nrow(dfMCL)) {
    currClusterId <- dfMCL[i, 2]
    if(currClusterId != prevClusterId)
      id <- id + 1
    dfMCL[i, 2] <- id
    prevClusterId <- currClusterId
  }
  return(dfMCL)
}
getMCL <- function(df) {
  adj <- getAdjMatrix(df)
  adj <- adj[which(rowSums(adj) > 0), which(colSums(adj) > 0)] # to prevent mcl errors, remove unconnected nodes
  g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE)
  mclOutput <- mcl(x = g, addLoops=FALSE)
  dfMCL <- data.frame(node = V(g)$name, cluster = mclOutput$Cluster)
  dfMCL <- dfMCL[order(dfMCL$cluster), ]
  dfMCL <- arrangeClusterIDs(dfMCL)
  return(dfMCL)
}

##########################################################

gosim <- readGosim(topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)
gosim <- removeNAs(gosim)
mets <- gosim[["inet"]][["mets"]][["Wang"]][["MF"]]
t2d <- gosim[["inet"]][["t2d"]][["Wang"]][["MF"]]
cad <- gosim[["inet"]][["cad"]][["Wang"]][["MF"]]

pGraph.mets <- getPerturbedGraphs(mets, 0.05)
pGraph.t2d <- getPerturbedGraphs(t2d, 0.05)
pGraph.cad <- getPerturbedGraphs(cad, 0.05)

### LINKCOMM

LC.mets <- getLinkcomm(mets)
LC.t2d <- getLinkcomm(t2d)
LC.cad <- getLinkcomm(cad)

LC.pert.mets <- getPerturbedClusters(pGraph.mets, "linkcomm")
LC.avgBMR.mets <- getAvgBestMatchRatio(LC.mets, LC.pert.mets)

LC.pert.t2d <- getPerturbedClusters(pGraph.t2d, "linkcomm")
LC.avgBMR.t2d <- getAvgBestMatchRatio(LC.t2d, LC.pert.t2d)

LC.pert.cad <- getPerturbedClusters(pGraph.cad, "linkcomm")
LC.avgBMR.cad <- getAvgBestMatchRatio(LC.cad, LC.pert.cad)

#### MCL

MCL.mets <- getMCL(mets)
MCL.t2d <- getMCL(t2d)
MCL.cad <- getMCL(cad)

MCL.pert.mets <- getPerturbedClusters(pGraph.mets, "mcl")
MCL.avgBMR.mets <- getAvgBestMatchRatio(MCL.mets, MCL.pert.mets)

MCL.pert.t2d <- getPerturbedClusters(pGraph.t2d, "mcl")
MCL.avgBMR.t2d <- getAvgBestMatchRatio(MCL.t2d, MCL.pert.t2d)

MCL.pert.cad <- getPerturbedClusters(pGraph.cad, "mcl")
MCL.avgBMR.cad <- getAvgBestMatchRatio(MCL.cad, MCL.pert.cad)


#### SPICi

writePerturbedGraphs <- function(perturbedGraphs) {
  pertDir <- "RES/SPICi/PERTURBED/"
  unlink(pertDir, recursive = TRUE)
  dir.create(pertDir)
  for(disease in names(perturbedGraphs)) {
    for(pg in 1:length(perturbedGraphs[[1]])) {
      fname <- paste(pertDir, disease, "_", pg, ".tsv", sep = "")
      df <- perturbedGraphs[[disease]][[pg]]
      write.table(df, fname, row.names = FALSE, col.names = FALSE, sep="\t")
    }
  }
}

writePerturbedGraphs(list("mets" = pGraph.mets, "t2d" = pGraph.t2d, "cad" = pGraph.cad))

# HERE RUN spiciPerturbedRunner.sh

readSpiciClustersHelper <- function(clusterFile) {
  fileContent <- readLines(clusterFile)
  clusterDF <- data.frame(node = character(), cluster = integer())
  if(length(fileContent) == 0) {
    return(clusterDF)
  }
  c <- 1
  for(i in 1:length(fileContent)) {
    for(n in strsplit(gsub('\"', '', fileContent[i]), split="\t")) {
      clusterDF <- rbind(clusterDF, data.frame(node = n, cluster = c,
                                               stringsAsFactors = FALSE))
    }
    c <- c + 1
  }
  return(clusterDF)
}

readSpiciPertClusters <- function(perturbedGraphs) {
  pertDir <- "RES/SPICi/LINE_BY_LINE_CLUSTERS/PERTURBED/"
  clusters <- list()
  for(disease in names(perturbedGraphs)) {
    for(pg in 1:length(perturbedGraphs[[1]])) {
      fname <- paste(pertDir, disease, "_", pg, ".tsv", sep = "")
          cat(fname, "...\n")
          clusters[[disease]][[pg]] <- readSpiciClustersHelper(fname)
    }
  }
  return(clusters)
}

spiciPertClusters <- readSpiciPertClusters(
  list("mets" = pGraph.mets, "t2d" = pGraph.t2d, "cad" = pGraph.cad))

readSpiciClusters <- function(perturbedGraphs) {
  clustDir <- "RES/SPICi/LINE_BY_LINE_CLUSTERS/byPROBEID/inet_"
  clusters <- list()
  for(disease in names(perturbedGraphs)) {
    fname <- paste(clustDir, disease, "_Wang_MF.tsv", sep = "")
    cat(fname, "...\n")
    clusters[[disease]] <- readSpiciClustersHelper(fname)
  }
  return(clusters)
}

spiciClusters <- readSpiciClusters(list("mets" = pGraph.mets, "t2d" = pGraph.t2d, "cad" = pGraph.cad))

SPICI.mets <- spiciClusters[["mets"]]
SPICI.t2d <- spiciClusters[["t2d"]]
SPICI.cad <- spiciClusters[["cad"]]

SPICI.avgBMR.mets <- getAvgBestMatchRatio(SPICI.mets, spiciPertClusters[["mets"]])
SPICI.avgBMR.t2d <- getAvgBestMatchRatio(SPICI.t2d, spiciPertClusters[["t2d"]])
SPICI.avgBMR.cad <- getAvgBestMatchRatio(SPICI.cad, spiciPertClusters[["cad"]])

