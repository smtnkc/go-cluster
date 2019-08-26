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

##########################################################

getAdjMatrices <- function(gosimObj) {
  cat("topology,subject,measure,ontology,dimension,time(s),memory(mb)\n")
  adjMatrices <- list()
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          t_start <- Sys.time()
          p <- profmem({
          df <- gosimObj[[t]][[s]][[m]][[o]]
          nodes <- union(df[, 1], df[, 2])
          N <- length(nodes)
          adj <- matrix(0, N*N, nrow=N, ncol=N, dimnames = list(nodes, nodes))
          for(r in 1:nrow(df)) {
            adj[df[r,1], df[r,2]] <- df[r, 3]
            adj[df[r,2], df[r,1]] <- df[r, 3]
          }
          adjMatrices[[t]][[s]][[m]][[o]] <- adj
          t_end <- Sys.time()
          })

          cat(paste(t,s,m,o,N,
                    round(difftime(t_end, t_start, units = "secs"),5),
                    total(p)/1000000, sep=","), "\n")
        }
      }
    }
  }
  return(adjMatrices)
}

#adjMatrices <- getAdjMatrices(gosim)

##########################################################

getMCL <- function(adjMatrices) {
  gosimClusters <- list()
  cat("topology,subject,measure,ontology,edge,node,cluster,time(s),memory(mb)\n")
  for(t in names(adjMatrices)) {
    for(s in names(adjMatrices[[t]])) {
      for(m in names(adjMatrices[[t]][[s]])) {
        for(o in names(adjMatrices[[t]][[s]][[m]])) {
          t_start <- Sys.time()
          p <- profmem({
          adj <- adjMatrices[[t]][[s]][[m]][[o]]
          adj <- adj[which(rowSums(adj) > 0), which(colSums(adj) > 0)] # to prevent mcl errors, remove unconnected nodes
          g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE)
          mclOutput <- mcl(x = g, addLoops=FALSE)
          })
          t_end <- Sys.time()
          cat(paste(t,s,m,o,
                    length(E(g)),
                    length(V(g)),
                    mclOutput$K,
                    round(difftime(t_end, t_start, units = "secs"),5),
                    total(p)/1000000, sep=","), "\n")
          df <- data.frame(node = V(g)$name, cluster = mclOutput$Cluster)
          gosimClusters[[t]][[s]][[m]][[o]] <- df[order(df$cluster), ]
        }
      }
    }
  }
  return(gosimClusters)
}

#gosimMCL <- getMCL(adjMatrices)

updateClusterIDs <- function(gosimMCL) {
  updatedGosimMCL <- list()
  for(t in names(gosimMCL)) {
    for(s in names(gosimMCL[[t]])) {
      for(m in names(gosimMCL[[t]][[s]])) {
        for(o in names(gosimMCL[[t]][[s]][[m]])) {
          df <- gosimMCL[[t]][[s]][[m]][[o]]
          id <- 1
          curr <- df[1, 2]
          prev <- curr

          for(i in 1:nrow(df)) {
            curr <- df[i, 2]
            if(curr != prev) id <- id + 1
            df[i, 2] <- id
            prev <- curr
          }
          updatedGosimMCL[[t]][[s]][[m]][[o]] <- df
        }
      }
    }
  }
  return(updatedGosimMCL)
}

#gosimMCL <- updateClusterIDs(gosimMCL)

##########################################################

writeMCL <- function(gosimMCL, naming) {
  for(t in names(gosimMCL)) {
    for(s in names(gosimMCL[[t]])) {
      for(m in names(gosimMCL[[t]][[s]])) {
        for(o in names(gosimMCL[[t]][[s]][[m]])) {
          df <- gosimMCL[[t]][[s]][[m]][[o]]
          fname <- paste("RES/MCL/CLUSTERS/by", naming, "/", t,
                         "_", s, "_", m, "_", o, ".csv", sep = "")
          cat(fname, "...\n")
          write.table(df, fname, row.names = FALSE, col.names = TRUE, sep=",")
        }
      }
    }
  }
}

#writeMCL(gosimMCL, naming = NAMING)

readMCL <- function(topologies, subjects, measures, ontTypes, includeComb, naming) {
  gosimMCL <- list()
  if(includeComb) measures <- c(measures, "Comb")

  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/MCL/CLUSTERS/by", naming, "/",
                         t, "_", s, "_", m, "_", o, ".csv", sep="")
          gosimMCL[[t]][[s]][[m]][[o]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
        }
      }
    }
  }
  return(gosimMCL)
}

gosimMCL <- readMCL(topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)

##########################################################

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

drawClusters <- function(gosimObj, gosimMCL, naming, BHIScores) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in c("Wang")) {
        #for(m in names(gosimObj[[t]][[s]])) {
        #for(o in c("BP")) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          nodes <- gosimMCL[[t]][[s]][[m]][[o]]
          edges <- gosimObj[[t]][[s]][[m]][[o]]
          edges <- edges[edges[,1] %in% nodes$node, ]
          edges <- edges[edges[,2] %in% nodes$node, ]

          if(nrow(nodes) != 0 && nrow(edges) != 0) {
            net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
            V(net)$color <- V(net)$cluster + 1

            fname <- paste("PLOTS/CLUSTERS/MCL/by", naming, "/",
                           t, "_", s, "_", m, "_", o, ".png", sep="")
            cat(fname, "...\n")
            png(filename=fname, width = 2280, height = 2280)
            lay <- layout_in_circle(net)

            info1 <- paste("Topology: ", toupper(t),
                           "     Subject: ", toupper(s),
                           "     Measure: ", toupper(m),
                           "     Ontology: ", toupper(o), sep="")

            info2 <- paste(
              "Total nodes: " , nrow(gosimMCL[[t]][[s]][[m]][[o]]),
              "     Clustered nodes: ", length(V(net)[V(net)$cluster != -1]),
              "     Clusters: ", length(unique(V(net)[V(net)$cluster != -1]$cluster)), sep="")

            if(naming == "PROBEID") {
              info2 <- paste(info2, "     BHI: ", BHIScores[[t]][[s]][[m]][[o]], sep="")
            }

            plot(net, layout=lay,
                 vertex.label=nodes$node, vertex.shape="circle",
                 vertex.size=4, #vertex.size=nchar(as.character(nodes$node))*4+5, vertex.size2=10,
                 vertex.label.color="black", vertex.label.cex=0.6,
                 vertex.frame.color="gray50", vertex.label.family="Helvetica",
                 edge.width=1, edge.label.color="black",
                 edge.color="gray50", edge.curved=0, edge.label.family="Helvetica")

            mtext(info1, side=3, line=0, cex=3)
            mtext(info2, side=1, line=0, cex=3)

            dev.off()
          }
        }
      }
    }
  }
}

drawClusters(gosim, gosimMCL, NAMING, BHIScoresMCL)
