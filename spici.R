source("vars.R")
library(igraph)
library(qgraph)
library(RColorBrewer)

### NOTE ###
### SPICi clustering is executed by spiciRunner.sh under RES/SPICi directory.
### In this file, only its results are read, formatted, and visualized.

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

######################################################

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

gosimNodes <- getGosimNodes(gosim)

######################################################

readSpiciHelper <- function(clusterFile) {
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

readSpici <- function(gosimObj, naming) {
  clusters <- list()
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          fname <- paste("RES/SPICi/LINE_BY_LINE_CLUSTERS/by", naming, "/",
                        t, "_", s, "_", m, "_", o, ".tsv", sep = "")
          cat(fname, "...\n")
          clusters[[t]][[s]][[m]][[o]] <- readSpiciHelper(fname)
        }
      }
    }
  }
  return(clusters)
}

gosimSpici <- readSpici(gosim, naming = NAMING)

addNonClusteredNodes <- function(gosimNodes, gosimSpici) {
  gosimClustersAll <- list()
  for(t in names(gosimSpici)) {
    for(s in names(gosimSpici[[t]])) {
      for(m in names(gosimSpici[[t]][[s]])) {
        for(o in names(gosimSpici[[t]][[s]][[m]])) {
          clusteredNodesDf <- gosimSpici[[t]][[s]][[m]][[o]]
          nonClusteredNodesDf <- data.frame(node = setdiff(gosimNodes[[t]][[s]][[m]][[o]],
                                                           gosimSpici[[t]][[s]][[m]][[o]]$node),
                                            cluster = 9999, stringsAsFactors = FALSE)
          allNodesDf <- rbind(clusteredNodesDf, nonClusteredNodesDf)
          gosimClustersAll[[t]][[s]][[m]][[o]] <- allNodesDf
        }
      }
    }
  }
  return(gosimClustersAll)
}

gosimSpiciExtended <- addNonClusteredNodes(gosimNodes, gosimSpici)

######################################################

writeFormattedSpici <- function(gosimSpiciExtended, naming) {
  for(t in names(gosimSpiciExtended)) {
    for(s in names(gosimSpiciExtended[[t]])) {
      for(m in names(gosimSpiciExtended[[t]][[s]])) {
        for(o in names(gosimSpiciExtended[[t]][[s]][[m]])) {
          df <- gosimSpiciExtended[[t]][[s]][[m]][[o]]
          fname <- paste("RES/SPICi/CLUSTERS/by", naming, "/", t,
                         "_", s, "_", m, "_", o, ".csv", sep = "")
          cat(fname, "...\n")
          write.table(df, fname, row.names = FALSE, col.names = TRUE, sep=",")
        }
      }
    }
  }
}

# writeFormattedSpici(gosimSpiciExtended, naming = NAMING)

######################################################

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

BHIScoresSpici <- readBHIScores("SPICi", topologies, subjects, measures, ontTypes,
                              includeComb, naming = NAMING)

######################################################

drawClusters <- function(gosimObj, gosimSpiciExtended, addNonClusteredNodes, naming, BHIScores) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in c("Wang")) {
      #for(m in names(gosimObj[[t]][[s]])) {
        #for(o in c("BP")) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          nodes <- gosimSpiciExtended[[t]][[s]][[m]][[o]]
          edges <- gosimObj[[t]][[s]][[m]][[o]]
          if(!addNonClusteredNodes) {
            nodes <- nodes[nodes$cluster != 0,]
            edges <- edges[edges[,1] %in% nodes$node, ]
            edges <- edges[edges[,2] %in% nodes$node, ]
          }
          if(nrow(nodes) != 0 && nrow(edges) != 0) {
            net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
            qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
            colors <- sample(col_vector, 95, replace = TRUE)
            V(net)[V(net)$cluster == 0]$color <- "gray90"
            V(net)[V(net)$cluster != 0]$color <- colors[V(net)$cluster]

            fname <- paste("PLOTS/CLUSTERS/SPICi/by", naming , "/", t,
                           "_", s, "_", m, "_", o, ".png", sep="")
            cat(fname, "...\n")
            png(filename=fname, width = 2280, height = 2280)
            lay <- layout_in_circle(net)

            info1 <- paste("Topology: ", toupper(t),
                           "     Subject: ", toupper(s),
                           "     Measure: ", toupper(m),
                           "     Ontology: ", toupper(o), sep="")

            info2 <- paste(
              "Total nodes: " , nrow(gosimSpiciExtended[[t]][[s]][[m]][[o]]),
              "     Clustered nodes: ", length(V(net)[V(net)$cluster != 9999]),
              "     Clusters: ", length(unique(V(net)[V(net)$cluster != 9999]$cluster)), sep="")

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

drawClusters(gosim, gosimSpiciExtended, addNonClusteredNodes = TRUE, NAMING, BHIScoresSpici)
