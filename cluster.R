source("vars.R")
library(igraph)
library(qgraph)
library(RColorBrewer)

gosim <- readGosim(topologies, subjects, measures, ontTypes, includeComb = FALSE)

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

readClustersHelper <- function(clusterFile) {
  fileContent <- readLines(clusterFile)
  clusterDF <- data.frame(node = character(), cluster = integer())
  if(length(fileContent) == 0) {
    return(clusterDF)
  }
  c <- 1
  for(i in 1:length(fileContent)) {
    for(n in strsplit(gsub('\"', '', fileContent[i]), split="\t")) {
      clusterDF <- rbind(clusterDF, data.frame(node = n, cluster = c, stringsAsFactors = FALSE))
    }
    c <- c + 1
  }
  return(clusterDF)
}

readClusters <- function(gosimObj) {
  clusters <- list()
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          fname <- paste("RES/SPICi/CLUSTERS/",
                        t, "_", s, "_", m, "_", o, ".tsv", sep = "")
          cat(fname, "...\n")
          clusters[[t]][[s]][[m]][[o]] <- readClustersHelper(fname)
        }
      }
    }
  }
  return(clusters)
}

gosimClusters <- readClusters(gosim)

addNonClusteredNodes <- function(gosimNodes, gosimClusters) {
  gosimClustersAll <- list()
  for(t in names(gosimClusters)) {
    for(s in names(gosimClusters[[t]])) {
      for(m in names(gosimClusters[[t]][[s]])) {
        for(o in names(gosimClusters[[t]][[s]][[m]])) {
          clusteredNodesDf <- gosimClusters[[t]][[s]][[m]][[o]]
          nonClusteredNodesDf <- data.frame(node = setdiff(gosimNodes[[t]][[s]][[m]][[o]],
                                                           gosimClusters[[t]][[s]][[m]][[o]]$node),
                                            cluster = 0, stringsAsFactors = FALSE)
          allNodesDf <- rbind(clusteredNodesDf, nonClusteredNodesDf)
          gosimClustersAll[[t]][[s]][[m]][[o]] <- allNodesDf
        }
      }
    }
  }
  return(gosimClustersAll)
}

gosimClustersAll <- addNonClusteredNodes(gosimNodes, gosimClusters)

######################################################

drawClusters <- function(gosimObj, gosimClustersAll, addNonClusteredNodes) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      #for(m in c("Wang")) {
      for(m in names(gosimObj[[t]][[s]])) {
        #for(o in c("BP")) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          edges <- gosimObj[[t]][[s]][[m]][[o]]
          if(addNonClusteredNodes && s != "cad") {
            nodes <- gosimClustersAll[[t]][[s]][[m]][[o]]
          } else {
            nodes <- gosimClusters[[t]][[s]][[m]][[o]]
            edges <- edges[edges$symbol1 %in% nodes$node, ]
            edges <- edges[edges$symbol2 %in% nodes$node, ]
          }
          if(nrow(nodes) != 0 && nrow(edges) != 0) {
            net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
            qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
            colors <- sample(col_vector, 95, replace = TRUE)
            V(net)[V(net)$cluster == 0]$color <- "gray90"
            V(net)[V(net)$cluster != 0]$color <- colors[V(net)$cluster]
            
            fname <- paste("PLOTS/CLUSTERS/", t, "_", s, "_", m, "_", o, ".png", sep="")
            cat(fname, "...\n")
            png(filename=fname, width = 960, height = 960)
            lay <- layout_in_circle(net)
            
            info1 <- paste("Topology: ", toupper(t),
                           "     Subject: ", toupper(s),
                           "     Measure: ", toupper(m),
                           "     Ontology: ", toupper(o), sep="")
            
            info2 <- paste(
              "Total nodes: " , nrow(gosimClustersAll[[t]][[s]][[m]][[o]]),
              "     Clustered nodes: ", length(V(net)[V(net)$cluster != 0]),
              "     Clusters: ", length(unique(V(net)[V(net)$cluster != 0]$cluster)), sep="")
            
            plot(net, main=info1, xlab=info2, layout=lay,
                 vertex.label=nodes$node, vertex.shape="circle",
                 vertex.size=4, #vertex.size=nchar(as.character(nodes$node))*4+5, vertex.size2=10,
                 vertex.label.color="black", vertex.label.cex=0.5,
                 vertex.frame.color="gray50", vertex.label.family="Helvetica",
                 edge.width=1, edge.label.color="black",
                 edge.color="gray50", edge.curved=0, edge.label.family="Helvetica")
            
            dev.off()
          }
        }
      }
    }
  }
}

drawClusters(gosim, gosimClustersAll, addNonClusteredNodes = TRUE)
