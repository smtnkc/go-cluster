source("vars.R")
library(igraph)
library(qgraph)
library(RColorBrewer)
library(MCL)

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

##########################################################

getAdjMatrices <- function(gosimObj) {
  adjMatrices <- list()
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in names(gosimObj[[t]][[s]])) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          df <- gosimObj[[t]][[s]][[m]][[o]]
          nodes <- union(df[, 1], df[, 2])
          N <- length(nodes)
          print(N)
          adj <- matrix(0, N*N, nrow=N, ncol=N, dimnames = list(nodes, nodes))
          for(r in 1:nrow(df)) {
            adj[df[r,1], df[r,2]] <- df[r, 3]
            adj[df[r,2], df[r,1]] <- df[r, 3]
          }
          adjMatrices[[t]][[s]][[m]][[o]] <- adj
        }
      }
    }
  }
  return(adjMatrices)
}

adjMatrices <- getAdjMatrices(gosim)

##########################################################

getMCL <- function(adjMatrices) {
  gosimClusters <- list()
  for(t in names(adjMatrices)) {
    for(s in names(adjMatrices[[t]])) {
      for(m in names(adjMatrices[[t]][[s]])) {
        for(o in names(adjMatrices[[t]][[s]][[m]])) {
          adj <- adjMatrices[[t]][[s]][[m]][[o]]
          adj <- adj[which(rowSums(adj) > 0), which(colSums(adj) > 0)] # to prevent mcl errors, remove unconnected nodes
          g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE)
          x <- paste("[[\"",t, "\"]]",
                     "[[\"",s, "\"]]",
                     "[[\"",m, "\"]]",
                     "[[\"",o, "\"]]", sep="")
          cat(x, "--> edge:", length(E(g)), "vertex:", length(V(g)))
          mclOutput <- mcl(x = g, addLoops=FALSE)
          cat(" cluster:", mclOutput$K, "\n")
          df <- data.frame(node = V(g)$name, cluster = mclOutput$Cluster)
          gosimClusters[[t]][[s]][[m]][[o]] <- df[order(df$cluster), ]
        }
      }
    }
  }
  return(gosimClusters)
}

gosimMCL <- getMCL(adjMatrices)

##########################################################

writeMCL <- function(gosimMCL) {
  for(t in names(gosimMCL)) {
    for(s in names(gosimMCL[[t]])) {
      for(m in names(gosimMCL[[t]][[s]])) {
        for(o in names(gosimMCL[[t]][[s]][[m]])) {
          df <- gosimMCL[[t]][[s]][[m]][[o]]
          fname <- paste("RES/MCL/CLUSTERS/", t, "_", s, "_", m, "_", o, ".csv", sep = "")
          cat(fname, "...\n")
          write.table(df, fname, row.names = FALSE, col.names = TRUE, sep=",")
        }
      }
    }
  }
}

writeMCL(gosimMCL)

readMCL <- function(topologies, subjects, measures, ontTypes, includeComb) {
  gosimMCL <- list()
  if(includeComb) measures <- c(measures, "Comb")
  
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/MCL/CLUSTERS/", t , "_", s, "_", m, "_", o, ".csv", sep="")
          gosimMCL[[t]][[s]][[m]][[o]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
        }
      }
    }
  }
  return(gosimMCL)
}

gosimMCL <- readMCL(topologies, subjects, measures, ontTypes, includeComb = FALSE)

##########################################################

drawClusters <- function(gosimObj, gosimMCL) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in c("Wang")) {
        #for(m in names(gosimObj[[t]][[s]])) {
        #for(o in c("BP")) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          nodes <- gosimMCL[[t]][[s]][[m]][[o]]
          edges <- gosimObj[[t]][[s]][[m]][[o]]
          
          if(nrow(nodes) != 0 && nrow(edges) != 0) {
            net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
            qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
            colors <- sample(col_vector, 95, replace = TRUE)
            V(net)[V(net)$cluster == 0]$color <- "gray90"
            V(net)[V(net)$cluster != 0]$color <- colors[V(net)$cluster]
            
            fname <- paste("PLOTS/CLUSTERS/MCL/", t, "_", s, "_", m, "_", o, ".png", sep="")
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

drawClusters(gosim, gosimMCL)
