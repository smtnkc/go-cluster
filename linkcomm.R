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

##########################################################

getLinkcomm <- function(gosim) {
  gosimClusters <- list()
  cat("topology,subject,measure,ontology,edge,node,cluster,time(s),memory(mb)\n")
  for(t in names(gosim)) {
    for(s in names(gosim[[t]])) {
      for(m in names(gosim[[t]][[s]])) {
        for(o in names(gosim[[t]][[s]][[m]])) {
          df <- gosim[[t]][[s]][[m]][[o]]
          t_start <- Sys.time()
          p <- profmem({
          lc <- getLinkCommunities(gosim[[t]][[s]][[m]][[o]], plot = FALSE, verbose = FALSE)
          })
          t_end <- Sys.time()
          dfLC <- lc[["nodeclusters"]]
          dfLC <- transform(dfLC, node = as.character(node))
          dfLC <- transform(dfLC, cluster = as.integer(cluster))
          gosimClusters[[t]][[s]][[m]][[o]] <- dfLC
          cat(paste(t,s,m,o,
                    lc[["numbers"]][1],
                    lc[["numbers"]][2],
                    lc[["numbers"]][3],
                    round(difftime(t_end, t_start, units = "secs"),5),
                    total(p)/1000000, sep=","), "\n")
        }
      }
    }
  }
  return(gosimClusters)
}

# gosimLinkcomm <- getLinkcomm(gosim)

##########################################################

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

# gosimLinkcommExtended <- addNonClusteredNodes(gosimNodes, gosimLinkcomm)

##########################################################

writeLinkcomm <- function(gosimLinkcommExtended, naming) {
  for(t in names(gosimLinkcommExtended)) {
    for(s in names(gosimLinkcommExtended[[t]])) {
      for(m in names(gosimLinkcommExtended[[t]][[s]])) {
        for(o in names(gosimLinkcommExtended[[t]][[s]][[m]])) {
          df <- gosimLinkcommExtended[[t]][[s]][[m]][[o]]
          fname <- paste("RES/LINKCOMM/CLUSTERS/by", naming, "/", t,
                         "_", s, "_", m, "_", o, ".csv", sep = "")
          cat(fname, "...\n")
          write.table(df, fname, row.names = FALSE, col.names = TRUE, sep=",")
        }
      }
    }
  }
}

# writeLinkcomm(gosimLinkcommExtended, naming = NAMING)

##########################################################

readLinkcomm <- function(topologies, subjects, measures, ontTypes, includeComb, naming) {
  gosimLinkcomm <- list()
  if(includeComb) measures <- c(measures, "Comb")

  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/LINKCOMM/CLUSTERS/by", naming, "/",
                         t, "_", s, "_", m, "_", o, ".csv", sep="")
          gosimLinkcomm[[t]][[s]][[m]][[o]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
        }
      }
    }
  }
  return(gosimLinkcomm)
}

gosimLinkcommExtended <- readLinkcomm(topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)

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

BHIScoresLinkcomm <- readBHIScores("LINKCOMM", topologies, subjects, measures, ontTypes,
                                   includeComb = FALSE, naming = NAMING)

######################################################

drawClusters <- function(gosimObj, gosimLinkcommExtended, addNonClusteredNodes, naming, BHIScores) {
  for(t in names(gosimObj)) {
    for(s in names(gosimObj[[t]])) {
      for(m in c("Wang")) {
        #for(m in names(gosimObj[[t]][[s]])) {
        #for(o in c("BP")) {
        for(o in names(gosimObj[[t]][[s]][[m]])) {
          nodes <- gosimLinkcommExtended[[t]][[s]][[m]][[o]]
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

            fname <- paste("PLOTS/CLUSTERS/LINKCOMM/by", naming , "/",
                           t, "_", s, "_", m, "_", o, ".png", sep="")
            cat(fname, "...\n")
            png(filename=fname, width = 2280, height = 2280)
            lay <- layout_in_circle(net)

            info1 <- paste("Topology: ", toupper(t),
                           "     Subject: ", toupper(s),
                           "     Measure: ", toupper(m),
                           "     Ontology: ", toupper(o), sep="")

            info2 <- paste(
              "Total nodes: " , nrow(gosimLinkcommExtended[[t]][[s]][[m]][[o]]),
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

# cannot draw plots because of duplicated nodes (Linkcomm has overlapping clusters)
# drawClusters(gosim, gosimLinkcommExtended, addNonClusteredNodes = TRUE, NAMING, BHIScoresLinkcomm)
