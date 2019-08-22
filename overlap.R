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

clusters_S <- list(
  "ms" = as.data.frame(fread(paste("RES/SPICi/CLUSTERS/by", NAMING, "/inet_mets_Wang_MF.csv", sep=""), header = TRUE, sep = ',')),
  "t2d" = as.data.frame(fread(paste("RES/SPICi/CLUSTERS/by", NAMING, "/inet_t2d_Wang_MF.csv", sep=""), header = TRUE, sep = ',')),
  "cad" = as.data.frame(fread(paste("RES/SPICi/CLUSTERS/by", NAMING, "/inet_cad_Wang_MF.csv", sep=""), header = TRUE, sep = ','))
)

clusters_L <- list(
  "ms" = as.data.frame(fread(paste("RES/LINKCOMM/CLUSTERS/by", NAMING, "/inet_mets_Wang_MF.csv", sep=""), header = TRUE, sep = ',')),
  "t2d" = as.data.frame(fread(paste("RES/LINKCOMM/CLUSTERS/by", NAMING, "/inet_t2d_Wang_MF.csv", sep=""), header = TRUE, sep = ',')),
  "cad" = as.data.frame(fread(paste("RES/LINKCOMM/CLUSTERS/by", NAMING, "/inet_cad_Wang_MF.csv", sep=""), header = TRUE, sep = ','))
)

removeUnclusteredNodes <- function(clusters_X) {
  clusters_Y <- list()
  for(m in names(clusters_X)) {
    clusters_Y[[m]] <- clusters_X[[m]][clusters_X[[m]]$cluster != 9999,]
  }
  return(clusters_Y)
}

# clusters_S <- removeUnclusteredNodes(clusters_S)
# clusters_L <- removeUnclusteredNodes(clusters_L)

getOverlaps <- function(clusters_X) {
  overlaps <- list()
  ms_t2d  <- sort(intersect(clusters_X[["ms"]][,1], clusters_X[["t2d"]][,1]))
  ms_cad  <- sort(intersect(clusters_X[["ms"]][,1], clusters_X[["cad"]][,1]))
  t2d_cad <- sort(intersect(clusters_X[["cad"]][,1], clusters_X[["t2d"]][,1]))
  # cat("ms_t2d -> ", ms_t2d, "(", length(ms_t2d),")\n")
  # cat("ms_cad -> ", ms_cad, "(", length(ms_cad),")\n")
  # cat("t2d_cad -> ", t2d_cad, "(", length(t2d_cad),")\n")
  overlaps[["ms_t2d"]] <- ms_t2d
  overlaps[["ms_cad"]] <- ms_cad
  overlaps[["t2d_cad"]] <- t2d_cad
  return(overlaps)
}

overlaps_S <- getOverlaps(clusters_S)
overlaps_L <- getOverlaps(clusters_L)

#####################################################################

### Linkcomm_MS_CAD
g <- gosim[["inet"]][["cad"]][["Wang"]][["MF"]]
g <- g[g$symbol1 %in% overlaps_L[["ms_cad"]] & g$symbol2 %in% overlaps_L[["ms_cad"]], ]
clust1 <- clusters_L[["ms"]][clusters_L[["ms"]]$node %in% overlaps_L[["ms_cad"]],]
clust2 <- clusters_L[["cad"]][clusters_L[["cad"]]$node %in% overlaps_L[["ms_cad"]],]
net <- graph_from_data_frame(d=g, vertices=overlaps_L[["ms_cad"]], directed=F)

### Linkcomm_T2D_CAD
# g <- gosim[["inet"]][["cad"]][["Wang"]][["MF"]]
# g <- g[g$symbol1 %in% overlaps_L[["t2d_cad"]] & g$symbol2 %in% overlaps_L[["t2d_cad"]], ]
# clust1 <- clusters_L[["t2d"]][clusters_L[["t2d"]]$node %in% overlaps_L[["t2d_cad"]],]
# clust2 <- clusters_L[["cad"]][clusters_L[["cad"]]$node %in% overlaps_L[["t2d_cad"]],]
# net <- graph_from_data_frame(d=g, vertices=overlaps_L[["t2d_cad"]], directed=F)

### SPICi_MS_CAD
# g <- gosim[["inet"]][["cad"]][["Wang"]][["MF"]]
# g <- g[g$symbol1 %in% overlaps_S[["ms_cad"]] & g$symbol2 %in% overlaps_S[["ms_cad"]], ]
# clust1 <- clusters_S[["ms"]][clusters_S[["ms"]]$node %in% overlaps_S[["ms_cad"]],]
# clust2 <- clusters_S[["cad"]][clusters_S[["cad"]]$node %in% overlaps_S[["ms_cad"]],]
# net <- graph_from_data_frame(d=g, vertices=overlaps_S[["ms_cad"]], directed=F)

### SPICi_T2D_CAD
# g <- gosim[["inet"]][["cad"]][["Wang"]][["MF"]]
# t2d_S <- clusters_S[[2]][clusters_S[[2]]$cluster != 9999,]
# t2d_L <- clusters_L[[2]][clusters_L[[2]]$cluster != 9999,]
# g <- g[g$symbol1 %in% t2d_S$node & g$symbol2 %in%  t2d_L$node, ]
# net <- graph_from_data_frame(d=g, vertices=t2d_L$node, directed=F)
# net <- graph_from_data_frame(d=g, vertices=t2d_L$node, directed=F)

# g <- g[g$symbol1 %in% overlaps_S[["t2d_cad"]] & g$symbol2 %in% overlaps_S[["t2d_cad"]], ]
# clust1 <- clusters_S[["t2d"]][clusters_S[["t2d"]]$node %in% overlaps_S[["t2d_cad"]],]
# clust2 <- clusters_S[["cad"]][clusters_S[["cad"]]$node %in% overlaps_S[["t2d_cad"]],]
# net <- graph_from_data_frame(d=g, vertices=overlaps_S[["t2d_cad"]], directed=F)

g
clust1
clust2

plot(net, vertex.size=20, layout=layout_nicely(net), edge.label="", edge.width=g$MF*10,
     vertex.label.color= "black", vertex.label.family="Helvetica", edge.label=g$MF)


########### CAD Cluster Tables for Supplementary Materials

cadL <- aggregate(node ~ cluster, data = clusters_L[[3]][clusters_L[[3]]$cluster != 9999,], c)
cadS <- aggregate(node ~ cluster, data = clusters_S[[3]][clusters_S[[3]]$cluster != 9999,], c)
for(i in 1:nrow(cadL)) { cadL[i, "node"] <- paste0(sort(cadL[i, "node"][[1]]), collapse=",") }
for(i in 1:nrow(cadS)) { cadS[i, "node"] <- paste0(sort(cadS[i, "node"][[1]]), collapse=",") }
cadL <- transform(cadL, node = as.character(node))
cadS <- transform(cadS, node = as.character(node))
write.table(cadL, "RES/cadL.txt", sep="\t", row.names = FALSE)
write.table(cadS, "RES/cadS.txt", sep="\t", row.names = FALSE)
