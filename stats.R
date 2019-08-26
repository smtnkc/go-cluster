source("vars.R")

##############################################

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
printGosimCounts <- function(gosim) {
  df <- data.frame(topology = as.character(),
                   subject = as.character(),
                   nCount = as.double(),
                   eCount = as.double())
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
      df <- rbind(df, data.frame(topology = t,
                                 subject = s,
                                 nCount = n,
                                 eCount = e,
                                 stringsAsFactors = FALSE))
    }
  }
  return(df)
}
gosimCounts <- printGosimCounts(gosim)

##############################################

NAMING <- "SYMBOL"
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
clusters <- list(
  "mcl" = readClusters("MCL", topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING),
  "spici" = readClusters("SPICi", topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING),
  "linkcomm" = readClusters("LINKCOMM", topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = NAMING)
)
printClusterStats <- function(clusters) {
  dfAvg <- data.frame(algorithm = as.character(), topology = as.character(),
                   subject = as.character(), nn_avg = as.double(), nc_avg = as.double(),
                   d_avg = as.double(), stringsAsFactors = FALSE)

  dfAll <- data.frame(algorithm = as.character(), topology = as.character(),
                      subject = as.character(), measure = as.character(), ontology = as.character(),
                      nn = as.double(), nc = as.double(), d = as.double(), stringsAsFactors = FALSE)

  for(c in names(clusters)) {
    for(t in names(clusters[[c]])) {
      for(s in names(clusters[[c]][[t]])) {
        nn_total<-0
        nc_total<-0
        n<-0
        #for(m in names(clusters[[c]][[t]][[s]])) {
          #for(o in names(clusters[[c]][[t]][[s]][[m]])) {
        for(m in c("Wang")) {
          for(o in c("MF")) {
            df <- clusters[[c]][[t]][[s]][[m]][[o]]
            nn_temp <- length(unique(df[df$cluster != 9999,]$node))
            nc_temp <- length(unique(df[df$cluster != 9999,]$cluster))
            if(nc_temp == 0) {
              d_temp <- 0
            }
            else {
              d_temp <- round(nn_temp/nc_temp,1)
            }

            dfAll <- rbind(dfAll, data.frame(algorithm = c, topology = t,
                                             subject = s, measure = m, ontology = o,
                                             nn = nn_temp, nc = nc_temp, d = d_temp, stringsAsFactors = FALSE))

            nn_total <- nn_total + nn_temp
            nc_total <- nc_total + nc_temp
            n <- n+1
          }
        }
        nn_avg <- round(nn_total/n,1)
        nc_avg <- round(nc_total/n,1)
        if(nc_total == 0) {
          d_avg <- 0
        }
        else {
          d_avg <- round(nn_total/nc_total,1)
        }
        print(paste(c,"_",t,"_",s,"      ",nn_avg,"/",nc_avg,"=",d_avg, sep=""))
        dfAvg <- rbind(dfAvg, data.frame(algorithm = c, topology = t, subject = s, nn_avg = nn_avg,
                                         nc_avg = nc_avg, d_avg = d_avg, stringsAsFactors = FALSE))
        colnames(dfAvg) <- c("algorithm", "topology", "subject", "nn_avg", "nc_avg", "d_avg")
        dfAvg <- transform(dfAvg, nn_avg = as.numeric(nn_avg), nc_avg = as.numeric(nc_avg), d_avg = as.numeric(d_avg))
        colnames(dfAll) <- c("algorithm", "topology", "subject", "measure", "ontology", "nn", "nc", "d")
        dfAll <- transform(dfAll, nn = as.numeric(nn), nc = as.numeric(nc), d = as.numeric(d))
      }
    }
  }
  return(dfAll)
  # return(dfAvg)
}
dfClusterStats <- printClusterStats(clusters)

##################### DRAW DOTPLOTS TO COMPARE MCL, SPICI, and LINKCOMM:
convertClusterStats <- function(dfClusterStats, gosimCounts) {
  dfSimple <- dfClusterStats[,c(1,2,3,6,7,8)]

  dfRes <- data.frame(algorithm = as.character(),
                      topology = as.character(),
                      subject = as.character(),
                      statType = as.character(),
                      statVal = as.double(),
                      stringsAsFactors = FALSE)

  for(a in c("mcl", "spici", "linkcomm")) {
    for(t in c("string", "inet")) {
      for(s in c("mets", "t2d", "cad")) {
        for(stty in c("nn", "nc", "d")) {
          vals <- dfSimple[dfSimple$algorithm == a &
                             dfSimple$topology == t &
                             dfSimple$subject == s, stty]

          if(stty== "nn") {
            vals <- 100*vals/gosimCounts[gosimCounts$topology==t & gosimCounts$subject==s, "nCount"]
          }
          else if(stty== "nc") {
            vals <- vals/gosimCounts[gosimCounts$topology==t & gosimCounts$subject==s, "nCount"]
          }
          dfRes <- rbind(dfRes, data.frame(algorithm = a,
                                             topology = t,
                                             subject = s,
                                             statType = stty,
                                             statVal = vals,
                                             stringsAsFactors = FALSE))
        }
      }
    }
  }
  return(dfRes)
}

dfcs <- convertClusterStats(dfClusterStats, gosimCounts)


generateDotPlots <- function(dfClusterStats, gosimCounts) {
  dfcs <- convertClusterStats(dfClusterStats, gosimCounts)
  dfcs[dfcs$subject=="mets", "subject"] <- "MS"
  dfcs[dfcs$subject=="cad", "subject"] <- "CAD"
  dfcs[dfcs$subject=="t2d", "subject"] <- "T2D"
  dfcs[dfcs$topology=="string", "topology"] <- "STRING"
  dfcs[dfcs$topology=="inet", "topology"] <- "INet"
  plotList <- list()
    for(stt in c("nn", "nc", "d")) {
      if(stt == "nn") { ttl <- "Coverage Rate (%)"}
      else if(stt == "nc") {ttl <- "Number of Clusters (per node)"}
      else { ttl <- "Cluster Size" }

      plotList[[stt]] <-
        ggplot(dfcs[dfcs$statType==stt,], aes(x = algorithm, y = statVal)) +
        geom_point(aes(fill=algorithm, size=3), colour="#666666", shape=21, stroke = 0.5) +
        scale_fill_manual(values=c(globalDarkColorPalette[3:5])) +
        facet_grid(subject ~ topology) + labs(title = ttl) +
        theme_bw() +
        theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(),
              plot.title = element_text(hjust = 0.5, size = 20, face="bold", margin=margin(0,0,10,0)),
              axis.text=element_text(size=16),
              axis.text.x=element_text(angle = 45, hjust = 1, face="bold"),
              strip.text.x = element_text(size = 20),
              strip.text.y = element_text(size = 20))

    }
  return(plotList)
}
pl <- generateDotPlots(dfClusterStats, gosimCounts)
svg(filename = paste("PLOTS/dot_clusters.svg", sep=""), width=14, height=9)
grid.arrange(grobs = pl, ncol = 3, nrow = 1)
dev.off()

#########################################################################################

bhi.MCL <- as.data.frame(fread("RES/MCL/SCORES.csv", header = TRUE, sep = ','))
bhi.SPICi <- as.data.frame(fread("RES/SPICi/SCORES.csv", header = TRUE, sep = ','))
bhi.Linkcomm <- as.data.frame(fread("RES/LINKCOMM/SCORES.csv", header = TRUE, sep = ','))

bhi.MCL$clustering <- c(rep("MCL", nrow(bhi.MCL)))
bhi.SPICi$clustering <- c(rep("SPICi", nrow(bhi.SPICi)))
bhi.Linkcomm$clustering <- c(rep("Linkcomm", nrow(bhi.Linkcomm)))
bhi.ALL <- rbind(bhi.MCL, bhi.SPICi, bhi.Linkcomm)

bhi.ALL[bhi.ALL$topology=="inet" & bhi.ALL$measure=="Wang" & bhi.ALL$ontology=="MF",] %>% group_by(subject,clustering) %>% summarise(BHI_mean=(mean(BHI)))
