library(ggplot2)
library(ggpubr)
bhi.MCL <- as.data.frame(fread("RES/MCL/SCORES.csv", header = TRUE, sep = ','))
bhi.SPICi <- as.data.frame(fread("RES/SPICi/SCORES.csv", header = TRUE, sep = ','))
bhi.Linkcomm <- as.data.frame(fread("RES/LINKCOMM/SCORES.csv", header = TRUE, sep = ','))

##### Means of ontologies

# o.means.ALL <- data.frame(rbind(
#   aggregate(bhi.MCL[, "BHI"], list(bhi.MCL$ontology), mean),
#   aggregate(bhi.SPICi[, "BHI"], list(bhi.SPICi$ontology), mean),
#   aggregate(bhi.Linkcomm[, "BHI"], list(bhi.Linkcomm$ontology), mean)))
# 
# o.means.ALL$clustering <- c(rep("MCL", 3), rep("SPICi", 3), rep("LINKCOMM", 3))
# colnames(o.means.ALL) <- c("ontology", "Mean_BHI", "clustering")
# o.means.ALL <- o.means.ALL[,c("clustering", "ontology", "Mean_BHI")]
# o.means.ALL$Mean_BHI <- round(o.means.ALL$Mean_BHI, 5)

##### Means of measures

# m.means.ALL <- data.frame(rbind(
#   aggregate(bhi.MCL[, "BHI"], list(bhi.MCL$measure), mean),
#   aggregate(bhi.SPICi[, "BHI"], list(bhi.SPICi$measure), mean),
#   aggregate(bhi.Linkcomm[, "BHI"], list(bhi.Linkcomm$measure), mean)))
# 
# m.means.ALL$clustering <- c(rep("MCL", 5), rep("SPICi", 5), rep("LINKCOMM", 5))
# colnames(m.means.ALL) <- c("measure", "Mean_BHI", "clustering")
# m.means.ALL <- m.means.ALL[,c("clustering", "measure", "Mean_BHI")]
# m.means.ALL$Mean_BHI <- round(m.means.ALL$Mean_BHI, 5)

###### Means of topologies

# t.means.ALL <- data.frame(rbind(
#   aggregate(bhi.MCL[, "BHI"], list(bhi.MCL$topology), mean),
#   aggregate(bhi.SPICi[, "BHI"], list(bhi.SPICi$topology), mean),
#   aggregate(bhi.Linkcomm[, "BHI"], list(bhi.Linkcomm$topology), mean)))
# 
# t.means.ALL$clustering <- c(rep("MCL", 2), rep("SPICi", 2), rep("LINKCOMM", 2))
# colnames(t.means.ALL) <- c("topology", "Mean_BHI", "clustering")
# t.means.ALL <- t.means.ALL[,c("clustering", "topology", "Mean_BHI")]
# t.means.ALL$Mean_BHI <- round(t.means.ALL$Mean_BHI, 5)

###### Write mean stats

# write.table(o.means.ALL, "RES/validation_stats/byOntology.csv",
#             row.names = FALSE, col.names = TRUE, sep=",")
# 
# write.table(m.means.ALL, "RES/validation_stats/byMeasure.csv",
#             row.names = FALSE, col.names = TRUE, sep=",")
# 
# write.table(t.means.ALL, "RES/validation_stats/byTopology.csv",
#             row.names = FALSE, col.names = TRUE, sep=",")

################################################################

# getBoxPlot <- function(bhiDf, aesX, aesY, aesFill, s) {
#   p <-  ggplot(bhiDf, aes(x=bhiDf[,aesX], y=bhiDf[,aesY], fill=bhiDf[,aesFill])) +
#         geom_boxplot() + 
#         theme_minimal() +
#         ylim(0.15, 0.5) +
#         labs(x = aesX, y = aesY, fill = aesFill) +
#         theme(legend.position="top",
#               axis.text=element_text(size=s-4),
#               axis.title.x=element_text(size=s),
#               axis.title.y=element_text(size=s),
#               legend.title=element_text(size=s), 
#               legend.text=element_text(size=s))
#   return(p)
# }
# 
# png(filename="PLOTS/VALIDATION/topology.png", width = 2400, height = 960)
# p1 <- getBoxPlot(bhi.MCL, "subject", "BHI", "topology", 18)
# p2 <- getBoxPlot(bhi.SPICi, "subject", "BHI", "topology", 18)
# p3 <- getBoxPlot(bhi.Linkcomm, "subject", "BHI", "topology", 18)
# ggarrange(p1, p2, p3, labels=c("MCL", "SPICi", "LinkComm"), font.label=list(size=22), ncol=3, nrow=1)
# dev.off()

##############################

# png(filename="PLOTS/VALIDATION/ontology.png", width = 2400, height = 960)
# p1 <- getBoxPlot(bhi.MCL, "subject", "BHI", "ontology", 18)
# p2 <- getBoxPlot(bhi.SPICi, "subject", "BHI", "ontology", 18)
# p3 <- getBoxPlot(bhi.Linkcomm, "subject", "BHI", "ontology", 18)
# ggarrange(p1, p2, p3, labels=c("MCL", "SPICi", "LinkComm"), font.label=list(size=22), ncol=3, nrow=1)
# dev.off()

##############################

# png(filename="PLOTS/VALIDATION/measure.png", width = 2400, height = 960)
# p1 <- getBoxPlot(bhi.MCL, "subject", "BHI", "measure", 18)
# p2 <- getBoxPlot(bhi.SPICi, "subject", "BHI", "measure", 18)
# p3 <- getBoxPlot(bhi.Linkcomm, "subject", "BHI", "measure", 18)
# ggarrange(p1, p2, p3, labels=c("MCL", "SPICi", "LinkComm"), font.label=list(size=22), ncol=3, nrow=1)
# dev.off()

##########################################################################################

bhi.MCL$clustering <- c(rep("MCL", nrow(bhi.MCL)))
bhi.SPICi$clustering <- c(rep("SPICi", nrow(bhi.SPICi)))
bhi.Linkcomm$clustering <- c(rep("Linkcomm", nrow(bhi.Linkcomm)))
bhi.ALL <- rbind(bhi.MCL, bhi.SPICi, bhi.Linkcomm)

MEAS <- "Wang"

bhi.ALL.MEAS.mets <- bhi.ALL[bhi.ALL$measure == MEAS & bhi.ALL$subject == "mets", ]
bhi.ALL.MEAS.t2d  <- bhi.ALL[bhi.ALL$measure == MEAS & bhi.ALL$subject == "t2d", ]
bhi.ALL.MEAS.cad  <- bhi.ALL[bhi.ALL$measure == MEAS & bhi.ALL$subject == "cad", ]

getScatterPlot <- function(bhiDf, aesX, aesY, aesShape, aesColor, s) {
  p <- ggplot(bhiDf, aes(x=bhiDf[,aesX], y=bhiDf[,aesY])) +
    geom_point(aes(shape=bhiDf[,aesShape], color=bhiDf[,aesColor]), size=5) +
    scale_shape_manual(values=c(16, 17))+
    scale_color_manual(values=c('#EA2027', '#0652DD', '#FFC312'))+
    labs(x = aesX, y = aesY, shape = aesShape, color = aesColor) +
    ylim(0.2, 0.5) +
    theme(legend.position="top",
          axis.text=element_text(size=s-4),
          axis.title.x=element_text(size=s),
          axis.title.y=element_text(size=s),
          legend.title=element_text(size=s-2), 
          legend.text=element_text(size=s-4))
  return(p)
}

p1 <- getScatterPlot(bhi.ALL.MEAS.mets, aesX="clustering", aesY="BHI", aesShape="topology", aesColor="ontology", s=18)
p2 <- getScatterPlot(bhi.ALL.MEAS.t2d, aesX="clustering", aesY="BHI", aesShape="topology", aesColor="ontology", s=18)
p3 <- getScatterPlot(bhi.ALL.MEAS.cad, aesX="clustering", aesY="BHI", aesShape="topology", aesColor="ontology", s=18)

png(filename=paste("PLOTS/VALIDATION/BHI_", MEAS, ".png", sep=""), width = 1800, height = 960)
ggarrange(p1, p2, p3, labels=c("METS", "T2D", "CAD"), font.label=list(size=22), ncol=3, nrow=1)
dev.off()
