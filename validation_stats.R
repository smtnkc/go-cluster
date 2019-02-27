bhi.MCL <- as.data.frame(fread("RES/MCL/SCORES.csv", header = TRUE, sep = ','))
bhi.SPICi <- as.data.frame(fread("RES/SPICi/SCORES.csv", header = TRUE, sep = ','))

##### Means of ontologies

o.means.ALL <- data.frame(rbind(
  aggregate(bhi.MCL[, "BHI"], list(bhi.MCL$ontology), mean),
  aggregate(bhi.SPICi[, "BHI"], list(bhi.SPICi$ontology), mean)))

o.means.ALL$clustering <- c(rep("MCL", 3), rep("SPICi", 3))
colnames(o.means.ALL) <- c("ontology", "Mean_BHI", "clustering")
o.means.ALL <- o.means.ALL[,c("clustering", "ontology", "Mean_BHI")]
o.means.ALL$Mean_BHI <- round(o.means.ALL$Mean_BHI, 3)

##### Means of measures

m.means.ALL <- data.frame(rbind(
  aggregate(bhi.MCL[, "BHI"], list(bhi.MCL$measure), mean),
  aggregate(bhi.SPICi[, "BHI"], list(bhi.SPICi$measure), mean)))

m.means.ALL$clustering <- c(rep("MCL", 5), rep("SPICi", 5))
colnames(m.means.ALL) <- c("measure", "Mean_BHI", "clustering")
m.means.ALL <- m.means.ALL[,c("clustering", "measure", "Mean_BHI")]
m.means.ALL$Mean_BHI <- round(m.means.ALL$Mean_BHI, 3)

###### Means of topologies

t.means.ALL <- data.frame(rbind(
  aggregate(bhi.MCL[, "BHI"], list(bhi.MCL$topology), mean),
  aggregate(bhi.SPICi[, "BHI"], list(bhi.SPICi$topology), mean)))

t.means.ALL$clustering <- c(rep("MCL", 2), rep("SPICi", 2))
colnames(t.means.ALL) <- c("topology", "Mean_BHI", "clustering")
t.means.ALL <- t.means.ALL[,c("clustering", "topology", "Mean_BHI")]
t.means.ALL$Mean_BHI <- round(t.means.ALL$Mean_BHI, 3)

###### Write mean stats

write.table(o.means.ALL, "RES/validation_stats/byOntology.csv",
            row.names = FALSE, col.names = TRUE, sep=",")

write.table(m.means.ALL, "RES/validation_stats/byMeasure.csv",
            row.names = FALSE, col.names = TRUE, sep=",")

write.table(t.means.ALL, "RES/validation_stats/byTopology.csv",
            row.names = FALSE, col.names = TRUE, sep=",")
