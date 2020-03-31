source("vars.R")

val_diab_lean <- as.data.frame(fread("RES/GSE133099_DIAB_LEAN.csv", header = TRUE, sep = ','))
val_diab_obese <- as.data.frame(fread("RES/GSE133099_DIAB_OBESE.csv", header = TRUE, sep = ','))

ms_cad <- list("ANAPC2","ARPC1B","COX5A","ENTPD2","FRG1",
          "GLUD2","GSPT2","HCFC1","ISYNA1","LHX2",
          "NFKBIB","PCGF6","POLR2L","RPS9","S100A8",
          "SIX3","TNFRSF13B","TNFSF13","VPS28")


t2d_cad <- list("ALAS2","HBD","SP1")

val_diab_lean_reduced <-
  val_diab_lean[which(abs(val_diab_lean$log2FoldChange) >= 0 &
                        val_diab_lean$pvalue <= 5),]

val_diab_obese_reduced <-
  val_diab_obese[which(abs(val_diab_obese$log2FoldChange) >= 0 &
                         val_diab_obese$pvalue <= 5),]

val_diab_all_reduces <- union(val_diab_lean_reduced[["gene_symbol"]],
           val_diab_obese_reduced[["gene_symbol"]])

Reduce(intersect, list(v1=ms_cad, v2=val_diab_all_reduces))

#####################################################

val_cad <- as.data.frame(fread("RES/GSE109048_CAD.csv", header = TRUE, sep = ','))

