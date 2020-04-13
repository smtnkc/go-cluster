if (!require('readr')) install.packages('readr')
if (!require('rstudioapi')) install.packages('rstudioapi')
if (!require('data.table')) install.packages('data.table')

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

t2d_cad <- list("ALAS2","HBD","SP1")
ms_cad  <- list("ANAPC2","ARPC1B","COX5A","ENTPD2","FRG1",
                "GLUD2","GSPT2","HCFC1","ISYNA1","LHX2",
                "NFKBIB","PCGF6","POLR2L","RPS9","S100A8",
                "SIX3","TNFRSF13B","TNFSF13","VPS28")

getSigDf <- function(df, pval, controlRange, caseRange) {
  sigGenes <- c()
  invalidRows <- c()
  for(r in 1:nrow(df)) {
    if(length(levels(factor(as.numeric(df[r, ])))) == 1) {
      #### eliminate the constant rows to overcome "t-test data are essentially constant" error
      invalidRows <- c(invalidRows, row.names(df)[r])
    }
    else if(t.test(df[r, controlRange], df[r, caseRange])$p.val <= pval) {
      sigGenes <- c(sigGenes, row.names(df)[r])
    }
  }
  cat("There are", length(invalidRows), "invalid rows:", invalidRows, "\n")
  sigDf <- subset(df, row.names(df) %in% sigGenes)
  return(sigDf)
}

getDegDf <- function(df, FC, controlRange, caseRange) {
  degGenes <- c()
  for(r in 1:nrow(df)) {
    controlMean <- apply(df[r, controlRange], 1, mean)
    caseMean <- apply(df[r, caseRange], 1, mean)
    if(abs(controlMean-caseMean) >= FC) {
      degGenes <- c(degGenes, row.names(df)[r])
    }
  }
  degDf <- subset(df, row.names(df) %in% degGenes)
  return(degDf)
}

getOverlap <- function(valDf, gpl, gplIDCol, parserRegex) {
  gene_symbols <- gpl[gpl[,gplIDCol] %in% row.names(valDf), 2]
  if(!is.null(parserRegex)) {
    gene_symbols_parsed <- c()
    for(e in gene_symbols) {
      for(x in strsplit(e, parserRegex)) {
        gene_symbols_parsed <- c(gene_symbols_parsed, x)
      }
    }
    overlap <- Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=gene_symbols_parsed))
  }
  else {
    overlap <- Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=gene_symbols))
  }
  return(overlap)
}


######################################## GSE133099 ###############################

gse_lean <- as.data.frame(fread("VALS/GSE133099_DIAB_LEAN.csv", header = TRUE, sep = ','))
gse_obese <- as.data.frame(fread("VALS/GSE133099_DIAB_OBESE.csv", header = TRUE, sep = ','))
gse_all <- unique(union(gse_lean[["gene_symbol"]], gse_obese[["gene_symbol"]]))
overlap <- Reduce(intersect, list(v1=union(ms_cad,t2d_cad), v2=gse_all))
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)


sigs_lean <- gse_lean[which(gse_lean$pvalue <= 0.05),]
sigs_obese <- gse_obese[which(gse_obese$pvalue <= 0.05),]
sigs_all <- union(sigs_lean[["gene_symbol"]], sigs_obese[["gene_symbol"]])
overlap <- Reduce(intersect, list(v1=union(ms_cad,t2d_cad), v2=sigs_all))
overlap


degs_lean <- gse_lean[which(abs(gse_lean$log2FoldChange) >= 1),]
degs_obese <- gse_obese[which(abs(gse_obese$log2FoldChange) >= 1),]
degs_all <- union(degs_lean[["gene_symbol"]], degs_obese[["gene_symbol"]])
overlap <- Reduce(intersect, list(v1=union(ms_cad,t2d_cad), v2=degs_all))
overlap


sigDegs_lean <- gse_lean[which(abs(gse_lean$log2FoldChange) >= 1 & gse_lean$pvalue <= 0.05),]
sigDegs_obese <- gse_obese[which(abs(gse_obese$log2FoldChange) >= 1 & gse_obese$pvalue <= 0.05),]
sigDegs_all <- union(sigDegs_lean[["gene_symbol"]], sigDegs_obese[["gene_symbol"]])
overlap <- Reduce(intersect, list(v1=union(ms_cad,t2d_cad), v2=sigDegs_all))
overlap





######################################## GSE109048 & GPL17586_V2 ###############################
gse_cad <- as.data.frame(fread("VALS/GSE109048.csv", header = TRUE, sep = ','))[,1:41]
gpl_cad <- as.data.frame(fread("VALS/GPL17586_V2.txt", header = TRUE, sep = '\t'))[,c(1,9)]
colnames(gpl_cad)[1] <- "ID"
colnames(gpl_cad)[2] <- "Symbol"
row.names(gse_cad) <- gse_cad[, 1]
gse_cad[,1] <- NULL
controls <- c(1:20)
cases <- c(21:40)

overlap <- getOverlap(gse_cad, gpl_cad, "ID", "\\|")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigs <- getSigDf(gse_cad, 0.05, controls, cases)
overlap <- getOverlap(sigs, gpl_cad, "ID", "\\|")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

degs <- getDegDf(gse_cad, 1, controls, cases)
overlap <- getOverlap(degs, gpl_cad, "ID", "\\|")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigDegs <- getDegDf(sigs, 1, controls, cases)
overlap <- getOverlap(sigDegs, gpl_cad, "ID", "\\|")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

original_markers <- c("S100A12","FKBP5","CLEC4E","SAMSN1","S100P")
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=original_markers))

################################### GSE12288 & GPL96 #############################

gse_cad <- as.data.frame(fread("VALS/GSE12288.csv", header = TRUE, sep = ','))
row.names(gse_cad) <- gse_cad[, 1]
gse_cad[,1] <- NULL
gpl_cad <- as.data.frame(fread("VALS/GPL96.csv", header = TRUE, sep = ','))
groups_cad <- as.data.frame(fread("VALS/GSE12288_GROUPS.csv", header = TRUE, sep = ';'))
controls <- groups_cad[groups_cad$group == "control",1]
cases <- groups_cad[groups_cad$group == "case",1]

overlap <- getOverlap(gse_cad, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigs <- getSigDf(gse_cad, 0.05, controls, cases)
overlap <- getOverlap(sigs, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

degs <- getDegDf(gse_cad, 1, controls, cases)
overlap <- getOverlap(degs, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigDegs <- getDegDf(sigs, 1, controls, cases)
overlap <- getOverlap(sigDegs, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)


original_markers <- scan("VALS/GSE12288_markers.txt", character(), quote = "")
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=original_markers))


################################### GSE42148 & GPL13607 #############################

gse_cad <- as.data.frame(fread("VALS/GSE42148.csv", header = TRUE, sep = ','))
row.names(gse_cad) <- gse_cad[, 1]
gse_cad[,1] <- NULL
gpl_cad <- as.data.frame(fread("VALS/GPL13607.txt", header = TRUE, sep = '\t'))[,c(1,6)]
controls <- c(1:11)
cases <- c(12:24)

overlap <- getOverlap(gse_cad, gpl_cad, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigs <- getSigDf(gse_cad, 0.05, controls, cases)
overlap <- getOverlap(sigs, gpl_cad, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

degs <- getDegDf(gse_cad, 1, controls, cases)
overlap <- getOverlap(degs, gpl_cad, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigDegs <- getDegDf(sigs, 1, controls, cases)
overlap <- getOverlap(sigDegs, gpl_cad, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

original_markers <- scan("VALS/GSE42148_markers.txt", character(), quote = "")
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=original_markers))

################################### GSE40234 & GPL6480 #############################

gse_t2d <- as.data.frame(fread("VALS/GSE40234.csv", header = TRUE, sep = ','))
row.names(gse_t2d) <- gse_t2d[, 1]
gse_t2d[,1] <- NULL
gpl_t2d <- as.data.frame(fread("VALS/GPL6480.txt", header = TRUE, sep = '\t'))[,c(1,7)]
groups_t2d <- as.data.frame(fread("VALS/GSE40234_GROUPS.csv", header = TRUE, sep = ';'))
controls <- groups_t2d[groups_t2d$Group == "Sensitive",1]
cases <- groups_t2d[groups_t2d$Group == "Resistant",1]

overlap <- getOverlap(gse_t2d, gpl_t2d, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigs <- getSigDf(gse_t2d, 0.05, controls, cases)
overlap <- getOverlap(sigs, gpl_t2d, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

degs <- getDegDf(gse_t2d, 1, controls, cases)
overlap <- getOverlap(degs, gpl_t2d, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigDegs <- getDegDf(sigs, 1, controls, cases)
overlap <- getOverlap(sigDegs, gpl_t2d, "ID", NULL)
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

original_markers <- scan("VALS/GSE40234_markers.txt", character(), quote = "")
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=original_markers))



################################### GSE3585 & GPL96 #############################

gse_cad <- as.data.frame(fread("VALS/GSE3585.txt", header = TRUE, sep = '\t'))
row.names(gse_cad) <- gse_cad[, 1]
gse_cad[,1] <- NULL
gpl_cad <- as.data.frame(fread("VALS/GPL96.csv", header = TRUE, sep = ','))
groups_cad <- as.data.frame(fread("VALS/GSE12288_GROUPS.csv", header = TRUE, sep = ';'))
controls <- c(1:5)
cases <- c(6:12)

overlap <- getOverlap(gse_cad, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigs <- getSigDf(gse_cad, 0.05, controls, cases)
overlap <- getOverlap(sigs, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

degs <- getDegDf(gse_cad, 1, controls, cases)
overlap <- getOverlap(degs, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)

sigDegs <- getDegDf(sigs, 1, controls, cases)
overlap <- getOverlap(sigDegs, gpl_cad, "ID", " \\/// |-")
overlap
# setdiff(union(ms_cad,t2d_cad), overlap)




