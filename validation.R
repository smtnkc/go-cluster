if (!require('readr')) install.packages('readr')
if (!require('rstudioapi')) install.packages('rstudioapi')
if (!require('data.table')) install.packages('data.table')

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)


rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

t2d_cad <- list("ALAS2","HBD","SP1")
ms_cad  <- list("ANAPC2","ARPC1B","COX5A","ENTPD2","FRG1",
                "GLUD2","GSPT2","HCFC1","ISYNA1","LHX2",
                "NFKBIB","PCGF6","POLR2L","RPS9","S100A8",
                "SIX3","TNFRSF13B","TNFSF13","VPS28")

getValSetDEGs <- function(valGSE, valGPL, valGSMS, P_VAL, FC, N, normalization) {
  gset <- getGEO(valGSE, GSEMatrix =TRUE, AnnotGPL=TRUE)
  if (length(gset) > 1) idx <- grep(valGPL, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # group names for all samples
  gsms <- valGSMS
  sml <- c()
  for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
  
  if(normalization) exprs(gset) <- log2(exprs(gset))
  
  # set up the data and proceed with analysis
  sml <- paste("G", sml, sep="")    # set group names
  fl <- as.factor(sml)
  gset$description <- fl
  design <- model.matrix(~ description + 0, gset)
  colnames(design) <- levels(fl)
  fit <- lmFit(gset, design)
  cont.matrix <- makeContrasts(G1-G0, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)

  tT <- topTable(fit2, adjust="none", sort.by="B", p.value=P_VAL, lfc=FC, number=N)
  
  degs_parsed <- c()
  if(nrow(tT) > 0) {
    tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
    for(e in tT$Gene.symbol) {
      for(x in strsplit(e, "\\///|-")) {
        degs_parsed <- c(degs_parsed, x)
      }
    }
  }
  return(degs_parsed)
  # return(tT)
}


###### CAD

valSetDEGs <- getValSetDEGs("GSE12288", "GPL96",
                            paste0("00000000000000000000000000000000001101111111110111",
                                   "11101011011101111010010111101111001011010111000010",
                                   "11111111111011010001101111111111011001110110110001",
                                   "01110001110100000010110001110101000010111001000000",
                                   "0000000101101011111000"), 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))


valSetDEGs <- getValSetDEGs("GSE3585", "GPL96", "000001111111", 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

# valSetDEGs <- getValSetDEGs("GSE109048", "GPL17586", "?", 0.05, 1, Inf, TRUE)
# Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

# valSetDEGs <- getValSetDEGs("GSE42148", "GPL13607", "?", 0.05, 1, Inf, TRUE)
# Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE18608", "GPL570", "1X1X1X1X1X1X1X1X1X1X0000", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE9820", "GPL6255",
                            paste0("1101010010011100000111101011100",
                                   "1100111100100111010010110100011",
                                   "1010011000110111010011001110110",
                                   "1011111001100110100111110101101",
                                   "01111010111110010110100110100"), 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE4172", "GPL570", "011111111000", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

icdc <- c("11111111011111101111111110010","1X1XXX1101XXX110111X1X1XX00X0","X1X111XX0X111XX0XXX1X1X110010")
valSetDEGs <- getValSetDEGs("GSE42955", "GPL6244", icdc[1], 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE48060", "GPL570", paste0(strrep("1",30), strrep("0",20), "10"),
                            0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE3586", "GPL3050", paste0(strrep("0",15), strrep("1",13)), 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE1145", "GPL8300", "1111111XXXXXX0000", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE9128", "GPL96", "00011111111", 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

###### T2D

valSetDEGs <- getValSetDEGs("GSE40234", "GPL6480",
                            "11111110010101001001101001000110111010101000110000110001111111", 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE16415", "GPL2986", "0000011111", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE26168", "GPL6883", "00000000XXXXXXX111111111", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE23343", "GPL570", "00000001111111111", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE25462", "GPL570", "00000000000000000000000000000000000011111111110000", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE15653", "GPL96", "00000000011111XXXX", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE20966", "GPL1352", "00000000001111111111", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE12643", "GPL8300", "11111111110000000000", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE38642", "GPL6244", "000000000011000000000000000000000000001100000111000000000000110",
                            0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE25724", "GPL96", "0000000111111", 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE13760", "GPL571", "101100010000101001111", 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

valSetDEGs <- getValSetDEGs("GSE9006", "GPL96",
                            paste0(strrep("0",20), strrep("X",81), strrep("0",4),
                                   strrep("1",12)), 0.05, 1, Inf, TRUE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

# valSetDEGs <- getValSetDEGs("GSE133099", "GPL16791", "?", 0.05, 1, Inf, TRUE)
# Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))

getValSetDEGsUsingMultipleGPLs <- function(valGSE, gplVector) {
  degs <- c()
  for(tempGPL in gplVector) {
    gset <- getGEO(valGSE, GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep(tempGPL, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    # make proper column names to match toptable 
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    
    # group names for all samples
    gsms <- "1111100000"
    sml <- c()
    for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
    
    # log2 transform
    # exprs(gset) <- log2(exprs(gset))
    
    # set up the data and proceed with analysis
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="none", sort.by="B", p.value=0.05, lfc=1, number=Inf)
    tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
    degs <- c(degs, tT$Gene.symbol)
  }
  degs <- unique(degs)
  degs_parsed <- c()
  for(e in degs) {
    for(x in strsplit(e, "\\///|-")) {
      degs_parsed <- c(degs_parsed, x)
    }
  }
  degs_parsed <- unique(degs_parsed)
  return(degs_parsed)
}

valSetDEGs <- getValSetDEGsUsingMultipleGPLs("GSE121",
                                             c("GPL80","GPL98","GPL99","GPL100","GPL101"))
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))


######### METS

valSetDEGs <- getValSetDEGs("GSE43760", "GPL6244", "1X1X1X1X1X1X0X0X0X0X0X0X", 0.05, 1, Inf, FALSE)
Reduce(intersect, list(v1=union(ms_cad, t2d_cad), v2=valSetDEGs))


##########

valDEGsMETS <- c("ALAS2", "S100A8")
valDEGsCAD <- c("SIX3","TNFSF13","S100A8","NFKBIB","SP1","FRG1","ALAS2")
valDEGsT2D <- c("POLR2L","ARPC1B","SIX3","NFKBIB","ALAS2","S100A8","SP1","HBD","TNFSF13","LHX2")
valDEGsALL <- union(union(valDEGsCAD, valDEGsT2D), valDEGsMETS)



