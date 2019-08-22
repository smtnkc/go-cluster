if (!require('readr')) install.packages('readr')
if (!require('rstudioapi')) install.packages('rstudioapi')
if (!require('data.table')) install.packages('data.table')
if (!require('igraph')) install.packages('igraph')
if (!require('dplyr')) install.packages('dplyr')
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('gridExtra')) install.packages('gridExtra')
if (!require('BiocManager')) install.packages('BiocManager')
if (!require('org.Hs.eg.db')) BiocManager::install('org.Hs.eg.db')
if (!require('GOSemSim')) BiocManager::install('GOSemSim')
if (!require('hgu133a.db')) BiocManager::install('hgu133a.db')
if (!require('annotate')) BiocManager::install('annotate')
if (!require('qgraph')) install.packages('qgraph')
if (!require("MCL")) install.packages("MCL")
if (!require("linkcomm")) install.packages("linkcomm")
if (!require("profmem")) install.packages("profmem")

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

raw <- as.data.frame(fread("RAW.csv", header = TRUE, sep = ','))

ontTypes <- c("MF", "BP", "CC")
measures <- c("Resnik", "Lin", "Rel", "Jiang", "Wang")
topologies <- c("string", "inet")
subjects <- c("mets", "t2d", "cad")

intervals_ctrl <- list( # control included intervals
  "mets" = c(1:10, 17:22),
  "t2d"  = c(1:10, 29:36),
  "cad"  = c(1:10, 23:28))

intervals <- list( # control excluded intervals
  "mets" = c(1, 17:22),
  "t2d"  = c(1, 29:36),
  "cad"  = c(1, 23:28))

cutoffs <- list(
  "10" = list(
    "string" = 407,
    "inet" = 0.175
  ))

FC <- 1
TOP_N = "10"
P_VAL = "0.05"

Stringify <- function(n) {
  # converts 0.05 to "005"
  s <- paste(substr(n, 1, regexpr("\\.", n)[1]-1), 
             substr(n, regexpr("\\.", n)[1]+1, nchar(n)), sep="")
  return(s)
}

FormatDf <- function(df) {
  df[is.na(df)] <- ""
  df[df == "-"] <- ""
  df[df == "."] <- ""
  df[df == "Null"] <- ""
  df[df == "--Null"] <- ""
  df[df == "EMPTY"] <- ""
  df[df == "--empty"] <- ""
  df[df == "--unknown"] <- ""
  df[df == "Dye Marker"] <- ""
  return (df)
}

readGosim <- function(topologies, subjects, measures, ontTypes, includeComb, naming) {
  gosim <- list()
  if(includeComb) measures <- c(measures, "Comb")
  
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/GOSIM/by", naming, "/", t , "_", s, "_", m, "_", o, ".tsv", sep="")
          gosim[[t]][[s]][[m]][[o]] <- as.data.frame(fread(fname, header = TRUE, sep = '\t'))
        }
      }
    }
  }
  return(gosim)
}

globalDarkColorPalette <- c("#528ECC", "#E8D08B", "#F19C99", "#6A6AB0", "#61AD85")
globalLightColorPalette<- c("#A6CFFF", "#CCFFCC", "#FFCCCC", "#FFE599", "#CDA2BE")