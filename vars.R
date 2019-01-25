library(readr)
library(rstudioapi)
library(data.table) # fread
library(org.Hs.eg.db) # to get entrez ids
library(GOSemSim)

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

mapTypes <- c("p2s", "p2e", "e2s")
idTypes <- c("s", "e")
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
    "s" = list(
      "string" = 407,
      "inet" = 0.168),
    "e" = list(
      "string" = 358,
      "inet" = 0.176)
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
