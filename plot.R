source("vars.R")

gosim <- readGosim(topologies, subjects, measures, ontTypes, includeComb = FALSE, naming = "SYMBOL")

####################################################################

generateBoxPlots <- function(gosim, subject) {
  boxCols <- c("#aaa69d", "#ccae62", "#40407a", "#cc8e35", "#218c74")
  onts <- list()
  for(t in topologies) {
    for(o in ontTypes) {
      df <- data.frame(measure = character(), score = double())
      for(m in measures) {
        temp <- data.frame(measure=m, score=gosim[[t]][[subject]][[m]][[o]][[o]])
        df <- rbind(df, temp)
      }
      onts[[t]][[o]] <- df
    }
  }
  svg(filename = paste("PLOTS/box_", subject,".svg", sep=""), width=10, height=8)
  par(mfrow=c(2,3), mar=c(4,4,2.3,0.0) + 0.5)
  for(t in topologies) {
    i <- 1
    for(o in ontTypes) {
      boxplot(score~measure,data=onts[[t]][[o]], main="", xlab="", ylab="", col=boxCols, cex.lab=1, cex.axis=1.3)
      mtext(side = 3, text = o, line = 0.8, cex=1, font=2)
      mtext(side = 1, text = "GO Similarity Measure", line = 3, cex=0.8)
      mtext(side = 2, text = "GO Similarity Score", line = 2.6, cex=0.8)
      i <- i+1
    }
  }
  dev.off()
}

generateBoxPlots(gosim, "mets")

####################################################################

generateDensityPlotsSeperated <- function(gosim) {
  lineColors <- c("orange", "red", "green", "purple", "blue")
  for(t in names(gosim)) {
    for(s in c("cad")) { #for(s in names(gosim[[t]])) {
      i <- 1
      for(m in rev(names(gosim[[t]][[s]]))) {
        for(o in c("BP")) { #for(o in names(gosim[[t]][[s]][[m]])) {
          df <- gosim[[t]][[s]][[m]][[o]]
          w <- df[!is.na(df[, 3]), 3]
          d <- density(w)
          title <- paste("Topology: ", toupper(t),
                         "   Subject: ", toupper(s),
                         "   Measure: ", toupper(m),
                         "   Ontology: ", toupper(o), sep="")
          fname <- paste("PLOTS/DENSITY/", t, "_", s, "_", m, "_", o, ".png", sep="")
          cat(fname, "...\n")
          png(filename=fname,
              width = 640, height = 640)
          xLabel <- paste("Weight (N = ", length(w), ")", sep="")
          plot(d, main=title, xlab=xLabel, col=lineColors[i], ylim=c(0,8),
               cex.lab=1.5, cex.main=1.5, cex.axis=1.5, lwd=2)
          dev.off()
        }
        i <- i + 1
      }
    }
  }
}

generateDensityPlotsSeperated(gosim)

generateDensityPlotsTogether <- function(gosim) {
  lineColors <- c("orange", "red", "green", "purple", "blue")
  for(t in names(gosim)) {
    for(s in c("cad")) { #for(s in names(gosim[[t]])) {
      fname <- paste("PLOTS/DENSITY/", t, "_", s, "_All.svg", sep="")
      cat(fname, "...\n")
      svg(filename = fname, width=9, height=3)
      par(mfrow=c(1,3), mar=c(2.8,2.8,1.0,0.2) + 0.2)
      for(o in c("BP", "CC", "MF")) { #for(o in names(gosim[[t]][[s]][[1]])) {
        i <- 1
        for(m in rev(names(gosim[[t]][[s]]))) {
          df <- gosim[[t]][[s]][[m]][[o]]
          w <- df[!is.na(df[, 3]), 3]
          d <- density(w)
          if(i == 1) {
            title <- paste("Topology: ", toupper(t),
                           "   Subject: ", toupper(s),
                           "   Measure: ", "All",
                           "   Ontology: ", toupper(o), sep="")
            plot(d, main=o, xlab="", ylab="", zero.line = FALSE,
                 lwd=2, col=lineColors[i], ylim=c(0,8), cex.lab=1, cex.axis=1)
            mtext(side = 1, text = "GO Semantic Similarity Score", line = 2, cex=0.7)
            mtext(side = 2, text = "Density", line = 2, cex=0.7)

          } else {
            lines(d, col=lineColors[i], lwd=2)
          }
          legend("topright", legend = names(rev(gosim[[t]][[s]])), col = lineColors, lty=1, lwd=2, cex=1)
          i <- i+1
        }
      }
      dev.off()
    }
  }
}

generateDensityPlotsTogether(gosim)

####################################################################

ctrl <- data.frame(subject="ctrl",exp=as.data.frame(apply(subsets[["mets"]][,1:9],1,median))[,1])
ms <- data.frame(subject="ms",exp=as.data.frame(apply(subsets[["mets"]][,10:15],1,median))[,1])
t2d  <- data.frame(subject="t2d",exp=as.data.frame(apply(subsets[["t2d"]][,10:17],1,median))[,1])
cad  <- data.frame(subject="cad",exp=as.data.frame(apply(subsets[["cad"]][,10:15],1,median))[,1])
df <- rbind(ctrl,ms,t2d,cad)

svg(filename = "PLOTS/dist.svg", width=6, height=3)
par(mfrow=c(1,1), mar=c(2.8,2.6,0.0,0.0) + 0.2)
plot(density(apply(subsets[["cad"]][,10:15],1,median)),
     main="", xlab="", ylab="", zero.line = FALSE,
     col = "#F19C99", lwd=2, cex.lab=0.7, cex.axis=0.7) #cad
mtext(side = 1, text = "Gene Expression", line = 2, cex=0.7)
mtext(side = 2, text = "Density", line = 2, cex=0.7)
legend("topright", legend=c("CTRL", "MS", "CAD", "T2D"),
       col=c("#000099","#9999FF","#F19C99", "#66B2FF"), lty=1, lwd=2, cex=0.7)

lines(density(apply(subsets[["mets"]][,10:15],1,median)), col = "#9999FF", lwd=2) #ms
lines(density(apply(subsets[["mets"]][,1:9],1,median)), col = "#000099", lwd=2) #ctrl
lines(density(apply(subsets[["t2d"]][,10:17],1,median)), col = "#66B2FF", lwd=2) #t2d
dev.off()
