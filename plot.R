source("vars.R")

readGosim <- function(topologies, subjects, measures, ontTypes, includeComb) {
  gosim <- list()
  if(includeComb) measures <- c(measures, "Comb")
  
  for(t in topologies) {
    for(s in subjects) {
      for(m in measures) {
        for(o in ontTypes) {
          fname <- paste("RES/GOSIM/", m , "/", t, "_", s, "_", o, ".tsv", sep="")
          if(o == "MF") {
            gosim[[t]][[s]][[m]] <- as.data.frame(fread(fname, header = TRUE, sep = '\t'))
          } else {
            gosim[[t]][[s]][[m]][,o] <- as.data.frame(fread(fname, header = TRUE, sep = '\t'))[,3]
          }
        }
      }
    }
  }
  return(gosim)
}

gosim <- readGosim(topologies, subjects, measures, ontTypes, includeComb = FALSE)

generateDensityPlotsSeperated <- function(gosim) {
  lineColors <- c("orange", "red", "green", "purple", "blue")
  for(t in names(gosim)) {
    for(s in c("cad")) {
      i <- 1
      for(m in rev(names(gosim[[t]][[s]]))) {
        for(o in c("BP")) {
          df <- gosim[[t]][[s]][[m]]
          w <- df[!is.na(df[,o]), o]
          d <- density(w)
          title <- paste(t,"_",s,"_",m,"_",o,sep="")
          png(filename=paste("PLOTS/SEP/", title, ".png", sep=""),
              width = 640, height = 640)
          xLabel <- paste("Weight (N = ", length(w), ")", sep="")
          plot(d, main=title, xlab=xLabel, col=lineColors[i], ylim=c(0,4),
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
    for(s in c("cad")) {
      i <- 1
      for(m in rev(names(gosim[[t]][[s]]))) {
        for(o in c("BP")) {
          df <- gosim[[t]][[s]][[m]]
          w <- df[!is.na(df[,o]), o]
          d <- density(w)
          if(i == 1) {
            title <- paste(t,"_", s, sep="")
            png(filename=paste("PLOTS/", title, ".png", sep=""),
                width = 960, height = 640)
            xLabel <- "Weight"
            plot(d, main=title, xlab=xLabel, col=lineColors[i], ylim=c(0,4),
                 cex.lab=1.5, cex.main=1.5, cex.axis=1.5, lwd=2)
          } else {
            lines(d, col=lineColors[i], lwd=2)
          }
        }
        legend("topright", legend = names(rev(gosim[[t]][[s]])), col = lineColors, lwd=2, cex=1.2)
        i <- i+1
      }
      dev.off()
    }
  }
}

generateDensityPlotsTogether(gosim)
