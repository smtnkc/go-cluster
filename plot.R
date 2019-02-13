source("vars.R")

gosim <- readGosim(topologies, subjects, measures, ontTypes, includeComb = FALSE)

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
      for(o in c("BP")) { #for(o in names(gosim[[t]][[s]][[1]])) {
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
            fname <- paste("PLOTS/DENSITY/", t, "_", s, "_All_", o,  ".png", sep="")
            cat(fname, "...\n")
            png(filename=fname, width = 960, height = 640)
            xLabel <- "Weight"
            plot(d, main=title, xlab=xLabel, col=lineColors[i], ylim=c(0,8),
                 cex.lab=1.5, cex.main=1.5, cex.axis=1.5, lwd=2)
          } else {
            lines(d, col=lineColors[i], lwd=2)
          }
          legend("topright", legend = names(rev(gosim[[t]][[s]])), col = lineColors, lwd=2, cex=1.2)
          i <- i+1
        }
        dev.off()
      }
    }
  }
}

generateDensityPlotsTogether(gosim)
