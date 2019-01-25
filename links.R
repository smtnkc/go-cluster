source("vars.R")

##### Create mapping list from mapping matrice

getMapLists <- function(topologies, mapTypes) {
  maps <- list()
  for(t in topologies) {
    fname <- paste("map_", t, "_unique.csv", sep="")
    maps[[t]] <- as.data.frame(fread(fname, header = TRUE, sep = ','))
    if(t == "inet") {
      maps[[t]] <- transform(maps[[t]], ENTREZID = as.character(ENTREZID))
    }
  }
  mapLists <- list()
  for(m in mapTypes) {
    
    mapLists[[m]] <- list()
    
    if(substr(m,1,1) == "p") t <- "string" else t <- "inet"
    if(substr(m,3,3) == "s") val <- 2 else val <- 3
    
    tmp <- maps[[t]]
    for(i in 1:nrow(tmp)) {
      mapLists[[m]][tmp[i, 1]] <- tmp[i, val]
    }
  }
  mapLists <- lapply(mapLists, function(x) x[!is.na(x)]) # remove NA values
  return(mapLists)
}

mapLists <- getMapLists(topologies, mapTypes)

##### Note that elinks_inet and plinks_string are downloaded files.
elinks_inet <- as.data.frame(fread("LINKS/elinks_inet.csv", header = TRUE, sep = ','))
plinks_string <- as.data.frame(fread("LINKS/plinks_string.csv", header = TRUE, sep = ' '))

createLinks <- function(df.links, mapLists, mapType) {
  portion <- 100000
  i <- 1
  n <- nrow(df.links)
  df.new_links <- df.links

  if(substr(mapType, 1, 1) == "p") {
    old_cols <- c("protein1", "protein2")
  } else {
    old_cols <- c("entrez1", "entrez2")
  }
  
  if(substr(mapType, 3, 3) == "s") {
    new_cols <- c("symbol1", "symbol2")
  } else {
    new_cols <- c("entrez1", "entrez2")
  }
  
  df.new_links[, new_cols] <- NA

  while(i <= n) {
    j <- i + portion - 1
    if(j > n) {
      j <- n
    }
    cat("range", i, "to", j, "...\n")

    df.new_links[i:j,new_cols[1]] <- 
      as.vector(as.character(mapLists[[mapType]][df.new_links[i:j,old_cols[1]]]))
    df.new_links[i:j,new_cols[2]] <- 
      as.vector(as.character(mapLists[[mapType]][df.new_links[i:j,old_cols[2]]]))

    i <- i+portion
  }

  df.new_links <- df.new_links[, c(new_cols, "combined_score")]
  n_old <- nrow(df.new_links)
  df.new_links <- df.new_links[df.new_links[, new_cols[1]] != "NULL", ] # remove NULL rows
  df.new_links <- df.new_links[df.new_links[, new_cols[2]] != "NULL", ] # remove NULL rows
  n_new <- nrow(df.new_links)
  
  cat("*********", n_old-n_new, "NULL rows are removed!\n")
  return(df.new_links)
}

getLinks <- function(df.links, mapLists, mapType) {
  if(substr(mapType, 1, 1) == "p") {
    # string mapping
    fname = paste("LINKS/", substr(mapType, 3, 3), "links_string.csv", sep="")
    
    if(file.exists(fname)) {
      print(paste("READING FROM -> ", fname, sep=""))
      return(as.data.frame(fread(fname, header = TRUE, sep = ',')))
    }
    print(paste("CREATING -> ", fname, sep=""))
    df.temp_links <- createLinks(df.links, mapLists, mapType)
    write.csv(df.temp_links, fname, row.names = FALSE)
    return(df.temp_links)
    
  } else {
    # inet mapping
    fname = "LINKS/slinks_inet.csv"
    if(file.exists(fname)) {
      print(paste("READING FROM -> ", fname, sep=""))
      return(as.data.frame(fread(fname, header = TRUE, sep = ',')))
    }
    print(paste("CREATING -> ", fname, sep=""))
    df.temp_links <- createLinks(df.links, mapLists, mapType)
    write.csv(df.temp_links, fname, row.names = FALSE)
    return(df.temp_links)
  }
}

elinks_string <- getLinks(plinks_string, mapLists, "p2e")
slinks_string <- getLinks(plinks_string, mapLists, "p2s")
slinks_inet   <- getLinks(elinks_inet, mapLists, "e2s")
