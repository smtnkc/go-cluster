source("vars.R")

##################################################

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

getMapRaw <- function() {
  fname <- "map_raw_unique.csv"
  if(file.exists(fname)) {
    cat("READING FROM", fname, "...\n")
    df.map_raw <- as.data.frame(fread(fname, header = TRUE, sep = ','))
  } else {
    cat("CREATING", fname, "...\n")
    s_raw <- as.data.frame(fread("s_RAW.csv", header = TRUE, sep = ','))
    GPL  <- as.data.frame(fread("GPL.csv", header = TRUE, sep = ','))
    GPL <- FormatDf(GPL)
    GPL <- GPL[GPL$`Symbol v12`!="",]
    
    symbols_initial <- sort(unique(GPL$`Symbol v12`)) # has duplicates because of casing
    symbols_uppered <- sort(s_raw$Symbol) # no duplicates but all uppercased
    symbols_uniqued <- sort(symbols_initial[! symbols_initial %in% symbols_initial[duplicated(toupper(symbols_initial))]])
    
    df.map_case <- data.frame(s_raw$Symbol, symbols_uniqued)
    df.map_case[,1] <- as.character(df.map_case[,1])
    df.map_case[,2] <- as.character(df.map_case[,2])
    df.map_raw <- select(org.Hs.eg.db, 
                         keys = df.map_case$symbols_uniqued,
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")
    
    #View(df.map_F635[df.map_F635$SYMBOL %in% (df.map_F635[duplicated(df.map_F635$SYMBOL),]$SYMBOL),]) # duplicated symbols of our dataset
    df.map_raw <- aggregate(df.map_raw, by = list(df.map_raw$SYMBOL), min)[,2:3]  # remove duplicated symbols by keeping min of entrez ids
    
    # if all symbols are same with F635 (when uppercased) then write to file
    if(nrow(df.map_raw[toupper(df.map_raw$SYMBOL) != df.map_raw$Symbol,]) == 0) {
      cat("WRITING", fname, "...\n")
      write.csv(df.map_raw, "map_raw_unique.csv", row.names = FALSE)  
    }
  }
  return(df.map_raw)
}

df.map_raw <- getMapRaw()

######################################################## INET

getMapInet <- function() {
  fname <- "map_inet_unique.csv"
  if(file.exists(fname)) {
    cat("READING FROM", fname, "...\n")
    df.map_inet <- as.data.frame(fread(fname, header = TRUE, sep = ','))
  } else {
    cat("CREATING", fname, "...\n")
    df.map_inet <- as.data.frame(fread("map_inet.csv", header = TRUE, sep = ','))
    #View(df.map_inet[df.map_inet$SYMBOL %in% (df.map_inet[duplicated(df.map_inet$SYMBOL),]$SYMBOL),]) # duplicated symbols of inet map
    df.map_inet <- aggregate(df.map_inet, by = list(df.map_inet$SYMBOL), min)[,2:3]  # remove duplicated symbols by keeping min of entrez ids
    df.map_inet[,1] <- as.character(df.map_inet[,1])
    cat("WRITING", fname, "...\n")
    write.csv(df.map_inet, fname, row.names = FALSE)
  }
  return(df.map_inet)
}

df.map_inet <- getMapInet()

######################################################## STRING

getMapString <- function() {
  fname <- "map_string_unique.csv"
  if(file.exists(fname)) {
    cat("READING FROM", fname, "...\n")
    df.map_string <- as.data.frame(fread("map_string_unique.csv",
                                         header = TRUE, sep = ','))
  } else {
    cat("CREATING", fname, "...\n")
    df.map_string <- as.data.frame(fread("map_string.csv",
                                         header = TRUE, sep = ','))
    df.map_tmp <- select(org.Hs.eg.db, 
                         keys = df.map_string$preferred_name,
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")
    
    #View(df.map[df.map$SYMBOL %in% (df.map[duplicated(df.map$SYMBOL),]$SYMBOL),]) # duplicated symbols of our dataset
    df.map_tmp <- aggregate(df.map_tmp, by = list(df.map_tmp$SYMBOL), min)[,2:3]  # remove duplicated symbols by keeping min of entrez ids
    
    df.map_string$ENTREZID <- df.map_tmp$ENTREZID
    colnames(df.map_string)[1] <- "PROTEIN"
    colnames(df.map_string)[2] <- "SYMBOL"
    print(nrow(df.map_string))
    #df.map_string <- na.omit(df.map_string)
    #print(nrow(df.map_string))
    cat("WRITING", fname, "...\n")
    write.csv(df.map_string, fname, row.names = FALSE)
  }
  return(df.map_string)
}

df.map_string <- getMapString()