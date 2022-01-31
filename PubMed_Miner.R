# Functions for Web-Scaping PubMed (and other NCBI databases)
# using the Entrez Programming Utilities API (https://www.ncbi.nlm.nih.gov/books/NBK25500/)

# init packages
library(stringr)    # for strings and regular expressions
library(xml2)       # for parsing data in XML (e.g. HTML)
library(rvest)      # for scraping XML and HTML content
library(tibble)     # for storing data frames
library(dplyr)       # for data frame processing

#functions
esearch <- function(search_term, mindate=2000, maxdate=2019, retmax=10000, db="pubmed"){
  entrez_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  
  params_esearch <- paste0(
    c(paste0("db=",db), #searched data base
      paste0("term=",search_term), #search term
      paste0("mindate=",mindate), #ealiest year of entries
      paste0("maxdate=",maxdate), #latest date of entries
      paste0("retmax=",retmax) #maximal number of returned entries
      ),
    collapse = "&")
  
  query <- paste0(entrez_url, "esearch.fcgi?", params_esearch)
  ret <- read_xml(query)
  
  nrecords <- length(html_text(html_nodes(ret, xpath = "//Id")))
  message(paste(nrecords,"records (PMIDs) extracted from",db))
  if(nrecords == retmax) warning("maximum record number reached. increase retmax argument")
  message(paste("search term:", search_term))
  message(paste("search filters: mindate=", mindate, ", maxdate=",maxdate))
  return(ret)
}


esummary <- function(esearch_entries, db="pubmed", maxquery=250){
  entrez_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  
  ids <- html_text(html_nodes(esearch_entries, xpath = "//Id")) # extract PMIDs
  ret <- tibble("PMID"=ids,
                last_author=as.character(NA),
                year=as.character(NA),
                title=as.character(NA),
                journal=as.character(NA),
                volume=as.character(NA),
                pages=as.character(NA),
                doi=as.character(NA))
  
  if(length(ids) <= maxquery){
    
    params_esummary <- paste0(
      paste0("db=",db,"&id="), 
      paste(ids, collapse = ",")
    )
    
    query <- paste0(entrez_url, "esummary.fcgi?", params_esummary)
    summaries <- read_xml(query)
    
    temp <- xml_children(summaries)
    pb <- utils::txtProgressBar(min = 0, max = length(temp), style = 3)
    for(i in 1:length(temp)){
      id <- which(ids == xml_text(xml_node(temp[i], xpath = 'Id')))
      ret[id, "last_author"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="LastAuthor"]'))
      ret[id, "year"] <- substr(xml_text(xml_node(temp[i], xpath = 'Item[@Name="PubDate"]')), start = 1, stop = 4)
      ret[id, "title"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Title"]'))
      ret[id, "journal"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Source"]'))
      ret[id, "volume"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Volume"]'))
      ret[id, "pages"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Pages"]'))
      ret[id, "doi"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="DOI"]'))
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
  
  }else{
    message("number of records > maxquery argument. query for article summaries is split. please be patient")
    
    ids_batch <- list()
    if(length(ids) %% maxquery == 0){
      nbatch <- length(ids) %/% maxquery
    }else{
      nbatch <- length(ids) %/% maxquery + 1
    }
    for(j in 1:nbatch) ids_batch[[j]] <- ids[((j-1)*maxquery+1):(j*maxquery)]
    
    pb <- utils::txtProgressBar(min = 0, max = nbatch, style = 3)
    for(j in 1:nbatch){
      params_esummary <- paste0(
        paste0("db=",db,"&id="), 
        paste(ids_batch[[j]], collapse = ",")
      )
      
      query <- paste0(entrez_url, "esummary.fcgi?", params_esummary)
      summaries <- read_xml(query)
      
      temp <- xml_children(summaries)
      
      for(i in 1:length(temp)){
        id <- which(ids == xml_text(xml_node(temp[i], xpath = 'Id')))
        ret[id, "last_author"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="LastAuthor"]'))
        ret[id, "year"] <- substr(xml_text(xml_node(temp[i], xpath = 'Item[@Name="PubDate"]')), start = 1, stop = 4)
        ret[id, "title"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Title"]'))
        ret[id, "journal"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Source"]'))
        ret[id, "volume"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Volume"]'))
        ret[id, "pages"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="Pages"]'))
        ret[id, "doi"] <- xml_text(xml_node(temp[i], xpath = 'Item[@Name="DOI"]'))

      }
      setTxtProgressBar(pb, j)
    }
    close(pb)
   
  }
  
  return(ret)
}


efetch <- function(summaries, db="pubmed", maxquery = 250){
  entrez_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  
  if(!("abstract" %in% names(summaries)) ) summaries$abstract <- as.character(NA)
  
  ids <- summaries$PMID[is.na(summaries$abstract)]
  
  if(length(ids) <= maxquery){
    
    params_efetch <- paste0(
      paste0("db=",db),
      "&id=", paste(ids, collapse = ","),
      "&retmode=abstract",
      "&rettype=txt"
    )
    
    query <- paste0(entrez_url, "efetch.fcgi?", params_efetch)
    txt_abstracts <- c(" "," ",readLines(query))
    eol_abstracts <- c(1,grep("PMID",substr(txt_abstracts, 1, 4)))
    
    pb <- utils::txtProgressBar(min = 1, max = length(eol_abstracts), style = 3)
    for(i in 2:length(eol_abstracts)){
      bol_abstract <- eol_abstracts[i-1]
      bol_abstract <- bol_abstract + grep(pattern = "[[:digit:]]", substr(txt_abstracts[(bol_abstract+1):(bol_abstract+25)],1,1))[1]
      if(is.na(bol_abstract)) bol_abstract <- eol_abstracts[i-1] + 1
      
      txt_abstract <- paste(txt_abstracts[bol_abstract:eol_abstracts[i]], collapse = " ")
      
      PMID <- unlist(strsplit(txt_abstracts[eol_abstracts[i]], " "))[2]
      if(PMID %in% summaries$PMID){
        summaries[which(summaries$PMID == PMID), "abstract"] <- txt_abstract
      }else{
        warning(paste0("The following abstract PMID ",PMID," could not be matched against PMIDs of summaries object:"))
        warning(txt_abstract)
      }
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
  }else{
    message("number of records > maxquery argument. query for article abstracts is split. please be patient")
    
    ids_batch <- list()
    if(length(ids) %% maxquery == 0){
      nbatch <- length(ids) %/% maxquery
    }else{
      nbatch <- length(ids) %/% maxquery + 1
    }
    for(j in 1:nbatch) ids_batch[[j]] <- ids[((j-1)*maxquery+1):(j*maxquery)]
    
    pb <- utils::txtProgressBar(min = 0, max = nbatch, style = 3)
    for(j in 1:nbatch){
      params_efetch <- paste0(
        paste0("db=",db),
        "&id=", paste(ids_batch[[j]], collapse = ","),
        "&retmode=abstract",
        "&rettype=txt"
      )
      
      query <- paste0(entrez_url, "efetch.fcgi?", params_efetch)
      txt_abstracts <- c(" "," ",readLines(query))
      eol_abstracts <- c(1,grep("PMID",substr(txt_abstracts, 1, 4)))
      
      for(i in 2:length(eol_abstracts)){
        bol_abstract <- eol_abstracts[i-1]
        bol_abstract <- bol_abstract + grep(pattern = "[[:digit:]]", substr(txt_abstracts[(bol_abstract+1):(bol_abstract+25)],1,1))[1]
        if(is.na(bol_abstract)) bol_abstract <- eol_abstracts[i-1] + 1
        
        txt_abstract <- paste(txt_abstracts[bol_abstract:eol_abstracts[i]], collapse = " ")
        
        PMID <- unlist(strsplit(txt_abstracts[eol_abstracts[i]], " "))[2]
        if(PMID %in% summaries$PMID){
          summaries[which(summaries$PMID == PMID), "abstract"] <- txt_abstract
        }else{
          warning(paste0("The following abstract PMID ",PMID," could not be matched against PMIDs of summaries object:"))
          warning(txt_abstract)
        }
        
      }
      setTxtProgressBar(pb, j)
    }
    close(pb)
  }
  
  return(summaries)
}

efilter <- function(summaries, filter_term, col="abstract", exclude = TRUE){
  
  if(col %in% names(summaries)){
    chartemp <- summaries %>% pull(col) %>% tolower
    efilter <- grepl(filter_term, chartemp)
    
    if(exclude){
      return(summaries[!efilter,])
    }else{
      return(summaries[efilter,])
    }
    
  }else{
    stop(paste0("filtering cannot be performed due to invalid col argument: '", col,"' is not available in summaries object"))
  }
}


#working example

esearch_entries <- esearch(search_term = "alcohol+stress", retmax = 1000) # search database using esearch 
summaries <- esummary(esearch_entries, maxquery = 100) # summarize esearch results using esummary 
summaries <- efetch(summaries, maxquery = 100) # append abstracts to summaries of esearch entries using efetch
summaries <- efilter(summaries, filter_term = "rodent|rat|mouse|mice", col = "title", exclude = TRUE) #remove entries with keywords in title column of summaries
  
View(summaries)
write.csv2(summaries, file = "pubmed_summaries.csv")
