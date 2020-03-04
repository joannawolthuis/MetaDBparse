install.packages("devtools")
install.packages("BiocManager")
install.packages("pacman")
library(pacman)

pacman::p_load(rcdk, rJava, parallel, pbapply, enviPat, data.table, 
               RSQLite, DBI, gsubfn, utils, RCurl, 
               XML, base, stringr, WikidataQueryServiceR, webchem, 
               openxlsx, jsonlite, R.utils, KEGGREST, zip, ChemmineR, 
               rvest, xml2, stringi, reshape2, Hmisc, httr, RJSONIO, 
               readxl, cmmr, progress, Rdisop, rlist)