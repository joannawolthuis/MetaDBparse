.onLoad <- function(libname, pkgname) {
  data("adduct_rules", "adducts", package=pkgname, envir=parent.env(environment()))
  data("isotopes", package="enviPat", envir=parent.env(environment()))
  }

#' @title Build WikiPathways
#' @description Parses WikiPathways chemical compound database, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[SPARQL]{SPARQL}}
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[rvest]{html_nodes}}
#'  \code{\link[xml2]{read_html}}
#' @rdname build.WIKIPATHWAYS
#' @export
#' @importFrom SPARQL SPARQL
#' @importFrom RCurl getURL
#' @importFrom RSQLite dbConnect dbWriteTable dbRemoveTable dbDisconnect
#' @importFrom data.table as.data.table
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes
#' @examples
#' \dontrun{build.WIKIPATHWAYS(outfolder=tempdir(), testMode=TRUE)}
build.WIKIPATHWAYS <- function(outfolder, testMode = FALSE) {
  . <- compoundname <- description <- NULL
  chebi.loc <- file.path(outfolder, "chebi.db")
  chebiExists <- file.exists(chebi.loc)
  if (!chebiExists) {
    print("Requires ChEBI. Building...")
    buildBaseDB(outfolder, "chebi")
  }
  chebi <- data.table::as.data.table(showAllBase(outfolder, "chebi"))
  chebi$identifier <- as.numeric(chebi$identifier)
  base.db <- SPARQL::SPARQL(url = "http://sparql.wikipathways.org/sparql", query = "prefix wp: <http://vocabularies.wikipathways.org/wp#>\n                                                   prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#>\n                                                   prefix dcterms: <http://purl.org/dc/terms/>\n                                                   prefix xsd:     <http://www.w3.org/2001/XMLSchema#>\n                                                   PREFIX wdt: <http://www.wikidata.org/prop/direct/>\n\n                                                   select  ?mb\n                                                   (group_concat(distinct str(?labelLit);separator=\", \") as ?label )\n                                                   ?idurl as ?csid\n                                                   (group_concat(distinct ?pwTitle;separator=\", \") as ?description)\n                                                   ?pathway\n                                                   where {\n                                                   ?mb a wp:Metabolite ;\n                                                   rdfs:label ?labelLit ;\n                                                   wp:bdbChEBI ?idurl ;\n                                                   dcterms:isPartOf ?pathway .\n                                                   ?pathway a wp:Pathway ;\n                                                   dc:title ?pwTitle .\n                                                   FILTER (BOUND(?idurl))\n                                                   }\n                                                   GROUP BY ?mb ?wp ?idurl ?pathway")
  chebi.join.table <- data.table::data.table(identifier = as.numeric(gsub(base.db$results$csid, pattern = ".*:|>", replacement = "")), description = base.db$results$description, widentifier = base.db$results$mb, pathway = base.db$results$pathway)
  base.db$identifier <- as.numeric(base.db$identifier)
  chebi <- chebi[, -"description", with = FALSE]
  db.formatted <- merge(chebi.join.table, chebi, by = "identifier")
  db.formatted <- data.table::as.data.table(db.formatted)
  db.formatted <- unique(db.formatted[, .(compoundname = paste0(unique(compoundname), collapse = ", "), description = paste0("Found in pathways: ", paste0(unique(description), collapse = ", "))), by = c("structure", "identifier", "charge")])
  page <- RCurl::getURL("https://www.wikipathways.org/index.php/WikiPathways")
  ver <- stringr::str_match(page, ">([\\w| ]*?) Release<")[, 2]
  list(db = db.formatted, version = ver)
}

#' @title Build MCDB
#' @description Parses the MCDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[XML]{readHTMLTable}},\code{\link[XML]{xmlValue}},\code{\link[XML]{xmlEventParse}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[pbapply]{pboptions}}
#'  \code{\link[base]{connections}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.MCDB
#' @export
#' @examples
#' \dontrun{build.MCDB(outfolder = tempdir(), testMode = TRUE)}
#' @importFrom utils download.file unzip
#' @importFrom RCurl getURL
#' @importFrom XML readHTMLTable xmlValue xmlEventParse
#' @importFrom data.table rbindlist
#' @importFrom pbapply startpb setpb
#' @importFrom base file
#' @importFrom stringr str_match
#' @examples
#' \dontrun{build.MCDB(outfolder=tempdir(), testMode=TRUE)}
build.MCDB <- function(outfolder, testMode = FALSE) {
  Description <- NULL
  oldpar <- options()
  options(stringsAsFactors = FALSE)
  on.exit(options(oldpar))
  file.url <- "http://www.mcdb.ca/system/downloads/current/milk_metabolites.zip"
  base.loc <- file.path(outfolder, "mcdb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "MCDB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc)
  input <- file.path(base.loc, "milk_metabolites.xml")
  header <- readLines(input, n = 10)
  version <- trimws(gsub(grep(pattern = "<version", header, value = TRUE), pattern = "<\\/?version>", replacement = ""))
  date <- trimws(gsub(grep(pattern = "update_date", header, value = TRUE), pattern = "<\\/?update_date>", replacement = ""))
  theurl <- RCurl::getURL("http://www.mcdb.ca/statistics", .opts = list(ssl.verifypeer = FALSE))
  tables <- XML::readHTMLTable(theurl)
  stats <- data.table::rbindlist(tables)
  n <- as.numeric(as.character(gsub(x = stats[Description == "Total Metabolites"]$Count, pattern = ",", replacement = "")))
  if (testMode) {
    n <- 10
  }
  envir <- environment()
  envir$db.formatted <- data.frame(compoundname = rep(NA, n), baseformula = rep(NA, n), identifier = rep(NA, n), structure = rep(NA, n), charge = rep(NA, n), description = rep("", n))
  envir$idx <- 1
  envir$maxn <- n
  envir$pb <- pbapply::startpb(min = envir$idx, max = envir$maxn)
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") {
      os <- "osx"
    }
  }
  else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)) {
      os <- "osx"
    }
    if (grepl("linux-gnu", R.version$os)) {
      os <- "linux"
    }
  }
  if (tolower(os) == "windows") {
    acc <- "primary"
    nm <- "primary"
    desc <- "primary"
    con <- base::file(input, "r")
    while (TRUE) {
      line <- readLines(con, n = 1, skipNul = TRUE)
      if (length(line) == 0 | envir$idx == envir$maxn) {
        break
      }
      if (line == "</metabolite>") {
        envir$idx <- envir$idx + 1
        pbapply::setpb(envir$pb, envir$idx)
        acc <- "primary"
        nm <- "primary"
        desc <- "primary"
      }
      tag <- stringr::str_match(line, pattern = "<(.*?)>")[, 2]
      switch(tag, accession = {
        if (acc == "primary") {
          envir$db.formatted[envir$idx, ]$identifier <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
          acc <- "secondary"
        }
      }, name = {
        if (nm == "primary") {
          envir$db.formatted[envir$idx, ]$compoundname <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
          nm <- "secondary"
        }
      }, smiles = {
        envir$db.formatted[envir$idx, ]$structure <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
      }, description = {
        if (desc == "primary") {
          envir$db.formatted[envir$idx, ]$description <- paste0(envir$db.formatted[envir$idx, ]$description, " HMDB: ", trimws(gsub(line, pattern = "(<.*?>)", replacement = "")))
          desc <- "secondary"
        }
      }, cs_description = {
        envir$db.formatted[envir$idx, ]$description <- paste0(envir$db.formatted[envir$idx, ]$description, "From ChemSpider: ", trimws(gsub(line, pattern = "(<.*?>)", replacement = "")))
      }, chemical_formula = {
        envir$db.formatted[envir$idx, ]$baseformula <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
      })
    }
    close(con)
  }
  else {
    metabolite <- function(currNode, currEnvir = envir) {
      if (currEnvir$idx %% 1000 == 0) {
        pbapply::setpb(currEnvir$pb, currEnvir$idx)
      }
      if (currEnvir$idx == currEnvir$maxn) {
        return()
      }
      currEnvir$db.formatted[currEnvir$idx, "compoundname"] <- XML::xmlValue(currNode[["name"]])
      currEnvir$db.formatted[currEnvir$idx, "identifier"] <- XML::xmlValue(currNode[["accession"]])
      currEnvir$db.formatted[currEnvir$idx, "baseformula"] <- XML::xmlValue(currNode[["chemical_formula"]])
      currEnvir$db.formatted[currEnvir$idx, "structure"] <- XML::xmlValue(currNode[["smiles"]])
      currEnvir$db.formatted[currEnvir$idx, "description"] <- paste(XML::xmlValue(currNode[["description"]]))
      x <- currNode[["predicted_properties"]]
      properties <- currNode[["predicted_properties"]]
      currEnvir$db.formatted[currEnvir$idx, "charge"] <- stringr::str_match(XML::xmlValue(properties), pattern = "formal_charge([+|\\-]\\d*|\\d*)")[, 2]
      currEnvir$idx <- currEnvir$idx + 1
    }
    XML::xmlEventParse(input, branches = list(metabolite = metabolite), replaceEntities = TRUE)
  }
  list(db = envir$db.formatted, version = version)
}

#' @title Build HMDB
#' @description Parses the HMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[XML]{readHTMLTable}},\code{\link[XML]{xmlValue}},\code{\link[XML]{xmlEventParse}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[pbapply]{pboptions}}
#'  \code{\link[base]{connections}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.HMDB
#' @export
#' @importFrom utils download.file unzip
#' @importFrom RCurl getURL
#' @importFrom XML readHTMLTable xmlValue xmlEventParse
#' @importFrom data.table rbindlist
#' @importFrom pbapply startpb setpb
#' @importFrom base file
#' @importFrom stringr str_match
#' @examples
#' \dontrun{build.HMDB(outfolder=tempdir(), testMode=TRUE)}
build.HMDB <- function(outfolder, testMode = FALSE) {
  Description <- DESCRIPTION <- NULL
  oldpar <- options()
  options(stringsAsFactors = FALSE)
  on.exit(options(oldpar))
  file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
  base.loc <- file.path(outfolder, "hmdb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "HMDB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc)
  input <- file.path(base.loc, "hmdb_metabolites.xml")
  header <- readLines(input, n = 10)
  version <- trimws(gsub(grep(pattern = "<version", header, value = TRUE), pattern = "<\\/?version>", replacement = ""))
  date <- trimws(gsub(grep(pattern = "update_date", header, value = TRUE), pattern = "<\\/?update_date>", replacement = ""))
  theurl <- RCurl::getURL("https://hmdb.ca/statistics", .opts = list(ssl.verifypeer = FALSE))
  tables <- XML::readHTMLTable(theurl)
  stats <- data.table::rbindlist(tables)
  n <- as.numeric(as.character(gsub(x = stats[Description == "Total Number of Metabolites"]$Count, pattern = ",", replacement = "")))
  if (testMode) {
    n <- 10
  }
  envir <- environment()
  envir$db.formatted <- data.frame(compoundname = rep(NA, n), baseformula = rep(NA, n), identifier = rep(NA, n), structure = rep(NA, n), charge = rep(NA, n), description = rep("", n))
  envir$idx <- 1
  envir$maxn <- n
  envir$pb <- pbapply::startpb(min = 1, max = envir$maxn)
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") {
      os <- "osx"
    }
  }
  else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)) {
      os <- "osx"
    }
    if (grepl("linux-gnu", R.version$os)) {
      os <- "linux"
    }
  }
  if (tolower(os) == "windows") {
    acc <- "primary"
    nm <- "primary"
    desc <- "primary"
    con <- base::file(input, "r")
    while (TRUE) {
      line <- readLines(con, n = 1, skipNul = TRUE)
      if (length(line) == 0 | envir$idx == envir$maxn) {
        break
      }
      if (line == "</metabolite>") {
        envir$idx <- envir$idx + 1
        if (envir$idx %% 100 == 0) {
          pbapply::setpb(envir$pb, envir$idx)
        }
        acc <- "primary"
        nm <- "primary"
        desc <- "primary"
      }
      if (envir$idx >= 2) {
        tag <- stringr::str_match(line, pattern = "<(.*?)>")[, 2]
        switch(tag, accession = {
          if (acc == "primary") {
            envir$db.formatted[envir$idx, ]$identifier <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
            acc <- "secondary"
          }
        }, name = {
          if (nm == "primary") {
            envir$db.formatted[envir$idx, ]$compoundname <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
            nm <- "secondary"
          }
        }, smiles = {
          envir$db.formatted[envir$idx, ]$structure <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
        }, description = {
          if (desc == "primary") {
            envir$db.formatted[envir$idx, ]$description <- paste0(envir$db.formatted[envir$idx, ]$description, " HMDB: ", trimws(gsub(line, pattern = "(<.*?>)", replacement = "")))
            desc <- "secondary"
          }
        }, cs_description = {
          envir$db.formatted[envir$idx, ]$description <- paste0(envir$db.formatted[envir$idx, ]$description, "From ChemSpider: ", trimws(gsub(line, pattern = "(<.*?>)", replacement = "")))
        }, chemical_formula = {
          envir$db.formatted[envir$idx, ]$baseformula <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
        })
      }
    }
    close(con)
  }
  else {
    metabolite <- function(currNode, currEnvir = envir) {
      if (currEnvir$idx %% 1000 == 0) {
        pbapply::setpb(currEnvir$pb, currEnvir$idx)
      }
      if (currEnvir$idx == currEnvir$maxn) {
        return()
      }
      currEnvir$db.formatted[currEnvir$idx, "compoundname"] <- XML::xmlValue(currNode[["name"]])
      currEnvir$db.formatted[currEnvir$idx, "identifier"] <- XML::xmlValue(currNode[["accession"]])
      currEnvir$db.formatted[currEnvir$idx, "baseformula"] <- XML::xmlValue(currNode[["chemical_formula"]])
      currEnvir$db.formatted[currEnvir$idx, "structure"] <- XML::xmlValue(currNode[["smiles"]])
      currEnvir$db.formatted[currEnvir$idx, "description"] <- paste("HMDB:", XML::xmlValue(currNode[["cs_description"]]), "CHEMSPIDER:", XML::xmlValue(currNode[["description"]]))
      x <- currNode[["predicted_properties"]]
      properties <- currNode[["predicted_properties"]]
      currEnvir$db.formatted[currEnvir$idx, "charge"] <- stringr::str_match(XML::xmlValue(properties), pattern = "formal_charge([+|\\-]\\d*|\\d*)")[, 2]
      currEnvir$idx <- currEnvir$idx + 1
    }
    XML::xmlEventParse(input, branches = list(metabolite = metabolite), replaceEntities = TRUE)
  }
  list(db = envir$db.formatted, version = version)
}

#' @title Build METACYC
#' @description Parses  METACYC, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @details Requires account creation! Then download SmartTable from 'https://trmetacyc.org/group?id=biocyc17-31223-3787684059' as 'All_compounds_of_MetaCyc.txt' and save in the databases/metacyc_source folder.
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[data.table]{fread}}
#'  \code{\link[pbapply]{pbapply}}
#' @rdname build.METACYC
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom data.table fread data.table
#' @importFrom pbapply pbsapply
build.METACYC <- function(outfolder) {
  base.loc <- file.path(outfolder, "metacyc_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  source.file <- file.path(base.loc, "All_compounds_of_MetaCyc.txt")
  if (!file.exists(source.file)) {
    msg <- "Please download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3787684059' as 'All_compounds_of_MetaCyc.txt' and save in the databases/metacyc_source folder."
    message(msg)
    return(NULL)
  }
  theurl <- RCurl::getURL("https://metacyc.org/group?id=biocyc17-31223-3787684059", .opts = list(ssl.verifypeer = FALSE))
  header <- stringr::str_match(theurl, pattern = "Created: (.*), Last Modified: (.*), Readable")
  date <- header[, 2]
  version <- header[, 3]
  metacyc.raw <- data.table::fread(source.file, fill = TRUE, quote = "", sep = "\t")
  colnames(metacyc.raw) <- gsub(x = as.character(colnames(metacyc.raw)), pattern = "\\\"", replacement = "")
  metacyc.raw[] <- lapply(metacyc.raw, gsub, pattern = "\\\"", replacement = "")
  compounds <- pbapply::pbsapply(metacyc.raw$`Common-Name`, FUN = function(pw) {
    pw <- iconv(pw, "latin1", "UTF-8", sub = "")
    pw <- pw[pw != " // "]
    pw <- gsub(pw, pattern = "&", replacement = "")
    pw <- gsub(pw, pattern = ";", replacement = "")
    res <- gsub(pw, pattern = "<(i|\\/i)>", replacement = "", perl = TRUE)
    res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "", perl = TRUE)
    paste0(res, collapse = " --- ")
  })
  db.formatted <- data.table::data.table(compoundname = compounds, description = sapply(1:nrow(metacyc.raw), function(i) {
    species <- metacyc.raw$Species[i]
    paste0(metacyc.raw$Summary[i], if (species != "") {
      paste0(" Found in: ", species, ".")
    } else {
      ""
    })
  }), baseformula = c(NA), identifier = metacyc.raw$`Object ID`, charge = c(0), structure = metacyc.raw$SMILES)
  list(db = db.formatted, version = version)
}

#' @title Build CHEBI
#' @description Parses CHEBI, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[utils]{download.file}}
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[ChemmineR]{datablock2ma}},\code{\link[ChemmineR]{datablock}}
#' @rdname build.CHEBI
#' @export
#' @importFrom RCurl getURL
#' @importFrom utils download.file
#' @importFrom data.table as.data.table
#' @importFrom ChemmineR datablock2ma datablock
#' @examples
#' \dontrun{build.CHEBI(outfolder=tempdir())}
build.CHEBI <- function(outfolder) {
  file.url <- "ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete.sdf.gz"
  base.loc <- file.path(outfolder, "chebi_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  zip.file <- file.path(base.loc, "chebi.sdf.gz")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  version <- Sys.Date()
  sdf.path <- zip.file
  desc <- function(sdfset) {
    mat <- NULL
    try(
      {
        db <- data.table::as.data.table(ChemmineR::datablock2ma(datablocklist = ChemmineR::datablock(sdfset)))
        mat <- data.table::data.table(identifier = gsub("CHEBI:", "", as.character(db$`ChEBI ID`)), compoundname = as.character(db$`ChEBI Name`), baseformula = as.character(db$Formulae), structure = sapply(1:nrow(db), function(i) {
          if (!is.empty(db$SMILES[i])) {
            db$SMILES[i]
          } else {
            db$InChI[i]
          }
        }), description = db$Definition)
        empty.smiles <- sapply(mat$structure, is.empty)
        mat <- mat[!empty.smiles, ]
        mat <- as.matrix(mat)
      },
      silent = TRUE
    )
    mat
  }
  parsedLoc <- file.path(base.loc, "chebi_parsed.csv")
  sdfStream.joanna(input = sdf.path, output = parsedLoc, append = FALSE, fct = desc, silent = TRUE)
  db.formatted <- data.table::fread(parsedLoc, fill = TRUE, header = TRUE)
  db.formatted$V1 <- NULL
  url <- "https://www.ebi.ac.uk/chebi/statisticsForward.do"
  version <- as.character(url %>% xml2::read_html() %>% rvest::html_nodes(xpath = "//*[@id=\"content\"]/h4"))
  version <- gsub("(.*Release)|(<.*?>)|\n", "", version)
  version <- trimws(version)
  list(db = db.formatted, version = version)
}

#' @title Build PharmGKB
#' @description Parses PharmGKB drugs and chemicals (drug metabolites, etc.), returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[data.table]{fread}}
#' @rdname build.PHARMGKB
#' @export
#' @importFrom utils download.file untar
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom data.table fread data.table
#' @examples
#' \dontrun{build.PHARMGKB(outfolder=tempdir())}
build.PHARMGKB <- function(outfolder) {
  SMILES <- `PharmGKB Accession Id` <- Name <- NULL
  file.url <- "https://s3.pgkb.org/data/drugs.zip"
  base.loc <- file.path(outfolder, "pharmgkb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "drugs.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  drug.db <- data.table::fread(file.path(base.loc, "drugs.tsv"), quote = "")
  drug.db <- drug.db[SMILES != ""]
  file.url <- "https://s3.pgkb.org/data/chemicals.zip"
  zip.file <- file.path(base.loc, "chemicals.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  chem.db <- data.table::fread(file.path(base.loc, "chemicals.tsv"), quote = "")
  chem.db <- chem.db[SMILES != ""]
  both.db <- data.table::rbindlist(list(drug.db, chem.db))
  cols <- c("Generic Names", "Trade Names", "Brand Mixtures", "Type", "Dosing Guideline", "Clinical Annotation Count", "Variant Annotation Count", "Pathway Count", "VIP Count", "Dosing Guideline Sources", "Top Clinical Annotation Level", "Top FDA Label Testing Level", "Top Any Drug Label Testing Level", "Label Has Dosing Info", "Has Rx Annotation")
  db.formatted <- unique(both.db[, list(compoundname = Name, description = do.call(paste, c(lapply(cols, function(x) paste(x, get(x), sep = ": ")), sep = ", ")), identifier = `PharmGKB Accession Id`, structure = SMILES)])
  db.formatted$description <- gsub("\",", " - ", db.formatted$description)
  db.formatted$description <- gsub("\"", "", db.formatted$description)
  db.formatted$description <- gsub(" (,*[\\w| ]+: [0|,|No]+)|([\\w| ]+: [0|No]*$)|^([\\w| ]+: [0|,|No]+)", "", db.formatted$description, perl = TRUE)
  db.formatted$description <- gsub("^ |,$", "", db.formatted$description, perl = TRUE)
  db.formatted$baseformula <- c(NA)
  db.formatted$charge <- c(0)
  date <- Sys.Date()
  version <- Sys.Date()
  list(db = db.formatted[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")], version = version)
}

#' @title Build FOODB
#' @description Parses the FOODB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{untar}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[data.table]{fread}}
#' @rdname build.FOODB
#' @export
#' @importFrom utils download.file untar
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom data.table fread data.table
#' @examples
#' \dontrun{build.FOODB(outfolder=tempdir())}
build.FOODB <- function(outfolder) {
  . <- orig_health_effect_name <- compound_id <- NULL
  file.url <- "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
  base.loc <- file.path(outfolder, "foodb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "foodb.tar.gz")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::untar(normalizePath(zip.file), exdir = normalizePath(base.loc))
  theurl <- RCurl::getURL("https://foodb.ca/downloads", .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(theurl, pattern = "FooDB Version <strong>(...)<")[, 2]
  date <- stringr::str_match(theurl, pattern = "FooDB CSV file<\\/td><td>(.{3,40})<\\/td>")[, 2]
  date <- as.Date(date, format = "%B%e%Y")
  subdate <- gsub(date, pattern = "-", replacement = "_")
  base.table <- data.table::fread(file = file.path(base.loc, paste0("foodb_", subdate, "_csv"), "Compound.csv"))
  health <- data.table::fread(file = file.path(base.loc, paste0("foodb_", subdate, "_csv"), "CompoundsHealthEffect.csv"))
  health_summ <- health[, .(health_effect = list(orig_health_effect_name)), by = compound_id]
  base.merged <- merge(base.table, health_summ, by.x = "id", by.y = "compound_id", all.x = TRUE)
  db.formatted <- data.table::data.table(compoundname = base.merged$name, description = paste0(base.table$annotation_quality, " ASSOCIATED HEALTH EFFECTS:", base.merged$health_effect), baseformula = c(NA), identifier = base.merged$public_id, charge = c(0), structure = base.merged$cas_number)
  db.formatted$description <- gsub("ASSOCIATED HEALTH EFFECTS:NULL", "", db.formatted$description)
  db.formatted <- unique(db.formatted)
  list(db = db.formatted, version = version)
}

#' @title Build Wikidata
#' @description Parses wikidata chemical compound database, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[WikidataQueryServiceR]{query_wikidata}}
#'  \code{\link[data.table]{as.data.table}}
#' @rdname build.WIKIDATA
#' @export
#' @importFrom WikidataQueryServiceR query_wikidata
#' @importFrom data.table as.data.table data.table
#' @examples
#' \dontrun{build.WIKIDATA(outfolder=tempdir())}
build.WIKIDATA <- function(outfolder) {
  date <- Sys.Date()
  version <- Sys.Date()
  sparql_query <- "PREFIX wd: <http://www.wikidata.org/entity/>
  PREFIX wds: <http://www.wikidata.org/entity/statement/>
  PREFIX wdv: <http://www.wikidata.org/value/>
  PREFIX wdt: <http://www.wikidata.org/prop/direct/>
  PREFIX wikibase: <http://wikiba.se/ontology#>
  PREFIX p: <http://www.wikidata.org/prop/>
  PREFIX ps: <http://www.wikidata.org/prop/statement/>
  PREFIX pq: <http://www.wikidata.org/prop/qualifier/>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX bd: <http://www.bigdata.com/rdf#>

  SELECT ?chemical_compound ?chemical_compoundLabel ?chemical_formula ?chemical_compoundDescription ?canonical_SMILES ?roleLabel WHERE {
  SERVICE wikibase:label { bd:serviceParam wikibase:language \"en, de\". }
  ?chemical_compound wdt:P31 wd:Q11173.
  ?chemical_compound wdt:P274 ?chemical_formula.
  ?chemical_compound wdt:P233 ?canonical_SMILES.
  OPTIONAL {?chemical_compound wdt:P2868 ?role.}}"
  db.1 <- WikidataQueryServiceR::query_wikidata(sparql_query, format = "simple")
  db.1$chemical_compound <- basename(db.1$chemical_compound)
  db.1$chemical_compoundDescription[db.1$chemical_compoundDescription == "chemical compound"] <- NA
  db.1$roleLabel[db.1$roleLabel == ""] <- NA
  db.2 <- data.table::as.data.table(aggregate(db.1, by = list(db.1$chemical_compoundLabel, db.1$chemical_formula), function(x) c(unique((x)))))
  db.2$description <- apply(db.2[, c("roleLabel", "chemical_compoundDescription")], 1, FUN = function(x) {
    x <- unlist(x)
    paste(x[!is.na(x)], collapse = ", ")
  })
  db.formatted <- data.table::data.table(compoundname = as.character(db.2$chemical_compoundLabel), description = as.character(db.2$description), baseformula = as.character(db.2$chemical_formula), identifier = as.character(db.2$chemical_compound), charge = c(), structure = as.character(db.2$canonical_SMILES))

  from = "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089"
  to = "0123456789"

  db.formatted$baseformula <- chartr(from, to, db.formatted$baseformula)
  db.formatted$description[db.formatted$description == ""] <- "Unknown"
  db.formatted$baseformula <- as.character(db.formatted$baseformula)
  db.formatted <- db.formatted[!is.na(db.formatted$baseformula), ]
  list(db = db.formatted, version = version)
}

#' @title Build RESPECT
#' @description Parses RESPECT, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname build.RESPECT
#' @export
#' @importFrom utils download.file unzip
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table rbindlist
#' @examples
#' \dontrun{build.RESPECT(outfolder=tempdir(), testMode=TRUE)}
build.RESPECT <- function(outfolder, testMode = FALSE) {
  baseformula <- NULL
  file.url <- "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"
  base.loc <- file.path(outfolder, "respect_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "respect.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = (base.loc))
  theurl <- RCurl::getURL("http://spectra.psc.riken.jp/menta.cgi/respect/download/download", .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(theurl, pattern = "<b>(.{1,30}) Update")[, 2]
  version <- gsub(version, pattern = "\\.", replacement = "-")
  date <- version
  cpd_files <- list.files(base.loc, full.names = TRUE)
  db_rows <- pbapply::pblapply(cpd_files, function(fn) {
    row <- data.table::data.table()
    try({
      lines <- readLines(fn, skipNul = TRUE, n = 100)
      split.lines <- sapply(lines, strsplit, ": ")
      names(split.lines) <- sapply(split.lines, function(x) x[1])
      split.lines <- lapply(split.lines, function(x) x[2:length(x)])
      row <- data.table::data.table(compoundname = split.lines$`CH$NAME`, description = split.lines$RECORD_TITLE, baseformula = split.lines$`CH$FORMULA`, identifier = split.lines$ACCESSION, charge = c(0), structure = split.lines$`CH$SMILES`)
      if (row$structure == "N/A") {
        row$structure <- split.lines$`CH$INCHI`
      }
    })
    row
  })
  db.formatted <- data.table::rbindlist(db_rows[!is.na(db_rows)], fill = TRUE)
  db.formatted <- unique(db.formatted[!is.na(baseformula), ])
  db.formatted <- aggregate(db.formatted, by = list(db.formatted$compoundname), FUN = function(x) paste0(unique(x), collapse = "/"))
  db.formatted <- db.formatted[, -1]
  list(db = db.formatted, version = version)
}

#' @title Build MACONDA
#' @description Parses MACONDA, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @param conn Connection to extended database (MaConDa writes directly to there due to anomalous adducts)
#' @param apikey ChemSpider API key
#' @return data table with parsed database
#' @seealso
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[data.table]{fread}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[DBI]{dbDisconnect}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname build.MACONDA
#' @export
#' @importFrom stringr str_match
#' @importFrom utils download.file unzip
#' @importFrom data.table fread data.table
#' @importFrom pbapply pbsapply
#' @importFrom RSQLite dbExecute dbConnect SQLite dbGetQuery dbWriteTable dbDisconnect
#' @importFrom DBI dbDisconnect
#' @importFrom gsubfn fn
#' @examples
#' \dontrun{build.MACONDA(outfolder=tempdir(), testMode=TRUE)}
build.MACONDA <- function(outfolder, testMode = FALSE, conn, apikey) {
  Name <- NULL
  file.url <- "https://www.maconda.bham.ac.uk/downloads/MaConDa__v1_0__csv.zip"
  base.loc <- file.path(outfolder, "maconda_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "maconda.zip")
  theurl <- paste0(readLines("https://www.maconda.bham.ac.uk/downloads.php"), collapse = "")
  version <- stringr::str_match(theurl, pattern = "Version (...)")[, 2]
  version <- gsub(version, pattern = "\\.", replacement = "-")
  date <- Sys.Date()
  utils::download.file(file.url, zip.file, mode = "wb", extra = "-k", method = "auto")
  utils::unzip(normalizePath(zip.file), files = "MaConDa__v1_0__extensive.csv", exdir = normalizePath(base.loc))
  base.table <- data.table::fread(file = file.path(base.loc, "MaConDa__v1_0__extensive.csv"))
  mysterious <- which(base.table$name == "unknown")
  base.table$formula[mysterious] <- paste0("IDK", 1:length(mysterious))
  has.inchi <- which(base.table$std_inchi != "")
  inchis <- base.table$std_inchi[has.inchi]
  smiles <- pbapply::pbsapply(inchis, function(x) {
    url <- gsubfn::fn$paste("https://cactus.nci.nih.gov/chemical/structure/$x/smiles")
    Sys.sleep(0.1)
    RCurl::getURL(url)
  })
  charges <- gsub(base.table$ion_form, pattern = ".*\\]", replacement = "")
  no.info <- which(charges == "")
  charges[no.info] <- sapply(base.table$ion_mode[no.info], function(x) ifelse(x == "POS", 1, -1))
  plus.minus <- which(charges == "+" | charges == "-")
  charges[plus.minus] <- sapply(charges[plus.minus], function(ch) switch(ch, `+` = 1, `-` = -1))
  db.base <- data.table::data.table(compoundname = base.table$name, description = paste(base.table$type_of_contaminant, base.table$ion_source_type, base.table$ion_mode), baseformula = base.table$formula, identifier = base.table$id, charge = charges, structure = paste0("[", base.table$formula, "]", charges))
  success <- !grepl("404", smiles)
  db.base$structure[has.inchi[success]] <- smiles[success]
  db.final <- cleanDB(db.base, cl = 0, silent = TRUE, blocksize = 400)
  dbpath <- file.path(outfolder, "maconda.db")
  if (file.exists(dbpath)) {
    file.remove(dbpath)
  }
  conn <- openBaseDB(outfolder, "maconda.db")
  writeDB(conn, data.table::data.table(date = Sys.Date(), version = version), "metadata")
  writeDB(conn, db.final, "base")
  RSQLite::dbExecute(conn, "CREATE INDEX b_idx1 ON base(structure)")
  DBI::dbDisconnect(conn)
  full.db <- file.path(outfolder, "extended.db")
  first.db <- !file.exists(full.db)
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA foreign_keys = ON"))
  RSQLite::dbExecute(full.conn, "CREATE TABLE IF NOT EXISTS structures(struct_id INT PRIMARY KEY,
                                                                       smiles TEXT,
                                                                       UNIQUE(struct_id, smiles))")
  RSQLite::dbExecute(full.conn, strwrap("CREATE TABLE IF NOT EXISTS extended(struct_id INT,
                                                                             fullmz decimal(30,13),
                                                                             adduct text,
                                                                             isoprevalence float)", width = 10000, simplify = TRUE))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA auto_vacuum = 1;"))
  base.db <- normalizePath(file.path(outfolder, "maconda.db"))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))
  if (first.db) {
    RSQLite::dbExecute(full.conn, "CREATE INDEX st_idx1 ON structures(smiles)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx1 on extended(struct_id)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx2 on extended(fullmz)")
    RSQLite::dbExecute(full.conn, "PRAGMA journal_mode=WAL;")
  }
  adducts_maconda <- unique(base.table[, c("ion_form", "ion_mode")])
  adducts_maconda <- adducts_maconda[adducts_maconda$ion_form != ""]
  adducts_maconda$ion_mode <- sapply(adducts_maconda$ion_mode, function(mode) switch(mode,
                                                                                     POS = "positive",
                                                                                     NEG = "negative"))
  adduct_table_maconda <- data.table::data.table(Name = adducts_maconda$ion_form, Ion_mode = adducts_maconda$ion_mode, Charge = c(NA), xM = c(NA), AddAt = c(NA), RemAt = c(NA), AddEx = c(NA), RemEx = c(NA), Nelec = c(NA), Rule = c(NA), Info = c("MACONDA ONLY"))
  if (!first.db) {
    new_adducts <- setdiff(adduct_table_maconda$Name, RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT Name FROM adducts")[, 1])
    if (length(new_adducts) == 0) {
      print("maconda is already in here!")
      return(NA)
    }
    RSQLite::dbWriteTable(full.conn, "adducts", adduct_table_maconda[Name %in% new_adducts], append = TRUE)
  }
  else {
    new_adducts <- c()
    RSQLite::dbWriteTable(full.conn, "adducts", adduct_table_maconda, overwrite = TRUE)
  }
  if (first.db | length(new_adducts) > 0) {
    to.do <- RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                             FROM tmp.base")
    to.do$struct_id <- c(NA)
  }
  else {
    to.do <- RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                             FROM tmp.base LEFT JOIN structures str
                                             ON base.structure = str.smiles
                                             WHERE str.smiles IS NULL")
  }
  if (nrow(to.do) == 0) {
    print("all already done")
  }
  done.structures <- if (first.db) {
    0
  } else {
    RSQLite::dbGetQuery(full.conn, "SELECT MAX(struct_id) FROM structures")[, 1]
  }
  start.id <- done.structures + 1
  structs <- unique(db.final[, c(structure)])
  mapper <- data.table::data.table(struct_id = seq(start.id, start.id + nrow(to.do) - 1, 1), smiles = to.do$structure)
  base.table$structure <- db.final$structure[match(base.table$id, table = db.final$identifier)]
  adduct.unknown <- which(base.table$ion_form == "")
  base.table$exact_adduct_mass[adduct.unknown] <- base.table$mz[adduct.unknown]
  base.table$ion_form[adduct.unknown] <- sapply(base.table$ion_mode[adduct.unknown], function(ion_mode) switch(ion_mode, POS = c("[M?]+?"), NEG = c("[M?]-?")))
  meta.table <- data.table::data.table(fullmz = base.table$exact_adduct_mass, adduct = base.table$ion_form, isoprevalence = c(100), structure = base.table$structure)
  ids <- mapper$struct_id[match(meta.table$structure, mapper$smiles)]
  meta.table$struct_id <- ids
  maconda.extended <- unique(meta.table[, -"structure"])
  RSQLite::dbWriteTable(full.conn, "extended", maconda.extended, append = TRUE)
  RSQLite::dbWriteTable(full.conn, "structures", mapper, append = TRUE)
  RSQLite::dbDisconnect(full.conn)
}

#' @title Build T3DB
#' @description Parses the T3DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[data.table]{fread}}
#' @rdname build.T3DB
#' @export
#' @importFrom utils download.file unzip
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom data.table fread data.table
#' @examples
#' \dontrun{build.T3DB(outfolder=tempdir(), testMode=TRUE)}
build.T3DB <- function(outfolder, testMode = FALSE) {
  file.url <- "http://www.t3db.ca/system/downloads/current/toxins.csv.zip"
  base.loc <- file.path(outfolder, "t3db_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "T3DB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  theurl <- RCurl::getURL("http://www.t3db.ca/downloads", .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(theurl, pattern = "T3DB Version <strong>(...)")[, 2]
  date <- Sys.Date()
  db.formatted <- data.table::fread(file.path(base.loc, "toxins.csv"), fill = TRUE)
  db.formatted <- data.table::data.table(compoundname = db.formatted$Name, baseformula = db.formatted$`Chemical Formula`, description = db.formatted$Description, charge = c(0), identifier = db.formatted$`T3DB ID`, structure = db.formatted$SMILES)
  list(db = db.formatted, version = version)
}

#' @title Build BLOOD EXPOSOME DB
#' @description Parses the BLOOD EXPOSOME DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}}
#'  \code{\link[openxlsx]{read.xlsx}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[jsonlite]{read_json}}
#' @rdname build.BLOODEXPOSOME
#' @export
#' @importFrom utils download.file
#' @importFrom openxlsx read.xlsx
#' @importFrom pbapply pbsapply
#' @importFrom jsonlite read_json
#' @importFrom data.table data.table
#' @examples
#' \dontrun{build.BLOODEXPOSOME(outfolder=tempdir(), testMode=TRUE)}
build.BLOODEXPOSOME <- function(outfolder, testMode = FALSE) {
  file.url <- "http://exposome1.fiehnlab.ucdavis.edu/download/BloodExpsomeDatabase_version_1.0.xlsx"
  base.loc <- file.path(outfolder, "bloodexposome_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  excel.file <- file.path(base.loc, "exposome.xlsx")
  utils::download.file(file.url, excel.file, mode = "wb", method = "auto")
  version <- "1.0"
  date <- Sys.Date()
  db.full <- openxlsx::read.xlsx(excel.file, sheet = 1, colNames = TRUE, startRow = 3)
  if (testMode) {
    db.full <- db.full[1:10, ]
  }
  print("getting iupac names from missing compound names...")
  new.names <- pbapply::pbsapply(db.full[which(sapply(db.full$Compound.Name, is.empty)), ]$PubChem.CID, function(cid) {
    try({
      url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%d/JSON", as.numeric(cid))
      Sys.sleep(0.1)
      json <- jsonlite::read_json(url)
      json$Record$RecordTitle
    })
  })
  db.full[which(sapply(db.full$Compound.Name, is.empty)), ]$Compound.Name <- new.names
  db.formatted <- unique(data.table::data.table(compoundname = db.full$Compound.Name, description = paste0("Found in ", db.full$BloodPaperCount, " papers related to blood exposome."), baseformula = db.full$Molecular.Formula, identifier = db.full$PubChem.CID, charge = db.full$Charge, structure = db.full$CanonicalSMILES))
  list(db = db.formatted, version = version)
}

#' @title Build EXPOSOME EXPLORER
#' @description Parses the EXPOSOME EXPLORER DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[data.table]{fread}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[R.utils]{capitalize}}
#' @rdname build.EXPOSOMEEXPLORER
#' @export
#' @importFrom utils download.file unzip
#' @importFrom stringr str_match
#' @importFrom RCurl getURL
#' @importFrom data.table fread data.table
#' @importFrom pbapply pbsapply
#' @importFrom R.utils decapitalize
#' @examples
#' \dontrun{build.EXPOSOMEEXPLORER(outfolder=tempdir(), testMode=TRUE)}
build.EXPOSOMEEXPLORER <- function(outfolder, testMode = FALSE) {
  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/biomarkers.csv.zip"
  base.loc <- file.path(outfolder, "expoexplorer_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "expoexpo_comp.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  version <- stringr::str_match(RCurl::getURL("http://exposome-explorer.iarc.fr/releases", .opts = list(ssl.verifypeer = FALSE)), pattern = "Release (...)")[, 2]
  date <- stringr::str_match(RCurl::getURL("http://exposome-explorer.iarc.fr/downloads", .opts = list(ssl.verifypeer = FALSE)), pattern = "(\\d{4}-\\d{2}-\\d{2})")[, 2]
  base.table <- data.table::fread(file = file.path(base.loc, "biomarkers.csv"))
  db.formatted <- data.table::data.table(compoundname = base.table$Name, description = base.table$Description, baseformula = base.table$Formula, identifier = base.table$ID, charge = c(NA), structure = base.table$SMILES)
  db.formatted <- unique(db.formatted)
  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/correlations.csv.zip"
  zip.file <- file.path(base.loc, "expoexpo_corr.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  corr.table <- data.table::fread(file = file.path(base.loc, "correlations.csv"))
  descriptions <- pbapply::pbsapply(1:nrow(corr.table), function(i) {
    row <- corr.table[i, ]
    desc <- paste("Found in", R.utils::decapitalize(row$Biospecimen), "of", R.utils::decapitalize(row$`Subject group`), "under", R.utils::decapitalize(row$Population), "in", row$Country, "after taking in", R.utils::decapitalize(row$Intake), paste0("(", row$`Analytical method`, ", p ", row$`Correlation p-value`, ")."))
  })
  corr.table$Pasted <- descriptions
  df <- corr.table[, c("Excretion ID", "Pasted")]
  aggr <- aggregate(Pasted ~ `Excretion ID`, df, function(x) toString(paste(unique(x), collapse = " ")))
  final.table <- merge(db.formatted, aggr, by.x = "identifier", by.y = "Excretion ID", all.x = TRUE)
  final.table$description <- pbapply::pbsapply(1:nrow(final.table), function(i) {
    row <- final.table[i, ]
    a <- if (!is.na(row$description) & row$description != "NA") {
      row$description
    } else {
      ""
    }
    b <- row$Pasted
    paste0(a, b)
  })
  db.formatted <- final.table[, -"Pasted"]
  list(db = db.formatted, version = version)
}

#' @title Build SMPDB
#' @description Parses the SMPDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{fread}},\code{\link[data.table]{rbindlist}}
#' @rdname build.SMPDB
#' @export
#' @importFrom utils download.file unzip
#' @importFrom RCurl getURL
#' @importFrom stringr str_match_all
#' @importFrom pbapply pblapply
#' @importFrom data.table fread rbindlist data.table
#' @examples
#' \dontrun{build.SMPDB(outfolder=tempdir(), testMode=TRUE)}
build.SMPDB <- function(outfolder, testMode = FALSE) {
  . <- description <- compoundname <- baseformula <- identifier <- NULL
  file.url <- "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
  base.loc <- file.path(outfolder, "smpdb_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  zip.file <- file.path(base.loc, "SMPDB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  theurl <- "http://smpdb.ca/release_notes"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  versions <- stringr::str_match_all(header, pattern = "Release (.{1,30}) -")[[1]]
  version <- gsub(pattern = "SMPDB ", replacement = "", x = versions[nrow(versions), 2])
  dates <- stringr::str_match_all(header, pattern = "Release .{1,30}- (.{1,30})<\\/h")[[1]]
  date <- dates[nrow(dates), 2]
  date <- as.Date(date, format = "%B%d,%Y")
  smpdb.paths <- list.files(path = base.loc, pattern = "\\.csv$", full.names = TRUE)
  if (testMode) {
    smpdb.paths <- smpdb.paths[1:10]
  }
  subtables <- pbapply::pblapply(smpdb.paths, data.table::fread)
  smpdb.tab <- data.table::rbindlist(subtables, fill = TRUE)
  db.formatted <- data.table::data.table(compoundname = smpdb.tab$`Metabolite Name`, description = paste0(smpdb.tab$`Pathway Name`, " (", smpdb.tab$`Pathway Subject`, " pathway)"), baseformula = smpdb.tab$Formula, identifier = smpdb.tab$`Metabolite ID`, structure = smpdb.tab$SMILES)
  db.formatted <- db.formatted[, .(description = paste0(unique(description), collapse = ", ")), by = list(compoundname, baseformula, identifier, structure)]
  db.formatted <- db.formatted[-1, ]
  list(db = db.formatted, version = version)
}

#' @title Build KEGG
#' @description Parses the KEGG DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[KEGGREST]{keggFind}},\code{\link[KEGGREST]{keggGet}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[rcdk]{load.molecules}},\code{\link[rcdk]{get.smiles}}
#' @rdname build.KEGG
#' @export
#' @importFrom pbapply pblapply
#' @importFrom KEGGREST keggFind keggGet
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom data.table data.table rbindlist
#' @importFrom rcdk load.molecules get.smiles get.total.formal.charge
#' @examples
#' \dontrun{build.KEGG(outfolder=tempdir(), testMode=TRUE)}
build.KEGG <- function(outfolder, testMode = FALSE) {
  batches <- split(0:2300, ceiling(seq_along(0:2300) / 100))
  cpds <- pbapply::pblapply(batches, FUN = function(batch) {
    names(KEGGREST::keggFind("compound", batch, "mol_weight"))
  })
  cpd.ids <- Reduce(c, cpds)
  if (testMode) {
    cpd.ids <- cpd.ids[1:10]
  }
  id.batches <- split(cpd.ids, ceiling(seq_along(cpd.ids) / 10))
  theurl <- "https://www.genome.jp/kegg/compound/"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Last updated: (.{1,30})<")[, 2]
  date <- as.Date(version, format = "%B%d, %Y")
  kegg.cpd.list <- pbapply::pblapply(id.batches, FUN = function(batch) {
    rest.result <- KEGGREST::keggGet(batch)
    base.list <- lapply(rest.result, FUN = function(cpd) {
      cpd$NAME_FILT <- gsub(cpd$NAME, pattern = ";", replacement = "")
      data.table::data.table(compoundname = c(paste(cpd$NAME_FILT, collapse = ", ")), description = paste0("Involved in pathways: ", paste0(cpd$PATHWAY, collapse = ", "), ". More specifically: ", paste0(cpd$MODULE, collapse = ", "), ". Also associated with compound classes: ", paste0(unique(trimws(gsub(cpd$BRITE, pattern = "\\[.*\\]|  D\\d* |\\(.*\\)|\\d*", replacement = ""))), collapse = ", ")), baseformula = c(cpd$FORMULA), identifier = c(cpd$ENTRY), charge = 0, structure = NA, pathway = if ("PATHWAY" %in%
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   names(cpd)) {
        names(cpd$PATHWAY)
      } else {
        NA
      })
    })
    res <- data.table::rbindlist(base.list)
    res
  })
  db.formatted <- data.table::rbindlist(kegg.cpd.list)
  base.loc <- file.path(outfolder, "kegg_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  kegg.mol.paths <- pbapply::pblapply(id.batches, FUN = function(batch) {
    bigmol <- KEGGREST::keggGet(batch, "mol")
    mols <- strsplit(x = paste0("\n \n \n", bigmol), split = "\\$\\$\\$\\$\n")[[1]]
    fps <- normalizePath(file.path(base.loc, paste0(stringr::str_match(mols, pattern = "<ENTRY>\ncpd:(.*)\n")[, 2], ".mol")), mustWork = FALSE)
    sapply(1:length(mols), function(i) writeLines(text = mols[[i]], con = fps[i]))
    fps
  })
  fns <- list.files(base.loc, full.names = TRUE)
  smiles.rows <- pbapply::pblapply(fns, function(fn) {
    smiles <- NA
    charge <- 0
    try({
      iatom <- rcdk::load.molecules(molfiles = fn)[[1]]
      smiles <- rcdk::get.smiles(molecule = iatom)
      charge <- rcdk::get.total.formal.charge(mol = iatom)
    })
    id <- gsub(basename(fn), pattern = "\\.mol", replacement = "")
    data.table::data.table(identifier = id, smiles = smiles, calcharge = charge)
  })
  smitable <- data.table::rbindlist(smiles.rows)
  db.merged <- merge(db.formatted, smitable, by = "identifier")
  db.formatted <- data.table::data.table(compoundname = db.merged$compoundname, description = db.merged$description, baseformula = db.merged$baseformula, identifier = db.merged$identifier, charge = db.merged$calcharge, structure = db.merged$smiles)
  list(db = db.formatted, version = version)
}

#' @title Build DRUGBANK
#' @description Parses the DRUGBANK DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @details Requires account creation! Then please download the full XML database from the website and place in databases/drugbank_source folder. Create it if it doesn't exist yet please!
#' @seealso
#'  \code{\link[utils]{unzip}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[XML]{readHTMLTable}},\code{\link[XML]{xmlToList}},\code{\link[XML]{xmlValue}},\code{\link[XML]{xmlEventParse}}
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[pbapply]{pboptions}}
#' @rdname build.DRUGBANK
#' @export
#' @importFrom utils unzip
#' @importFrom stringr str_match
#' @importFrom RCurl getURL
#' @importFrom XML readHTMLTable xmlToList xmlValue xmlEventParse
#' @importFrom data.table as.data.table
#' @importFrom pbapply startpb setpb
build.DRUGBANK <- function(outfolder) {
  Description <- NULL
  base.loc <- file.path(outfolder, "drugbank_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  zip.file <- file.path(base.loc, "drugbank.zip")
  if (length(list.files(base.loc)) == 0) {
    msg <- "Please create a DrugBank account to download the database: https://www.drugbank.ca/releases/latest. Save the .xml.zip file in the databases/drugbank_source folder."
    message(msg)
    return(NULL)
  }
  if (!file.exists(zip.file)) {
    file.rename(file.path(base.loc, "drugbank_all_full_database.xml.zip"), zip.file)
  }
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  input <- file.path(base.loc, "full database.xml")
  header <- readLines(input, n = 10)
  hasInfo <- grep(x = header, pattern = "version", value = TRUE, perl = TRUE)[2]
  version <- stringr::str_match(string = hasInfo, pattern = "version=\"(.*)\" exported")[, 2]
  theurl <- RCurl::getURL("https://www.drugbank.ca/stats", .opts = list(ssl.verifypeer = FALSE))
  tables <- XML::readHTMLTable(theurl, header = FALSE)
  stats <- data.table::as.data.table(tables[[1]])
  colnames(stats) <- c("Description", "Count")
  n <- as.numeric(as.character(gsub(x = stats[Description == "Total Number of Drugs"]$Count, pattern = ",", replacement = "")))
  envir <- environment()
  envir$db.formatted <- data.frame(compoundname = rep(NA, n), baseformula = rep(NA, n), identifier = rep(NA, n), structure = rep(NA, n), charge = rep(NA, n), description = rep(NA, n))
  envir$pb <- pbapply::startpb(min = 0, max = n)
  envir$idx <- 0
  metabolite <- function(currNode, currEnvir = envir) {
    if (currEnvir$idx %% 10 == 0) {
      pbapply::setpb(currEnvir$pb, currEnvir$idx)
    }
    currEnvir$idx <- currEnvir$idx + 1
    properties <- currNode[["calculated-properties"]]
    if (is.null(properties)) {
      properties <- currNode[["experimental-properties"]]
    }
    proplist <- XML::xmlToList(properties)
    if (length(proplist) == 0) {
      return(NULL)
    }
    which.form <- which(sapply(proplist, function(x) {
      if ("kind" %in% names(x)) {
        res <- x[["kind"]] == "Molecular Formula"
      }
      else {
        res <- FALSE
      }
      res
    }))
    which.struc <- which(sapply(proplist, function(x) {
      if ("kind" %in% names(x)) {
        res <- x[["kind"]] == "SMILES"
      }
      else {
        res <- FALSE
      }
      res
    }))
    which.charge <- which(sapply(proplist, function(x) {
      if ("kind" %in% names(x)) {
        res <- x[["kind"]] == "Physiological Charge"
      }
      else {
        res <- FALSE
      }
      res
    }))
    if (length(which.form) == 0 & length(which.struc) == 0) {
      return(NULL)
    }
    currEnvir$db.formatted[currEnvir$idx, "compoundname"] <- XML::xmlValue(currNode[["name"]])
    currEnvir$db.formatted[currEnvir$idx, "identifier"] <- XML::xmlValue(currNode[["drugbank-id"]])
    currEnvir$db.formatted[currEnvir$idx, "baseformula"] <- proplist[[which.form]][["value"]]
    currEnvir$db.formatted[currEnvir$idx, "structure"] <- if (length(which.struc) > 0) {
      proplist[[which.struc]][["value"]]
    }
    else {
      ""
    }
    currEnvir$db.formatted[currEnvir$idx, "description"] <- XML::xmlValue(currNode[["description"]])
    currEnvir$db.formatted[currEnvir$idx, "charge"] <- if (length(which.charge) > 0) {
      proplist[[which.charge]][["value"]]
    }
    else {
      0
    }
  }
  res <- XML::xmlEventParse(file = input, branches = list(drug = metabolite, `drugbank-metabolite-id-value` = print))
  envir$db.formatted <- envir$db.formatted[-1, ]
  list(db = envir$db.formatted, version = version)
}

#' @title Build LIPID MAPS
#' @description Parses the LIPID MAPS DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @param apikey ChemSpider API key
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}}
#'  \code{\link[zip]{unzip}}
#'  \code{\link[data.table]{as.data.table}},\code{\link[data.table]{fread}},\code{\link[data.table]{rbindlist}}
#'  \code{\link[ChemmineR]{datablock2ma}},\code{\link[ChemmineR]{datablock}}
#'  \code{\link[xml2]{read_xml}}
#'  \code{\link[rvest]{html_nodes}},\code{\link[rvest]{html_text}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[stringi]{stri_detect}}
#' @rdname build.LIPIDMAPS
#' @export
#' @importFrom utils download.file
#' @importFrom zip unzip
#' @importFrom data.table as.data.table data.table fread rbindlist
#' @importFrom ChemmineR datablock2ma datablock
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_text
#' @importFrom stringr str_match_all str_match
#' @importFrom pbapply pblapply pbsapply
#' @importFrom stringi stri_detect_fixed
#' @examples
#' \dontrun{build.LIPIDMAPS(outfolder=tempdir(), testMode=TRUE)}
build.LIPIDMAPS <- function(outfolder, testMode = FALSE, apikey) {
  file.url <- "https://www.lipidmaps.org/files/?file=LMSD_20191002&ext=sdf.zip"
  base.loc <- file.path(outfolder, "lipidmaps_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  zip.file <- file.path(base.loc, "lipidmaps.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  zip::unzip(zipfile = normalizePath(zip.file), exdir = normalizePath(base.loc))
  version <- Sys.Date()
  sdf.path <- list.files(base.loc, pattern = "sdf", full.names = TRUE, recursive = TRUE)
  desc <- function(sdfset) {
    mat <- NULL
    db <- data.table::as.data.table(ChemmineR::datablock2ma(datablocklist = ChemmineR::datablock(sdfset)))
    mat <- data.table::data.table(identifier = as.character(db$LM_ID), compoundname = sapply(1:nrow(db), function(i) {
      if (!is.empty(db$NAME[i])) {
        db$NAME[i]
      } else {
        db$SYSTEMATIC_NAME[i]
      }
    }), baseformula = as.character(db$FORMULA), structure = sapply(1:nrow(db), function(i) {
      if (!is.empty(db$SMILES[i])) {
        db$SMILES[i]
      } else {
        db$INCHI[i]
      }
    }), description = sapply(1:nrow(db), function(i) {
      string <- ""
      db$SYSTEMATIC_NAME
      if (!is.empty(db$ABBREVIATION[i])) {
        string <- paste0(string, "Often abbreviated as ", db$ABBREVIATION[i], ".")
      }
      if (!is.empty(db$NAME[i])) {
        string <- paste(string, "Systematic name:", db$SYSTEMATIC_NAME[i])
      }
      if (is.empty(string)) {
        string <- "Unknown"
      }
      trimws(string)
    }))
    empty.smiles <- which(sapply(mat$structure, is.empty))
    if (length(empty.smiles) > 0) {
      print(mat$structure[empty.smiles])
      mat$structure[empty.smiles] <- pbapply::pbsapply(mat$structure[empty.smiles], function(x) {
        url <- gsubfn::fn$paste("https://cactus.nci.nih.gov/chemical/structure/$x/smiles")
        Sys.sleep(0.1)
        RCurl::getURL(url)
      })
    }
    as.matrix(mat)
  }
  sdfStream.joanna(input = sdf.path, output = file.path(base.loc, "lipidmaps_parsed.csv"), append = FALSE, fct = desc, silent = TRUE)
  db.base <- data.table::fread(file.path(base.loc, "lipidmaps_parsed.csv"), fill = TRUE, header = TRUE)
  db.base$charge <- c(NA)
  doc <- xml2::read_html("https://www.lipidmaps.org/data/classification/LM_classification_exp.php")
  categories <- doc %>%
    rvest::html_nodes("div:nth-child(2)") %>%
    rvest::html_text()
  categories <- gsub(x = categories, pattern = "\n|\t", replacement = "")
  mainlist <- categories[5]
  filt_cats <- stringr::str_match_all(mainlist, pattern = "(\\[.*?) \\[")[[1]][, 2]
  tbl.rows <- pbapply::pblapply(filt_cats, function(cat) {
    data.table::data.table(catcode = stringr::str_match(cat, pattern = "\\[(.*?)\\]")[, 2], catdesc = gsub(cat, pattern = "\\[.*\\]", replacement = ""))
  })
  conv_tbl <- data.table::rbindlist(tbl.rows)
  db.base$description <- pbapply::pbsapply(1:nrow(db.base), function(i) {
    row <- db.base[i, ]
    matching <- stringi::stri_detect_fixed(row$identifier, conv_tbl$catcode, fixed = TRUE)
    paste0(row$description, ". ", "In class(es):", tolower(paste(conv_tbl$catdesc[matching], collapse = ", ")))
  })
  db.formatted <- db.base[, -1, with = FALSE]
  list(db = db.formatted, version = version)
}

#' @title Build METABOLIGHTS DB
#' @description Parses the METABOLIGHTS DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}}
#'  \code{\link[XML]{xmlToList}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[jsonlite]{read_json}}
#' @rdname build.METABOLIGHTS
#' @export
#' @importFrom utils download.file
#' @importFrom XML xmlToList
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table rbindlist
#' @importFrom jsonlite read_json
#' @examples
#' \dontrun{build.METABOLIGHTS(outfolder=tempdir(), testMode=TRUE)}
build.METABOLIGHTS <- function(outfolder, testMode = FALSE) {
  identifier <- study <- NULL
  file.url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/eb-eye_metabolights_complete.xml"
  base.loc <- file.path(outfolder, "metabolights_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  xml.file <- file.path(base.loc, "metabolights.xml")
  utils::download.file(file.url, xml.file, mode = "wb", method = "auto")
  version <- Sys.Date()
  mtbls <- XML::xmlToList(file.path(base.loc, "metabolights.xml"))
  compendium <- pbapply::pblapply(1:length(mtbls$entries), function(i) {
    study <- mtbls$entries[[i]]
    cpd.ids <- unlist(sapply(study$cross_references, function(ref) {
      if (ref[2] == "MetaboLights") {
        return(ref[1])
      }
      else {
        return(NULL)
      }
    }))
    if (length(cpd.ids) > 0) {
      data.table::data.table(identifier = cpd.ids, study = c(paste0("(", study$.attrs, "): ", study$name)))
    }
    else {
      data.table::data.table()
    }
  })
  overview <- data.table::rbindlist(compendium)
  ids <- unique(overview$identifier)
  if (testMode) {
    ids <- ids[1:10]
  }
  db.rows <- pbapply::pblapply(ids, cl = 0, FUN = function(id) {
    url <- paste0("https://www.ebi.ac.uk/metabolights/webservice/beta/compound/", id)
    res <- data.table::data.table()
    try(
      {
        info <- jsonlite::read_json(url)
        res <- data.table::data.table(compoundname = info$name, description = info$definition, baseformula = info$formula, identifier = id, charge = info$charge, structure = info$smiles)
      },
      silent = TRUE
    )
    res
  })
  db.rows.final <- pbapply::pblapply(db.rows, function(row) {
    id <- row$identifier
    if (is.null(id)) {
      return(data.table::data.table())
    }
    else {
      studies <- overview[identifier == id, study]
      study.summary <- paste0(unlist(studies), collapse = " ")
      row$description <- paste0(row$description, " Mentioned in the following studies --> ", study.summary)
      return(row)
    }
  })
  db.formatted <- data.table::rbindlist(db.rows.final, use.names = TRUE)
  list(db = db.formatted, version = version)
}

#' @title Build DIMEDB
#' @description Parses the DIMEDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[data.table]{fread}}
#'  \code{\link[reshape2]{cast}}
#'  \code{\link[Hmisc]{capitalize}}
#' @rdname build.DIMEDB
#' @export
#' @importFrom pbapply pbsapply
#' @importFrom utils download.file unzip
#' @importFrom data.table fread data.table
#' @importFrom reshape2 dcast
#' @importFrom Hmisc capitalize
#' @examples
#' \dontrun{build.DIMEDB(outfolder=tempdir(), testMode=TRUE)}
build.DIMEDB <- function(outfolder, testMode = FALSE) {
  files <- c("dimedb_pathways.zip", "dimedb_sources.zip", "dimedb_pc_info.zip", "dimedb_id_info.zip")
  file.url <- "https://dimedb.ibers.aber.ac.uk/help/downloads/"
  file.urls <- paste0(file.url, files)
  base.loc <- file.path(outfolder, "dimedb_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  pbapply::pbsapply(file.urls, function(url) {
    zip.file <- file.path(base.loc, basename(url))
    utils::download.file(url, zip.file, mode = "wb", method = "auto")
    utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  })
  version <- Sys.Date()
  atom <- data.table::fread(file.path(base.loc, "dimedb_pc_info.tsv"))
  ids <- data.table::fread(file.path(base.loc, "dimedb_id_info.tsv"))
  source <- data.table::fread(file.path(base.loc, "dimedb_sources.tsv"))
  pathway <- data.table::fread(file.path(base.loc, "dimedb_pathways.tsv"))
  unique.inchi <- unique(ids$InChIKey)
  joined <- rbind(ids, atom)
  casted <- reshape2::dcast(joined, InChIKey ~ Property, function(vec) paste0(vec, collapse = ","))
  db.formatted <- data.table::data.table(compoundname = Hmisc::capitalize(tolower(casted$Name)), description = do.call(paste0, c(casted[, c("IUPAC Name", "Synonym")], col = "-")), baseformula = casted$`Molecular Formula`, identifier = casted$InChIKey, charge = casted$`Formal Charge`, structure = casted$SMILES)
  db.formatted[which(db.formatted$description == "-")] <- c("Unknown")
  db.formatted$description <- gsub(db.formatted$description, pattern = "^-|-$", replacement = "")
  rmv <- which(db.formatted$baseformula == "Unknown")
  db.formatted <- db.formatted[-rmv, ]
  list(db = db.formatted, version = version)
}

#' @title Build VMH
#' @description Parses the VMH DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[httr]{GET}},\code{\link[httr]{content_type}},\code{\link[httr]{content}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname build.VMH
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom pbapply pblapply
#' @importFrom httr GET accept content
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist data.table
#' @examples
#' \dontrun{build.VMH(outfolder=tempdir(), testMode=TRUE)}
build.VMH <- function(outfolder, testMode = FALSE) {
  api_url <- "https://vmh.uni.lu/_api/metabolites/"
  pagerange <- 150
  theurl <- "https://www.vmh.life/files/release-notes/release.html"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "(\\w+ 20\\d\\d)")[, 2]
  if (testMode) {
    pagerange <- 1
  }
  table_list <- pbapply::pblapply(1:pagerange, function(i) {
    tbl <- NA
    try({
      url <- paste0("http://vmh.uni.lu/_api/metabolites/?page=", i)
      r <- httr::GET(url, httr::accept(".json"))
      lst <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
      tbl <- lst[[4]]
      Sys.sleep(0.1)
    })
    tbl
  })
  table_main <- data.table::rbindlist(table_list[!is.na(table_list)])
  db.formatted <- data.table::data.table(compoundname = table_main$fullName, description = table_main$description, baseformula = table_main$chargedFormula, identifier = table_main$abbreviation, charge = table_main$charge, structure = table_main$smile, isHuman = table_main$isHuman, isMicrobe = table_main$isMicrobe)
  missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "" | is.na(db.formatted$description))
  replacements <- table_main$synonyms
  db.formatted$description[missing.desc] <- replacements[missing.desc]
  missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "" | is.na(db.formatted$description))
  db.formatted$description[missing.desc] <- c("Unknown")
  descriptions <- sapply(1:nrow(db.formatted), function(i) {
    row <- db.formatted[i, ]
    if (row$isHuman & row$isMicrobe) {
      suffix <- "Found in humans and microbes."
    }
    else if (row$isHuman & !row$isMicrobe) {
      suffix <- "Found in humans."
    }
    else {
      suffix <- "Found in microbes"
    }
    if (length(row$description) > 1) {
      if (substring(row$description, nchar(row$description)) == ".") {
        paste0(row$description, " -- ", suffix, " -- ")
      }
      else {
        paste0(row$description, ". -- ", suffix, " -- ")
      }
    }
    else {
      paste0(row$description, " -- ", suffix, " -- ")
    }
  })
  db.formatted$description <- descriptions
  db.formatted <- db.formatted[, -c("isHuman", "isMicrobe")]
  list(db = db.formatted, version = version)
}

#' @title Build PHENOL EXPLORER DB
#' @description Parses the PHENOL EXPLORER DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[openxlsx]{read.xlsx}}
#'  \code{\link[data.table]{fread}}
#' @rdname build.PHENOLEXPLORER
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom utils download.file unzip
#' @importFrom openxlsx read.xlsx
#' @importFrom data.table fread data.table
#' @examples
#' \dontrun{build.PHENOLEXPLORER(outfolder=tempdir(), testMode=TRUE)}
build.PHENOLEXPLORER <- function(outfolder, testMode = FALSE) {
  . <- description <- NULL
  file.urls <- paste0("http://phenol-explorer.eu/system/downloads/current/", c("composition-data.xlsx.zip", "compounds.csv.zip", "compounds-structures.csv.zip", "metabolites.csv.zip", "metabolites-structures.csv.zip"))
  theurl <- "http://phenol-explorer.eu/"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Welcome to Phenol-Explorer (\\d\\.\\d)")[, 2]
  base.loc <- file.path(outfolder, "phenolexplorer_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  for (url in file.urls) {
    zip.file <- file.path(base.loc, basename(url))
    utils::download.file(url, destfile = zip.file, mode = "wb", method = "auto")
    utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  }
  compo_tbl <- openxlsx::read.xlsx(file.path(base.loc, "composition-data.xlsx"), sheet = 1)
  compo_tbl$description <- paste0("Present in ", tolower(compo_tbl$food), "(", tolower(compo_tbl$food_group), ", ", tolower(compo_tbl$food_sub_group), "). ", "Belongs to the compound class of ", tolower(compo_tbl$compound_group), " (", tolower(compo_tbl$compound_sub_group), "). ", "PMIDS: ", compo_tbl$pubmed_ids)
  db.base <- unique(compo_tbl[, c("compound", "description")])
  struct_tbl <- data.table::fread(file.path(base.loc, "compounds-structures.csv"))
  struct_tbl_keep <- struct_tbl[, c("id", "smiles", "name", "formula")]
  merged.cpds <- merge(db.base, struct_tbl_keep, by.x = "compound", by.y = "name")
  met_tbl <- data.table::fread(file.path(base.loc, "metabolites.csv"))
  met_struct <- data.table::fread(file.path(base.loc, "metabolites-structures.csv"))
  met_struct_keep <- met_struct[, c("id", "smiles", "formula", "name")]
  missing <- (!(met_tbl$name %in% compo_tbl$compound))
  mis_mets <- met_tbl[missing, ]
  merged.mets <- merge(mis_mets, met_struct_keep, by.x = "name", by.y = "name")
  merged.mets <- merged.mets[, c("name", "synonyms", "id.x", "formula.x", "smiles")]
  colnames(merged.mets) <- c("compound", "description", "id", "formula", "smiles")
  merged.mets$description <- paste0("Metabolite of compound in food.", merged.mets$description)
  merged <- unique(rbind(merged.cpds, merged.mets))
  db.formatted <- data.table::data.table(compoundname = merged$compound, description = merged$description, baseformula = merged$formula, identifier = merged$id, charge = c(NA), structure = merged$smiles)
  db.formatted <- db.formatted[, .(description = paste(description, collapse = ". ")), by = c("compoundname", "baseformula", "identifier", "charge", "structure")]
  list(db = db.formatted, version = version)
}

#' @title Build MASSBANK DB
#' @description Parses MASSBANK, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname build.MASSBANK
#' @export
#' @importFrom stringr str_match
#' @importFrom utils download.file unzip
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table rbindlist
#' @examples
#' \dontrun{build.MASSBANK(outfolder=tempdir(), testMode=TRUE)}
build.MASSBANK <- function(outfolder, testMode = FALSE) {
  baseformula <- NULL
  theurl <- "https://massbank.eu/MassBank/"
  header <- paste0(readLines(theurl), collapse = " ")
  version <- stringr::str_match(header, pattern = "Update (.*?):")[, 2]
  file.url <- "https://github.com/MassBank/MassBank-data/archive/2020.06.zip"
  base.loc <- file.path(outfolder, "massbank_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "massbank.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  cpd_files <- list.files(base.loc, pattern = ".txt$", full.names = TRUE, recursive = TRUE)
  if (testMode) {
    cpd_files <- cpd_files[1:10]
  }
  db_rows <- pbapply::pblapply(cpd_files, function(fn) {
    row <- NA
    try({
      lines <- readLines(fn)
      split.lines <- sapply(lines, strsplit, ": ")
      names(split.lines) <- sapply(split.lines, function(x) x[1])
      split.lines <- lapply(split.lines, function(x) x[2:length(x)])
      row <- data.table::data.table(compoundname = split.lines$`CH$NAME`, description = split.lines$RECORD_TITLE, baseformula = split.lines$`CH$FORMULA`, identifier = split.lines$ACCESSION, charge = NA, structure = {
        struct <- "N/A"
        try({
          struct <- split.lines$`CH$SMILES`
        })
        struct
      })
      if (row$structure[[1]] == "N/A") {
        row$structure <- split.lines$`CH$INCHI`
      }
    })
    row
  })
  db.formatted <- data.table::rbindlist(db_rows[!is.na(db_rows)], fill = TRUE)
  db.formatted <- db.formatted[!is.na(baseformula), ]
  db.formatted <- aggregate(db.formatted, by = list(db.formatted$structure), FUN = function(x) paste0(unique(x), collapse = "/"))
  db.formatted <- db.formatted[, -1]
  list(db = db.formatted, version = version)
}

#' @title Build BMDB
#' @description Parses the BMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[pbapply]{pboptions}}
#'  \code{\link[base]{connections}}
#'  \code{\link[XML]{xmlValue}},\code{\link[XML]{xmlEventParse}}
#' @rdname build.BMDB
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match str_match_all
#' @importFrom utils download.file unzip
#' @importFrom pbapply startpb setpb
#' @importFrom base file
#' @importFrom XML xmlValue xmlEventParse
#' @examples
#' \dontrun{build.BMDB(outfolder=tempdir(), testMode=TRUE)}
build.BMDB <- function(outfolder, testMode = FALSE) {
  oldpar <- options()
  options(stringsAsFactors = FALSE)
  on.exit(options(oldpar))
  theurl <- "http://www.bovinedb.ca/about"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "BMDB Version <strong>(\\d.\\d)")[, 2]
  file.url <- "https://www.bovinedb.ca/system/downloads/current/bmdb_metabolites.zip"
  base.loc <- file.path(outfolder, "bmdb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "BMDB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc)
  input <- file.path(base.loc, "bmdb_metabolites.xml")
  header <- readLines(input, n = 10)
  version <- trimws(gsub(grep(pattern = "<version", header, value = TRUE), pattern = "<\\/?version>", replacement = ""))
  date <- trimws(gsub(grep(pattern = "update_date", header, value = TRUE), pattern = "<\\/?update_date>", replacement = ""))
  theurl <- RCurl::getURL("https://bovinedb.ca/metabolites")
  n <- as.numeric(stringr::str_match_all(theurl, "of <.+?>(.*?)<.+?>")[[1]][1, 2])
  if (testMode) {
    n <- 10
  }
  envir <- environment()
  envir$db.formatted <- data.frame(compoundname = rep(NA, n), baseformula = rep(NA, n), identifier = rep(NA, n), structure = rep(NA, n), charge = rep(NA, n), description = rep("", n))
  envir$idx <- 1
  envir$maxn <- n
  envir$pb <- pbapply::startpb(min = 1, max = envir$maxn)
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") {
      os <- "osx"
    }
  }
  else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)) {
      os <- "osx"
    }
    if (grepl("linux-gnu", R.version$os)) {
      os <- "linux"
    }
  }
  if (tolower(os) == "windows") {
    acc <- "primary"
    nm <- "primary"
    desc <- "primary"
    con <- base::file(input, "r")
    while (TRUE) {
      line <- readLines(con, n = 1, skipNul = TRUE)
      if (length(line) == 0 | envir$idx == envir$maxn) {
        break
      }
      if (line == "</metabolite>") {
        envir$idx <- envir$idx + 1
        pbapply::setpb(envir$pb, envir$idx)
        acc <- "primary"
        nm <- "primary"
        desc <- "primary"
      }
      tag <- stringr::str_match(line, pattern = "<(.*?)>")[, 2]
      switch(tag, accession = {
        if (acc == "primary") {
          envir$db.formatted[envir$idx, ]$identifier <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
          acc <- "secondary"
        }
      }, name = {
        if (nm == "primary") {
          envir$db.formatted[envir$idx, ]$compoundname <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
          nm <- "secondary"
        }
      }, smiles = {
        envir$db.formatted[envir$idx, ]$structure <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
      }, description = {
        if (desc == "primary") {
          envir$db.formatted[envir$idx, ]$description <- paste0(envir$db.formatted[envir$idx, ]$description, " HMDB: ", trimws(gsub(line, pattern = "(<.*?>)", replacement = "")))
          desc <- "secondary"
        }
      }, cs_description = {
        envir$db.formatted[envir$idx, ]$description <- paste0(envir$db.formatted[envir$idx, ]$description, "From ChemSpider: ", trimws(gsub(line, pattern = "(<.*?>)", replacement = "")))
      }, chemical_formula = {
        envir$db.formatted[envir$idx, ]$baseformula <- trimws(gsub(line, pattern = "(<.*?>)", replacement = ""))
      })
    }
    close(con)
  }
  else {
    metabolite <- function(currNode, currEnvir = envir) {
      if (currEnvir$idx %% 1000 == 0) {
        pbapply::setpb(currEnvir$pb, currEnvir$idx)
      }
      if (currEnvir$idx == currEnvir$maxn) {
        return()
      }
      currEnvir$db.formatted[currEnvir$idx, "compoundname"] <- XML::xmlValue(currNode[["name"]])
      currEnvir$db.formatted[currEnvir$idx, "identifier"] <- XML::xmlValue(currNode[["accession"]])
      currEnvir$db.formatted[currEnvir$idx, "baseformula"] <- XML::xmlValue(currNode[["chemical_formula"]])
      currEnvir$db.formatted[currEnvir$idx, "structure"] <- XML::xmlValue(currNode[["smiles"]])
      currEnvir$db.formatted[currEnvir$idx, "description"] <- paste(XML::xmlValue(currNode[["description"]]))
      x <- currNode[["predicted_properties"]]
      properties <- currNode[["predicted_properties"]]
      currEnvir$db.formatted[currEnvir$idx, "charge"] <- stringr::str_match(XML::xmlValue(properties), pattern = "formal_charge([+|\\-]\\d*|\\d*)")[, 2]
      currEnvir$idx <- currEnvir$idx + 1
    }
    XML::xmlEventParse(input, branches = list(metabolite = metabolite), replaceEntities = TRUE)
  }
  list(db = envir$db.formatted, version = version)
}

#' @title Build RMDB
#' @description Parses the RMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}}
#'  \code{\link[stringr]{str_split}},\code{\link[stringr]{str_extract}},\code{\link[stringr]{str_match}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname build.RMDB
#' @export
#' @importFrom utils download.file
#' @importFrom stringr str_split str_extract str_match
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table rbindlist
#' @examples
#' \dontrun{build.RMDB(outfolder=tempdir(), testMode=TRUE)}
build.RMDB <- function(outfolder, testMode = FALSE) {
  file.url <- "http://www.rumendb.ca/public/downloads/current/metabocards.gz"
  base.loc <- file.path(outfolder, "rmdb_source")
  version <- Sys.Date()
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  gz.file <- file.path(base.loc, "rmdb.gz")
  utils::download.file(file.url, gz.file, mode = "wb", cacheOK = TRUE, method = "auto")
  zz <- gzfile(gz.file, "rt")
  dat <- readLines(zz)
  dat.pasted <- paste0(dat, collapse = ";")
  n <- sum(grepl("#BEGIN_METABOCARD", dat))
  split <- stringr::str_split(dat.pasted, pattern = "END_METABOCARD")[[1]]
  if (testMode) {
    split <- split[1:10]
  }
  db.rows <- pbapply::pblapply(split, function(l) {
    data.table::data.table(identifier = stringr::str_extract(string = l, pattern = "RMDB\\d+"), compoundname = stringr::str_match(string = l, pattern = "name:;(.*?);;")[, 2], baseformula = stringr::str_match(string = l, pattern = "chemical_formula:;(.*?);;")[, 2], description = paste0(stringr::str_match(string = l, pattern = "description:;(.*?);")[, 2], " Found in ", tolower(stringr::str_match(string = l, pattern = "biofluid_location:;(.*?);")[, 2]), "."), structure = stringr::str_match(
      string = l,
      pattern = "smiles_canonical:;(.*?);;"
    )[, 2])
  })
  db.formatted <- data.table::rbindlist(db.rows)
  db.formatted$charge <- c(0)
  list(db = db.formatted, version = version)
}

#' @title Build ECMDB
#' @description Parses the ECMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[RJSONIO]{fromJSON}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname build.ECMDB
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom utils download.file unzip
#' @importFrom RJSONIO fromJSON
#' @importFrom data.table rbindlist data.table
#' @examples
#' \dontrun{build.ECMDB(outfolder=tempdir(), testMode=TRUE)}
build.ECMDB <- function(outfolder, testMode = FALSE) {
  theurl <- "http://ecmdb.ca/downloads"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Version <strong>(\\d.\\d)")[, 2]
  file.url <- "http://ecmdb.ca/download/ecmdb.json.zip"
  base.loc <- file.path(outfolder, "ecmdb_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  zip.file <- file.path(base.loc, "ecmdb.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc)
  json <- file.path(base.loc, "ecmdb.json")
  json.rows <- RJSONIO::fromJSON(json)
  if (testMode) {
    json.rows <- json.rows[1:10]
  }
  db.base <- data.table::rbindlist(json.rows)
  db.formatted <- data.table::data.table(compoundname = db.base$name, description = db.base$description, baseformula = db.base$moldb_formula, identifier = db.base$met_id, charge = db.base$moldb_formal_charge, structure = db.base$moldb_smiles)
  list(db = db.formatted, version = version)
}

#' @title LMDB
#' @description LMDB database, included with permission from the DB creators.
"lmdb"

#' @title Adduct table
#' @description Table with all adducts included by default in MetaDBparse.
"adducts"

#' @title Adduct rule table
#' @description Table with all adduct rules included by default in MetaDBparse.
"adduct_rules"

#' @title Build LMDB
#' @description Parses the LMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.LMDB
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @examples
#' \dontrun{build.LMDB(outfolder=tempdir())}
build.LMDB <- function(outfolder) {
  lmdb <- NULL
  theurl <- "http://lmdb.ca/"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Version <strong>(\\d.\\d)")[, 2]
  data(lmdb, package = "MetaDBparse", envir = environment())
  db.formatted <- lmdb
  list(db = db.formatted, version = version)
}

#' @title Build YMDB
#' @description Parses the YMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[data.table]{as.data.table}},\code{\link[data.table]{fread}},\code{\link[data.table]{rbindlist}}
#'  \code{\link[ChemmineR]{datablock2ma}},\code{\link[ChemmineR]{datablock}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.YMDB
#' @export
#' @importFrom utils download.file unzip
#' @importFrom jsonlite fromJSON
#' @importFrom data.table data.table as.data.table fread rbindlist
#' @importFrom ChemmineR datablock2ma datablock
#' @importFrom pbapply pblapply
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @examples
#' \dontrun{build.YMDB(outfolder=tempdir())}
build.YMDB <- function(outfolder) {
  file.url <- "http://www.ymdb.ca/system/downloads/current/ymdb.json.zip"
  base.loc <- file.path(outfolder, "ymdb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "ymdb.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc)
  json <- file.path(base.loc, "ymdb.json")
  line <- readLines(json)[[1]]
  fixed_lines <- gsub(line, pattern = "\\]\\}\\{", replacement = "]},{")
  jsonParsed <- jsonlite::fromJSON(fixed_lines)
  db.partial <- data.table::data.table(compoundname = jsonParsed$name, description = paste0(jsonParsed$description, " Synonyms: ", paste0(jsonParsed$synonyms[[1]], collapse = ",")), baseformula = "", identifier = jsonParsed$ymdb_id, charge = jsonParsed$physiological_charge, structure = "")
  sdf.url <- "http://www.ymdb.ca/system/downloads/current/ymdb.sdf.zip"
  zip.file <- file.path(base.loc, "ymdb_sdf.zip")
  utils::download.file(sdf.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc)
  sdf.path <- list.files(base.loc, pattern = "sdf$", full.names = TRUE, recursive = TRUE)
  desc <- function(sdfset) {
    mat <- NULL
    db <- data.table::as.data.table(ChemmineR::datablock2ma(datablocklist = ChemmineR::datablock(sdfset)))
    info <- data.table::data.table(identifier = db$DATABASE_ID, compoundname = db$GENERIC_NAME, structure = db$SMILES, baseformula = db$JCHEM_FORMULA, description = paste0("Synonyms: ", db$SYNONYMS))
    info
  }
  out.csv <- file.path(base.loc, "ymdb_parsed.csv")
  if (file.exists(out.csv)) {
    file.remove(out.csv)
  }
  sdfStream.joanna(input = sdf.path, output = out.csv, append = FALSE, fct = desc, silent = TRUE)
  db.struct <- data.table::fread(file.path(base.loc, "ymdb_parsed.csv"), fill = TRUE, header = TRUE)
  db.merged <- merge(db.partial, db.struct, by = "identifier", all.y = TRUE)
  db.rows <- pbapply::pblapply(1:nrow(db.merged), function(i) {
    row <- db.merged[i, ]
    data.table::data.table(identifier = row$identifier, compoundname = row$compoundname.y, structure = row$structure.y, baseformula = row$baseformula.y, description = if (is.na(row$description.x)) {
      row$description.y
    } else {
      row$description.x
    })
  })
  db.formatted <- data.table::rbindlist(db.rows)
  theurl <- "http://www.ymdb.ca/about"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Version <strong>(\\d.\\d)")[, 2]
  list(db = db.formatted, version = version)
}

#' @title Build PAMDB
#' @description Parses the PAMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}}
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[readxl]{read_excel}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.PAMDB
#' @export
#' @importFrom utils download.file
#' @importFrom data.table as.data.table data.table
#' @importFrom readxl read_excel
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @examples
#' \dontrun{build.PAMDB(outfolder=tempdir(), testMode=TRUE)}
build.PAMDB <- function(outfolder, testMode = FALSE) {
  file.url <- "http://pseudomonas.umaryland.edu/PaDl/PaMet.xlsx"
  base.loc <- file.path(outfolder, "pamdb_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  xlsx.file <- file.path(base.loc, "pamdb.xlsx")
  utils::download.file(file.url, xlsx.file, mode = "wb", cacheOK = TRUE, method = "auto")
  db.base <- data.table::as.data.table(readxl::read_excel(xlsx.file, sheet = 1))
  db.formatted <- data.table::data.table(identifier = db.base$MetID, compoundname = db.base$Name, structure = db.base$SMILES, baseformula = c(NA), description = gsub(gsub(db.base$Reactions, pattern = "\\r", replacement = ", "), pattern = "\\n", replacement = ""), charge = db.base$Charge)
  theurl <- "http://pseudomonas.umaryland.edu/"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Version  <STRONG>(\\d.\\d)")[, 2]
  list(db = db.formatted, version = version)
}

#' @title Build mVOC db
#' @description Parses the mVOC db, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[XML]{getNodeSet}},\code{\link[XML]{xmlAttrs}},\code{\link[XML]{readHTMLTable}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.mVOC
#' @export
#' @importFrom XML xpathSApply xmlAttrs readHTMLTable
#' @importFrom pbapply pbsapply pblapply
#' @importFrom data.table data.table rbindlist
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @examples
#' \dontrun{build.mVOC(outfolder=tempdir(), testMode=TRUE)}
build.mVOC <- function(outfolder, testMode = FALSE) {
  categories <- c("\\(", "[", "$", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "H", "I", "M", "N", "O", "P", "Q", "S", "T", "U")
  hrefFun <- function(x) {
    XML::xpathSApply(x, "./a", XML::xmlAttrs)
  }
  urlbase <- "http://bioinformatics.charite.de/mvoc/"
  if (testMode) {
    categories <- c("A")
  }
  search_urls <- pbapply::pbsapply(categories, function(categ) {
    theurl <- paste0("http://bioinformatics.charite.de/mvoc/index.php?site=browse&char=", categ)
    tables <- XML::readHTMLTable(theurl, elFun = hrefFun)
    tables.nonull <- tables[sapply(tables, function(x) !is.null(x))]
    tables.wrows <- tables.nonull[sapply(tables.nonull, function(x) nrow(x) > 0)]
    urls <- gsub(unlist(tables.wrows), pattern = "\\./", replacement = urlbase)
    urls[!is.na(urls)]
  })
  db_rows <- pbapply::pblapply(unlist(search_urls), function(theurl) {
    try({
      tables <- XML::readHTMLTable(theurl)
      tables.nonull <- tables[sapply(tables, function(x) !is.null(x))]
      end.nameblock <- min(which(!is.na(tables.nonull[[1]][, 2])))
      names <- tables.nonull[[1]][1:end.nameblock - 1, 1]
      restblock <- tables.nonull[[1]][-c(1:end.nameblock), ]
      trans_restblock <- as.data.frame(t(restblock))
      colnames(trans_restblock) <- trans_restblock[1, ]
      trans_restblock <- trans_restblock[2, ]
      has.desc <- which(sapply(tables.nonull, function(x) "Biological Function" %in% colnames(x)))
      tbl_desc <- tables.nonull[[has.desc]]
      desc <- paste0(sapply(1:nrow(tbl_desc), function(i) {
        row <- tbl_desc[i, 1:2]
        paste0(row[2], "(", row[1], ")")
      }), collapse = ". ")
      db.formatted <- data.table::data.table(identifier = trans_restblock$`PubChem ID`, compoundname = names[1], structure = trans_restblock$SMILES, baseformula = trans_restblock$Formula, description = paste0("Microbes producing this compound:", desc, ". Other names:", paste0(names[-1], collapse = ",")), charge = c(0))
      Sys.sleep(1)
      db.formatted
    })
  })
  db.formatted <- data.table::rbindlist(db_rows[sapply(db_rows, function(x) !is.null(nrow(x)))])
  theurl <- "http://bioinformatics.charite.de/mvoc/index.php?site=home"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "mVOC (\\d.\\d)")[, 2]
  list(db = db.formatted, version = version)
}

#' @title Build NANPDB
#' @description Parses the NANPDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}}
#'  \code{\link[data.table]{fread}}
#' @rdname build.NANPDB
#' @export
#' @importFrom utils download.file
#' @importFrom data.table fread
#' @examples
#' \dontrun{build.NANPDB(outfolder=tempdir(), testMode=TRUE)}
build.NANPDB <- function(outfolder, testMode = FALSE) {
  . <- V2 <- V1 <- V3 <- NULL
  n <- 13141
  if (testMode) {
    n <- 10
  }
  base.url <- "http://african-compounds.org/nanpdb/get_compound_card/"
  db.rows <- pbapply::pblapply(1:n, function(i) {
    url <- paste0(base.url, i, "/")
    db.formatted <- data.table::data.table()
    try({
      htmlfile <- RCurl::getURL(url, .opts = list(ssl.verifypeer = FALSE))
      cpdname <- stringr::str_match(htmlfile, "Compound Card of <b>(.*?)<\\/b>")[, 2]
      htmlfile <- gsub("\t|\n", "", htmlfile)
      htmlfile <- gsub("<\\/li>", ". ", htmlfile)
      htmlfile <- gsub(":", ": ", htmlfile)
      tbl <- xml2::read_html(htmlfile) %>%
        rvest::html_node("table") %>%
        rvest::html_table(fill = TRUE)
      tbl <- t(tbl)
      colnames(tbl) <- tbl[1, ]
      tbl <- data.table::as.data.table(tbl)
      db.formatted <- data.table::data.table(identifier = i, compoundname = cpdname, structure = tbl$`SMILES:`[2], baseformula = tbl$`Molecular Formula:`[2], description = paste0(tbl$`Source Species Information`[2], tbl$`Predicted toxicity from pkCSM`[2]), charge = c(0))
    })
    db.formatted
  })
  db.formatted <- data.table::rbindlist(db.rows, fill = TRUE)
  version <- Sys.Date()
  list(db = db.formatted, version = version)
}

#' @title Build STOFF db
#' @description Parses the STOFF db, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[readxl]{read_excel}}
#' @rdname build.STOFF
#' @export
#' @importFrom utils download.file unzip
#' @importFrom data.table as.data.table
#' @importFrom readxl read_excel
#' @examples
#' \dontrun{build.STOFF(outfolder=tempdir())}
build.STOFF <- function(outfolder) {
  . <- Name <- Formula <- SMILES <- `Additional Names` <- Index <- compoundname <- NULL
  file.url <- "http://www.lfu.bayern.de/stoffident/stoffident-static-content/html/download/SI_Content.zip"
  base.loc <- file.path(outfolder, "stoff_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  zip.file <- file.path(base.loc, "stoff.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  excel.file <- file.path(base.loc, "SI_Content.xls")
  xlsx.file <- gsub(excel.file, pattern = "xls", replacement = "xlsx")
  file.copy(excel.file, xlsx.file)
  db.base <- data.table::as.data.table(readxl::read_excel(xlsx.file, sheet = 1))
  db.base.aggr <- db.base[, .(compoundname = unique(Name), baseformula = unique(Formula), structure = unique(SMILES), description = paste0("Synonyms: ", paste(`Additional Names`, collapse = ","))), by = Index]
  db.base.aggr$charge <- c(0)
  db.formatted <- db.base.aggr[!is.na(compoundname)]
  colnames(db.formatted)[1] <- "identifier"
  version <- Sys.Date()
  list(db = db.formatted, version = version)
}

#' @title Build REACTOME db
#' @description Parses the REACTOME db, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[rvest]{html_nodes}}
#'  \code{\link[xml2]{read_html}}
#' @rdname build.REACTOME
#' @export
#' @importFrom rvest html_nodes
#' @importFrom xml2 read_html
#' @importFrom data.table as.data.table fread
#' @examples
#' \dontrun{build.REACTOME(outfolder=tempdir())}
build.REACTOME <- function(outfolder) {
  . <- description <- organism <- identifier <- NULL
  theurl <- "https://reactome.org/download/current/ChEBI2Reactome.txt"
  chebiExists <- length(list.files(outfolder, pattern = "chebi\\.db")) > 0
  if (!chebiExists) {
    print("Requires ChEBI. Building...")
    buildBaseDB(outfolder, "chebi")
  }
  chebi <- data.table::as.data.table(showAllBase(outfolder, "chebi"))
  chebi$identifier <- as.numeric(chebi$identifier)
  base.db <- data.table::fread(theurl)
  colnames(base.db) <- c("identifier", "pathway_id", "pathway_url", "description", "annotation_type", "organism")
  base.db <- base.db[, c("identifier", "description", "organism")]
  base.db <- unique(base.db[, .(description = paste0(unique(description), collapse = ", "), organism = paste0(unique(organism), collapse = ", ")), by = c("identifier")])
  base.db$description <- paste0("Found in pathways: ", base.db$description, ".", " Found in organisms: ", base.db$organism, ".")
  base.db <- unique(base.db[, c("identifier", "description")])
  base.db$identifier <- as.character(base.db$identifier)
  chebi <- chebi[, -"description", with = FALSE]
  base.db$identifier <- as.numeric(base.db$identifier)
  db.formatted <- merge(chebi, base.db)
  colnames(db.formatted)[colnames(db.formatted) == "formula"] <- "baseformula"
  db.formatted <- db.formatted[, c("identifier", "compoundname", "structure", "baseformula", "description", "charge")]
  url <- "https://reactome.org/"
  ver <- as.character(url %>% xml2::read_html() %>% rvest::html_nodes(xpath = "//*[@id=\"fav-portfolio1\"]/div/h3/text()"))
  ver <- gsub(" released on ", " - ", ver)
  ver <- gsub("Version ", "", ver)
  return(list(db = db.formatted, version = ver))
}
