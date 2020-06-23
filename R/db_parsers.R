# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# should return data table!
build.MCDB <- function(outfolder){ # WORKS

  options(stringsAsFactors = F)
  file.url <- "http://www.mcdb.ca/system/downloads/current/milk_metabolites.zip"
  base.loc <- file.path(outfolder, "mcdb_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "MCDB.zip")
  utils::download.file(file.url, zip.file,mode = "wb",cacheOK = T)
  utils::unzip(zip.file, exdir = base.loc)

  input = file.path(base.loc, "milk_metabolites.xml")
  header = readLines(input,n = 10)
  version = trimws(gsub(grep(pattern = "<version", header, value = T),
                        pattern = "<\\/?version>",
                        replacement = ""))
  date = trimws(gsub(grep(pattern = "update_date", header, value = T),
                     pattern = "<\\/?update_date>",
                     replacement = ""))

  theurl <- RCurl::getURL("http://www.mcdb.ca/statistics",.opts = list(ssl.verifypeer = FALSE) )
  tables <- XML::readHTMLTable(theurl)
  stats = data.table::rbindlist(tables)
  n = as.numeric(as.character(gsub(x = stats[Description == "Total Metabolites"]$Count,
                                   pattern = ",",
                                   replacement="")))

  db.formatted <- data.frame(
    compoundname = rep(NA, n),
    baseformula = rep(NA, n),
    identifier = rep(NA, n),
    structure = rep(NA, n),
    charge = rep(NA, n),
    description = rep("", n)
  )

  idx = 1 # which metabolite are we on
  pb <- pbapply::startpb(min = idx, max = n)

  # FOR WINDOWS
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }

  if(tolower(os) == "windows"){
    acc = "primary"
    nm = "primary"
    desc = "primary"
    con = base::file(input, "r")
    while (TRUE) {

      line = readLines(con, n = 1,skipNul = T)

      if (length(line) == 0){
        break
      }

      if(line == "</metabolite>"){
        m = m+1
        pbapply::setpb(pb, idx)
        acc = "primary"
        nm = "primary"
        desc = "primary"
      }

      tag = stringr::str_match(line, pattern = "<(.*?)>")[,2]

      switch(tag,
             accession = {
               if(acc == "primary"){
                 db.formatted[idx,]$identifier <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
                 acc <- "secondary"
               }
             },
             name = {
               if(nm == "primary"){
                 db.formatted[idx,]$compoundname <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
                 nm = "secondary"
               }
             },
             smiles = {
               db.formatted[idx,]$structure <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
             },
             description = {
               if(desc == "primary"){
                 db.formatted[idx,]$description <- paste0(db.formatted[idx,]$description,
                                                          " HMDB: ",
                                                          trimws(gsub(line, pattern = "(<.*?>)", replacement="")))
                 desc = "secondary"
               }
             },
             cs_description = {
               db.formatted[idx,]$description <- paste0(db.formatted[idx,]$description,
                                                        "From ChemSpider: ",
                                                        trimws(gsub(line, pattern = "(<.*?>)", replacement="")))
             },
             chemical_formula = {
               db.formatted[idx,]$baseformula <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
             })
    }
    close(con)
  }else{
    metabolite = function(currNode){
      if(idx %% 1000 == 0){
        pbapply::setpb(pb, idx)
      }

      currNode <<- currNode

      db.formatted[idx, "compoundname"] <<- XML::xmlValue(currNode[['name']])
      db.formatted[idx, "identifier"] <<- XML::xmlValue(currNode[['accession']])
      db.formatted[idx, "baseformula"] <<- XML::xmlValue(currNode[['chemical_formula']])
      db.formatted[idx, "structure"] <<- XML::xmlValue(currNode[['smiles']])
      db.formatted[idx, "description"] <<- paste(XML::xmlValue(currNode[['description']])
      )
      x <- currNode[['predicted_properties']]
      properties <- currNode[['predicted_properties']]
      db.formatted[idx, "charge"] <<- stringr::str_match(XML::xmlValue(properties),
                                                         pattern = "formal_charge([+|\\-]\\d*|\\d*)")[,2]

      idx <<- idx + 1
    }

    XML::xmlEventParse(input,
                       branches = list(metabolite = metabolite),
                       replaceEntities=T)
  }
  list(db = db.formatted, version = version)
}

# should return data table!
build.HMDB <- function(outfolder){ # WORKS

  options(stringsAsFactors = F)

  file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
  base.loc <- file.path(outfolder, "hmdb_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "HMDB.zip")
  utils::download.file(file.url, zip.file,mode = "wb",cacheOK = T)
  utils::unzip(zip.file, exdir = base.loc)
  input = file.path(base.loc, "hmdb_metabolites.xml")

  header = readLines(input,n = 10)
  version = trimws(gsub(grep(pattern = "<version", header, value = T),
                     pattern = "<\\/?version>",
                     replacement = ""))
  date = trimws(gsub(grep(pattern = "update_date", header, value = T),
                        pattern = "<\\/?update_date>",
                        replacement = ""))

  theurl <- RCurl::getURL("http://www.hmdb.ca/statistics",.opts = list(ssl.verifypeer = FALSE) )
  tables <- XML::readHTMLTable(theurl)
  stats = data.table::rbindlist(tables)
  n = as.numeric(as.character(gsub(x = stats[Description == "Total Number of Metabolites"]$Count,
                                   pattern = ",",
                                   replacement="")))

  db.formatted <- data.frame(
    compoundname = rep(NA, n),
    baseformula = rep(NA, n),
    identifier = rep(NA, n),
    structure = rep(NA, n),
    charge = rep(NA, n),
    description = rep("", n)
  )

  idx = 1 # which metabolite are we on
  pb <- pbapply::startpb(min = idx, max = n)

  # FOR WINDOWS
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }

  if(tolower(os) == "windows"){
    acc = "primary"
    nm = "primary"
    desc = "primary"
    con = base::file(input, "r")
    while (TRUE) {

      line = readLines(con, n = 1,skipNul = T)

      if (length(line) == 0){
        break
      }

      if(line == "</metabolite>"){
        m = m+1
        pbapply::setpb(pb, idx)
        acc = "primary"
        nm = "primary"
        desc = "primary"
      }

      tag = stringr::str_match(line, pattern = "<(.*?)>")[,2]

      switch(tag,
             accession = {
               if(acc == "primary"){
                 db.formatted[idx,]$identifier <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
                 acc <- "secondary"
               }
             },
             name = {
               if(nm == "primary"){
                 db.formatted[idx,]$compoundname <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
                 nm = "secondary"
               }
             },
             smiles = {
               db.formatted[idx,]$structure <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
             },
             description = {
               if(desc == "primary"){
                 db.formatted[idx,]$description <- paste0(db.formatted[idx,]$description,
                                                          " HMDB: ",
                                                          trimws(gsub(line, pattern = "(<.*?>)", replacement="")))
                 desc = "secondary"
               }
             },
             cs_description = {
               db.formatted[idx,]$description <- paste0(db.formatted[idx,]$description,
                                                        "From ChemSpider: ",
                                                        trimws(gsub(line, pattern = "(<.*?>)", replacement="")))
             },
             chemical_formula = {
               db.formatted[idx,]$baseformula <- trimws(gsub(line, pattern = "(<.*?>)", replacement=""))
             })
    }
    close(con)
  }else{
    metabolite = function(currNode){

      if(idx %% 1000 == 0){
        pbapply::setpb(pb, idx)
      }

      currNode <<- currNode

      db.formatted[idx, "compoundname"] <<- XML::xmlValue(currNode[['name']])
      db.formatted[idx, "identifier"] <<- XML::xmlValue(currNode[['accession']])
      db.formatted[idx, "baseformula"] <<- XML::xmlValue(currNode[['chemical_formula']])
      db.formatted[idx, "structure"] <<- XML::xmlValue(currNode[['smiles']])
      db.formatted[idx, "description"] <<- paste("HMDB:",
                                                 XML::xmlValue(currNode[['cs_description']]),
                                                 "CHEMSPIDER:",
                                                 XML::xmlValue(currNode[['description']])
      )
      x <- currNode[['predicted_properties']]
      properties <- currNode[['predicted_properties']]
      db.formatted[idx, "charge"] <<- stringr::str_match(XML::xmlValue(properties),
                                                         pattern = "formal_charge([+|\\-]\\d*|\\d*)")[,2]

      idx <<- idx + 1
    }

    XML::xmlEventParse(input,
                       branches = list(metabolite = metabolite),
                       replaceEntities=T)
  }
  list(db = db.formatted, version = version)
}

build.METACYC <- function(outfolder){ # WORKS
  # May need to remake smartTable if anything on the website changes unfortunately
  # TODO: download file directly from link, will need a javascript. Maybe Rselenium??

  base.loc <- file.path(outfolder, "metacyc_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  source.file = file.path(base.loc, "All_compounds_of_MetaCyc.txt")
  if(!file.exists(source.file)){
    message("Please download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3787684059' as 'All_compounds_of_MetaCyc.txt' and save in the outfolder/metacyc_source folder.")
    return(NULL)
  }

  theurl <- RCurl::getURL("https://metacyc.org/group?id=biocyc17-31223-3787684059",.opts = list(ssl.verifypeer = FALSE))
  header = stringr::str_match(theurl, pattern = "Created: (.*), Last Modified: (.*), Readable")
  date = header[,2]
  version = header[,3]

  #tables <- XML::readHTMLTable(theurl)
  #stats = data.table::rbindlist(tables)

  metacyc.raw = data.table::fread(source.file,
                                  fill=TRUE,
                                  quote="",
                                  sep ="\t")

  colnames(metacyc.raw) <- gsub(x=as.character(colnames(metacyc.raw)), pattern = '\\"', replacement="")
  metacyc.raw[] <- lapply(metacyc.raw, gsub, pattern = '\\"', replacement = "")

  compounds <- pbapply::pbsapply(metacyc.raw$`Common-Name`, FUN=function(pw){
    pw <- iconv(pw, "latin1", "UTF-8",sub='')
    pw <- pw[pw != " // "]
    pw <- gsub(pw, pattern = "&", replacement="")
    pw <- gsub(pw, pattern = ";", replacement="")
    res <- gsub(pw, pattern = "<(i|\\/i)>", replacement = "",perl = T)
    res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "",perl = T)
    paste0(res, collapse=" --- ")
  })

  db.formatted <- data.table::data.table(compoundname = compounds,
                                         description = sapply(1:nrow(metacyc.raw), function(i){
                                           species = metacyc.raw$Species[i]
                                           paste0(metacyc.raw$Summary[i], if(species != "") paste0(" Found in: ", species, ".") else "")
                                         }),
                                         baseformula = c(NA),
                                         identifier = metacyc.raw$`Object ID`,
                                         charge = c(0),
                                         structure = metacyc.raw$SMILES)
  list(db = db.formatted, version = version)
}

build.CHEBI <- function(outfolder){ # WORKS
  db.full <- {
    release = "latest"
    woAssociations = FALSE
    chebi_download <- tempdir()
    url = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/"
    filenames = RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)

    releases <- strsplit(filenames, split="\n")[[1]]
    releases <- as.numeric(gsub(x = releases, pattern = "rel", replacement = ""))
    message("Validating ChEBI release number ... ", appendLF = FALSE)

    if (release == "latest") {
      release <- max(releases)
    } else {
      release <- releases[match(release, releases)]
    }

    version = release

    message("OK")
    ftp <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel",
                  release, "/Flat_file_tab_delimited/")
    message("Downloading compounds ... ", appendLF = FALSE)
    utils::download.file(paste0(ftp, "compounds.tsv.gz"), paste0(chebi_download,
                                                                 "compounds.tsv"),
                         quiet = TRUE,
                         mode="wb")
    compounds <- as.data.frame.array(read.delim2(paste0(chebi_download,
                                                        "compounds.tsv")))
    message("DONE", appendLF = TRUE)
    message("Downloading synonyms ... ", appendLF = FALSE)
    utils::download.file(paste0(ftp, "names.tsv.gz"), paste0(chebi_download,
                                                             "names.tsv"), quiet = TRUE, mode="wb")
    names <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                     "names.tsv"))))
    message("DONE", appendLF = TRUE)
    message("Downloading formulas ... ", appendLF = FALSE)
    utils::download.file(paste0(ftp, "chemical_data.tsv"), paste0(chebi_download,
                                                                  "formulas.tsv"), quiet = TRUE, mode="wb")
    formulas <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                        "formulas.tsv"))))

    message("Downloading structures ... ", appendLF = FALSE)

    utils::download.file(paste0(ftp, "structures.csv.gz"), paste0(chebi_download,
                                                                  "structures.csv"), quiet = TRUE, mode="wb")

    message("DONE", appendLF = TRUE)

    structures <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                          "structures.csv"), sep=",")))
    # mol.rows <- which(structures$TYPE == "mol")
    # inchi.rows <- which(structures$TYPE == "InChI")
    # inchikey.rows <- which(structures$TYPE == "InChIKey")

    smile.rows <- which(structures$TYPE == "SMILES")

    structures <- structures[smile.rows,]

    message("DONE", appendLF = TRUE)
    message("Building ChEBI ... ", appendLF = TRUE)

    #compounds <- compounds[compounds[, "STAR"] >= 3, ]
    latest <- compounds[, c("ID", "NAME", "DEFINITION")]
    old <- compounds[, c("ID", "PARENT_ID")]
    old <- merge(x = old, y = latest, by.x = "PARENT_ID", by.y = "ID")
    compounds <- rbind(latest, old[, c("ID", "NAME", "DEFINITION")])
    compounds[compounds[, "NAME"] == "null", "NAME"] <- NA
    compounds <- compounds[complete.cases(compounds), ]
    DB <- suppressWarnings((merge(compounds[, c("ID", "NAME", "DEFINITION")],
                                  names[, c("COMPOUND_ID", "SOURCE", "NAME")], by.x = "ID",
                                  by.y = "COMPOUND_ID", all.x = TRUE)))
    ChEBI <- unique(DB[, c("ID", "NAME.x", "DEFINITION")])
    colnames(ChEBI) <- c("ID", "ChEBI", "DEFINITION")

    message(" KEGG Associations ... ", appendLF = FALSE)
    KEGG <- unique(DB[DB[, "SOURCE"] == "KEGG COMPOUND", c("ID",
                                                           "NAME.y")])
    KEGG <- KEGG[complete.cases(KEGG), ]
    colnames(KEGG) <- c("ID", "KEGG")
    message("DONE", appendLF = TRUE)
    message(" IUPAC Associations ... ", appendLF = FALSE)
    IUPAC <- unique(DB[DB[, "SOURCE"] == "IUPAC", c("ID", "NAME.y")])
    IUPAC <- IUPAC[complete.cases(IUPAC), ]
    colnames(IUPAC) <- c("ID", "IUPAC")
    message("DONE", appendLF = TRUE)
    message(" MetaCyc Associations ... ", appendLF = FALSE)
    MetaCyc <- unique(DB[DB[, "SOURCE"] == "MetaCyc", c("ID",
                                                        "NAME.y")])
    MetaCyc <- MetaCyc[complete.cases(MetaCyc), ]
    colnames(MetaCyc) <- c("ID", "MetaCyc")
    message("DONE", appendLF = TRUE)
    message(" ChEMBL Associations ... ", appendLF = FALSE)
    ChEMBL <- unique(DB[DB[, "SOURCE"] == "ChEMBL", c("ID",
                                                      "NAME.y")])
    ChEMBL <- ChEMBL[complete.cases(ChEMBL), ]
    colnames(ChEMBL) <- c("ID", "ChEMBL")
    message("DONE", appendLF = TRUE)
    DB <- unique(merge(DB["ID"], ChEBI, by = "ID", all.x = TRUE))
    DB <- unique(merge(DB, KEGG, by = "ID", all.x = TRUE))
    DB <- unique(merge(DB, IUPAC, by = "ID", all.x = TRUE))
    DB <- unique(merge(DB, MetaCyc, by = "ID", all.x = TRUE))
    DB <- unique(merge(DB, ChEMBL, by = "ID", all.x = TRUE))
    rm(ChEBI, ChEMBL, compounds, IUPAC, KEGG, latest, MetaCyc,
       names, old)
    if ("FORMULA" %in% unique(formulas[, "TYPE"])) {
      message(" Formula Associations ... ", appendLF = FALSE)
      formula <- formulas[formulas[, "TYPE"] == "FORMULA",
                          c("COMPOUND_ID", "CHEMICAL_DATA")]
      colnames(formula) <- c("ID", "FORMULA")
      DB <- merge(DB, formula, by = "ID", all.x = TRUE)
      DB <- merge(DB, DB[, c("ChEBI", "FORMULA")], by = "ChEBI",
                  all.x = TRUE)
      DB[is.na(DB[, "FORMULA.x"]), "FORMULA.x"] <- "null"
      DB[is.na(DB[, "FORMULA.y"]), "FORMULA.y"] <- "null"
      DB[DB[, "FORMULA.x"] != "null" & DB[, "FORMULA.y"] ==
           "null", "FORMULA.y"] <- DB[DB[, "FORMULA.x"] !=
                                        "null" & DB[, "FORMULA.y"] == "null", "FORMULA.x"]
      DB[DB[, "FORMULA.y"] != "null" & DB[, "FORMULA.x"] ==
           "null", "FORMULA.x"] <- DB[DB[, "FORMULA.y"] !=
                                        "null" & DB[, "FORMULA.x"] == "null", "FORMULA.y"]
      DB <- unique(DB[DB[, "FORMULA.x"] != "null" & DB[, "FORMULA.y"] !=
                        "null", c("ID", "DEFINITION","ChEBI", "KEGG", "IUPAC", "MetaCyc",
                                  "ChEMBL", "FORMULA.x")])
      rm(formula)
      message("DONE", appendLF = TRUE)
    }
    else {
      message("NOT AVAILABLE FOR THIS RELEASE")
    }
    message("Downloading molecular weights ... ", appendLF = FALSE)
    if ("MASS" %in% unique(formulas[, "TYPE"])) {
      mass <- formulas[formulas[, "TYPE"] == "MASS", c("COMPOUND_ID",
                                                       "CHEMICAL_DATA")]
      colnames(mass) <- c("ID", "MASS")
      DB <- merge(DB, mass, by = "ID", all.x = TRUE)
      DB <- merge(DB, DB[, c("ChEBI", "MASS")], by = "ChEBI",
                  all.x = TRUE)
      DB[is.na(DB[, "MASS.x"]), "MASS.x"] <- "null"
      DB[is.na(DB[, "MASS.y"]), "MASS.y"] <- "null"
      DB[DB[, "MASS.x"] != "null" & DB[, "MASS.y"] == "null",
         "MASS.y"] <- DB[DB[, "MASS.x"] != "null" & DB[,
                                                       "MASS.y"] == "null", "MASS.x"]
      DB[DB[, "MASS.y"] != "null" & DB[, "MASS.x"] == "null",
         "MASS.x"] <- DB[DB[, "MASS.y"] != "null" & DB[,
                                                       "MASS.x"] == "null", "MASS.y"]
      DB <- unique(DB[, c("ID", "DEFINITION", "ChEBI", "KEGG", "IUPAC",
                          "MetaCyc", "ChEMBL", "FORMULA.x", "MASS.x")])
      rm(mass)
      message("DONE", appendLF = TRUE)
    }
    else {
      message("NOT AVAILABLE FOR THIS RELEASE")
    }
    message("Downloading monoisotopic molecular weights ... ",
            appendLF = FALSE)
    if ("MONOISOTOPIC MASS" %in% unique(formulas[, "TYPE"])) {
      mmass <- formulas[formulas[, "TYPE"] == "MONOISOTOPIC MASS",
                        c("COMPOUND_ID", "CHEMICAL_DATA")]
      colnames(mmass) <- c("ID", "MONOISOTOPIC")
      DB <- merge(DB, mmass, by = "ID", all.x = TRUE)
      DB <- merge(DB, DB[, c("ChEBI", "MONOISOTOPIC")], by = "ChEBI",
                  all.x = TRUE)
      DB[is.na(DB[, "MONOISOTOPIC.x"]), "MONOISOTOPIC.x"] <- "null"
      DB[is.na(DB[, "MONOISOTOPIC.y"]), "MONOISOTOPIC.y"] <- "null"
      DB[DB[, "MONOISOTOPIC.x"] != "null" & DB[, "MONOISOTOPIC.y"] ==
           "null", "MONOISOTOPIC.y"] <- DB[DB[, "MONOISOTOPIC.x"] !=
                                             "null" & DB[, "MONOISOTOPIC.y"] == "null", "MONOISOTOPIC.x"]
      DB[DB[, "MONOISOTOPIC.y"] != "null" & DB[, "MONOISOTOPIC.x"] ==
           "null", "MONOISOTOPIC.x"] <- DB[DB[, "MONOISOTOPIC.y"] !=
                                             "null" & DB[, "MONOISOTOPIC.x"] == "null", "MONOISOTOPIC.y"]
      DB <- unique(DB[, c("ID","DEFINITION", "ChEBI", "KEGG", "IUPAC",
                          "MetaCyc", "ChEMBL", "FORMULA.x", "MASS.x", "MONOISOTOPIC.x")])
      rm(mmass)
      message("DONE", appendLF = TRUE)
    }
    else {
      message("NOT AVAILABLE FOR THIS RELEASE")
    }
    message("Downloading molecular charges ... ", appendLF = FALSE)
    if ("CHARGE" %in% unique(formulas[, "TYPE"])) {
      charge <- formulas[formulas[, "TYPE"] == "CHARGE", c("COMPOUND_ID",
                                                           "CHEMICAL_DATA")]
      colnames(charge) <- c("ID", "CHARGE")
      DB <- merge(DB, charge, by = "ID", all.x = TRUE)
      DB <- merge(DB, DB[, c("ChEBI", "CHARGE")], by = "ChEBI",
                  all.x = TRUE)
      DB[is.na(DB[, "CHARGE.x"]), "CHARGE.x"] <- "null"
      DB[is.na(DB[, "CHARGE.y"]), "CHARGE.y"] <- "null"
      DB[DB[, "CHARGE.x"] != "null" & DB[, "CHARGE.y"] ==
           "null", "CHARGE.y"] <- DB[DB[, "CHARGE.x"] != "null" &
                                       DB[, "CHARGE.y"] == "null", "CHARGE.x"]
      DB[DB[, "CHARGE.y"] != "null" & DB[, "CHARGE.x"] ==
           "null", "CHARGE.x"] <- DB[DB[, "CHARGE.y"] != "null" &
                                       DB[, "CHARGE.x"] == "null", "CHARGE.y"]
      DB <- unique(DB[, c("ID", "DEFINITION", "ChEBI", "KEGG", "IUPAC",
                          "MetaCyc", "ChEMBL", "FORMULA.x", "MASS.x", "MONOISOTOPIC.x",
                          "CHARGE.x")])
      message("DONE", appendLF = TRUE)
    }
    else {
      message("NOT AVAILABLE FOR THIS RELEASE")
    }
    DB[DB == "null"] <- NA

    DB = merge(DB, structures[,-1], by.x = "ID",  by.y = "COMPOUND_ID", all.x = TRUE, incomparables = "unknown")

    DB <- unique(DB[complete.cases(DB[, c("ID", "DEFINITION","ChEBI", "FORMULA.x",
                                          "MASS.x", "MONOISOTOPIC.x", "CHARGE.x", "STRUCTURE")]), ])
    colnames(DB) <- c("ID", "DEFINITION","ChEBI", "KEGG", "IUPAC", "MetaCyc",
                      "ChEMBL", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE", "STRUCTURE")

    if (woAssociations == TRUE) {
      compounds <- unique(rbind(setNames(DB[, c("ChEBI", "DEFINITION","FORMULA",
                                                "MASS", "MONOISOTOPIC", "CHARGE")], c("NAME","DEFINITION","FORMULA",
                                                                                      "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")), setNames(DB[,
                                                                                                                                                  c("KEGG","DEFINITION", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")],
                                                                                                                                               c("NAME", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")),
                                setNames(DB[, c("IUPAC", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC",
                                                "CHARGE")], c("NAME", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC",
                                                              "CHARGE")), setNames(DB[, c("MetaCyc", "DEFINITION","FORMULA",
                                                                                          "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")], c("NAME", "DEFINITION",
                                                                                                                                            "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")),
                                setNames(DB[, c("ChEMBL","DEFINITION", "FORMULA", "MASS", "MONOISOTOPIC",
                                                "CHARGE","STRUCTURE")], c("NAME", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC",
                                                                          "CHARGE","STRUCTURE"))))
      compounds <- compounds[complete.cases(compounds), ]
      compounds
    }
    else {
      DB
    }} # uses an altered version from minval package
  db.full <- data.table::as.data.table(db.full)
  db.formatted <- unique(db.full[, list(compoundname = ChEBI,
                                        description = DEFINITION,
                                        baseformula = FORMULA,
                                        identifier = ID,
                                        charge = gsub(CHARGE,pattern = "$\\+\\d", replacement = ""),
                                        structure = toupper(STRUCTURE)
  )])
  list(db = db.formatted, version = version)
}

build.FOODB <- function(outfolder){ # WORKS
  file.url <- "http://foodb.ca/public/system/downloads/foodb_2017_06_29_csv.tar.gz"
  base.loc <- file.path(outfolder, "foodb_source")
  if(dir.exists(base.loc))(unlink(base.loc,recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "foodb.tar.gz")
  utils::download.file(file.url, zip.file, mode = 'wb', method = 'libcurl')
  utils::untar(normalizePath(zip.file), exdir = normalizePath(base.loc))

  theurl <- RCurl::getURL("http://foodb.ca/downloads",.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(theurl, pattern = "FooDB Version <strong>(...)<")[,2]
  date = stringr::str_match(theurl, pattern = "FooDB CSV file<\\/td><td>(.{3,40})<\\/td>")[,2]
  date = as.Date(date, format = "%B%e%Y")
  subdate = gsub(date, pattern = "-", replacement = "_")

  base.table <- data.table::fread(file = file.path(base.loc, gsubfn::fn$paste("foodb_$subdate_csv"), "compounds.csv"))

  db.formatted <- data.table::data.table(compoundname = base.table$name,
                                         description = base.table$description,
                                         baseformula = base.table$moldb_formula,
                                         identifier= base.table$id,
                                         charge= base.table$charge,
                                         structure= base.table$moldb_smiles)

  db.formatted <- unique(db.formatted)

  list(db = db.formatted, version = version)

}

build.WIKIDATA <- function(outfolder){ # WORKS

  date = Sys.Date()
  version = Sys.Date()

  sparql_query <- 'PREFIX wd: <http://www.wikidata.org/entity/>
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
                           SERVICE wikibase:label { bd:serviceParam wikibase:language "en, de". }
                           ?chemical_compound wdt:P31 wd:Q11173.
                           ?chemical_compound wdt:P274 ?chemical_formula.
                           ?chemical_compound wdt:P233 ?canonical_SMILES.
                           OPTIONAL {?chemical_compound wdt:P2868 ?role.}}'

  db.1 <- WikidataQueryServiceR::query_wikidata(sparql_query,
                                                format = "simple")

  db.1$chemical_compound <- basename(db.1$chemical_compound)
  db.1$chemical_compoundDescription[db.1$chemical_compoundDescription == "chemical compound"] <- NA
  db.1$roleLabel[db.1$roleLabel == ""] <- NA

  # https://spark.apache.org/ for speed increases? is it useful locally? more an HPC thing?

  # NOTE: tool - myFAIR, SOAR, EUDAT, JUNIPER (ML, can use R, jupyter notebooks?), GALAXY'S GUI/API CAN BE CHANGED??

  db.2 <-  data.table::as.data.table(aggregate(db.1, by=list(db.1$chemical_compoundLabel, db.1$chemical_formula), function(x) c(unique((x)))))

  db.2$description = apply(db.2[,c("roleLabel", "chemical_compoundDescription")], 1, FUN=function(x){
    x <- unlist(x)
    paste(x[!is.na(x)], collapse=", ")
  })

  db.formatted <- data.table::data.table(compoundname = as.character(db.2$chemical_compoundLabel),
                                         description = as.character(db.2$description),
                                         baseformula = as.character(db.2$chemical_formula),
                                         identifier= as.character(db.2$chemical_compound),
                                         charge= c(),
                                         structure = as.character(db.2$canonical_SMILES))

  from = "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089"
  to = "0123456789"
  db.formatted$baseformula <- chartr(from, to, db.formatted$baseformula)

  db.formatted$description[db.formatted$description == ""] <- "Unknown"
  db.formatted$baseformula <- as.character(db.formatted$baseformula)
  db.formatted <- db.formatted[!is.na(db.formatted$baseformula),]

  # - - write - -

  list(db = db.formatted, version = version)

}

# RE-ENABLE AFTER TALKING TO EGON
# build.WIKIPATHWAYS <- function(outfolder){
#   chebi.loc <- file.path(outfolder, "chebi.db")
#   # ---------------------------------------------------
#   chebi <- SPARQL::SPARQL(url="http://sparql.wikipathways.org/",
#                           query='prefix wp:      <http://vocabularies.wikipathways.org/wp#>
#                                                    prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
#                                                    prefix dcterms: <http://purl.org/dc/terms/>
#                                                    prefix xsd:     <http://www.w3.org/2001/XMLSchema#>
#                                                    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
#
#                                                    select  ?mb
#                                                    (group_concat(distinct str(?labelLit);separator=", ") as ?label )
#                                                    ?idurl as ?csid
#                                                    (group_concat(distinct ?pwTitle;separator=", ") as ?description)
#                                                    ?pathway
#                                                    where {
#                                                    ?mb a wp:Metabolite ;
#                                                    rdfs:label ?labelLit ;
#                                                    wp:bdbChEBI ?idurl ;
#                                                    dcterms:isPartOf ?pathway .
#                                                    ?pathway a wp:Pathway ;
#                                                    dc:title ?pwTitle .
#                                                    FILTER (BOUND(?idurl))
#                                                    }
#                                                    GROUP BY ?mb ?wp ?idurl ?pathway')
#   chebi.ids <- gsub(chebi$results$csid, pattern = ".*:|>", replacement = "")
#   conn.chebi <- RSQLite::dbConnect(RSQLite::SQLite(), chebi.loc)
#   chebi.join.table <- data.table::data.table(identifier = chebi.ids,
#                                              description = chebi$results$description,
#                                              widentifier = chebi$results$mb,
#                                              pathway = chebi$results$pathway)
#   RSQLite::dbWriteTable(conn.chebi, "wikipathways", chebi.join.table, overwrite=TRUE)
#   db.formatted <- RSQLite::dbGetQuery(conn.chebi, "SELECT DISTINCT  b.compoundname,
#                                                                b.identifier,
#                                                                w.description,
#                                                                b.baseformula,
#                                                                b.charge,
#                                                                b.structure
#                                                                FROM base b
#                                                                JOIN wikipathways w
#                                                                ON b.identifier = w.identifier")
#
#   db.formatted$identifier <- gsub(db.formatted$identifier, pattern = "<|>", replacement = "")
#   RSQLite::dbRemoveTable(conn.chebi, "wikipathways")
#   RSQLite::dbDisconnect(conn.chebi)
#   db.formatted
# }

build.RESPECT <- function(outfolder){ # WORKS
  file.url <- "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"

  base.loc <- file.path(outfolder, "respect_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "respect.zip")
  utils::download.file(file.url, zip.file,mode = "wb")
  utils::unzip(normalizePath(zip.file), exdir = (base.loc))

  theurl <- RCurl::getURL("http://spectra.psc.riken.jp/menta.cgi/respect/download/download",.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(theurl, pattern = "<b>(.{1,30}) Update")[,2]
  version = gsub(version, pattern="\\.", replacement = "-")
  date = version

  cpd_files <- list.files(base.loc,
                          full.names = T)

  db_rows <- pbapply::pblapply(cpd_files, function(fn){
    row = data.table::data.table()
    try({
      lines <- readLines(fn,skipNul = T, n = 100)
      split.lines <- sapply(lines, strsplit, ": ")
      names(split.lines) <- sapply(split.lines, function(x) x[1])
      split.lines <- lapply(split.lines, function(x) x[2:length(x)])
      row <- data.table::data.table(
        compoundname = split.lines$`CH$NAME`,
        description = split.lines$RECORD_TITLE,
        baseformula = split.lines$`CH$FORMULA`,
        identifier = split.lines$ACCESSION,
        charge = c(0),
        structure = split.lines$`CH$SMILES`
      )
      if(row$structure == "N/A") row$structure <- split.lines$`CH$INCHI`
    })
    row
  })

  db.formatted <- data.table::rbindlist(db_rows[!is.na(db_rows)],fill = T)
  db.formatted <- unique(db.formatted[!is.na(baseformula),])
  db.formatted <- aggregate(db.formatted, by = list(db.formatted$compoundname), FUN = function(x) paste0(unique(x), collapse="/"))
  db.formatted <- db.formatted[,-1]

  # - - -

  list(db = db.formatted, version = version)

}

build.MACONDA <- function(outfolder, conn){ # NEEDS SPECIAL FUNCTIONALITY

  file.url = "https://www.maconda.bham.ac.uk/downloads/MaConDa__v1_0__csv.zip"

  base.loc <- file.path(outfolder, "maconda_source")

  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "maconda.zip")

  theurl <- RCurl::getURL("https://www.maconda.bham.ac.uk/downloads.php",.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(theurl, pattern = "Version (...)")[,2]
  version = gsub(version, pattern="\\.", replacement = "-")
  date = Sys.Date()

  utils::download.file(file.url, zip.file,mode = "wb",extra = "-k",method = "curl")

  utils::unzip(normalizePath(zip.file),files = "MaConDa__v1_0__extensive.csv",exdir = normalizePath(base.loc))

  base.table <- data.table::fread(file = file.path(base.loc, "MaConDa__v1_0__extensive.csv"))

  mysterious = which(base.table$name == "unknown")

  base.table$formula[mysterious] <- paste0("IDK", 1:length(mysterious))

  has.inchi <- which(base.table$std_inchi != "")
  inchis <- base.table$std_inchi[has.inchi]
  smiles = pbapply::pbsapply(inchis, webchem::cs_inchi_smiles)

  charges <- gsub(base.table$ion_form, pattern = ".*\\]", replacement = "")
  no.info <- which(charges == "")
  charges[no.info] <- sapply(base.table$ion_mode[no.info], function(x) ifelse(x=="POS", 1, -1))
  plus.minus <- which(charges == "+" | charges == "-")
  charges[plus.minus] <- sapply(charges[plus.minus], function(ch) switch(ch, "+" = 1, "-" = -1))

  # - - can stop here, am gonna re-calculate charges etc! - -

  db.base <- data.table::data.table(compoundname = base.table$name,
                                    description = paste(base.table$type_of_contaminant,
                                                        base.table$ion_source_type,
                                                        base.table$ion_mode),
                                    baseformula = base.table$formula,
                                    identifier=base.table$id,
                                    charge=charges,
                                    structure=c(NA))

  db.base$structure[has.inchi] <- smiles

  db.final <- cleanDB(db.base, cl=0, silent=F, blocksize=400)

  #write to base db
  writeDB(conn, data.table::data.table(date = Sys.Date(),
                                       version = version),
          "metadata")
  writeDB(conn, db.final, "base")
  RSQLite::dbExecute(conn, "CREATE INDEX b_idx1 ON base(structure)")
  DBI::dbDisconnect(conn)

  # WRITE TO EXTENDED DB TOO
  full.db <- file.path(outfolder, "extended.db")
  first.db <- !file.exists(full.db)

  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA foreign_keys = ON"))
  RSQLite::dbExecute(full.conn, "CREATE TABLE IF NOT EXISTS structures(
                                 struct_id INT PRIMARY KEY,
                                 smiles TEXT,
                                 UNIQUE(struct_id, smiles))")

  RSQLite::dbExecute(full.conn, strwrap("CREATE TABLE IF NOT EXISTS extended(
                                         struct_id INT,
                                         fullmz decimal(30,13),
                                         adduct text,
                                         isoprevalence float)", width=10000, simplify=TRUE))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA auto_vacuum = 1;"))

  # B
  base.db <- normalizePath(file.path(outfolder, "maconda.db"))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))

  if(first.db){
    RSQLite::dbExecute(full.conn, "CREATE INDEX st_idx1 ON structures(smiles)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx1 on extended(struct_id)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx2 on extended(fullmz)")
    RSQLite::dbExecute(full.conn, "PRAGMA journal_mode=WAL;")
  }

  adducts_maconda = unique(base.table[,c("ion_form","ion_mode")])
  adducts_maconda = adducts_maconda[adducts_maconda$ion_form != ""]
  adducts_maconda$ion_mode <- sapply(adducts_maconda$ion_mode, function(mode) switch(mode, POS="positive", NEG="negative"))

  adduct_table_maconda <- data.table::data.table(
    Name = adducts_maconda$ion_form,
    Ion_mode = adducts_maconda$ion_mode,
    Charge = c(NA),
    xM = c(NA),
    AddAt = c(NA),
    RemAt = c(NA),
    AddEx = c(NA),
    RemEx = c(NA),
    Nelec = c(NA),
    Rule = c(NA),
    Info = c("MACONDA ONLY")
  )

  # custom adduct table???
  if(!first.db){
    new_adducts <- setdiff(adduct_table_maconda$Name,
                           RSQLite::dbGetQuery(full.conn,
                                               "SELECT DISTINCT Name FROM adducts")[,1])

    if(length(new_adducts)==0){
      print("maconda is already in here!")
      return(NA)
    }
    RSQLite::dbWriteTable(full.conn, "adducts", adduct_table_maconda[Name %in% new_adducts], append=T)
  }else{
    RSQLite::dbWriteTable(full.conn, "adducts", adduct_table_maconda, overwrite=T)
  }

  if(first.db | length(new_adducts) > 0){
    to.do = RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                            FROM tmp.base")
    to.do$struct_id <- c(NA)
  }else{
    to.do = RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                            FROM tmp.base LEFT JOIN structures str
                                            ON base.structure = str.smiles
                                            WHERE str.smiles IS NULL")
  }
  if(nrow(to.do) == 0){
    print("all already done")
  }

  done.structures <- if(first.db) 0 else RSQLite::dbGetQuery(full.conn, "SELECT MAX(struct_id) FROM structures")[,1]
  start.id = done.structures + 1
  structs = unique(db.final[,c(structure)])
  # new will be written

  mapper = data.table::data.table(struct_id = seq(start.id, start.id + nrow(to.do) - 1, 1),
                                  smiles = to.do$structure)
  # - - - - - -
  base.table$structure <- db.final$structure[match(base.table$id, table = db.final$identifier)]

  adduct.unknown <- which(base.table$ion_form == "")
  base.table$exact_adduct_mass[adduct.unknown] <- base.table$mz[adduct.unknown]
  base.table$ion_form[adduct.unknown] <- sapply(base.table$ion_mode[adduct.unknown], function(ion_mode) switch(ion_mode, POS=c("[M?]+?"), NEG=c("[M?]-?")))

  meta.table <- data.table::data.table(fullmz = base.table$exact_adduct_mass,
                                       adduct = base.table$ion_form,
                                       isoprevalence = c(100),
                                       structure = base.table$structure)

  # map SMILES to smile_id
  ids <- mapper$struct_id[match(meta.table$structure,
                                mapper$smiles)]
  meta.table$struct_id <- ids

  # === RETURN ===

  maconda.extended = unique(meta.table[,-"structure"])

  # map SMILES to smile_id
  RSQLite::dbWriteTable(full.conn, "extended", maconda.extended, append=T)
  RSQLite::dbWriteTable(full.conn, "structures", mapper, append=T)

  RSQLite::dbDisconnect(full.conn)
}

build.T3DB <- function(outfolder){ # WORKS
  # t3db
  file.url <- "http://www.t3db.ca/system/downloads/current/toxins.csv.zip"
  # ----
  base.loc <- file.path(outfolder, "t3db_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "T3DB.zip")
  utils::download.file(file.url, zip.file,mode = "wb")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  theurl <- RCurl::getURL("http://www.t3db.ca/downloads",.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(theurl, pattern = "T3DB Version <strong>(...)")[,2]
  date = Sys.Date()

  db.formatted <- data.table::fread(file.path(base.loc, "toxins.csv"), fill=T)

  db.formatted <- data.table::data.table(
    compoundname = db.formatted$Name,
    baseformula = db.formatted$`Chemical Formula`,
    description = db.formatted$Description,
    charge = c(0),
    identifier = db.formatted$`T3DB ID`,
    structure = db.formatted$SMILES
  )

  list(db = db.formatted, version = version)

}

build.HSDB <- function(outfolder){ # NEEDS WORK
  file.url = "ftp://ftp.nlm.nih.gov/nlmdata/.hsdblease/hsdb.xml.20190528.zip"
  base.loc <- file.path(outfolder, "hsdb_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "HSDB.zip")
  utils::download.file(file.url, zip.file,mode = "wb")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  input = list.files(base.loc, pattern = "\\.xml",full.names = T)
  theurl <- getURL("https://toxnet.nlm.nih.gov/help/hsdbcasrn.html",.opts = list(ssl.verifypeer = FALSE) )

  n = str_count(as.character(theurl), pattern = "cgi-bin")

  db.formatted <- data.frame(
    compoundname = rep(NA, n),
    baseformula = rep(NA, n),
    identifier = rep(NA, n),
    structure = rep(NA, n),
    charge = rep(NA, n),
    description = rep(NA, n)
  )

  pb <- pbapply::startpb(min = 0, max = n)

  idx <<- 0

  parsed <- XML::xmlTreeParse(input)
  i=1

  db.formatted <- data.table::data.table(
    compoundname = xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['NameOfSubstance'][[1]])),
    baseformula = xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['mf'][[1]])),
    identifier = xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['DOCNO'][[1]]))
  )
  parsed[[1]]$children$hsdb[[i]]
  parsed[[1]]$children$hsdb[[i]]['mf'][1]
  parsed[[1]]$children$hsdb[[i]]['ocpp']
  parsed[[1]]$children$hsdb[[i]]['sy']
  parsed[[1]]$children$hsdb[[i]]['mf']

  identifier = xmlValue(parsed[[1]]$children$hsdb[[i]]['CASRegistryNumber']$CASRegistryNumber)

  cas_ids <- xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['CASRegistryNumber'][[1]]))

  smiles = pbapply::pbsapply(cas_ids, function(id) webchem::cir_query(id, representation = "smiles", resolver = NULL,
                                                                      first = FALSE)[[1]])
  # TBA
  # - - - - - - - - - - - -
  db.formatted
}

build.BLOODEXPOSOME <- function(outfolder){ # WORKS
  file.url = "https://exposome1.fiehnlab.ucdavis.edu/download/BloodExpsomeDatabase_version_1.0.xlsx"

  base.loc <- file.path(outfolder, "bloodexposome_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)
  excel.file <- file.path(base.loc, "exposome.xlsx")
  utils::download.file(file.url, excel.file, mode="wb")

  version = "1.0"
  date = Sys.Date()

  db.full <- openxlsx::read.xlsx(excel.file, sheet = 1, colNames=T, startRow = 3)

  print("getting iupac names from missing compound names...")
  new.names <- pbapply::pbsapply(db.full[which(sapply(db.full$Compound.Name, is.empty)),]$PubChem.CID, function(cid){
    try({
      url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%d/JSON",
                     as.numeric(cid))
      Sys.sleep(0.1)
      json = jsonlite::read_json(url)
      json$Record$RecordTitle
    })
  })

  db.full[which(sapply(db.full$Compound.Name, is.empty)),]$Compound.Name <- new.names

  db.formatted <- unique(data.table::data.table(compoundname = db.full$Compound.Name,
                                                description = paste0("Found in ", db.full$BloodPaperCount, " papers related to blood exposome."),
                                                baseformula = db.full$Molecular.Formula,
                                                identifier = db.full$PubChem.CID,
                                                charge = db.full$Charge,
                                                structure = db.full$CanonicalSMILES))

  list(db = db.formatted, version = version)

}

build.EXPOSOMEEXPLORER <- function(outfolder){ # WORKS
  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/biomarkers.csv.zip"

  base.loc <- file.path(outfolder, "expoexplorer_source")

  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "expoexpo_comp.zip")
  utils::download.file(file.url, zip.file,mode = "wb")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  version = stringr::str_match(RCurl::getURL("http://exposome-explorer.iarc.fr/releases",.opts = list(ssl.verifypeer = FALSE)),
                               pattern = "Release (...)")[,2]
  date = stringr::str_match(RCurl::getURL("http://exposome-explorer.iarc.fr/downloads",.opts = list(ssl.verifypeer = FALSE)),
                               pattern = "(\\d{4}-\\d{2}-\\d{2})")[,2]


  base.table <- data.table::fread(file = file.path(base.loc, "biomarkers.csv"))

  db.formatted <- data.table::data.table(compoundname = base.table$Name,
                                         description = base.table$Description,
                                         baseformula = base.table$Formula,
                                         identifier= base.table$ID,
                                         charge= c(NA),
                                         structure= base.table$SMILES)

  db.formatted <- unique(db.formatted)

  # - - use correlations to get some custom descriptions :) - -

  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/correlations.csv.zip"

  zip.file <- file.path(base.loc, "expoexpo_corr.zip")
  utils::download.file(file.url, zip.file,mode = "wb")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  corr.table <- data.table::fread(file = file.path(base.loc, "correlations.csv"))

  descriptions <- pbapply::pbsapply(1:nrow(corr.table), function(i){
    row = corr.table[i,]
    desc <- paste("Found in", R.utils::decapitalize(row$Biospecimen),
                  "of", R.utils::decapitalize(row$`Subject group`),
                  "under", R.utils::decapitalize(row$Population),
                  "in", row$Country,
                  "after taking in", R.utils::decapitalize(row$Intake),
                  paste0("(", row$`Analytical method`, ", p ", row$`Correlation p-value`, ")."))
  })

  corr.table$`Pasted` <- descriptions

  df <- corr.table[,c("Excretion ID", "Pasted")]
  aggr = aggregate( Pasted ~ `Excretion ID`, df, function(x) toString(paste(unique(x),collapse = " ")))

  final.table <- merge(db.formatted, aggr, by.x = "identifier", by.y = "Excretion ID", all.x=T)

  final.table$description <- pbapply::pbsapply(1:nrow(final.table), function(i){
    row = final.table[i,]
    a = if(!is.na(row$description) & row$description != "NA") row$description else ""
    b = row$Pasted
    paste0(a,b)
  })

  db.formatted <- final.table[,-"Pasted"]

  list(db = db.formatted, version = version)

}

build.SMPDB <- function(outfolder){ # OK I THINK
  file.url <- "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
  # ----
  base.loc <- file.path(outfolder, "smpdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)
  zip.file <- file.path(base.loc, "SMPDB.zip")
  utils::download.file(file.url, zip.file, mode="wb")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  # -------------------------------

  theurl = "http://smpdb.ca/release_notes"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  versions = stringr::str_match_all(header,
                               pattern = "Release (.{1,30}) -")[[1]]
  version = gsub(pattern = "SMPDB ", replacement = "", x = versions[nrow(versions),2])
  dates = stringr::str_match_all(header,
                                  pattern = "Release .{1,30}- (.{1,30})<\\/h")[[1]]
  date = dates[nrow(dates),2]
  date = as.Date(date, format = "%B%d,%Y")

  smpdb.paths <- list.files(path = base.loc, pattern = "\\.csv$", full.names = T)

  subtables <- pbapply::pblapply(smpdb.paths, data.table::fread)
  smpdb.tab <- data.table::rbindlist(subtables, fill = T)

  db.formatted <- data.table::data.table(compoundname = smpdb.tab$`Metabolite Name`,
                                         description = paste0(smpdb.tab$`Pathway Name`,
                                                              " (",
                                                              smpdb.tab$`Pathway Subject`,
                                                              " pathway)"),
                                         baseformula = smpdb.tab$Formula,
                                         identifier = smpdb.tab$`Metabolite ID`,
                                         structure = smpdb.tab$SMILES)


  db.formatted <- db.formatted[,.(description=paste0(unique(description),collapse=", ")),by=list(compoundname, baseformula, identifier, structure)]
  db.formatted <- db.formatted[-1,]

  list(db = db.formatted, version = version)

}

build.KEGG <- function(outfolder){ # WORKS
  batches <- split(0:2300, ceiling(seq_along(0:2300)/100))
  cpds <- pbapply::pblapply(batches, FUN=function(batch){
    names(KEGGREST::keggFind("compound", batch, "mol_weight"))
  })
  cpd.ids <- Reduce(c, cpds)
  id.batches <- split(cpd.ids, ceiling(seq_along(cpd.ids)/10))

  theurl = "https://www.genome.jp/kegg/compound/"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                                    pattern = "Last updated: (.{1,30})<")[,2]
  date = as.Date(date, format = "%B%d,%Y")
  #version = date

  # --- GET COMPOUNDS ---

  kegg.cpd.list <- pbapply::pblapply(id.batches, FUN=function(batch){
    rest.result <- KEGGREST::keggGet(batch)
    # ---------------------------
    base.list <- lapply(rest.result, FUN=function(cpd){
      cpd$NAME_FILT <- gsub(cpd$NAME, pattern = ";", replacement = "")
      data.table::data.table(compoundname = c(paste(cpd$NAME_FILT, collapse=", ")),
                             description = paste0("Involved in pathways: ",
                                                  paste0(cpd$PATHWAY, collapse = ", "),
                                                  ". More specifically: ",
                                                  paste0(cpd$MODULE, collapse = ", "),
                                                  ". Also associated with compound classes: ",
                                                  paste0(
                                                    unique(trimws(
                                                      gsub(cpd$BRITE, pattern = "\\[.*\\]|  D\\d* |\\(.*\\)|\\d*", replacement= "")
                                                    )
                                                    ), collapse = ", ")
                             ),
                             baseformula = c(cpd$FORMULA),
                             identifier = c(cpd$ENTRY),
                             charge = 0,
                             structure = NA
                             ,pathway = if("PATHWAY" %in% names(cpd)) names(cpd$PATHWAY) else{NA}
      )
    })
    res = data.table::rbindlist(base.list)
    res
  })

  # - - -
  db.formatted <- data.table::rbindlist(kegg.cpd.list)

  base.loc <- file.path(outfolder, "kegg_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  kegg.mol.paths = pbapply::pblapply(id.batches, FUN=function(batch){
    #mols = Rcpi::getMolFromKEGG(batch, parallel = 1)
    bigmol = KEGGREST::keggGet(batch, "mol")
    mols = strsplit(x = paste0("\n \n \n",bigmol), split = "\\$\\$\\$\\$\n")[[1]]
    fps = normalizePath(file.path(base.loc, paste0(stringr::str_match(mols, pattern = "<ENTRY>\ncpd:(.*)\n")[,2], ".mol")),mustWork = F)
    sapply(1:length(mols), function(i) writeLines(text = mols[[i]],
                                                  con = fps[i]))
    fps
  })

  fns = list.files(base.loc,full.names = T)

  smiles.rows = pbapply::pblapply(fns, function(fn){
    smiles=NA
    charge = 0
    try({
      iatom = rcdk::load.molecules(molfiles = fn)[[1]]
      smiles = rcdk::get.smiles(molecule = iatom)
      charge = rcdk::get.total.formal.charge(mol = iatom)
    })
    id = gsub(basename(fn), pattern = "\\.mol", replacement="")
    data.table::data.table(identifier = id, smiles = smiles, calcharge = charge)
  })

  smitable <- data.table::rbindlist(smiles.rows)

  db.merged <- merge(db.formatted, smitable, by = "identifier")

  db.formatted <- data.table::data.table(compoundname = db.merged$compoundname,
                                         description = db.merged$description,
                                         baseformula = db.merged$baseformula,
                                         identifier = db.merged$identifier,
                                         charge = db.merged$calcharge,
                                         structure = db.merged$smiles)
  # - - - - - -

  list(db = db.formatted, version = version)

}

build.DRUGBANK <- function(outfolder){  # WORKS
  base.loc <- file.path(outfolder, "drugbank_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  zip.file <- file.path(base.loc, "drugbank.zip")
  if(!file.exists(zip.file)){
    file.rename(file.path(base.loc, "drugbank_all_full_database.xml.zip"), zip.file)
  }
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  input = file.path(base.loc, "full database.xml")
  header = readLines(input, n = 10)
  hasInfo = grep(x = header, pattern = "version", value = T, perl = T)[2]
  version = stringr::str_match(string = hasInfo, pattern = "version=\"(.*)\" exported")[,2]

  # theurl = "https://www.drugbank.ca/releases/latest"
  # header = XML::readHTMLTable(RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE)))[[1]]
  # date = header$`Released on`
  # version = header$Version

  theurl <- RCurl::getURL("https://www.drugbank.ca/stats",.opts = list(ssl.verifypeer = FALSE) )
  tables <- XML::readHTMLTable(theurl,header = F)
  stats =  data.table::as.data.table(tables[[1]])
  colnames(stats) <- c("Description", "Count")
  n = as.numeric(as.character(gsub(x = stats[Description == "Total Number of Drugs"]$Count, pattern = ",", replacement="")))

  db.formatted <- data.frame(
    compoundname = rep(NA, n),
    baseformula = rep(NA, n),
    identifier = rep(NA, n),
    structure = rep(NA, n),
    charge = rep(NA, n),
    description = rep(NA, n)
  )

  pb <- pbapply::startpb(min = 0, max = n)

  idx <<- 0

  metabolite = function(currNode){

    if(idx %% 10 == 0){
      pbapply::setpb(pb, idx)
    }

    idx <<- idx + 1

    currNode <<- currNode

    properties <- currNode[['calculated-properties']]

    if(is.null(properties)){
      properties <- currNode[['experimental-properties']]
    }

    proplist <- XML::xmlToList(properties)
    if(length(proplist) == 0){
      return(NULL)
    }

    # find formula
    which.form <- which(sapply(proplist, function(x){
      if("kind" %in% names(x)){
        res = x[['kind']] == "Molecular Formula"
      }else{
        res = FALSE
      }
      res
    }))

    which.struc <- which(sapply(proplist, function(x){
      if("kind" %in% names(x)){
        res = x[['kind']] == "SMILES"
      }else{
        res = FALSE
      }
      res
    }))

    which.charge <- which(sapply(proplist, function(x){
      if("kind" %in% names(x)){
        res = x[['kind']] == "Physiological Charge"
      }else{
        res = FALSE
      }
      res
    }))

    # find structure

    if(length(which.form) == 0 & length(which.struc) == 0){
      return(NULL)
    }

    db.formatted[idx, "compoundname"] <<- XML::xmlValue(currNode[['name']])
    db.formatted[idx, "identifier"] <<- XML::xmlValue(currNode[['drugbank-id']])
    db.formatted[idx, "baseformula"] <<- proplist[[which.form]][['value']]
    db.formatted[idx, "structure"] <<- if(length(which.struc) > 0){
      proplist[[which.struc]][['value']]
    }else{
      ""
    }
    db.formatted[idx, "description"] <<- XML::xmlValue(currNode[['description']])
    db.formatted[idx, "charge"] <<- if(length(which.charge) > 0){
      proplist[[which.charge]][['value']]
    }else{
      0
    }
  }

  res = XML::xmlEventParse(file = input, branches =
                             list("drug" = metabolite, "drugbank-metabolite-id-value" = print))

  db.formatted <- db.formatted[-1,]

  # - - -

  list(db = db.formatted, version = version)

}

build.LIPIDMAPS <- function(outfolder){ # WORKS (description needs some tweaking)

  file.url = "https://www.lipidmaps.org/resources/downloads/LMSD/LMSD_20190711.sdf.zip"

  # ----
  base.loc <- file.path(outfolder, "lipidmaps_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  zip.file <- file.path(base.loc, "lipidmaps.zip")
  utils::download.file(file.url, zip.file, mode="wb")
  zip::unzip(zipfile = normalizePath(zip.file), exdir = normalizePath(base.loc))
  # -------------------------------

  version = Sys.Date()
  sdf.path <- list.files(base.loc,
                         pattern = "sdf",
                         full.names = T,
                         recursive = T)

  desc <- function(sdfset){
    mat <- NULL
    #last_sdf <<- sdfset
    db <- data.table::as.data.table(ChemmineR::datablock2ma(datablocklist=ChemmineR::datablock(sdfset)))
    #last_db <<- db
    mat = as.matrix(data.table::data.table(
      identifier = as.character(db$LM_ID),
      compoundname = sapply(1:nrow(db), function(i) if(!is.empty(db$NAME[i])) db$NAME[i] else db$SYSTEMATIC_NAME[i]),
      baseformula = as.character(db$FORMULA),
      structure = sapply(1:nrow(db), function(i) if(!is.empty(db$SMILES[i])) db$SMILES[i] else webchem::cs_convert(db$INCHI[i],
                                                                                                                   from="inchi", to="smiles")),
      description = sapply(1:nrow(db), function(i){
        string=""
        db$SYSTEMATIC_NAME
        if(!is.empty(db$ABBREVIATION[i])) string <- paste0(string, "Often abbreviated as ", db$ABBREVIATION[i],".")
        if(!is.empty(db$NAME[i])) string <- paste(string, "Systematic name:", db$SYSTEMATIC_NAME[i])
        if(is.empty(string)) string <- "Unknown"
        trimws(string)
      })
    ))
    mat
  }

  sdfStream.joanna(input=sdf.path, output=file.path(base.loc,
                                                    "lipidmaps_parsed.csv"),
                   append=FALSE,
                   fct=desc,
                   silent = T)

  db.base <- data.table::fread(file.path(base.loc, "lipidmaps_parsed.csv"), fill = T, header=T)

  db.base$charge <- c(NA)

  # - - - add classification - - -
  require(rvest)
  doc <- xml2::read_html("https://www.lipidmaps.org/data/classification/LM_classification_exp.php")
  categories = doc %>%
    rvest::html_nodes("div:nth-child(2)") %>%
    rvest::html_text()

  categories = gsub(x=categories, pattern = "\n|\t", replacement="")
  mainlist = categories[5]
  filt_cats = stringr::str_match_all(mainlist, pattern = "(\\[.*?) \\[")[[1]][,2]
  tbl.rows <- pbapply::pblapply(filt_cats, function(cat){
    data.table::data.table(catcode = stringr::str_match(cat, pattern="\\[(.*?)\\]")[,2],
                           catdesc = gsub(cat, pattern = "\\[.*\\]", replacement=""))
  })
  conv_tbl = data.table::rbindlist(tbl.rows)

  db.base$description <- pbapply::pbsapply(1:nrow(db.base), function(i){
    row = db.base[i,]
    matching = stringi::stri_detect_fixed(row$identifier,conv_tbl$catcode,fixed=T)
    paste0(row$description, ". ", "In class(es):", tolower(paste(conv_tbl$catdesc[matching], collapse=", ")))
  })

  # - - - - - - - - - - - - - - - -
  db.formatted <- db.base[,-1,with=F]

  list(db = db.formatted, version = version)

}

build.METABOLIGHTS <- function(outfolder){

  file.url = "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/eb-eye_metabolights_complete.xml"

  # ----
  base.loc <- file.path(outfolder, "metabolights_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  xml.file <- file.path(base.loc, "metabolights.xml")
  utils::download.file(file.url, xml.file, mode="wb")

  # ----
  version = Sys.Date()

  mtbls = XML::xmlToList(file.path(base.loc, "metabolights.xml"))

  compendium = pbapply::pblapply(1:length(mtbls$entries), function(i){

    study = mtbls$entries[[i]]

    cpd.ids <- unlist(sapply(study$cross_references, function(ref){
      if(ref[2] == "MetaboLights"){
        return(ref[1])
      }else{
        return(NULL)
      }
    }))

    # cpd.names <- unlist(sapply(study$additional_fields, function(field){
    #   if(".attrs" %in% names(field)){
    #     if(field[[".attrs"]] == "metabolite_name"){
    #       return(field[['text']])
    #     }
    #   }else{NULL}
    # }))
    #    if(length(cpd.ids) != length(cpd.names)) cpd.names <- cpd.names[1:length(cpd.ids)]

    if(length(cpd.ids)>0){
      data.table::data.table(identifier = cpd.ids,
                             study = c(paste0("(",study$.attrs,"): ", study$name)))
    }else{
      data.table::data.table()
    }
  })

  overview = data.table::rbindlist(compendium)

  ids = unique(overview$identifier)

  db.rows <- pbapply::pblapply(ids, cl=cl, FUN=function(id){
    url <- paste0("https://www.ebi.ac.uk/metabolights/webservice/beta/compound/", id)
    res = data.table::data.table()
    try({
      info = jsonlite::read_json(url)

      res = data.table::data.table(compoundname = info$name,
                                   description = info$definition,
                                   baseformula = info$formula,
                                   identifier =  id,
                                   charge = info$charge,
                                   structure = info$smiles
      )
    })
    res
  })

  db.rows.final <- pbapply::pblapply(db.rows, function(row){
    id = row$identifier
    if(is.null(id)){
      return(data.table::data.table())
    }else{
      studies = overview[identifier == id, study]
      study.summary = paste0(unlist(studies), collapse=" ")
      row$description <- paste0(row$description, " Mentioned in the following studies --> ", study.summary)
      return(row)
    }
  })

  db.formatted <- data.table::rbindlist(db.rows.final,
                                        use.names = T)

  list(db = db.formatted, version = version)

}

build.DIMEDB <- function(outfolder){ # WORKS
  files = c(#"structures.zip",
    "dimedb_pathways.zip",
    "dimedb_sources.zip",
    "dimedb_pc_info.zip",
    "dimedb_id_info.zip")
  file.url <- "https://dimedb.ibers.aber.ac.uk/help/downloads/"
  file.urls <- paste0(file.url, files)
  # ----
  base.loc <- file.path(outfolder, "dimedb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)
  pbapply::pbsapply(file.urls, function(url){
    zip.file <- file.path(base.loc, basename(url))
    utils::download.file(url, zip.file, mode="wb")
    utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  })
  version = Sys.Date()
  atom <- data.table::fread(file.path(base.loc,"dimedb_pc_info.tsv"))
  ids <- data.table::fread(file.path(base.loc,"dimedb_id_info.tsv"))
  source <- data.table::fread(file.path(base.loc,"dimedb_sources.tsv"))
  pathway <- data.table::fread(file.path(base.loc, "dimedb_pathways.tsv"))

  unique.inchi <- unique(ids$InChIKey)

  joined <- rbind(ids, atom)
  casted <- reshape2::dcast(joined, InChIKey ~ Property, function(vec) paste0(vec, collapse=","))

  db.formatted <- data.table::data.table(compoundname = Hmisc::capitalize(tolower(casted$Name)),
                                         description = do.call(paste0, c(casted[,c("IUPAC Name", "Synonym")], col="-")),
                                         baseformula = casted$`Molecular Formula`,
                                         identifier =  casted$InChIKey,
                                         charge = casted$`Formal Charge`,
                                         structure = casted$SMILES
                                         #pathway = c(NA)
  )
  db.formatted[which(db.formatted$description == "-")] <- c("Unknown")
  db.formatted$description <- gsub(db.formatted$description, pattern = "^-|-$", replacement = "")

  rmv = which(db.formatted$baseformula == "Unknown")

  db.formatted = db.formatted[-rmv,]

  list(db = db.formatted, version = version)
}

build.VMH <- function(outfolder){ # WORKS
  api_url <- "https://vmh.uni.lu/_api/metabolites/"

  pagerange = 150
  # get the first page

  theurl = "https://www.vmh.life/files/release-notes/release.html"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "(\\w+ 20\\d\\d)")[,2]

  table_list <- pbapply::pblapply(1:pagerange, function(i){
    tbl = NA
    try({
      url = gsubfn::fn$paste("http://vmh.uni.lu/_api/metabolites/?page=$i")
      r <- httr::GET(url, httr::accept(".json"))
      lst <- jsonlite::fromJSON(httr::content(r, "text"))
      tbl <- lst[[4]]
      Sys.sleep(.1)
    })
    # - - return - -
    tbl
  })

  table_main <- data.table::rbindlist(table_list[!is.na(table_list)])

  db.formatted <- data.table::data.table(compoundname = table_main$fullName,
                                         description = table_main$description,
                                         baseformula = table_main$chargedFormula,
                                         identifier= table_main$abbreviation,
                                         charge= table_main$charge,
                                         structure= table_main$smile,
                                         isHuman = table_main$isHuman,
                                         isMicrobe = table_main$isMicrobe)

  #filter for weights that are way too much w/ weird formulas

  missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "" | is.na(db.formatted$description))
  replacements <- table_main$synonyms # use synonum instead
  db.formatted$description[missing.desc] <- replacements[missing.desc]
  missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "" | is.na(db.formatted$description))
  db.formatted$description[missing.desc] <- c("Unknown")

  descriptions <- sapply(1:nrow(db.formatted), function(i){

    row = db.formatted[i,]

    if(row$isHuman & row$isMicrobe){
      suffix = "Found in humans and microbes."
    }else if(row$isHuman & !row$isMicrobe){
      suffix = "Found in humans."
    }else{
      suffix = "Found in microbes"
    }

    if(length(row$description) > 1){
      if(substring(row$description, nchar(row$description)) == "."){
        paste0(row$description," -- ", suffix, " -- ")
      }else{
        paste0(row$description, ". -- ", suffix, " -- ")
      }
    }else{
      paste0(row$description," -- ", suffix, " -- ")
    }
  })

  db.formatted$description <- descriptions

  db.formatted <- db.formatted[,-c("isHuman", "isMicrobe")]

  list(db = db.formatted, version = version)

}

build.PHENOLEXPLORER <- function(outfolder){ # WORKS
  file.urls = paste0(
    "http://phenol-explorer.eu/system/downloads/current/",
    c("composition-data.xlsx.zip",
      "compounds.csv.zip",
      "compounds-structures.csv.zip",
      "metabolites.csv.zip",
      "metabolites-structures.csv.zip"))

  theurl="http://phenol-explorer.eu/"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "Welcome to Phenol-Explorer (\\d\\.\\d)")[,2]

  base.loc <- file.path(outfolder, "phenolexplorer_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  for(url in file.urls){
    zip.file <- file.path(base.loc, basename(url))
    utils::download.file(url, destfile = zip.file, mode="wb")
    utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  }

  # start from the composition table
  compo_tbl <- openxlsx::read.xlsx(file.path(base.loc, "composition-data.xlsx"), sheet = 1)
  compo_tbl$description <- paste0("Present in ", tolower(compo_tbl$food),
                                  "(", tolower(compo_tbl$food_group), ", ",
                                  tolower(compo_tbl$food_sub_group), "). ",
                                  "Belongs to the compound class of ",
                                  tolower(compo_tbl$compound_group),
                                  " (", tolower(compo_tbl$compound_sub_group), "). ",
                                  "PMIDS: ", compo_tbl$pubmed_ids)

  db.base <- unique(compo_tbl[, c("compound", "description")])


  # unique structure info
  struct_tbl <- data.table::fread(file.path(base.loc,"compounds-structures.csv"))
  struct_tbl_keep <- struct_tbl[, c("id", "smiles", "name", "formula")]

  merged.cpds <- merge(db.base, struct_tbl_keep, by.x = "compound", by.y = "name")

  # metabolites
  met_tbl <- data.table::fread(file.path(base.loc,"metabolites.csv"))
  met_struct <- data.table::fread(file.path(base.loc,"metabolites-structures.csv"))
  met_struct_keep <- met_struct[, c("id", "smiles", "formula", "name")]

  missing = (!(met_tbl$name %in% compo_tbl$compound))
  mis_mets <- met_tbl[missing,]

  merged.mets <- merge(mis_mets, met_struct_keep, by.x = "name", by.y = "name")
  merged.mets <- merged.mets[, c("name", "synonyms", "id.x", "formula.x", "smiles")]
  colnames(merged.mets) <- c("compound", "description", "id", "formula", "smiles")
  merged.mets$description <- paste0("Metabolite of compound in food.", merged.mets$description)

  merged <- unique(rbind(merged.cpds, merged.mets))

  db.formatted <- data.table::data.table(compoundname = merged$compound,
                                         description = merged$description,
                                         baseformula = merged$formula,
                                         identifier= merged$id,
                                         charge= c(NA),
                                         structure= merged$smiles)

  db.formatted <- db.formatted[ , .(description = paste(description, collapse=". ")),
                                by = c("compoundname", "baseformula", "identifier", "charge", "structure")]

  list(db = db.formatted, version = version)

}

build.MASSBANK <- function(outfolder){ # WORKS

  theurl = "https://massbank.eu/MassBank/"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "Update (.* 20\\d\\d):")[,2]

  file.url <- "https://github.com/MassBank/MassBank-data/archive/master.zip"
  base.loc <- file.path(outfolder, "massbank_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "massbank.zip")
  utils::download.file(file.url, zip.file,mode = "wb", method='libcurl')
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  cpd_files <- list.files(base.loc,
                          pattern = ".txt$",
                          full.names = T,
                          recursive = T)

  # - - - - - - - - - - - - - - -

  db_rows <- pbapply::pblapply(cpd_files, function(fn){
    row = NA

    try({
      lines <- readLines(fn)
      split.lines <- sapply(lines, strsplit, ": ")
      names(split.lines) <- sapply(split.lines, function(x) x[1])
      split.lines <- lapply(split.lines, function(x) x[2:length(x)])
      row <- data.table::data.table(
        compoundname = split.lines$`CH$NAME`,
        description = split.lines$RECORD_TITLE,
        baseformula = split.lines$`CH$FORMULA`,
        identifier = split.lines$ACCESSION,
        charge = NA,
        structure = {
          struct = "N/A"
          try({
            struct = split.lines$`CH$SMILES`
          })
          struct
        }
      )
      if(row$structure[[1]] == "N/A"){
        row$structure <- split.lines$`CH$INCHI`
      }
    })
    row
  })

  db.formatted <- data.table::rbindlist(db_rows[!is.na(db_rows)], fill=TRUE)
  db.formatted <- db.formatted[!is.na(baseformula),]

  db.formatted <- aggregate(db.formatted, by = list(db.formatted$structure),
                            FUN = function(x) paste0(unique(x), collapse="/"))

  db.formatted <- db.formatted[,-1]

  list(db = db.formatted, version = version)
}

# ======= HERE WITH VERSION NUMBERS =======

build.BMDB <- function(outfolder){

  theurl = "http://bmdb.wishartlab.com/about"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "BMDB Version <strong>(\\d.\\d)")[,2]

  file.url = "http://www.cowmetdb.ca/public/downloads/current/metabocards.gz"
  base.loc <- file.path(outfolder, "bmdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  gz.file <- file.path(base.loc, "bmdb.gz")
  utils::download.file(file.url, gz.file,mode = "wb", cacheOK = T)
  zz = gzfile(gz.file, 'rt')
  dat = readLines(zz)
  dat.pasted = paste0(dat, collapse = ";")
  n = sum(grepl("END_METABOCARD",dat))
  split = stringr::str_split(dat.pasted, pattern = "END_METABOCARD")[[1]]
  db.rows = pbapply::pblapply(split, function(l){
    data.table::data.table(identifier = stringr::str_extract(string = l,
                                                             pattern = "BMDB\\d+"),
                           compoundname = stringr::str_match(string = l,
                                                             pattern = "name:;(.*?);;")[,2],
                           baseformula = stringr::str_match(string = l,
                                                            pattern = "chemical_formula:;(.*?);;")[,2],
                           description = paste0(stringr::str_match(string = l,
                                                                   pattern = "description:;(.*?);")[,2],
                                                " Found in ",
                                                tolower(stringr::str_match(string = l,
                                                                           pattern = "biofluid_location:;(.*?);")[,2]),
                                                "."),
                           structure = stringr::str_match(string = l,
                                                          pattern = "smiles_canonical:;(.*?);;")[,2]
    )
  })
  db.formatted = data.table::rbindlist(db.rows)
  db.formatted$charge <- c(0)

  list(db = db.formatted, version = version)

}

build.RMDB <- function(outfolder){
  file.url = "http://www.rumendb.ca/public/downloads/current/metabocards.gz"
  base.loc <- file.path(outfolder, "rmdb_source")
  version = Sys.Date()
  if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  gz.file <- file.path(base.loc, "rmdb.gz")
  utils::download.file(file.url, gz.file,mode = "wb", cacheOK = T)
  zz = gzfile(gz.file,'rt')
  dat = readLines(zz)
  dat.pasted = paste0(dat, collapse = ";")
  n = sum(grepl("#BEGIN_METABOCARD",dat))
  split = stringr::str_split(dat.pasted, pattern = "END_METABOCARD")[[1]]
  db.rows = pbapply::pblapply(split, function(l){
    data.table::data.table(identifier = stringr::str_extract(string = l,
                                                             pattern = "RMDB\\d+"),
                           compoundname = stringr::str_match(string = l,
                                                             pattern = "name:;(.*?);;")[,2],
                           baseformula = stringr::str_match(string = l,
                                                            pattern = "chemical_formula:;(.*?);;")[,2],
                           description = paste0(stringr::str_match(string = l,
                                                                   pattern = "description:;(.*?);")[,2],
                                                " Found in ",
                                                tolower(stringr::str_match(string = l,
                                                                           pattern = "biofluid_location:;(.*?);")[,2]),
                                                "."),
                           structure = stringr::str_match(string = l,
                                                          pattern = "smiles_canonical:;(.*?);;")[,2]
    )
  })
  db.formatted = data.table::rbindlist(db.rows)
  db.formatted$charge = c(0)

  list(db = db.formatted, version = version)

}

build.ECMDB <- function(outfolder){

  theurl = "http://ecmdb.ca/downloads"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "Version <strong>(\\d.\\d)")[,2]

  file.url = "http://ecmdb.ca/download/ecmdb.json.zip"
  base.loc <- file.path(outfolder, "ecmdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  zip.file <- file.path(base.loc, "ecmdb.zip")
  utils::download.file(file.url, zip.file,mode = "wb", cacheOK = T)
  utils::unzip(zip.file, exdir = base.loc)
  json = file.path(base.loc, "ecmdb.json")
  json.rows = RJSONIO::fromJSON(json)
  db.base = data.table::rbindlist(json.rows)
  db.formatted = data.table::data.table(
    compoundname = db.base$name,
    description = db.base$description,
    baseformula = db.base$moldb_formula,
    identifier = db.base$met_id,
    charge = db.base$moldb_formal_charge,
    structure = db.base$moldb_smiles
  )

  list(db = db.formatted, version = version)

}

build.LMDB <- function(outfolder){
  # file.url = "http://lmdb.ca/system/downloads/current/structures.zip"
  # base.loc <- file.path(outfolder, "lmdb_source")
  # if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  # zip.file <- file.path(base.loc, "lmdb.zip")
  # utils::download.file(file.url, zip.file,mode = "wb", cacheOK = T)
  # utils::unzip(zip.file, exdir = base.loc)
  # sdf.path <- list.files(base.loc,
  #                        pattern = "sdf$",
  #                        full.names = T,
  #                        recursive = T)
  #
  # desc <- function(sdfset){
  #   mat <- NULL
  #   db <- data.table::as.data.table(ChemmineR::datablock2ma(datablocklist=ChemmineR::datablock(sdfset)))
  #   info = data.table::data.table(identifier = db$DATABASE_ID,
  #                                 compoundname = db$GENERIC_NAME,
  #                                 structure = db$SMILES,
  #                                 baseformula = db$JCHEM_FORMULA,
  #                                 description = paste0("Synonyms: ", db$SYNONYMS))
  #   info
  # }
  #
  # out.csv = file.path(base.loc,
  #                     "lmdb_parsed.csv")
  # if(file.exists(out.csv)){
  #   file.remove(out.csv)
  # }
  #
  # require(ChemmineR)
  # sdfStream.joanna(input=sdf.path, output=out.csv,
  #                  append=FALSE,
  #                  fct=desc,
  #                  silent = T)
  #db.struct <- data.table::fread(file.path(base.loc, "lmdb_parsed.csv"), fill = T, header=T)
  # - lmdb is only available as seperate file here -
  #
  #   lmdb <- XML::xmlParse("~/MetaDBparse/inst/files/lmdb.xml")
  #   xmldf <- data.table::as.data.table(XML::xmlToDataFrame(nodes = XML::getNodeSet(lmdb, "//row")))
  # db.extra.info <- data.table::as.data.table(xmldf[, .(identifier = lmdb_id, description)])
  # lmdb <- xmldf[, .(identifier = lmdb_id,
  #                   compoundname = name,
  #                   baseformula = moldb_formula,
  #                   structure = moldb_smiles,
  #                   description = description,
  #                   charge = c(0))]
  theurl = "http://lmdb.ca/"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "Version <strong>(\\d.\\d)")[,2]

  data(lmdb, package = "MetaDBparse")
  # descs = data.table::fread("~/Downloads/lmdb_descriptions.csv", header=T)
  # descs <- data.table::data.table(identifier = c(colnames(descs)[1], descs[,1][[1]]),
  #                                  description = c(colnames(descs)[2], descs[,2][[1]]))
  # db.a = db.formatted[ , -c("description")]
  # lmdb = merge(db.a, descs)
  #lmdb <- db.final
  db.formatted <- lmdb

  list(db = db.formatted, version = version)

}

build.YMDB <- function(outfolder){
  file.url = "http://www.ymdb.ca/system/downloads/current/ymdb.json.zip"
  base.loc <- file.path(outfolder, "ymdb_source")
  if(dir.exists(base.loc))(unlink(base.loc, recursive = T)); dir.create(base.loc, recursive = T);
  zip.file <- file.path(base.loc, "ymdb.zip")
  utils::download.file(file.url, zip.file,mode = "wb", cacheOK = T)
  utils::unzip(zip.file, exdir = base.loc)
  json = file.path(base.loc, "ymdb.json")
  line = readLines(json)[[1]]
  #stringr::str_match_all(line, "ymdb_id")
  fixed_lines = gsub(line, pattern = "\\]\\}\\{", replacement = "]},{")
  jsonParsed = jsonlite::fromJSON(fixed_lines)
  db.partial =
    data.table::data.table(
      compoundname = jsonParsed$name,
      description = paste0(jsonParsed$description, " Synonyms: ", paste0(jsonParsed$synonyms[[1]], collapse = ",")),
      baseformula = "",
      identifier = jsonParsed$ymdb_id,
      charge = jsonParsed$physiological_charge,
      structure = ""
    )
  require(ChemmineR)
  sdf.url = "http://www.ymdb.ca/system/downloads/current/ymdb.sdf.zip"
  zip.file <- file.path(base.loc, "ymdb_sdf.zip")
  utils::download.file(sdf.url, zip.file,mode = "wb",cacheOK = T)
  utils::unzip(zip.file, exdir = base.loc)

  sdf.path <- list.files(base.loc,
                         pattern = "sdf$",
                         full.names = T,
                         recursive = T)

  desc <- function(sdfset){
    mat <- NULL
    db <- data.table::as.data.table(ChemmineR::datablock2ma(datablocklist=ChemmineR::datablock(sdfset)))
    info = data.table::data.table(identifier = db$DATABASE_ID,
                                  compoundname = db$GENERIC_NAME,
                                  structure = db$SMILES,
                                  baseformula = db$JCHEM_FORMULA,
                                  description = paste0("Synonyms: ", db$SYNONYMS))
    info
  }

  out.csv = file.path(base.loc,
                      "ymdb_parsed.csv")
  if(file.exists(out.csv)){
    file.remove(out.csv)
  }

  sdfStream.joanna(input=sdf.path, output=out.csv,
                   append=FALSE,
                   fct=desc,
                   silent = T)

  db.struct <- data.table::fread(file.path(base.loc, "ymdb_parsed.csv"), fill = T, header=T)
  db.merged <- merge(db.partial, db.struct, by = "identifier", all.y = T)
  db.rows <- pbapply::pblapply(1:nrow(db.merged), function(i){
    row = db.merged[i,]
    data.table::data.table(identifier = row$identifier,
                           compoundname = row$compoundname.y,
                           structure = row$structure.y,
                           baseformula = row$baseformula.y,
                           description = if(is.na(row$description.x)) row$description.y else row$description.x)
  })
  db.formatted = data.table::rbindlist(db.rows)

  theurl = "http://www.ymdb.ca/about"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "Version <strong>(\\d.\\d)")[,2]

  list(db = db.formatted, version = version)

}

build.PAMDB <- function(outfolder){
  file.url = "http://pseudomonas.umaryland.edu/PaDl/PaMet.xlsx"
  base.loc <- file.path(outfolder, "pamdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  xlsx.file <- file.path(base.loc, "pamdb.xlsx")
  utils::download.file(file.url, xlsx.file, mode = "wb", cacheOK = T)
  db.base = data.table::as.data.table(readxl::read_excel(xlsx.file,sheet = 1))
  db.formatted <- data.table::data.table(identifier = db.base$MetID,
                                         compoundname = db.base$Name,
                                         structure = db.base$SMILES,
                                         baseformula = c(NA),
                                         description = gsub(gsub(db.base$Reactions, pattern="\\r", replacement=", "), pattern="\\n", replacement=""),
                                         charge=db.base$Charge)
  theurl = "http://pseudomonas.umaryland.edu/"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "Version  <STRONG>(\\d.\\d)")[,2]

  list(db = db.formatted, version = version)
}

build.mVOC <- function(outfolder){

  categories = c("\\(", "[",	"$",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",
                 "A",	"B",	"C",	"D",	"E",	"F",	"H",	"I",	"M",	"N",	"O",	"P",	"Q",	"S",	"T",	"U")
  hrefFun <- function(x){
    XML::xpathSApply(x,'./a',XML::xmlAttrs)
  }

  urlbase = "http://bioinformatics.charite.de/mvoc/"

  search_urls <- pbapply::pbsapply(categories, function(categ){
    theurl=paste0("http://bioinformatics.charite.de/mvoc/index.php?site=browse&char=", categ)
    print(theurl)
    tables <- XML::readHTMLTable(theurl,elFun = hrefFun)
    tables.nonull <- tables[sapply(tables, function(x) !is.null(x))]
    tables.wrows <- tables.nonull[sapply(tables.nonull, function(x) nrow(x)>0)]
    urls = gsub(unlist(tables.wrows), pattern = "\\./", replacement = urlbase)
    urls[!is.na(urls)]
  })

  db_rows <- pbapply::pblapply(unlist(search_urls), function(theurl){
    try({
      tables <- XML::readHTMLTable(theurl)
      tables.nonull <- tables[sapply(tables, function(x) !is.null(x))]
      end.nameblock <- min(which(!is.na(tables.nonull[[1]][,2])))
      names = tables.nonull[[1]][1:end.nameblock-1,1]
      restblock = tables.nonull[[1]][-c(1:end.nameblock),]
      trans_restblock = as.data.frame(t(restblock))
      colnames(trans_restblock) <- trans_restblock[1,]
      trans_restblock <- trans_restblock[2,]
      # - - - -
      has.desc <- which(sapply(tables.nonull, function(x) "Biological Function" %in% colnames(x)))
      tbl_desc <- tables.nonull[[has.desc]]
      desc <- paste0(sapply(1:nrow(tbl_desc), function(i){
        row = tbl_desc[i, 1:2]
        paste0(row[2], "(",row[1],")")
        #paste0(row, collapse=",")
      }), collapse=". ")
      db.formatted <- data.table::data.table(identifier = trans_restblock$`PubChem ID`,
                                             compoundname = names[1],
                                             structure = trans_restblock$SMILES,
                                             baseformula = trans_restblock$Formula,
                                             description = paste0("Microbes producing this compound:",desc, ". Other names:", paste0(names[-1], collapse=",")),
                                             charge=c(0))
      Sys.sleep(1)
      db.formatted
    })
  })

  db.formatted <- data.table::rbindlist(db_rows[sapply(db_rows, function(x) !is.null(nrow(x)))])

  theurl = "http://bioinformatics.charite.de/mvoc/index.php?site=home"
  header = RCurl::getURL(theurl,.opts = list(ssl.verifypeer = FALSE))
  version = stringr::str_match(header,
                               pattern = "mVOC (\\d.\\d)")[,2]

  list(db = db.formatted, version = version)

}

build.NANPDB <- function(outfolder){
  file.url = "http://african-compounds.org/nanpdb/downloads/smiles/"
  base.loc <- file.path(outfolder, "nanpdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  smi.file <- file.path(base.loc, "nanpdb.smi")
  utils::download.file(file.url, smi.file, mode = "wb", cacheOK = T)
  db.base = data.table::fread(smi.file, sep = "\t")
  db.formatted = db.base[, .(identifier = V2, structure = V1, compoundname = V3)]
  db.formatted$charge = c(0)
  db.formatted$description = c("Found in North Africa. For more info,check the NANPDB website.")

  version = Sys.Date()

  list(db = db.formatted, version = version)
}

build.STOFF <- function(outfolder){
  file.url = "http://www.lfu.bayern.de/stoffident/stoffident-static-content/html/download/SI_Content.zip"
  base.loc <- file.path(outfolder, "stoff_source")
  if(!dir.exists(base.loc)) dir.create(base.loc, recursive = T)
  zip.file <- file.path(base.loc, "stoff.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = T)
  # - - - - - - - - - - - - - - - - - - - - - - - - -
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  excel.file = file.path(base.loc, "SI_Content.xls")
  xlsx.file = gsub(excel.file, pattern = "xls", replacement = "xlsx")
  file.copy(excel.file, xlsx.file)
  db.base = data.table::as.data.table(readxl::read_excel(xlsx.file,sheet = 1))
  db.base.aggr = db.base[ , .(compoundname = unique(Name),
                              baseformula = unique(Formula),
                              structure = unique(SMILES),
                              description = paste0("Synonyms: ", paste(`Additional Names`, collapse=","))), by = Index]
  db.base.aggr$charge = c(0)
  db.formatted = db.base.aggr[!is.na(compoundname)]
  colnames(db.formatted)[1] <- "identifier"

  version = Sys.Date()
  list(db = db.formatted, version = version)
}

