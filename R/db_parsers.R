# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# should return data table!
build.HMDB <- function(outfolder){ # WORKS
  file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
  # ----
  base.loc <- file.path(outfolder, "hmdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "HMDB.zip")
  utils::download.file(file.url, zip.file,mode = "w")
  utils::unzip(zip.file, exdir = base.loc)

  input = file.path(base.loc,"hmdb_metabolites.xml")

  library(XML)
  library(RCurl)
  library(rlist)
  theurl <- getURL("http://www.hmdb.ca/statistics",.opts = list(ssl.verifypeer = FALSE) )
  tables <- readHTMLTable(theurl)
  stats = data.table::rbindlist(tables)
  n = as.numeric(as.character(gsub(x = stats[Description == "Total Number of Metabolites"]$Count, pattern = ",", replacement="")))

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

    if(idx %% 1000 == 0){
      pbapply::setpb(pb, idx)
    }

    idx <<- idx + 1

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
    db.formatted[idx, "charge"] <<- stringr::str_match(xmlValue(properties),
                                              pattern = "formal_charge([+|\\-]\\d*|\\d*)")[,2]
  }

  XML::xmlEventParse(input, branches =
                  list(metabolite = metabolite))

  # - - - - - - - - - - - -
  db.formatted
}

build.METACYC <- function(outfolder){ # WORKS
  # May need to remake smartTable if anything on the website changes unfortunately
  # TODO: download file directly from link, will need a javascript. Maybe Rselenium??

  base.loc <- file.path(outfolder, "metacyc_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  source.file = file.path(base.loc, "All_compounds_of_MetaCyc.txt")
  if(!file.exists(source.file)){
    message("Please download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3729417004' as 'All_compounds_of_MetaCyc.txt' and save in the outfolder/metacyc_source folder.")
    return(NULL)
  }
  metacyc.raw = read.table(source.file,header = T, sep = "\t",
                           check.names = F,fill=T,quote="")

  colnames(metacyc.raw) <- gsub(x=as.character(colnames(metacyc.raw)), pattern = '\\"', replacement="")
  metacyc.raw[] <- lapply(metacyc.raw, gsub, pattern = '\\"', replacement = "")

  compounds <- pbapply::pbsapply(metacyc.raw$Compound, FUN=function(pw){
    pw <- iconv(pw, "latin1", "UTF-8",sub='')
    pw <- pw[pw != " // "]
    pw <- gsub(pw, pattern = "&", replacement="")
    pw <- gsub(pw, pattern = ";", replacement="")
    res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "",perl = T)
    paste0(res, collapse=" --- ")
  })

  db.formatted <- data.table::data.table(compoundname = compounds,
                             description = metacyc.raw$Summary,
                             baseformula = metacyc.raw$`Chemical Formula`,
                             identifier = paste0("METACYC_CP_", 1:nrow(metacyc.raw)),
                             charge = c(0),
                             structure = metacyc.raw$SMILES)

  db.formatted
}

build.CHEBI <- function(outfolder){ # WORKS

  db.full <- {
    release = "latest"
    woAssociations = FALSE
    chebi_download <- tempdir()
  utils::download.file("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/",
                       paste0(chebi_download, "releases.txt"), quiet = TRUE)
  releases <- gsub("rel", "", read.table(paste0(chebi_download,
                                                "releases.txt"), quote = "\"", comment.char = "")[,
                                                                                                  9])
  message("Validating ChEBI release number ... ", appendLF = FALSE)
  if (release == "latest") {
    release <- max(releases)
  }
  else {
    release <- releases[match(release, releases)]
  }
  message("OK")
  ftp <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel",
                release, "/Flat_file_tab_delimited/")
  message("Downloading compounds ... ", appendLF = FALSE)
  utils::download.file(paste0(ftp, "compounds.tsv.gz"), paste0(chebi_download,
                                                               "compounds.tsv"), quiet = TRUE)
  compounds <- as.data.frame.array(read.delim2(paste0(chebi_download,
                                                      "compounds.tsv")))
  message("DONE", appendLF = TRUE)
  message("Downloading synonyms ... ", appendLF = FALSE)
  utils::download.file(paste0(ftp, "names.tsv.gz"), paste0(chebi_download,
                                                           "names.tsv"), quiet = TRUE)
  names <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                   "names.tsv"))))
  message("DONE", appendLF = TRUE)
  message("Downloading formulas ... ", appendLF = FALSE)
  utils::download.file(paste0(ftp, "chemical_data.tsv"), paste0(chebi_download,
                                                                "formulas.tsv"), quiet = TRUE)
  formulas <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                      "formulas.tsv"))))

  message("Downloading structures ... ", appendLF = FALSE)

  utils::download.file(paste0(ftp, "structures.csv.gz"), paste0(chebi_download,
                                                                "structures.csv"), quiet = TRUE)

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

  db.formatted
}

build.FOODB <- function(outfolder){ # WORKS
  file.url <- "http://www.foodb.ca/system/foodb_2017_06_29_csv.tar.gz"
  base.loc <- file.path(outfolder, "foodb_source")

  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "foodb.tar.gz")
  utils::download.file(file.url, zip.file, mode = 'wb', method = 'libcurl')
  utils::untar(normalizePath(zip.file), exdir = normalizePath(base.loc))
  base.table <- data.table::fread(file = file.path(base.loc, "foodb_2017_06_29_csv", "compounds.csv"))

  db.formatted <- data.table::data.table(compoundname = base.table$name,
                                         description = base.table$description,
                                         baseformula = base.table$moldb_formula,
                                         identifier= base.table$id,
                                         charge= base.table$charge,
                                         structure= base.table$moldb_smiles)

  db.formatted <- unique(db.formatted)

  db.formatted
}

build.WIKIDATA <- function(outfolder){ # WORKS
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

  db.formatted
}

# RE-ENABLE AFTER TALKING TO EGON
# build.WIKIPATHWAYS <- function(outfolder){
#   chebi.loc <- file.path(outfolder, "chebi.db")
#   require(SPARQL)
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
  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "respect.zip")
  utils::download.file(file.url, zip.file,mode = "w")
  utils::unzip(normalizePath(zip.file), exdir = (base.loc))

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
        charge = c(NA),
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

  db.formatted
}

build.MACONDA <- function(outfolder){ # NEEDS SPECIAL FUNCTIONALITY

  file.url = "https://www.maconda.bham.ac.uk/downloads/MaConDa__v1_0__csv.zip"

  base.loc <- file.path(outfolder, "maconda_source")

  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "maconda.zip")
  utils::download.file(file.url, zip.file,mode = "w",extra = "-k")
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
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return(db.base)
}

build.T3DB <- function(outfolder){ # WORKS
  # t3db
  file.url <- "http://www.t3db.ca/system/downloads/current/toxins.csv.zip"
  # ----
  base.loc <- file.path(outfolder, "t3db_source")
  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "T3DB.zip")
  utils::download.file(file.url, zip.file,mode = "w")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  db.formatted <- data.table::fread(file.path(base.loc, "toxins.csv"), fill=T)

  db.formatted <- data.table::data.table(
    compoundname = db.formatted$Name,
    baseformula = db.formatted$`Chemical Formula`,
    description = db.formatted$Description,
    charge = c(0),
    identifier = db.formatted$`T3DB ID`,
    structure = db.formatted$SMILES
  )
  return(db.formatted)
}

build.HSDB <- function(outfolder){ # NEEDS WORK
  file.url = "ftp://ftp.nlm.nih.gov/nlmdata/.hsdblease/hsdb.xml.20190528.zip"
  base.loc <- file.path(outfolder, "hsdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "HSDB.zip")
  utils::download.file(file.url, zip.file,mode = "w")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  library(XML)
  library(RCurl)
  library(rlist)

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

  require(webchem)
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
  utils::download.file(file.url, excel.file)

  db.full <- openxlsx::read.xlsx(excel.file, sheet = 1, colNames=T, startRow = 3)

  db.formatted <- unique(data.table::data.table(compoundname = db.full$Compound.Name,
                                                description = paste0("Found in ", db.full$BloodPaperCount, " papers related to blood exposome."),
                                                baseformula = db.full$Molecular.Formula,
                                                identifier = db.full$PubChem.CID,
                                                charge = db.full$Charge,
                                                structure = db.full$CanonicalSMILES))

  db.formatted
}

build.EXPOSOMEEXPLORER <- function(outfolder){ # WORKS
  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/biomarkers.csv.zip"

  base.loc <- file.path(outfolder, "exex_source")

  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "expoexpo_comp.zip")
  utils::download.file(file.url, zip.file,mode = "w")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  base.table <- data.table::fread(file = file.path(base.loc, "biomarkers.csv"))

  db.formatted <- data.table::data.table(compoundname = base.table$Name,
                                         description = base.table$Description,
                                         baseformula = base.table$Formula,
                                         identifier= base.table$ID,
                                         charge= c(NA),
                                         structure= base.table$SMILES)

  db.formatted <- unique(db.formatted)

  # - - use correlations to get some custom descriptions :) - -

  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/correlation_values.csv.zip"

  zip.file <- file.path(base.loc, "expoexpo_corr.zip")
  utils::download.file(file.url, zip.file,mode = "w")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

  corr.table <- data.table::fread(file = file.path(base.loc, "correlation_values.csv"))

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

  db.formatted
}

build.SMPDB <- function(outfolder){ # OK I THINK
  file.url <- "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
  # ----
  base.loc <- file.path(outfolder, "smpdb_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)
  zip.file <- file.path(base.loc, "SMPDB.zip")
  utils::download.file(file.url, zip.file)
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  # -------------------------------

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

  db.formatted
}

build.KEGG <- function(outfolder){ # WORKS
  batches <- split(0:2300, ceiling(seq_along(0:2300)/100))
  cpds <- pbapply::pblapply(batches, FUN=function(batch){
    names(KEGGREST::keggFind("compound", batch, "mol_weight"))
  })
  cpd.ids <- Reduce(c, cpds)
  id.batches <- split(cpd.ids, ceiling(seq_along(cpd.ids)/10))

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

  fns = unlist(kegg.mol.paths)

  rJava::.jcall("java/lang/System","V","gc")
  gc()

  smiles.rows = pbapply::pblapply(fns, function(fn){
    smiles=NA
    try({
      iatom = rcdk::load.molecules(molfiles = fn)[[1]]
      smiles = rcdk::get.smiles(molecule = iatom)
      charge = rcdk::get.total.formal.charge(molecule = iatom)
    })
    id = gsub(basename(fn), pattern = "\\.mol", replacement="")
    data.table::data.table(identifier = id, smiles = smiles, calcharge = charge)
  })

  smitable <- data.table::rbindlist(smiles.rows)

  db.merged <- merge(db.formatted, smitable, by = "identifier")

  db.formatted <- data.table(compoundname = db.merged$compoundname,
                             description = db.merged$description,
                             baseformula = db.merged$baseformula,
                             identifier = db.merged$identifier,
                             charge = db.merged$calcharge,
                             structure = db.merged$smiles)
  # - - - - - -
  db.formatted
}

build.DRUGBANK <- function(outfolder){  # WORKS
  base.loc <- file.path(outfolder, "drugbank_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  zip.file <- file.path(base.loc, "drugbank.zip")
  if(!file.exists(zip.file)){
    file.rename(file.path(base.loc, "drugbank_all_full_database.xml.zip"), zip.file)
  }
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))

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

  res = XML::xmlEventParse(file = file.path(base.loc, "full database.xml"), branches =
                        list("drug" = metabolite, "drugbank-metabolite-id-value" = print))

  db.formatted <- db.formatted[-1,]

  # - - -

  db.formatted
}

build.LIPIDMAPS <- function(outfolder){ # WORKS (description needs some tweaking)
  file.url = "https://www.lipidmaps.org/resources/downloads/LMSD/LMSD_20190711.sdf.zip"

  # ----
  base.loc <- file.path(outfolder, "lipidmaps_source")
  zip.file <- file.path(base.loc, "lipidmaps.zip")
  utils::download.file(file.url, zip.file)
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  # -------------------------------

  sdf.path <- list.files(base.loc,
                         pattern = "sdf",
                         full.names = T,
                         recursive = T)

  desc <- function(sdfset){
    mat <- NULL
    try({
      datablock = data.table::as.data.table(ChemmineR::datablock2ma(datablocklist=ChemmineR::datablock(sdfset)))

      last_cn <<- colnames(datablock)

      if(!("FORMULA" %in% last_cn)){
        mat = as.matrix(data.table::data.table(
          identifier=datablock$PUBCHEM_COMPOUND_CID,
          compoundname = datablock$PUBCHEM_IUPAC_NAME,
          baseformula = datablock$PUBCHEM_MOLECULAR_FORMULA,
          structure = datablock$PUBCHEM_OPENEYE_CAN_SMILES,
          description = paste0("Alternative names: ",
                               apply( datablock[ , grep(colnames(datablock), pattern="NAME"),with=F ] , 1 , paste , collapse = "-" ))
        ))
      }else{
        mat = as.matrix(data.table::data.table(
          identifier = as.character(datablock$LM_ID),
          compoundname = as.character(if("NAME" %in% colnames(datablock)) datablock$NAME else datablock$SYSTEMATIC_NAME),
          baseformula = as.character(datablock$FORMULA),
          structure = as.character(if("SMILES" %in% last_cn) datablock$SMILES else datablock$INCHI),
          description = as.character(paste0("Also known as ", paste0(datablock$ABBREVIATION, " or ", datablock$SYSTEMATIC_NAME)))
        ))
      }
    })
    mat
  }

  sdfStream.joanna(input=sdf.path, output=file.path(base.loc, "lipidmaps_parsed.csv"),
                       append=FALSE,
                       fct=desc,
                       silent = T)

  db.base <- data.table::fread(file.path(base.loc, "lipidmaps_parsed.csv"), fill = T, header=T)

  db.base$charge <- c(NA)

  db.formatted <- db.base[,-1,with=F]

  db.formatted
}

build.METABOLIGHTS <- function(outfolder){
  require(webchem)
  require(rcdk)
  require(curl)

  all_ids <- read.table("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/unichem.tsv", sep="\t")
  colnames(all_ids) <- c("identifier", "inchi", "inchikey")

  metabs <- all_ids$identifier

  # require(RCurl)
  #
  uris <- pbapply::pblapply(metabs, FUN=function(id){
    met_info = NA
    url <- paste0("https://www.ebi.ac.uk/metabolights/webservice/beta/compound/", id)
  })
  #
  # all_info <- jsonlite::read_json("~/Downloads/mapping.json")
  # View(all_info$studyMapping$MTBLS59)
  # View(all_info$compoundMapping$`CSID 16569894`)
  #
  # for(item in all_info$compoundMapping){
  #   for(info in item){
  #     print(info$study)
  #   }
  # }
  #
  # metab_rows <- pbapply::pblapply(all_info$compoundMapping, cl = F, FUN=function(cpd){
  #   met_info <- cpd[[1]]$mafEntry
  #   data.table::data.table(compoundname = met_info$metaboliteIdentification,
  #              description = if(!is.null(met_info$description)) met_info$description else "unknown",
  #              baseformula = met_info$chemicalFormula,
  #              identifier = if(!is.null(met_info$identifier)) met_info$identifier else "unknown",
  #              structure = met_info$smiles,
  #              charge = met_info$charge)
  #
  # })
  #metab_tbl <- rbindlist(metab_rows)


  db_rows <- pbapply::pblapply(metabs[!is.na(metabs)],FUN=function(met_info){
    formula <- NA
    charge <- NA
    try({
      formula <- met_info$formula
      charge = met_info$charge
    })
    smi <- toupper(gsub(" ", "", met_info$smiles))
    # - - - - - - - - - - - - - - --
    if(is.na(formula) | is.na(charge)){
      # inchi or smiles
      worked = FALSE
      try({
        smi <- met_info$smiles
        iatom <- rcdk::parse.smiles(smi)[[1]]
        charge = rcdk::get.total.charge(iatom)
        formula = rcdk::get.mol2formula(iatom, charge = charge)@string
        worked = TRUE
      })
      if(!worked){
        smi <- toupper(webchem::cs_inchi_smiles(inchi = met_info$inchi,verbose = TRUE))
        iatom <- rcdk::parse.smiles(smi)[[1]]
        charge <- rcdk::get.total.charge(iatom)
        formula <- rcdk::get.mol2formula(iatom, charge = charge)@string
      }
    }
    met_dt <- data.table::data.table(compoundname = met_info$name,
                         description = if(!is.null(met_info$definition)) met_info$definition else "Unknown",
                         baseformula = formula,
                         identifier = met_info$id,
                         structure = smi,
                         charge = charge)
    # -------
    list(metabs=met_dt, pathway=met_info$pathways)
  })

  db_rows_a <- lapply(db_rows, function(x) x$metabs)
  db.formatted <- data.table::rbindlist(db_rows_a[!is.na(db_rows_a)])

  db.formatted
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
    utils::download.file(url, zip.file)
    utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
  })
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

  # - - - -

  db.formatted
}

build.VMH <- function(outfolder){ # WORKS
  api_url <- "https://vmh.uni.lu/_api/metabolites/"

  pagerange = 150
  # get the first page

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

  db.formatted
}

build.PHENOLEXPLORER <- function(outfolder){ # WORKS
  file.urls = paste0(
    "http://phenol-explorer.eu/system/downloads/current/",
    c("composition-data.xlsx.zip",
      "compounds.csv.zip",
      "compounds-structures.csv.zip",
      "metabolites.csv.zip",
      "metabolites-structures.csv.zip"))

  base.loc <- file.path(outfolder, "phenolexplorer_source")
  if(!dir.exists(base.loc)) dir.create(base.loc)

  for(url in file.urls){
    zip.file <- file.path(base.loc, basename(url))
    utils::download.file(url, destfile = zip.file)
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

  db.formatted
}

build.MASSBANK <- function(outfolder){ # WORKS
  file.url <- "https://github.com/MassBank/MassBank-data/archive/master.zip"
  base.loc <- file.path(outfolder, "massbank_source")
  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
  zip.file <- file.path(base.loc, "massbank.zip")
  utils::download.file(file.url, zip.file,mode = "w", method='libcurl')
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
  db.formatted
}

build.SUPERNATURAL <- function(outfolder){ # NEEDS WORK, REALLY SLOW TOO INFEASIBLE
  base.loc <- file.path(outfolder, "supernatural_source")
  if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)

  # http://bioinf-applied.charite.de/supernatural_new/src/download_mols.php?sn_id=SN00270303,SN00332107,SN00371241,SN00399119,SN00429427,SN00332221,SN00228244,SN00355833,SN00295142,SN00295267,SN00228769,SN00268639,SN00288836,SN00372987,SN00297293

  library(XML)
  library(RCurl)
  library(rlist)
  theurl <- getURL("http://bioinf-applied.charite.de/supernatural_new/",.opts = list(ssl.verifypeer = FALSE) )
  tables <- readHTMLTable(theurl)
  stats = data.table::rbindlist(tables)
  n = as.numeric(
    gsub(stringr::str_match(theurl, pattern="contains (.*) natural compounds")[,2], pattern = ",", replacement = '')
  )

  # http://bioinf-applied.charite.de/supernatural_new/src/download_mol.php?sn_id=SN00000001
  base.url = "http://bioinf-applied.charite.de/supernatural_new/src/download_mol.php?sn_id="
  id.nr = stringr::str_pad(i, 8, pad = "0")
  id.str = paste0("SN", id.nr)
  file.url = paste0(base.url, id.str)

  all.ids = paste0("SN", stringr::str_pad(1:n, 8, pad = "0"))

  # break into blocks


  # this might be too big... mail the ppl if they want to upload the whole thing..

}
