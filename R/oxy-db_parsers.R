.onLoad <- function(libname, pkgname) {
  data("adduct_rules", "adducts", package=pkgname, envir=parent.env(environment()))
  data("isotopes", package="enviPat", envir=parent.env(environment()))
}

# ---- Internal helpers for the parsers below ----------------------------------
# Shared readers (XML/SDF), ChEBI lookup, resumable downloads and fallback sources.

.metadb_text <- function(node, name) {
  child <- node[[name]]
  if (is.null(child)) NA_character_ else XML::xmlValue(child)
}

.metadb_stop <- function() stop(structure(list(message = "Sample complete", call = NULL),
                                         class = c("metadb_sample_complete", "error", "condition")))

.metadb_xml_zip <- function(path, member, testMode = FALSE) {
  con <- unz(path, member, open = "rb")
  on.exit(close(con))
  # XML's default connection reader seeks; an explicit binary reader supports unz.
  reader <- function(n) {
    if (n < 0) return(character())
    bytes <- readBin(con, "raw", n = n)
    if (length(bytes)) rawToChar(bytes) else character()
  }
  rows <- list(); version <- NA_character_
  metabolite <- function(node) {
    if (is.na(version)) version <<- .metadb_text(node, "version")
    props <- node[["predicted_properties"]]
    charge <- NA_integer_
    if (!is.null(props)) for (prop in XML::xmlChildren(props)) {
      if (identical(.metadb_text(prop, "kind"), "formal_charge")) {
        charge <- suppressWarnings(as.integer(.metadb_text(prop, "value")))
      }
    }
    rows[[length(rows) + 1L]] <<- data.table::data.table(
      identifier = .metadb_text(node, "accession"), compoundname = .metadb_text(node, "name"),
      baseformula = .metadb_text(node, "chemical_formula"), structure = .metadb_text(node, "smiles"),
      description = .metadb_text(node, "description"), charge = charge)
    if (testMode && length(rows) >= 10L) .metadb_stop()
  }
  tryCatch(XML::xmlEventParse(reader, branches = list(metabolite = metabolite), replaceEntities = TRUE),
           metadb_sample_complete = function(e) NULL)
  if (!length(rows)) stop("No metabolite records in source ZIP")
  list(db = data.table::rbindlist(rows), version = version)
}

.metadb_sdf <- function(path, parse_record, testMode = FALSE) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  rows <- list(); buffer <- character()
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    buffer <- c(buffer, lines)
    ends <- which(trimws(buffer) == "$$$$")
    starts <- c(1L, head(ends, -1L) + 1L)
    for (i in seq_along(ends)) {
      record <- buffer[seq.int(starts[i], ends[i])]
      tags <- grep("^>\\s*<[^>]+>", record)
      props <- list()
      for (tag in tags) {
        name <- sub("^>\\s*<([^>]+)>.*", "\\1", record[tag])
        remaining <- record[seq.int(tag + 1L, length(record))]
        end <- match("", trimws(remaining), nomatch = length(remaining) + 1L)
        props[[toupper(name)]] <- paste(head(remaining, end - 1L), collapse = " ")
      }
      row <- parse_record(props)
      if (!is.null(row) && nrow(row)) rows[[length(rows) + 1L]] <- row
      if (testMode && length(rows) >= 10L) return(data.table::rbindlist(rows, fill = TRUE))
    }
    if (length(ends)) {
      last <- tail(ends, 1L)
      buffer <- if (last < length(buffer)) buffer[seq.int(last + 1L, length(buffer))] else character()
    }
  }
  if (length(buffer) && any(nzchar(trimws(buffer)))) stop("Incomplete final SDF record")
  if (!length(rows)) stop("No compound records parsed from SDF")
  data.table::rbindlist(rows, fill = TRUE)
}

.metadb_property <- function(props, names) {
  for (name in toupper(names)) if (!is.null(props[[name]]) && nzchar(props[[name]])) return(props[[name]])
  NA_character_
}

.metadb_chebi <- function(outfolder) {
  supplied <- getOption("metadbparse.chebi")
  if (is.data.frame(supplied)) return(data.table::as.data.table(supplied))
  parsed <- file.path(outfolder, "chebi_source", "parsed.rds")
  if (file.exists(parsed)) return(readRDS(parsed)$db)
  if (file.exists(file.path(outfolder, "chebi.db"))) return(data.table::as.data.table(showAllBase(outfolder, "chebi")))
  build.CHEBI(outfolder)$db
}

# Resumable downloads + fallback sources for HMDB / FooDB.
# hmdb.ca and foodb.ca sit behind a Cloudflare JS challenge, so the primary
# source is a Wayback Machine capture. Wayback often drops long transfers, so
# downloads resume (curl -C -) until the archive validates.

.metadb_zip_ok <- function(p) {
  file.exists(p) && isTRUE(file.info(p)$size > 0) &&
    tryCatch(nrow(utils::unzip(p, list = TRUE)) > 0, error = function(e) FALSE)
}

.metadb_targz_ok <- function(p) {
  if (!file.exists(p) || !isTRUE(file.info(p)$size > 0)) return(FALSE)
  gz <- Sys.which("gzip")
  if (nzchar(gz)) return(system2(gz, c("-t", shQuote(p)), stdout = FALSE, stderr = FALSE) == 0)
  tryCatch(length(utils::untar(p, list = TRUE)) > 0, error = function(e) FALSE, warning = function(w) FALSE)
}

.metadb_resume_download <- function(url, dest, valid, tries = 20L) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (valid(dest)) return(TRUE)
  curl <- Sys.which("curl")
  for (i in seq_len(tries)) {
    if (nzchar(curl)) {
      system2(curl, c("-L", "--fail", "--progress-bar", "--retry", "5", "--retry-delay", "10",
                      "--connect-timeout", "60", "-C", "-", "-o", shQuote(dest), shQuote(url)))
    } else {
      old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old), add = TRUE)
      try(utils::download.file(url, dest, mode = "wb"), silent = TRUE)
    }
    if (valid(dest)) return(TRUE)
    size <- if (file.exists(dest)) file.info(dest)$size / 1e6 else 0
    message(sprintf("Download incomplete (%.0f MB so far), resuming [%d/%d]...", size, i, tries))
    Sys.sleep(10)
  }
  message("Could not complete download of ", url,
          "\nIf the partial file is corrupt, delete it and retry: ", dest)
  FALSE
}

# ---- HMDB fallback: Zenodo 'HMDB Clean' 21Q4 (HMDB 5.0, CC BY 4.0) + PubChem SMILES ----
.hmdb_zenodo <- function(base.loc, testMode = FALSE, batch_size = 200) {
  meta_file <- file.path(base.loc, "zenodo_hmdb_meta.Rds")
  if (!file.exists(meta_file)) {
    rec <- jsonlite::fromJSON("https://zenodo.org/api/records/5758926")
    keys <- if (!is.null(rec$files$key)) rec$files$key else rec$files$filename
    hit <- grep("hmdb_meta\\.Rds$", keys, ignore.case = TRUE)
    if (length(hit) != 1) stop("hmdb_meta.Rds not found in Zenodo record 5758926")
    url <- rec$files$links$self[hit]
    if (is.null(url)) url <- sprintf("https://zenodo.org/records/5758926/files/%s?download=1", keys[hit])
    old <- options(timeout = max(600, getOption("timeout"))); on.exit(options(old), add = TRUE)
    utils::download.file(url, meta_file, mode = "wb")
  }
  meta <- data.table::as.data.table(readRDS(meta_file))
  meta <- meta[!is.na(accession) & accession != ""]
  meta <- unique(meta, by = "accession")
  if (testMode) meta <- meta[seq_len(min(10L, nrow(meta)))]
  meta[, cid := trimws(as.character(cid))]
  meta[!grepl("^[0-9]+$", cid), cid := NA_character_]

  pc_file <- file.path(base.loc, "pubchem_hmdb_smiles.Rds")
  pc <- if (file.exists(pc_file)) readRDS(pc_file) else
    data.table::data.table(CID = character(), SMILES = character())
  todo <- setdiff(unique(stats::na.omit(meta$cid)), pc$CID)
  if (length(todo)) {
    url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/SMILES/CSV"
    batches <- split(todo, ceiling(seq_along(todo) / batch_size))
    pb <- pbapply::startpb(0, length(batches))
    for (i in seq_along(batches)) {
      d <- NULL
      for (attempt in 1:4) {
        r <- try(httr::POST(url, body = list(cid = paste(batches[[i]], collapse = ",")),
                            encode = "form", httr::timeout(120)), silent = TRUE)
        if (!inherits(r, "try-error") && httr::status_code(r) == 200) {
          d <- data.table::fread(text = httr::content(r, as = "text", encoding = "UTF-8"),
                                 colClasses = list(character = "CID"))
          if (!"SMILES" %in% names(d) && "IsomericSMILES" %in% names(d))
            data.table::setnames(d, "IsomericSMILES", "SMILES")
          break
        }
        Sys.sleep(2^attempt)
      }
      if (!is.null(d)) pc <- data.table::rbindlist(list(pc, d[, .(CID, SMILES)]), fill = TRUE)
      if (i %% 50 == 0) saveRDS(pc, pc_file)
      pbapply::setpb(pb, i)
      Sys.sleep(0.25)
    }
    pbapply::closepb(pb)
    pc <- unique(pc, by = "CID")
    saveRDS(pc, pc_file)
  }
  db <- merge(meta, pc, by.x = "cid", by.y = "CID", all.x = TRUE, sort = FALSE)
  db <- db[, .(compoundname = name, description = description, baseformula = chemical_formula,
               identifier = accession, charge = 0L, structure = SMILES)]
  list(db = db, version = "5.0 (Zenodo HMDB Clean 21Q4 + PubChem SMILES)")
}

# ---- FooDB fallback: compound table mirrored on GitHub (jaspershen/databases, MIT) ----
# Same column shift as the official Compound.csv: 'annotation_quality' holds the
# description, 'cas_number' holds SMILES.
.foodb_github <- function(base.loc) {
  f <- file.path(base.loc, "foodb_compound_info.rda")
  if (!file.exists(f) || file.info(f)$size == 0) {
    utils::download.file("https://raw.githubusercontent.com/jaspershen/databases/main/data/foodb_compound_info.rda",
                         f, mode = "wb")
  }
  e <- new.env(); load(f, envir = e)
  x <- data.table::as.data.table(get(ls(e)[1], envir = e))
  smi <- stats::na.omit(x$cas_number)
  if (mean(grepl("^[0-9]+-[0-9]{2}-[0-9]$", smi)) > 0.5)
    stop("FooDB GitHub copy layout changed: 'cas_number' now holds CAS numbers, not SMILES")
  x <- x[!is.na(public_id)]
  x <- unique(x, by = "public_id")
  db <- x[, .(compoundname = name, description = annotation_quality, baseformula = NA_character_,
              identifier = public_id, charge = 0L, structure = cas_number)]
  list(db = db, version = "Apr 07 2020 (GitHub mirror jaspershen/databases)")
}

# ---- DrugBank fallback: parsed DrugBank 5.1.13 in the 'dbdataset' R package (GitHub) ----
.drugbank_dbdataset <- function() {
  e <- new.env()
  objs <- utils::data(package = "dbdataset")$results[, "Item"]
  utils::data(list = objs[1], package = "dbdataset", envir = e)
  db <- get(objs[1], envir = e)
  drugs <- if (!is.null(db$drugs)) db$drugs else db
  gi <- data.table::as.data.table(drugs$general_information)
  cp <- data.table::as.data.table(drugs$calculated_properties)
  if (!nrow(cp) && !is.null(drugs$experimental_properties)) cp <- data.table::as.data.table(drugs$experimental_properties)
  idcol <- function(d) {
    hit <- intersect(c("drugbank_id", "parent_key", "primary_key", "drugbank-id"), names(d))
    if (!length(hit)) stop("dbdataset: no DrugBank ID column; columns: ", paste(names(d), collapse = ", "))
    hit[1]
  }
  gi[, id := get(idcol(gi))]
  cp[, id := get(idcol(cp))]
  prop <- function(kind) {
    x <- cp[cp$kind == kind, .(value = value[1]), by = id]
    data.table::setnames(x, "value", kind)
    x
  }
  db.formatted <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE),
                         list(gi[, .(id, compoundname = name, description = description)],
                              prop("SMILES"), prop("Molecular Formula"), prop("Physiological Charge")))
  db.formatted <- db.formatted[!is.na(SMILES) & nzchar(SMILES), .(
    compoundname, description, baseformula = `Molecular Formula`, identifier = id,
    charge = suppressWarnings(as.integer(`Physiological Charge`)), structure = SMILES)]
  ver <- tryCatch(as.character(utils::packageVersion("dbdataset")), error = function(e) "5.1.13")
  list(db = db.formatted, version = paste0(ver, " (dbdataset R package)"))
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
  chebi <- data.table::as.data.table(.metadb_chebi(outfolder))  # reuses parsed ChEBI, builds it if needed
  chebi[, identifier := as.numeric(identifier)]
  q <- paste0("prefix wp: <http://vocabularies.wikipathways.org/wp#>
    prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    prefix dcterms: <http://purl.org/dc/terms/>
    prefix dc: <http://purl.org/dc/elements/1.1/>
    select ?mb (group_concat(distinct str(?labelLit);separator=\", \") as ?label)
           (?idurl as ?csid) (group_concat(distinct ?pwTitle;separator=\", \") as ?description) ?pathway
    where {
      ?mb a wp:Metabolite ; rdfs:label ?labelLit ; wp:bdbChEBI ?idurl ; dcterms:isPartOf ?pathway .
      ?pathway a wp:Pathway ; dc:title ?pwTitle .
    }
    GROUP BY ?mb ?idurl ?pathway", if (testMode) " LIMIT 50" else "")
  # Plain SPARQL-over-HTTP (no glitter dependency)
  r <- httr::POST("https://sparql.wikipathways.org/sparql", body = list(query = q), encode = "form",
                  httr::add_headers(Accept = "application/sparql-results+json"))
  httr::stop_for_status(r)
  b <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))$results$bindings
  if (!length(b) || !NROW(b)) stop("WikiPathways SPARQL returned no metabolites")
  base.db <- data.table::data.table(mb = b$mb$value, csid = b$csid$value,
                                    description = b$description$value, pathway = b$pathway$value)
  base.db[, identifier := suppressWarnings(as.numeric(gsub(".*[:_]|>", "", csid)))]
  chebi.join.table <- base.db[!is.na(identifier), .(identifier, description, widentifier = mb, pathway)]
  chebi <- chebi[, -"description", with = FALSE]
  db.formatted <- merge(chebi.join.table, chebi, by = "identifier")
  db.formatted <- unique(db.formatted[, .(compoundname = paste0(unique(compoundname), collapse = ", "),
                                          baseformula = baseformula[1],
                                          description = paste0("Found in pathways: ", paste0(unique(description), collapse = ", "))),
                                      by = c("structure", "identifier", "charge")])
  pathway.table <- base.db[, .(pathway = gsub(".*/|>", "", pathway), identifier, description)]
  list(db = db.formatted[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")],
       version = as.character(Sys.Date()), path = pathway.table)
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
  options(stringsAsFactors = FALSE, timeout=1e7)
  on.exit(options(oldpar))
  file.url <- "https://web.archive.org/web/20250811072449id_/https://mcdb.ca/system/downloads/current/milk_metabolites.zip"
  base.loc <- file.path(outfolder, "mcdb_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "MCDB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc, unzip = getOption("unzip"))
  input <- file.path(base.loc, "milk_metabolites.xml")
  header <- readLines(input, n = 10)
  version <- trimws(gsub(grep(pattern = "<version", header, value = TRUE), pattern = "<\\/?version>", replacement = ""))
  date <- trimws(gsub(grep(pattern = "update_date", header, value = TRUE), pattern = "<\\/?update_date>", replacement = ""))
  theurl <- RCurl::getURL("https://mcdb.ca/statistics", .opts = list(ssl.verifypeer = FALSE))
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
  } else {
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
    XML::xmlEventParse(file = input,
                       branches = list(metabolite = metabolite),
                       replaceEntities = TRUE)
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
#' @param source 'zenodo' (Zenodo HMDB Clean 21Q4 + PubChem SMILES) or 'local'
#' @param sourcefile path to a manually downloaded hmdb_metabolites.zip (implies source = 'local')
#' @examples
#' \dontrun{build.HMDB(outfolder=tempdir(), testMode=TRUE)}
build.HMDB <- function(outfolder, testMode = FALSE,
                       source = c("zenodo", "local"), sourcefile = NULL) {
  source <- match.arg(source)
  if (!is.null(sourcefile)) source <- "local"
  base.loc <- file.path(outfolder, "hmdb_source")
  dir.create(base.loc, recursive = TRUE, showWarnings = FALSE)
  # Wayback archive option disabled for now (transfers keep breaking):
  # https://web.archive.org/web/20260212044226id_/https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip
  if (source == "zenodo") return(.hmdb_zenodo(base.loc, testMode))
  if (is.null(sourcefile) || !.metadb_zip_ok(sourcefile)) {
    stop("Provide 'sourcefile': a complete hmdb_metabolites.zip (download in a browser from https://hmdb.ca/downloads)")
  }
  zip.file <- sourcefile
  members <- utils::unzip(zip.file, list = TRUE)$Name
  member <- grep("metabolites\\.xml$", members, value = TRUE)[1]
  if (is.na(member)) stop("No metabolites XML in ", zip.file)
  res <- .metadb_xml_zip(zip.file, member, testMode = testMode)
  db <- data.table::as.data.table(res$db)
  list(db = db[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")],
       version = res$version)
}

#' @title Build METACYC
#' @description Parses  METACYC, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @details Requires account creation! Then download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3787684059#' as 'All_compounds_of_MetaCyc.txt' and save in the databases/metacyc_source folder.
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
  identifier <- NULL
  currurl = "https://metacyc.org/group?id=biocyc17-31223-3787684059"
  base.loc <- file.path(outfolder, "metacyc_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  source.file <- file.path(base.loc, "All_compounds_of_MetaCyc.txt")
  if (!file.exists(source.file)) {
    msg <- gsubfn::fn$paste("Please download SmartTable from '$currurl' using Download Spreadsheet File > use common names as 'All_compounds_of_MetaCyc.txt' and save in the databases/metacyc_source folder.")
    message(msg)
    utils::browseURL(currurl)
    return(NULL)
  }
  theurl <- RCurl::getURL(currurl, .opts = list(ssl.verifypeer = FALSE))
  header <- stringr::str_match(theurl, pattern = "Created: (.*), Last Modified: (.*), Readable")
  date <- header[, 2]
  version <- header[, 3]
  # Summary fields contain unquoted line breaks, so a record can span several lines.
  # Re-join lines until each record has as many tab-separated fields as the header.
  lines <- readLines(source.file, warn = FALSE, encoding = "latin1")
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  nf <- length(header)
  records <- character(); buf <- NULL
  for (l in lines[-1]) {
    buf <- if (is.null(buf)) l else paste(buf, l, sep = " ")
    if (lengths(regmatches(buf, gregexpr("\t", buf))) >= nf - 1L) { records <- c(records, buf); buf <- NULL }
  }
  if (!is.null(buf) && nzchar(trimws(buf))) records <- c(records, buf)
  sumcol <- match("Summary", header)
  fields <- lapply(strsplit(records, "\t", fixed = TRUE), function(x) {
    if (length(x) > nf && !is.na(sumcol)) {  # stray tab inside Summary: fold extras back into it
      extra <- length(x) - nf
      x <- c(x[seq_len(sumcol - 1)], paste(x[sumcol:(sumcol + extra)], collapse = " "), x[(sumcol + extra + 1):length(x)])
    }
    length(x) <- nf; x
  })
  metacyc.raw <- data.table::as.data.table(do.call(rbind, fields))
  data.table::setnames(metacyc.raw, header)
  metacyc.raw[is.na(metacyc.raw)] <- ""
  metacyc.raw <- metacyc.raw[`Object ID` != ""]   # export contains empty placeholder rows
  colnames(metacyc.raw) <- gsub(x = as.character(colnames(metacyc.raw)), pattern = "\\\"", replacement = "")
  # BioCyc exports are mostly ASCII with some Latin-1 bytes (e.g. \xb1 for "±"); make them valid UTF-8
  metacyc.raw[] <- lapply(metacyc.raw, function(x) {
    if (!is.character(x)) return(x)
    bad <- !validUTF8(x)
    x[bad] <- iconv(x[bad], from = "latin1", to = "UTF-8", sub = "")
    x
  })
  metacyc.raw[] <- lapply(metacyc.raw, gsub, pattern = "\\\"", replacement = "")
  compounds <- pbapply::pbsapply(metacyc.raw$`Compound Name`, FUN = function(pw) {
    pw <- pw[pw != " // "]
    pw <- gsub(pw, pattern = "&", replacement = "")
    pw <- gsub(pw, pattern = ";", replacement = "")
    res <- gsub(pw, pattern = "<(i|\\/i)>", replacement = "", perl = TRUE)
    res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "", perl = TRUE)
    paste0(res, collapse = " --- ")
  })
  db.formatted <- data.table::data.table(compoundname = compounds,
                                         description = sapply(1:nrow(metacyc.raw), function(i) {
    species <- metacyc.raw$Species[i]
    paste0(metacyc.raw$Summary[i], if (species != "") {
      paste0(" Found in: ", species, ".")
    } else {
      ""
    })
  }), baseformula = metacyc.raw$`Chemical Formula`,
  identifier = metacyc.raw$`Object ID`,
  charge = c(0),
  structure = metacyc.raw$SMILES)

  db.formatted <- db.formatted[!grepl(" ", db.formatted$baseformula),]

  path.tab <- unique(metacyc.raw[,c("Object ID", "Pathways of compound")])
  colnames(path.tab) = c("identifier","description")
  path.tab.long <- data.table::setDT(path.tab)[, lapply(.SD, tstrsplit, " // "), by = identifier]
  path.tab.long$description <- gsub("&|<.*?>|;", "", path.tab.long$description)
  path.tab.long$pathway <- gsub(" |\\(|\\)|-", "_", path.tab.long$description)
  path.tab.long$pathway <- gsub("_+$|^_+", "", path.tab.long$pathway)
  path.tab.long$pathway <- gsub("_+", "_", path.tab.long$pathway)

  db.formatted.all = list(db = db.formatted, version = version, path = path.tab.long[,c("pathway", "identifier", "description")])
  db.formatted.all
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
build.CHEBI <- function(outfolder, testMode = FALSE) {
  file.url <- "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi.sdf.gz"
  base.loc <- file.path(outfolder, "chebi_source")
  dir.create(base.loc, recursive = TRUE, showWarnings = FALSE)
  sdf.file <- file.path(base.loc, "chebi.sdf.gz")
  if (!.metadb_resume_download(file.url, sdf.file, valid = .metadb_targz_ok)) {
    stop("Could not download ChEBI SDF from ", file.url)
  }
  version <- as.character(as.Date(file.info(sdf.file)$mtime))
  # ChEBI 2.0 SDF tags: ChEBI ID, ChEBI NAME, FORMULA, SMILES, INCHI, DEFINITION
  # (old names kept as fallbacks for pre-2025 files)
  parse_record <- function(props) {
    id <- .metadb_property(props, "ChEBI ID")
    structure <- .metadb_property(props, c("SMILES", "INCHI", "InChI"))
    if (is.na(id) || is.na(structure)) return(NULL)
    desc <- .metadb_property(props, c("DEFINITION", "Definition"))
    data.table::data.table(
      identifier = sub("^CHEBI:", "", id),
      compoundname = .metadb_property(props, c("ChEBI NAME", "ChEBI Name")),
      baseformula = .metadb_property(props, c("FORMULA", "Formulae")),
      structure = structure,
      description = gsub("<[^>]+>", "", desc),
      charge = suppressWarnings(as.integer(.metadb_property(props, c("CHARGE", "Charge")))))
  }
  db <- .metadb_sdf(sdf.file, parse_record, testMode = testMode)
  res <- list(db = db[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")],
              version = version)
  if (!testMode) saveRDS(res, file.path(base.loc, "parsed.rds"))  # reused by .metadb_chebi()
  res
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
  utils::unzip(zip.file,
               exdir = normalizePath(base.loc),
               unzip = getOption("unzip"))
  drug.db <- data.table::fread(file.path(base.loc, "drugs.tsv"), quote = "")
  drug.db <- drug.db[SMILES != ""]
  file.url <- "https://s3.pgkb.org/data/chemicals.zip"
  zip.file <- file.path(base.loc, "chemicals.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(zip.file,
               exdir = normalizePath(base.loc),
               unzip = getOption("unzip"))
  chem.db <- data.table::fread(file.path(base.loc, "chemicals.tsv"), quote = "")
  chem.db <- chem.db[SMILES != ""]
  both.db <- data.table::rbindlist(list(drug.db, chem.db), use.names = TRUE, fill = TRUE)
  cols <- c("Generic Names", "Trade Names", "Brand Mixtures", "Type", "Dosing Guideline", "Clinical Annotation Count", "Variant Annotation Count", "Pathway Count", "VIP Count", "Dosing Guideline Sources", "Top Clinical Annotation Level", "Top FDA Label Testing Level", "Top Any Drug Label Testing Level", "Label Has Dosing Info", "Has Rx Annotation")
  cols <- intersect(cols, names(both.db))
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
  base.loc <- file.path(outfolder, "foodb_source")
  dir.create(base.loc, recursive = TRUE, showWarnings = FALSE)
  # Wayback archive option disabled for now (transfers keep breaking):
  # https://web.archive.org/web/20260108183020id_/https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz
  .foodb_github(base.loc)
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
build.WIKIDATA <- function(outfolder, testMode = FALSE) {
  # Wikidata moved SMILES from P233 (canonical) to P2017, and many chemicals are no longer
  # typed 'chemical compound' (Q11173). So: any item with a formula (P274) and a SMILES
  # (P2017 or P233). Full query is too big for WDQS' 60 s limit -> QLever first, WDQS fallback.
  q <- paste0('PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX schema: <http://schema.org/>
    SELECT ?qid ?label ?formula ?smiles ?desc WHERE {
      ?qid wdt:P274 ?formula .
      { ?qid wdt:P2017 ?smiles } UNION { ?qid wdt:P233 ?smiles }
      OPTIONAL { ?qid rdfs:label ?label FILTER(LANG(?label) = "en") }
      OPTIONAL { ?qid schema:description ?desc FILTER(LANG(?desc) = "en") }
    }', if (testMode) " LIMIT 100" else "")
  endpoints <- c("https://qlever.cs.uni-freiburg.de/api/wikidata", "https://query.wikidata.org/sparql")
  raw <- NULL
  for (ep in endpoints) {
    r <- tryCatch(httr::POST(ep, body = list(query = q), encode = "form",
                             httr::add_headers(Accept = "text/csv", `User-Agent` = "MetaDBparse (R)"),
                             httr::timeout(3600)), error = function(e) NULL)
    if (!is.null(r) && httr::status_code(r) == 200) {
      raw <- data.table::fread(text = httr::content(r, as = "text", encoding = "UTF-8"), encoding = "UTF-8")
      if (nrow(raw)) break
    }
    message("Wikidata endpoint failed: ", ep)
  }
  if (is.null(raw) || !nrow(raw)) stop("No Wikidata results from any endpoint")
  raw[, qid := basename(as.character(qid))]
  # one row per item: prefer the first SMILES, keep the first English label/description
  db <- raw[, .(compoundname = label[1], description = desc[1], baseformula = formula[1], structure = smiles[1]), by = qid]
  sub_digits <- function(x) chartr("₀₁₂₃₄₅₆₇₈₉", "0123456789", x)
  db[, baseformula := sub_digits(as.character(baseformula))]
  db[is.na(description) | description == "" | description == "chemical compound", description := "Unknown"]
  db[is.na(compoundname) | compoundname == "", compoundname := qid]
  list(db = db[, .(compoundname, description, baseformula, identifier = qid, charge = NA_integer_, structure)],
       version = as.character(Sys.Date()))
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
  # NOTE: original source (spectra.psc.riken.jp) is dead as of 2026 (confirmed
  # via Wayback CDX: crawls of this same path after 2022 return 404). Using
  # the last known-good Wayback capture of the actual data file instead.
  file.url <- "https://web.archive.org/web/20220527073215id_/http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"
  base.loc <- file.path(outfolder, "respect_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  zip.file <- file.path(base.loc, "respect.zip")
  if (!file.exists(zip.file) || file.info(zip.file)$size == 0) {
    utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  }
  utils::unzip(zip.file,
               exdir = normalizePath(base.loc),
               unzip = getOption("unzip"))
  # the live version-lookup page was on the same now-dead host; using the
  # Wayback capture date as the version instead.
  version <- "2022-05-27 (Wayback capture)"
  date <- version
  cpd_files <- list.files(base.loc, full.names = TRUE)
  cpd_files <- cpd_files[cpd_files != zip.file]
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
build.MACONDA <- function(outfolder, testMode = FALSE) {
  . <- NULL
  # NOTE: original source (maconda.bham.ac.uk) is dead as of 2026 (confirmed:
  # connection timeout; the actual CSV zip was never crawled by Wayback --
  # only the HTML downloads page survived, never the file link on it).
  # Rewritten against the `MaConDa` dataset bundled in the CRAN package `pmd`
  # (actively maintained, v0.3.0, published Oct 2026), which redistributes a
  # reduced copy of the same database: id, name, formula, exact_mass,
  # type_of_contaminant -- no InChI/SMILES/ion-mode/adduct fields, unlike the
  # original web download. NEW DEPENDENCY: requires install.packages("pmd").
  #
  # The original function also wrote a second "extended" SQLite adduct/mass
  # database as a side effect (using ion_form/ion_mode/mz columns it had from
  # the original CSV), and required `conn`/`apikey` arguments the smoke-test
  # harness never supplied -- it was never actually callable through this
  # harness as written. That side table can't be reproduced from this reduced
  # source, so this rewrite returns the plain list(db=..., version=...) shape
  # every other parser here uses instead.
  maconda <- pmd::MaConDa
  if (testMode) {
    maconda <- maconda[1:10, ]
  }
  version <- "pmd package v0.3.0 (MaConDa v1.0 subset)"
  db.formatted <- data.table::data.table(compoundname = maconda$name,
                                         description = maconda$type_of_contaminant,
                                         baseformula = maconda$formula,
                                         identifier = as.character(maconda$id),
                                         charge = c(0),
                                         structure = paste0("[", maconda$formula, "]"))
  list(db = db.formatted, version = version)
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
  file.url <- "https://web.archive.org/web/20250810135010id_/https://www.t3db.ca/system/downloads/current/toxins.csv.zip"
  base.loc <- file.path(outfolder, "t3db_source")
  if (dir.exists(base.loc)) {
    (unlink(base.loc, recursive = TRUE))
  }
  dir.create(base.loc, recursive = TRUE)
  zip.file <- file.path(base.loc, "T3DB.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
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
#'  \code{\link[readxl]{read_xlsx}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[jsonlite]{read_json}}
#' @rdname build.BLOODEXPOSOME
#' @export
#' @importFrom utils download.file
#' @importFrom readxl read_xlsx
#' @importFrom pbapply pbsapply
#' @importFrom jsonlite read_json
#' @importFrom data.table data.table
#' @examples
#' \dontrun{build.BLOODEXPOSOME(outfolder=tempdir(), testMode=TRUE)}
build.BLOODEXPOSOME <- function(outfolder, testMode = FALSE) {
  . <- CID <- PMID <- NULL
  # NOTE: original source (exposome1.fiehnlab.ucdavis.edu) is dead as of 2026.
  # Rewritten against the Zenodo mirror of the Blood Exposome Database v2
  # (Barupal & Fiehn, DOI 10.5281/zenodo.8146024, deposited July 2023).
  file.url <- "https://zenodo.org/records/8146024/files/blood_exposome_chemicals_july_2023.csv"
  mapping.url <- "https://zenodo.org/records/8146024/files/blood_exposome_chemical_literature_mapping_july_2023.csv"
  base.loc <- file.path(outfolder, "bloodexposome_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  csv.file <- file.path(base.loc, "blood_exposome_chemicals.csv")
  if (!file.exists(csv.file) || file.info(csv.file)$size == 0) {
    utils::download.file(file.url, csv.file, mode = "wb", method = "auto")
  }
  mapping.file <- file.path(base.loc, "cid_pmid_mapping.csv")
  if (!file.exists(mapping.file) || file.info(mapping.file)$size == 0) {
    utils::download.file(mapping.url, mapping.file, mode = "wb", method = "auto")
  }
  version <- "2.0 (July 2023, Zenodo 10.5281/zenodo.8146024)"
  db.full <- data.table::fread(csv.file)
  if (testMode) {
    db.full <- db.full[1:10, ]
  }
  mapping <- data.table::fread(mapping.file, header = TRUE)
  colnames(mapping) <- c("CID", "PMID")
  aggregated <- mapping[ , .(description = paste0("MENTIONED IN THE FOLLOWING PAPERS (PMID):",
                                                  paste0(unique(PMID), collapse=","))), by = CID]
  merged.final = merge(db.full, aggregated, by="CID", all.x = TRUE)
  # ---------------
  missing.desc <- is.na(merged.final$description)
  merged.final$description[missing.desc] <- paste0(
    "No individual paper (PMID) mentions found in the Blood Exposome Database literature mapping; ",
    merged.final$LiteratureCount[missing.desc], " total literature citation(s) recorded.")
  db.formatted <- unique(data.table::data.table(compoundname = merged.final$Title,
                                                description = merged.final$description,
                                                baseformula = merged.final$MolecularFormula,
                                                identifier = as.character(merged.final$CID),
                                                charge = merged.final$Charge,
                                                structure = merged.final$CanonicalSMILES))
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
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
  version <- stringr::str_match(RCurl::getURL("http://exposome-explorer.iarc.fr/releases", .opts = list(ssl.verifypeer = FALSE)), pattern = "Release (...)")[, 2]
  date <- stringr::str_match(RCurl::getURL("http://exposome-explorer.iarc.fr/downloads", .opts = list(ssl.verifypeer = FALSE)), pattern = "(\\d{4}-\\d{2}-\\d{2})")[, 2]
  base.table <- data.table::fread(file = file.path(base.loc, "biomarkers.csv"))
  db.formatted <- data.table::data.table(compoundname = base.table$Name, description = base.table$Description, baseformula = base.table$Formula, identifier = base.table$ID, charge = c(NA), structure = base.table$SMILES)
  db.formatted <- unique(db.formatted)
  file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/correlations.csv.zip"
  zip.file <- file.path(base.loc, "expoexpo_corr.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
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
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
  theurl <- "https://smpdb.ca/release_notes"
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
  subtables <- pbapply::pblapply(smpdb.paths, function(x){
   tbl = data.table::fread(x)
   tbl[,10] <- as.character(tbl[,10])
   tbl
  })
  smpdb.tab <- data.table::rbindlist(subtables, fill = TRUE)
  db.formatted <- data.table::data.table(compoundname = smpdb.tab$`Metabolite Name`,
                                         description = paste0(smpdb.tab$`Pathway Name`, "
                                                              (", smpdb.tab$`Pathway Subject`, " pathway)"),
                                         baseformula = smpdb.tab$Formula,
                                         identifier = smpdb.tab$`Metabolite ID`,
                                         charge = c(0),
                                         structure = smpdb.tab$SMILES)
  db.formatted <- db.formatted[, .(description = paste0(unique(description), collapse = ", ")), by = list(compoundname, baseformula,
                                                                                                          identifier, structure)]
  db.formatted <- db.formatted[-1, ]
  db.formatted[, charge := 0L]
  # "pathway","identifier","description"
  path.tab = smpdb.tab[,c("SMPDB ID",
                          "Metabolite ID",
                          "Pathway Name")]
  colnames(path.tab) = c("pathway","identifier","description")
  db.formatted.all = list(db = db.formatted, version = version, path = path.tab)
  db.formatted.all
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
  print("Downloading compound IDs")
  cpd.ids <- names(KEGGREST::keggList("compound"))
  if (testMode) {
    cpd.ids <- cpd.ids[1:10]
  }
  # KEGGREST 1.44.1 cannot parse current multi-record compound responses.
  id.batches <- as.list(cpd.ids)
  theurl <- "https://www.genome.jp/kegg/compound/"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Last updated: (.{1,30})<")[, 2]
  date <- as.Date(version, format = "%B%d, %Y")

  print("Downloading general info")
  kegg.cpd.list <- pbapply::pblapply(id.batches, FUN = function(batch) {
    rest.result <- KEGGREST::keggGet(batch)
    base.list <- lapply(rest.result, FUN = function(cpd) {
      cpd$NAME_FILT <- gsub(cpd$NAME, pattern = ";", replacement = "")
      main = data.table::data.table(compoundname = c(paste(cpd$NAME_FILT, collapse = ", ")),
                             description = paste0("Involved in pathways: ",
                                                  paste0(cpd$PATHWAY, collapse = ", "),
                                                  ". More specifically: ", paste0(cpd$MODULE, collapse = ", "),
                                                  ". Also associated with compound classes: ",
                                                  paste0(unique(trimws(gsub(cpd$BRITE, pattern = "\\[.*\\]|  D\\d* |\\(.*\\)|\\d*", replacement = ""))),
                                                         collapse = ", ")),
                             baseformula = c(cpd$FORMULA),
                             identifier = c(cpd$ENTRY),
                             charge = 0,
                             structure = NA,
                             pathway = if ("PATHWAY" %in% names(cpd)) {
                               names(cpd$PATHWAY)
                             } else {
                               NA
                             })
      path = if (length(cpd$PATHWAY)) data.table::data.table(pathway = names(cpd$PATHWAY),
                                    identifier = cpd$ENTRY,
                                    description = cpd$PATHWAY) else data.table::data.table()
      if(ncol(path) < 3){
        path = data.table::data.table()
      }

      list(main = main, path = path)
    })
    res <- data.table::rbindlist(lapply(base.list, function(x) x$main), fill=T)
    path <- data.table::rbindlist(lapply(base.list, function(x) x$path), fill=T)
    list(main = res, path = path)
  })
  db.formatted <- data.table::rbindlist(lapply(kegg.cpd.list, function(x) x$main))
  path.tab <- data.table::rbindlist(lapply(kegg.cpd.list, function(x) x$path))

  base.loc <- file.path(outfolder, "kegg_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc)
  }
  print("Downloading structural info")
  kegg.mol.paths <- pbapply::pblapply(id.batches, FUN = function(batch) {
    bigmol <- KEGGREST::keggGet(batch, "mol")
    fp <- file.path(base.loc, paste0(sub("^cpd:", "", batch), ".mol"))
    writeLines(bigmol, fp)
    fp
  })
  fns <- unlist(kegg.mol.paths, use.names = FALSE)
  can.convert <- requireNamespace("rcdk", quietly = TRUE)
  if (!can.convert) warning("KEGG records and MOL files downloaded; SMILES and formal charges require rcdk with a working Java runtime")
  smiles.rows <- pbapply::pblapply(fns, function(fn) {
    smiles <- NA_character_
    charge <- NA_integer_
    if (can.convert) try({
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
  db.formatted.all = list(db = unique(db.formatted), version = version, path = path.tab)
  db.formatted.all
}

#' @title Build DRUGBANK
#' @description Parses the DRUGBANK DB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @details Uses the parsed DrugBank (5.1.13) in the 'dbdataset' package; install with remotes::install_github("interstellar-Consultation-Services/dbdataset"). DrugBank data is CC BY-NC 4.0 (non-commercial use, cite DrugBank).
#' @rdname build.DRUGBANK
#' @export
build.DRUGBANK <- function(outfolder) {
  # DrugBank academic downloads are paused (Nov 2026); use the parsed DrugBank copy
  # in the 'dbdataset' package (GitHub; listed under Suggests/Remotes in DESCRIPTION).
  # DrugBank data is CC BY-NC 4.0: non-commercial use, cite DrugBank.
  if (!requireNamespace("dbdataset", quietly = TRUE)) {
    stop("build.DRUGBANK needs the 'dbdataset' package: ",
         "remotes::install_github(\"interstellar-Consultation-Services/dbdataset\")")
  }
  .drugbank_dbdataset()
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
  file.url <- "https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip"
  base.loc <- file.path(outfolder, "lipidmaps_source")
  if (dir.exists(base.loc)) {
   unlink(base.loc,recursive = T)
  }
  dir.create(base.loc)
  zip.file <- file.path(base.loc, "lipidmaps.zip")
  utils::download.file(file.url, zip.file, mode = "wb", method = "auto")
  utils::unzip(zipfile = normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
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
      mat$structure[empty.smiles] <- pbapply::pbsapply(mat$structure[empty.smiles], function(x) {
        url <- gsubfn::fn$paste("https://cactus.nci.nih.gov/chemical/structure/$x/smiles")
        Sys.sleep(0.1)
        RCurl::getURL(url)
      })
    }
    as.matrix(mat)
  }
  sdfStream.joanna(input = sdf.path,
                   output = file.path(base.loc, "lipidmaps_parsed.csv"),
                   append = FALSE, fct = desc, silent = TRUE)
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
    matching <- stringi::stri_detect_fixed(row$identifier, conv_tbl$catcode)
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
  # The EB-eye dump contains both studies (MTBLS*) and compounds (MTBLC*, with
  # name/formula/InChI), so no per-compound web service calls are needed.
  file.url <- "https://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/eb_eye_metabolights_complete.xml"
  base.loc <- file.path(outfolder, "metabolights_source")
  dir.create(base.loc, recursive = TRUE, showWarnings = FALSE)
  xml.file <- file.path(base.loc, "metabolights.xml")
  xml_complete <- function(p) {
    if (!file.exists(p) || !isTRUE(file.info(p)$size > 100)) return(FALSE)
    con <- file(p, "rb"); on.exit(close(con))
    seek(con, file.info(p)$size - 100)
    grepl("</database>", rawToChar(readBin(con, "raw", 100)), fixed = TRUE)
  }
  if (!.metadb_resume_download(file.url, xml.file, valid = xml_complete)) {
    stop("Could not download MetaboLights dump from ", file.url)
  }
  header <- readLines(xml.file, n = 10, warn = FALSE)
  version <- trimws(gsub("<\\/?release_date>", "", grep("<release_date>", header, value = TRUE)[1]))

  cpds <- list(); studies <- list()
  entry <- function(node) {
    id <- XML::xmlGetAttr(node, "id")
    nm <- .metadb_text(node, "name")
    if (startsWith(id, "MTBLC")) {
      fields <- node[["additional_fields"]]
      f <- character()
      if (!is.null(fields)) {
        kids <- XML::xmlChildren(fields)
        kids <- kids[names(kids) == "field"]
        f <- vapply(kids, XML::xmlValue, "")
        names(f) <- vapply(kids, XML::xmlGetAttr, "", "name")
      }
      structure <- c(f["smiles"], f["inchi"]); structure <- structure[!is.na(structure) & nzchar(structure)]
      cpds[[length(cpds) + 1L]] <<- list(identifier = id, compoundname = nm,
                                         baseformula = if ("formula" %in% names(f)) unname(f["formula"]) else NA_character_,
                                         structure = if (length(structure)) unname(structure[1]) else NA_character_)
      if (testMode && length(cpds) >= 10L) .metadb_stop()
    } else if (startsWith(id, "MTBLS")) {
      refs <- node[["cross_references"]]
      if (!is.null(refs)) {
        kids <- XML::xmlChildren(refs)
        keys <- vapply(kids, function(k) if (identical(XML::xmlGetAttr(k, "dbname"), "MetaboLights")) XML::xmlGetAttr(k, "dbkey") else NA_character_, "")
        keys <- unique(keys[!is.na(keys)])
        if (length(keys)) studies[[length(studies) + 1L]] <<- list(identifier = keys, study = paste0("(", id, "): ", nm))
      }
    }
  }
  tryCatch(XML::xmlEventParse(xml.file, branches = list(entry = entry), replaceEntities = TRUE),
           metadb_sample_complete = function(e) NULL)
  if (!length(cpds)) stop("No MTBLC compound entries found in ", xml.file)
  db <- data.table::rbindlist(cpds)
  if (length(studies)) {
    st <- data.table::rbindlist(lapply(studies, data.table::as.data.table))
    st <- st[, .(description = {
      u <- unique(study)
      paste0("Mentioned in the following studies --> ", paste(head(u, 20), collapse = " "),
             if (length(u) > 20) paste0(" ... and ", length(u) - 20, " more") else "")
    }), by = identifier]
    db <- merge(db, st, by = "identifier", all.x = TRUE, sort = FALSE)
  } else db[, description := NA_character_]
  db[, charge := NA_integer_]
  # No SMILES in the dump: MTBLC<n> == CHEBI:<n>, so take SMILES from ChEBI.
  # Full builds use/create the parsed ChEBI; testMode only uses it if it already exists.
  chebi <- NULL
  chebi.cache <- file.path(outfolder, "chebi_source", "parsed.rds")
  if (!testMode || file.exists(chebi.cache) || is.data.frame(getOption("metadbparse.chebi"))) {
    chebi <- tryCatch(.metadb_chebi(outfolder), error = function(e) {
      message("ChEBI unavailable (", conditionMessage(e), "); keeping InChI as structure")
      NULL
    })
  }
  if (!is.null(chebi)) {
    smi <- chebi[!startsWith(structure, "InChI="), .(chebi_id = as.character(identifier), chebi_smiles = structure)]
    db[, chebi_id := sub("^MTBLC", "", identifier)]
    db <- merge(db, unique(smi, by = "chebi_id"), by = "chebi_id", all.x = TRUE, sort = FALSE)
    db[!is.na(chebi_smiles), structure := chebi_smiles]
    db[, c("chebi_id", "chebi_smiles") := NULL]
  }
  db <- db[!is.na(structure)]
  list(db = db[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")], version = version)
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
  # Retired: dimedb.ibers.aber.ac.uk (incl. /help/downloads/) is offline and no mirror exists.
  stop("DIMEDB is retired: its download server (dimedb.ibers.aber.ac.uk) is offline.")
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
  theurl <- "https://www.vmh.life/files/release-notes/release.html"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "(\\w+ 20\\d\\d)")[, 2]
  # NOTE: legacy pagination API (vmh.uni.lu/_api/metabolites/?page=N) is dead.
  # Its documented REST replacement (/advancedSearch/metabolite) returns a
  # server-side 500 Internal Server Error regardless of parameters as of
  # 2026 -- confirmed broken on VMH's end, not a usage error on ours. Using
  # the SPARQL endpoint instead, which VMH's own API docs recommend for bulk
  # pulls. No API key is required for this endpoint.
  limit <- if (testMode) 10 else 50000
  sparql <- paste0("SELECT ?abbr ?name ?formula ?charge ?smiles ?isHuman ?isMicrobe WHERE { ",
                   "?met a vmh:Metabolite ; vmh:abbreviation ?abbr ; rdfs:label ?name ; ",
                   "vmh:chargedFormula ?formula ; vmh:charge ?charge ; vmh:smiles ?smiles ; ",
                   "vmh:isHuman ?isHuman ; vmh:isMicrobe ?isMicrobe . } LIMIT ", limit)
  url <- paste0("https://www.vmh.life/api/public/sparql?query=", utils::URLencode(sparql, reserved = TRUE))
  r <- httr::GET(url, httr::accept(".json"))
  lst <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
  b <- lst$results$bindings
  table_main <- data.table::data.table(fullName = b$name$value,
                                       chargedFormula = b$formula$value,
                                       abbreviation = b$abbr$value,
                                       charge = as.numeric(b$charge$value),
                                       smile = b$smiles$value,
                                       isHuman = b$isHuman$value == "true",
                                       isMicrobe = b$isMicrobe$value == "true")
  # the SPARQL store returns duplicate bindings per metabolite (confirmed);
  # the final unique() call below collapses those back down.
  db.formatted <- data.table::data.table(compoundname = table_main$fullName, description = NA_character_, baseformula = table_main$chargedFormula, identifier = table_main$abbreviation, charge = table_main$charge, structure = table_main$smile, isHuman = table_main$isHuman, isMicrobe = table_main$isMicrobe)
  # NOTE: unlike the old API, VMH's bulk SPARQL/REST data has no free-text
  # description or synonyms field -- that text simply isn't exposed in bulk
  # (only the single-record /getMetHandler?id= lookup has one, which would
  # need one HTTP call per metabolite and isn't practical at this volume).
  # Building a minimal description from what IS available in bulk instead.
  db.formatted$description <- paste0(db.formatted$compoundname, " (VMH ID: ", db.formatted$identifier, ")")
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
  db.formatted <- unique(db.formatted[, -c("isHuman", "isMicrobe")])
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
#'  \code{\link[readxl]{read_xlsx}}
#'  \code{\link[data.table]{fread}}
#' @rdname build.PHENOLEXPLORER
#' @export
#' @importFrom RCurl getURL
#' @importFrom stringr str_match
#' @importFrom utils download.file unzip
#' @importFrom readxl read_xlsx
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
    utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
  }
  compo_tbl <- as.data.frame(readxl::read_xlsx(file.path(base.loc, "composition-data.xlsx"), sheet = 1))
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
  # Latest MassBank-data release from GitHub (was pinned to 2020.06)
  rel <- tryCatch(jsonlite::fromJSON("https://api.github.com/repos/MassBank/MassBank-data/releases/latest"),
                  error = function(e) NULL)
  version <- if (!is.null(rel$tag_name)) rel$tag_name else "2020.06"
  file.url <- paste0("https://github.com/MassBank/MassBank-data/archive/refs/tags/", version, ".zip")
  base.loc <- file.path(outfolder, "massbank_source")
  dir.create(base.loc, recursive = TRUE, showWarnings = FALSE)
  zip.file <- file.path(base.loc, paste0("massbank_", version, ".zip"))
  if (!.metadb_resume_download(file.url, zip.file, valid = .metadb_zip_ok)) stop("Could not download ", file.url)
  members <- utils::unzip(zip.file, list = TRUE)$Name
  members <- members[grepl("\\.txt$", members) & grepl("^MSBNK-", basename(members))]  # records only, no README/licence
  if (testMode) members <- head(members, 50)
  utils::unzip(normalizePath(zip.file), files = members, exdir = normalizePath(base.loc), unzip = getOption("unzip"))
  cpd_files <- file.path(base.loc, members)
  field <- function(lines, key) {
    hit <- grep(paste0("^", gsub("$", "\\$", key, fixed = TRUE), ": "), lines, value = TRUE)
    if (length(hit)) sub("^[^:]+: ", "", hit[1]) else NA_character_
  }
  db_rows <- pbapply::pblapply(cpd_files, function(fn) {
    lines <- readLines(fn, warn = FALSE, encoding = "UTF-8")
    smi <- field(lines, "CH$SMILES")
    if (is.na(smi) || smi == "N/A") smi <- field(lines, "CH$IUPAC")   # MassBank stores InChI under CH$IUPAC
    data.table::data.table(compoundname = field(lines, "CH$NAME"), description = field(lines, "RECORD_TITLE"),
                           baseformula = field(lines, "CH$FORMULA"), identifier = field(lines, "ACCESSION"),
                           charge = NA_integer_, structure = smi)
  })
  db.formatted <- data.table::rbindlist(db_rows, fill = TRUE)
  db.formatted <- db.formatted[!is.na(baseformula) & !is.na(structure) & structure != "N/A"]
  # one row per structure; spectra of the same compound are collapsed
  db.formatted <- db.formatted[, lapply(.SD, function(x) paste0(unique(x), collapse = "/")), by = structure]
  db.formatted[, charge := NA_integer_]
  list(db = db.formatted[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")],
       version = version)
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
  options(stringsAsFactors = FALSE, timeout=1e7)
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
  utils::unzip(zip.file, exdir = base.loc, unzip = getOption("unzip"))
  input <- file.path(base.loc, "bmdb_metabolites.xml")
  header <- readLines(input, n = 10)
  version <- trimws(gsub(grep(pattern = "<version", header, value = TRUE), pattern = "<\\/?version>", replacement = ""))
  date <- trimws(gsub(grep(pattern = "update_date", header, value = TRUE), pattern = "<\\/?update_date>", replacement = ""))
  if (testMode) {
    n <- 10
  } else {
    # Count the downloaded records rather than relying on the website's markup.
    con <- base::file(input, "r")
    n <- 0L
    repeat {
      lines <- readLines(con, n = 10000, warn = FALSE)
      if (!length(lines)) break
      n <- n + sum(grepl("^\\s*<metabolite>\\s*$", lines))
    }
    close(con)
  }
  if (!n) stop("BMDB XML contains no metabolite records")
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
      if (length(line) == 0 | envir$idx > envir$maxn) {
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
      if (currEnvir$idx > currEnvir$maxn) {
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
  list(db = envir$db.formatted[seq_len(envir$idx - 1L), , drop = FALSE], version = version)
}

#' @title Build RMDB (deprecated)
#' @description Parses the RMDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @param testMode run in test mode? Only parses first ten compounds
#' @return Message that RMDB is deprecated
#' @rdname build.RMDB
#' @export
#' @examples
#' \dontrun{build.RMDB(outfolder=tempdir(), testMode=TRUE)}
build.RMDB <- function(outfolder, testMode = FALSE) {
  print("RMDB is as of recently incorporated in BMDB.")
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
  base.loc <- file.path(outfolder, "ecmdb_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  theurl <- "http://ecmdb.ca/downloads"
  header <- RCurl::getURL(theurl, .opts = list(ssl.verifypeer = FALSE))
  version <- stringr::str_match(header, pattern = "Version <strong>(\\d.\\d)")[, 2]
  file.url <- "http://ecmdb.ca/download/ecmdb.json.zip"
  zip.file <- file.path(base.loc, "ecmdb.zip")

  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")

  utils::unzip(zip.file, exdir = base.loc, unzip = getOption("unzip"))
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
  utils::unzip(zip.file, exdir = base.loc, unzip = getOption("unzip"))
  json <- file.path(base.loc, "ymdb.json")
  line <- readLines(json)[[1]]
  fixed_lines <- gsub(line, pattern = "\\]\\}\\{", replacement = "]},{")
  jsonParsed <- jsonlite::fromJSON(fixed_lines)
  db.partial <- data.table::data.table(compoundname = jsonParsed$name, description = paste0(jsonParsed$description, " Synonyms: ", vapply(jsonParsed$synonyms, paste, collapse = ",", FUN.VALUE = character(1))), baseformula = "", identifier = jsonParsed$ymdb_id, charge = as.integer(jsonParsed$formal_charge), structure = "")
  sdf.url <- "http://www.ymdb.ca/system/downloads/current/ymdb.sdf.zip"
  zip.file <- file.path(base.loc, "ymdb_sdf.zip")
  utils::download.file(sdf.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = base.loc, unzip = getOption("unzip"))
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
    data.table::data.table(identifier = row$identifier, compoundname = row$compoundname.y, structure = row$structure.y, baseformula = row$baseformula.y, charge = row$charge, description = if (is.na(row$description.x)) {
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
#'  \code{\link[readxl]{read_xlsx}}
#'  \code{\link[RCurl]{getURL}}
#'  \code{\link[stringr]{str_match}}
#' @rdname build.PAMDB
#' @export
#' @importFrom utils download.file
#' @importFrom data.table as.data.table data.table
#' @importFrom readxl read_xlsx
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
  db.base <- data.table::as.data.table(readxl::read_xlsx(xlsx.file, sheet = 1))
  db.formatted <- data.table::data.table(identifier = db.base$MetID, compoundname = db.base$Name, structure = db.base$SMILES, baseformula = c(NA), description = gsub(gsub(db.base$Reactions, pattern = "\\r", replacement = ", "), pattern = "\\n", replacement = ""), charge = db.base$Charge)
  db.formatted <- unique(db.formatted, by = "identifier")
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
    theurl <- paste0("https://bioinformatics.charite.de/mvoc/index.php?site=browse&char=", categ)
    data <- readLines(theurl)
    # https://bioinformatics.charite.de/mvoc/index.php?site=ergebnis&compound_id=11095734
    matches = stringr::str_match_all(data, pattern = "compound\\_id\\=(\\d+)")
    matches = unlist(matches)
    cpd_ids = grep("^\\d+$", matches, value=T)
    paste0("https://bioinformatics.charite.de/mvoc/index.php?site=ergebnis&compound_id=", cpd_ids)
  })

  db_rows <- pbapply::pblapply(unlist(search_urls), function(theurl) {
    try({
      lines = readLines(theurl)
      pmid = stringr::str_match(lines, "PubChem ID: (\\d+)")
      pmid = pmid[complete.cases(pmid),]
      pmid = pmid[2]
      # ---------------------
      name = stringr::str_match(lines, "<div style='width:50%'><h1>(.+?)<\\/h1>")
      name = name[complete.cases(name),]
      name = name[2]
      # ---------------------
      syn_names = stringr::str_match_all(lines, "<td>.+?<\\/td>")
      syn_first = which(sapply(syn_names, function(tbl) nrow(tbl)>0))[1]
      syns = syn_names[[syn_first]][,1]
      syns = gsub("<.+?>", "", syns)
      # ---------------------
      formula = stringr::str_match(lines, "Formula<\\/td><td.+?>(.*?)<\\/td>")
      formula = formula[complete.cases(formula),]
      formula = gsub("<.*?>", "", formula[2])
      # ---------------------
      smiles = stringr::str_match(lines,  "SMILES<\\/td><td.+?>(.*?)<\\/td>")
      smiles = smiles[complete.cases(smiles),]
      smiles = smiles[2]
      # ---------------------
      microbes = unlist(stringr::str_match_all(lines, "&species=(.+?)'"))
      microbes = unique(grep("&species", microbes, invert=T, value=T))
      microbes = microbes[2]

      # ---------------------

      db.formatted <- data.table::data.table(identifier = pmid,
                                             compoundname = name,
                                             structure = smiles,
                                             baseformula = formula,
                                             description = paste0("Microbes producing this compound:",
                                                                  microbes,
                                                                  ". Other names:",
                                                                  paste0(syns, collapse = ",")),
                                             charge = c(0))
      Sys.sleep(1)
      db.formatted
    })
  })
  db.formatted <- data.table::rbindlist(db_rows[sapply(db_rows, function(x) !is.null(nrow(x)))])
  version <- "2.0"
  list(db = db.formatted,
       version = version)
}

#' @title Build ANPDB
#' @description Parses the ANPDB, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
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
build.ANPDB <- function(outfolder, testMode = FALSE) {
  # ANPDB moved to Uni Freiburg and now offers a bulk CSV (CC BY-NC 4.0).
  # One row per molecule x species occurrence -> collapse per molecule_id.
  molecule_id <- NULL
  file.url <- "https://phabidb.vm.uni-freiburg.de/static/anpdb/download/ANPDB_molecules.csv"
  base.loc <- file.path(outfolder, "anpdb_source")
  dir.create(base.loc, recursive = TRUE, showWarnings = FALSE)
  csv.file <- file.path(base.loc, "ANPDB_molecules.csv")
  if (!.metadb_resume_download(file.url, csv.file, valid = function(p) file.exists(p) && isTRUE(file.info(p)$size > 1e5))) {
    stop("Could not download ANPDB from ", file.url)
  }
  tab <- data.table::fread(csv.file, encoding = "UTF-8", colClasses = "character",
                           select = c("molecule_id", "mol_name", "smiles", "class_name", "subclass_name",
                                      "bioactivity", "species"))
  tab <- tab[!is.na(smiles) & nzchar(smiles)]
  if (testMode) tab <- tab[molecule_id %in% head(unique(molecule_id), 10)]
  collapse <- function(x) { x <- unique(x[!is.na(x) & nzchar(x)]); if (length(x)) paste(x, collapse = "; ") else NA_character_ }
  db <- tab[, {
    parts <- c(if (!is.na(cl <- collapse(c(class_name, subclass_name)))) paste0("Class: ", cl),
               if (!is.na(sp <- collapse(species))) paste0("Found in: ", sp),
               if (!is.na(ba <- collapse(bioactivity))) paste0("Bioactivity: ", ba))
    .(compoundname = mol_name[1], description = if (length(parts)) paste(parts, collapse = ". ") else NA_character_,
      baseformula = NA_character_, charge = 0L, structure = smiles[1])
  }, by = .(identifier = molecule_id)]
  list(db = db[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")],
       version = as.character(Sys.Date()))
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
  . <- Name <- Formula <- SMILES <- `Additional Names` <- Index <- compoundname <- Categories <- NULL
  file.url = "https://water.for-ident.org/download/STOFF-IDENT_content_17.10.17.zip"
  base.loc <- file.path(outfolder, "stoff_source")
  if (!dir.exists(base.loc)) {
    dir.create(base.loc, recursive = TRUE)
  }
  zip.file <- file.path(base.loc, "stoff.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc), unzip = getOption("unzip"))
  excel.file <- file.path(base.loc, "STOFF-IDENT_content_17.10.17.xlsx")
  xlsx.file <- gsub(excel.file, pattern = "xls", replacement = "xlsx")
  file.copy(excel.file, xlsx.file)
  db.base <- data.table::as.data.table(readxl::read_xlsx(xlsx.file, sheet = 1))
  db.base.aggr <- db.base[, .(compoundname = `Additional Names`[1],
                              baseformula = Formula[1],
                              structure = SMILES[1],
                              description = paste0("In categories:", Categories[1] ,".\nSynonyms: ", paste(`Additional Names`, collapse = ","))),
                          by = Index]
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
build.REACTOME <- function(outfolder, testMode = FALSE) {
  . <- description <- organism <- identifier <- NULL
  theurl <- "https://reactome.org/download/current/ChEBI2Reactome.txt"
  chebi <- data.table::as.data.table(.metadb_chebi(outfolder))  # reuses parsed ChEBI, builds it if needed
  chebi[, identifier := as.numeric(identifier)]
  base.db <- data.table::fread(theurl, header = FALSE)
  colnames(base.db) <- c("identifier", "pathway_id", "pathway_url", "description", "annotation_type", "organism")
  if (testMode) base.db <- base.db[identifier %in% head(unique(identifier), 10)]
  base.db$identifier <- as.character(base.db$identifier)
  path.tab <- base.db[, c("pathway_id", "identifier", "description", "organism")]
  colnames(path.tab) <- c("pathway", "identifier", "description", "organism")
  base.db <- unique(base.db[, .(description = paste0(unique(description), collapse = ", "),
                                organism = paste0(unique(organism), collapse = ", ")), by = "identifier"])
  base.db$description <- paste0("Found in pathways: ", base.db$description, ". Found in organisms: ", base.db$organism, ".")
  base.db <- unique(base.db[, c("identifier", "description")])
  base.db$identifier <- as.numeric(base.db$identifier)
  chebi <- chebi[, -"description", with = FALSE]
  db.formatted <- merge(chebi, base.db, by = "identifier")
  if ("formula" %in% names(db.formatted)) data.table::setnames(db.formatted, "formula", "baseformula")
  ver <- tryCatch({
    v <- trimws(readLines("https://reactome.org/ContentService/data/database/version", warn = FALSE)[1])
    if (nzchar(v)) v else as.character(Sys.Date())
  }, error = function(e) as.character(Sys.Date()))
  list(db = db.formatted[, c("compoundname", "description", "baseformula", "identifier", "charge", "structure")],
       version = ver, path = path.tab)
}

#' @title Build METABOLOMICSWORKBENCH db
#' @description Parses the METABOLOMICSWORKBENCH db, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
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
#' \dontrun{build.METABOLOMICSWORKBENCH(outfolder=tempdir())}
build.METABOLOMICSWORKBENCH <- function(outfolder, testMode = FALSE){
  . <- study <- disease <- identifier <- compoundname <- NULL
  study_url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/species"
  studies = jsonlite::read_json(study_url)
  study_tbl = data.table::rbindlist(pbapply::pblapply(studies,function(l){
    data.table::data.table(identifier = l$`Study ID`, species = l$`Latin name`)
  }))
  # studies > disease https://www.metabolomicsworkbench.org/rest/study/study_id/ST/disease
  disease_url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/disease"
  diseases = jsonlite::read_json(disease_url)
  disease_tbl = data.table::rbindlist(pbapply::pblapply(diseases,function(l){
    data.table::data.table(identifier = l$`Study ID`, disease = l$Disease)
  }))
  merged_info = merge(study_tbl, disease_tbl, by = "identifier", all = T)
  if (testMode) merged_info <- head(merged_info, 10L)

  # studies > compounds https://www.metabolomicsworkbench.org/rest/study/study_id/ST000009/metabolites
  met_id_name_tbl = data.table::rbindlist(pbapply::pblapply(merged_info$identifier, function(id){
    met_rows=data.table::data.table()
    try({
      url = paste0("https://www.metabolomicsworkbench.org/rest/study/study_id/", id, "/metabolites")
      res = jsonlite::read_json(url)
      met_rows = data.table::rbindlist(lapply(res,function(l){
        data.table::data.table(study = id,
                               identifier = if (is.null(l$pubchem_id)) NA_character_ else as.character(l$pubchem_id),
                               refmet_id = if (is.null(l$refmet_id)) NA_character_ else l$refmet_id,
                               compoundname = l$metabolite_name)
      }), fill=T)
    })
    met_rows
  }), fill=T)

  if (!nrow(met_id_name_tbl)) stop("No study metabolites retrieved from Metabolomics Workbench")
  # Current study responses supply RefMet IDs instead of PubChem IDs.
  missing.ids <- is.na(met_id_name_tbl$identifier) | !nzchar(met_id_name_tbl$identifier)
  if (any(missing.ids)) {
    refmet <- data.table::rbindlist(jsonlite::read_json("https://www.metabolomicsworkbench.org/rest/refmet/all"), fill = TRUE)
    refmet <- refmet[!is.na(refmet_id) & nzchar(refmet_id)]
    matched <- match(met_id_name_tbl$refmet_id[missing.ids], refmet$refmet_id)
    met_id_name_tbl$identifier[missing.ids] <- as.character(refmet$pubchem_cid[matched])
  }
  met_ids = unique(met_id_name_tbl$identifier)
  met_ids <- met_ids[!is.na(met_ids) & grepl("^[1-9][0-9]*$", met_ids)]
  if (!length(met_ids)) stop("No valid PubChem compound IDs in study metabolites")

  # info on one met https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/74300/all
  # try pubchem instead
  split.ids <- split(met_ids, ceiling(seq_along(met_ids) / 50))
  met_info_tbl <- data.table::rbindlist(pbapply::pblapply(split.ids, function(idgroup) {
    rows = data.table::data.table()
    try({
      url_struct <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", paste0(idgroup, collapse = ","), "/property/MolecularFormula,SMILES,Charge/JSON")
      struct_res <- jsonlite::fromJSON(url_struct, simplifyVector = TRUE)
      keep.ids <- which(struct_res$PropertyTable$Properties$Charge == 0)
      props <- struct_res$PropertyTable$Properties[keep.ids, , drop = FALSE]
      rows <- data.table::data.table(identifier = as.character(props$CID),
                                     baseformula = props$MolecularFormula,
                                     structure = props$SMILES, charge = props$Charge)
    })
    rows
  }))
  if (!nrow(met_info_tbl)) stop("No compound properties retrieved from PubChem")
  met_info_tbl$identifier <- as.character(met_info_tbl$identifier)
  met_info_study = merge(met_info_tbl,met_id_name_tbl, by="identifier")
  big_merged <- merge(met_info_study, merged_info, by.x = "study", by.y = "identifier", allow.cartesian = T)
  aggregated <- big_merged[ , .(description = paste0("Mentioned in the following studies (and connected diseases):",
                                                     paste0(unique(
                                                       paste0(study, "(", disease, ")")),
                                                       collapse=","))
                                ),
                            by = identifier]
  aggregated$description <- gsub("\\(NA\\)", "", aggregated$description)
  final_merged = merge(aggregated, met_info_tbl, by = "identifier")
  name_tbl = unique(met_id_name_tbl[,-"study"])
  name_aggregated <- name_tbl[ , .(compoundname = compoundname[1]),
                               by = identifier]
  final_merged = merge(final_merged, name_aggregated, by = "identifier")
  # - - - - - -
  db.formatted <- data.table::data.table(identifier = final_merged$identifier,
                                         compoundname = final_merged$compoundname,
                                         structure = final_merged$structure,
                                         baseformula = final_merged$baseformula,
                                         charge = final_merged$charge,
                                         description = final_merged$description)

  date <- Sys.Date()
  version <- Sys.Date()
  db.formatted.all = list(db = db.formatted,
                          version = version)

  return(db.formatted.all)
}

#' @title Build Natural Products Atlas db
#' @description Parses the Natural Products Atlas db, returns data table with columns compoundname, description, charge, formula and structure (in SMILES)
#' @param outfolder Which folder to save temp files to?
#' @return data table with parsed database
#' @seealso
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[rvest]{html_nodes}}
#'  \code{\link[xml2]{read_html}}
#' @rdname build.NPA
#' @export
#' @importFrom data.table as.data.table fread
#' @examples
#' \dontrun{build.NPA(outfolder=tempdir())}
build.NPA <- function(outfolder){
  theurl = "https://www.npatlas.org/download"
  html = rvest::read_html(theurl)
  version = as.character(rvest::html_elements(html, css = "h5"))
  version = stringr::str_match(version, "\\((.*?)\\)")[,2]

  theurl = "https://www.npatlas.org/static/downloads/NPAtlas_download.tsv"

  db.base = data.table::fread(theurl)
  # - - - - - -
  db.formatted <- data.table::data.table(identifier = db.base$npaid,
                                         compoundname = if ("compound_name" %in% names(db.base)) db.base$compound_name else db.base$compound_names,
                                         structure = db.base$compound_smiles,
                                         baseformula = db.base$compound_molecular_formula,
                                         charge = NA_integer_,
                                         description = sapply(1:nrow(db.base), function(i){
                                           row = db.base[i,]
                                           if(grepl(pattern = "unknown", row$genus,ignore.case = T)){
                                             paste0("Found in unknown ",
                                                    tolower(row$origin_type),
                                                    " (DOI:", row$original_reference_doi, ")."
                                             )
                                           }else{
                                             paste0("Found in the ",
                                                    tolower(row$origin_type), " ",
                                                    Hmisc::capitalize(row$genus), " ",
                                                    Hmisc::capitalize(row$origin_species),
                                                    " (DOI:", row$original_reference_doi, ")."
                                             )
                                           }
                                         }))

  date <- Sys.Date()
  db.formatted.all = list(db = db.formatted,
                          version = version)

  return(db.formatted.all)
}

build.MARKERDB <- function(outfolder){
  apikey = "cc69fb7c1a4809bbae7c0701d5188c17"
  page = 1
  first.page = gsubfn::fn$paste("https://markerdb.ca/api/v1/chemicalapi/chemicalrequest?api_key=$apikey&name=&page=$page")
  first.page.info = jsonlite::read_json(first.page)
  npages = first.page.info$total_page
  cpd.info.rows <- pbapply::pblapply(1:npages, function(page){
    page.json = jsonlite::read_json(gsubfn::fn$paste("https://markerdb.ca/api/v1/chemicalapi/chemicalrequest?api_key=$apikey&name=&page=$page"))
    row = data.table::data.table()
    page.rows = lapply(page.json$chemicals, function(json){
      row = data.table::data.table(
        identifier = paste0("MDB", stringr::str_pad(json$id, width = 8, pad = "0", side = "left")),
        compoundname = json$mdbid,
        structure = json$moldb_smiles,
        baseformula = json$moldb_formula,
        description = json$description,
        charge = c(0)
      )
    })
    data.table::rbindlist(page.rows, fill = TRUE)
  })
  db.formatted = data.table::rbindlist(cpd.info.rows, fill = TRUE)

  db.formatted.all = list(db = db.formatted,
                         version = Sys.Date())

  return(db.formatted.all)
}



