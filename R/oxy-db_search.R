#' @title Show all compounds in base db
#' @description Shows all compounds in this base database as data table.
#' @param outfolder Which folder is your database in?
#' @param base.dbname Which base database do you want to explore? (exclude .db suffix)
#' @return Data table with whole database.
#' @details This may be quite memory consuming for larger databases!!
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#' @rdname showAllBase
#' @export
#' @importFrom RSQLite dbConnect SQLite dbGetQuery dbDisconnect
#' @examples
#'  \dontrun{myFolder = tempdir()}
#'  \dontrun{buildBaseDB(outfolder = myFolder, "lmdb", test = TRUE)}
#'  \dontrun{showAllBase(outfolder = myFolder, "lmdb")}
showAllBase <- function(outfolder, base.dbname) {
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(base.dbname, ".db")))
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT identifier,compoundname, baseformula as formula, structure as structure, description as description, charge as charge FROM base")
  RSQLite::dbDisconnect(conn)
  result
}

#' @title Find matches for m/z value in given database
#' @description This function takes user m/z, ppm error, base database and the extended database to return hits of interest.
#' @param mzs Vector of m/z values
#' @param ionmodes Vector of pos/negative mode for each m/z value
#' @param outfolder Which folder are your databases stored in?
#' @param base.dbname Which base database do you want to retrieve info from? (without .db suffix)
#' @param ppm Parts per million accepted error range
#' @param ext.dbname Name of extended database (without .db suffix), Default: 'extended'
#' @param append Use this when searching muiltiple base databases, so only one result table is created, Default: FALSE
#' @return Data table with match results
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[DBI]{dbExecute}},\code{\link[DBI]{dbGetQuery}},\code{\link[DBI]{dbDisconnect}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname searchMZ
#' @export
#' @importFrom RSQLite dbConnect SQLite dbExecute dbWriteTable
#' @importFrom data.table data.table rbindlist
#' @importFrom DBI dbExecute dbGetQuery dbDisconnect
#' @importFrom gsubfn fn
#' @examples
#'  \dontrun{myFolder = tempdir()}
#'  \dontrun{buildBaseDB(outfolder = myFolder, "lmdb", test = TRUE)}
#'  \dontrun{file.remove(file.path(myFolder, "extended.db"))}
#'  \dontrun{data(adducts)}
#'  \dontrun{data(adduct_rules)}
#'  \dontrun{buildExtDB(outfolder = myFolder, base.dbname = "lmdb",
#'  silent=FALSE, adduct_table = adducts, adduct_rules = adduct_rules)}
#'  \dontrun{searchMZ(c("104.3519421"), "positive", outfolder = myFolder, "lmdb", ppm = 3)}
searchMZ <- function(mzs, ionmodes, outfolder, base.dbname, ppm, ext.dbname = "extended", append = FALSE) {
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(ext.dbname, ".db")))
  mzvals <- data.table::data.table(mzmed = as.numeric(mzs), foundinmode = ionmodes)
  eachPPM <- length(ppm) > 1
  mzranges <- data.table::rbindlist(lapply(1:length(mzs), function(i) {
    mz <- as.numeric(mzs[i])
    if (eachPPM) {
      ppm <- as.numeric(ppm[i])
    }
    data.table::data.table(mzmin = mz - (mz * (ppm * 1e-06)), mzmax = mz + (mz * (ppm * 1e-06)))
  }))
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS mzvals")
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width = 10000, simplify = TRUE)
  RSQLite::dbExecute(conn, sql.make.meta)
  RSQLite::dbExecute(conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  RSQLite::dbWriteTable(conn, "mzvals", mzvals, append = TRUE)
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS mzranges")
  sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges USING rtree(ID INTEGER PRIMARY KEY AUTOINCREMENT,
                                                                       mzmin decimal(30,13),
                                                                       mzmax decimal(30,13));", width = 10000, simplify = TRUE)
  RSQLite::dbExecute(conn, sql.make.rtree)
  RSQLite::dbWriteTable(conn, "mzranges", mzranges, append = TRUE)
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS unfiltered")
  query <- "CREATE TABLE unfiltered AS SELECT DISTINCT cpd.adduct as adduct,
                                                       cpd.isoprevalence as isoprevalence,
                                                       cpd.fullformula,
                                                       cpd.finalcharge,
                                                       struc.smiles as structure,
                                                       mz.mzmed as query_mz,
                                                       (1e6*ABS(mz.mzmed - cpd.fullmz)/cpd.fullmz) AS dppm
                                                       FROM mzvals mz
                                                       JOIN mzranges rng ON rng.ID = mz.ID
                                                       JOIN extended cpd indexed by e_idx2
                                                       ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
                                                       JOIN adducts
                                                       ON cpd.adduct = adducts.Name
                                                       AND mz.foundinmode = adducts.Ion_Mode
                                                       JOIN structures struc
                                                       ON cpd.struct_id = struc.struct_id"

  RSQLite::dbExecute(conn, query)
  table.per.db <- lapply(base.dbname, function(db) {
    dbpath <- file.path(outfolder, paste0(db, ".db"))
    try(
      {
        DBI::dbExecute(conn, gsubfn::fn$paste("DETACH base"))
      },
      silent = TRUE
    )
    query <- gsubfn::fn$paste("ATTACH '$dbpath' AS base")
    RSQLite::dbExecute(conn, query)
    query <- strwrap("SELECT b.compoundname,
                             b.baseformula,
                             u.adduct,
                             u.isoprevalence as perciso,
                             u.fullformula,
                             u.finalcharge,
                             u.dppm,
                             b.identifier,
                             b.description,
                             u.structure,
                             u.query_mz
                             FROM unfiltered u
                             JOIN base.base b
                             ON u.structure = b.structure",
                     width = 10000, simplify = TRUE
    )
    results <- DBI::dbGetQuery(conn, query)
    if (nrow(results) > 0) {
      results$perciso <- round(results$perciso, 2)
      results$dppm <- signif(results$dppm, 2)
      results$source <- c(db)
      colnames(results)[which(colnames(results) == "perciso")] <- "%iso"
    }
    else {
      results <- data.table::data.table()
    }
    results
  })
  merged.results <- data.table::rbindlist(table.per.db)
  DBI::dbDisconnect(conn)
  return(merged.results)
}

#' @title Find matches based on molecular formula
#' @description Goes through database of choice (base database) and retrieves hits that have the  molecular formula of interest.
#' @param formula Molecular formula (should be checked by enviPat::check_chemform first!)
#' @param charge Charge of formula
#' @param outfolder Which folder are your databases stored in?
#' @param base.dbname Base database name (without .db suffix)
#' @return Data table with compounds with this molecular formula and the other available information
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname searchFormula
#' @export
#' @importFrom RSQLite dbConnect SQLite dbExecute dbWriteTable dbGetQuery
#' @importFrom data.table data.table rbindlist
#' @examples
#'  \dontrun{myFolder = tempdir()}
#'  \dontrun{buildBaseDB(outfolder = myFolder, "lmdb", test = TRUE)}
#'  \dontrun{searchFormula(formula = c("C7H11N3O2"), charge = 0,
#'  outfolder = myFolder, base.dbname = c("lmdb"))}
searchFormula <- function(formula, charge, outfolder, base.dbname) {
  table.per.db <- lapply(base.dbname, function(db) {
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(db, ".db")))
    RSQLite::dbExecute(conn, "CREATE TEMP TABLE formulas(formula TEXT,
                                                         charge INTEGER)")
    formTbl = data.table::data.table(formula = formula,
                                     charge = charge)
    print(formTbl)
    RSQLite::dbWriteTable(conn = conn,
                          name = "formulas",
                          value = formTbl,
                          append = TRUE)
    res <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname,
                                      baseformula as formula,
                                      identifier,
                                      description,
                                      structure
                                      FROM base
                                      JOIN formulas
                                      ON base.baseformula = formulas.formula
                                      AND base.charge = formulas.charge")
    res
  })
  merged.results <- data.table::rbindlist(table.per.db)
  return(merged.results)
}

#' @title Reverse searching
#' @description Takes a SMILES structure and finds m/z values for all adducts and isotopes matching that structure.
#' @param structure SMILES structure string
#' @param ext.dbname Name of extended database (without .db), Default: 'extended'
#' @param outfolder Which folder are your databases in?
#' @return Data table with m/z values, additionally molecular formula, charge, adduct, isotope %.
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#' @rdname searchRev
#' @export
#' @importFrom RSQLite dbConnect SQLite dbSendStatement dbBind dbFetch dbClearResult dbDisconnect
#' @examples
#'  \dontrun{myFolder = tempdir()}
#'  \dontrun{buildBaseDB(outfolder = myFolder, "lmdb", test = TRUE)}
#'  \dontrun{file.remove(file.path(myFolder, "extended.db"))}
#'  \dontrun{data(adducts)}
#'  \dontrun{data(adduct_rules)}
#'  \dontrun{buildExtDB(outfolder = myFolder, base.dbname = "lmdb",
#'  silent=FALSE, adduct_table = adducts, adduct_rules = adduct_rules)}
#'  \dontrun{searchRev("O=C(O)C(N)CC=1N=CN(C1)C", outfolder = myFolder)}
searchRev <- function(structure, ext.dbname = "extended", outfolder) {
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(ext.dbname, ".db")))
  result <- RSQLite::dbSendStatement(conn, "SELECT DISTINCT fullmz, fullformula, finalcharge, adduct, isoprevalence
                                            FROM extended ext
                                            JOIN structures struct
                                            ON ext.struct_id = struct.struct_id
                                            WHERE struct.smiles = $structure")
  RSQLite::dbBind(result, list(structure = structure))
  res <- RSQLite::dbFetch(result)
  RSQLite::dbClearResult(result)
  RSQLite::dbDisconnect(conn)
  res
}

#' @title Search CMMR
#' @description Queries Ceu Mass Mediator through their API
#' @param cmm_url API base url, Default: 'http://ceumass.eps.uspceu.es/mediator/api/v3/batch'
#' @param metabolites_type Which type of metabolites to consider?, Default: 'all-except-peptides'
#' @param databases Which databases to consider?, Default: '["all-except-mine"]'
#' @param masses_mode Format of input compound, Default: 'mz'
#' @param ion_mode Which ion mode was the compound found in?, Default: 'positive'
#' @param adducts Adducts to be considered, Default: switch(ion_mode, positive = "[\"M+H\", \"M+2H\", \"M+Na\", \"M+K\", \"M+NH4\", \"M+H-H2O\"]",
#'    negative = "[\"M-H\", \"M+Cl\", \"M+FA-H\", \"M-H-H2O\"]")
#' @param tolerance Error margin, units of 'tolerance_mode', Default: 10
#' @param tolerance_mode Mode of error margin, Default: 'ppm'
#' @param unique_mz M/z(s) to use in query
#' @return Data table with match results
#' @seealso
#'  \code{\link[cmmr]{create_batch_body}}
#'  \code{\link[httr]{POST}},\code{\link[httr]{content_type}},\code{\link[httr]{content}}
#'  \code{\link[RJSONIO]{fromJSON}}
#'  \code{\link[progress]{progress_bar}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname searchCMMR
#' @export
#' @importFrom cmmr create_batch_body
#' @importFrom httr POST content_type content
#' @importFrom RJSONIO fromJSON
#' @importFrom progress progress_bar
#' @importFrom data.table rbindlist
#' @examples
#' searchCMMR(unique_mz = "170.09240307", ion_mode = "positive")
searchCMMR <- function(cmm_url = "http://ceumass.eps.uspceu.es/mediator/api/v3/batch", metabolites_type = "all-except-peptides", databases = "[\"all-except-mine\"]", masses_mode = "mz", ion_mode = "positive", adducts = switch(ion_mode, positive = "[\"M+H\", \"M+2H\", \"M+Na\", \"M+K\", \"M+NH4\", \"M+H-H2O\"]", negative = "[\"M-H\", \"M+Cl\", \"M+FA-H\", \"M-H-H2O\"]"), tolerance = 10, tolerance_mode = "ppm", unique_mz) {
  body <- cmmr::create_batch_body(metabolites_type, databases, masses_mode, ion_mode, adducts, tolerance, tolerance_mode, unique_mz)
  if (cmm_url == "http://ceumass.eps.uspceu.es/mediator/api/v3/batch") {
    cat("Using the CEU Mass Mediator server API.\n")
  }
  else {
    cat("Using the local/3rd party server API.\n")
    cat(paste0(cmm_url, "\n"))
  }
  r <- httr::POST(url = cmm_url, body = body, httr::content_type("application/json"))
  if (r$status_code == 200) {
    cat(paste0("Status: ", r$status_code, ", Success!\n"))
    cat(paste0("Date: ", r$date, "\n"))
    json_file <- RJSONIO::fromJSON(httr::content(r, "text", encoding = "UTF-8"))$results
    if (length(json_file) == 0) {
      cat("No compounds found in the database search.\n")
      return("No compounds found in the database search.")
    }
    pb <- progress::progress_bar$new(format = "  Parsing database search results [:bar] :percent in :elapsed", total = length(json_file) - 1, clear = FALSE, width = 100)
    df <- data.table::rbindlist(json_file, fill = TRUE)
    return(df)
  }
  else {
    cat(paste0("Status: ", r$status_code, ", Fail to connect the API service!\n"))
    cat(paste0("Date: ", r$date, "\n"))
  }
}

#' @title Find m/z matches with CMMR, ChemSpider or PubChem
#' @description Wrapper function for all online searches.
#' @param mz M/z of interest, Default: 178.1219
#' @param mode Is m/z positive or negative mode?, Default: 'positive'
#' @param adducts Which adducts will you consider (for cmmr only)
#' @param ppm Allowed error margin in parts per million, Default: 2
#' @param which_db Which online database do you want to search?, Default: 'cmmr'
#' @param apikey ChemSpider API key. Only required if searching in ChemSpider.
#' @return Table with match information
#' @seealso
#'  \code{\link[pbapply]{pbapply}}
#' @rdname searchMZonline
#' @export
#' @importFrom data.table data.table
#' @importFrom pbapply pbsapply
#' @importFrom webchem cs_inchikey_inchi cs_inchi_smiles
#' @examples
#' \dontrun{searchMZonline(mz = 170.09240307, mode = "positive", which_db = "cmmr")}
searchMZonline <- function(mz = 178.1219, mode = "positive", adducts, ppm = 2, which_db = "cmmr", apikey = "") {
  switch(which_db, cmmr = {
    results <- searchCMMR("http://ceumass.eps.uspceu.es/mediator/api/v3/batch", masses_mode = "mz", ion_mode = "positive", tolerance = ppm, tolerance_mode = "ppm", unique_mz = mz)
    if (typeof(results) == "character") {
      data.table::data.table()
    } else {
      base.table <- data.table::data.table(compoundname = results$name, baseformula = results$formula, adduct = results$adduct, `%iso` = c(100), fullformula = c(NA), finalcharge = c(NA), dppm = results$error_ppm, identifier = results$identifier, description = c(""), structure = results$inChiKey, query_mz = results$EM)
      has.inchi <- which(!is.na(base.table$structure))
      if (length(has.inchi) > 0) {
        inchis <- base.table$structure[has.inchi]
        smiles <- pbapply::pbsapply(inchis, function(x) {
          url <- gsubfn::fn$paste("https://cactus.nci.nih.gov/chemical/structure/$x/smiles")
          Sys.sleep(0.1)
          RCurl::getURL(url)
        })
        valid.smiles <- !grepl("Page not found", smiles)
        has.newline <- grep("\n", smiles)
        smiles[has.newline] <- sapply(stringr::str_split(smiles[has.newline], "\n"), function(x) x[[1]])
        base.table$structure[!has.inchi] <- c("")
        base.table$structure[has.inchi[valid.smiles]] <- smiles[valid.smiles]
      } else {
        base.table$structure <- c("")
      }
      return(base.table)
    }
  })
}
