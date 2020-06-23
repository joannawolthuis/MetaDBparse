# === MOVE TO db_search.R ===

showAllBase <- function(outfolder, base.dbname){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(base.dbname,".db"))) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname as name, baseformula as formula, structure as structure, description as description, charge as charge FROM base")
  RSQLite::dbDisconnect(conn)
  # --- result ---
  result
}


searchMZ <- function(mzs, ionmodes, outfolder,
                     base.dbname, ppm,
                     ext.dbname="extended",
                     append=F){

  # connect to extended DB
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(ext.dbname, ".db")))

  # create temp range tables for chosen mz (mimic mzvals, mzrangea
  mzvals = data.table::data.table(mzmed = as.numeric(mzs),
                                  foundinmode=ionmodes)
  eachPPM <- length(ppm) > 1

  mzranges = data.table::rbindlist(lapply(1:length(mzs), function(i){
    mz = as.numeric(mzs[i])
    if(eachPPM) ppm <- as.numeric(ppm[i])
    data.table::data.table(mzmin = mz - (mz*(ppm*1e-6)),
                           mzmax = mz + (mz*(ppm*1e-6)))
  }))

  # create tables
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS mzvals")
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, sql.make.meta)
  RSQLite::dbExecute(conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  RSQLite::dbWriteTable(conn, "mzvals", mzvals, append=TRUE) # insert into

  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS mzranges")
  sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges USING rtree(
                            ID INTEGER PRIMARY KEY AUTOINCREMENT,
                            mzmin decimal(30,13),
                            mzmax decimal(30,13));",width=10000, simplify=TRUE)

  RSQLite::dbExecute(conn, sql.make.rtree)
  RSQLite::dbWriteTable(conn, "mzranges", mzranges, append=TRUE) # insert into

  # query
  RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS unfiltered')

  query = "CREATE TABLE unfiltered AS
                            SELECT DISTINCT
                            cpd.adduct as adduct,
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

  table.per.db <- lapply(base.dbname, function(db){
    dbpath = file.path(outfolder, paste0(db, ".db"))
    try({
      DBI::dbExecute(conn, gsubfn::fn$paste("DETACH base"))
    },silent=T)
    query = gsubfn::fn$paste("ATTACH '$dbpath' AS base")
    RSQLite::dbExecute(conn, query)

    query = strwrap("SELECT
                  b.compoundname as name,
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
                  ON u.structure = b.structure",width=10000, simplify=TRUE)

    results <- DBI::dbGetQuery(conn, query)
    if(nrow(results)>0){
      results$perciso <- round(results$perciso, 2)
      results$dppm <- signif(results$dppm, 2)
      results$source <- c(db)
      colnames(results)[which(colnames(results) == "perciso")] <- "%iso"
    }else{
      results = data.table::data.table()
    }
    # return
    results
  })

  merged.results <- data.table::rbindlist(table.per.db)

  DBI::dbDisconnect(conn)

  return(merged.results)
}

searchFormula <- function(formula, outfolder, base.dbname){
  table.per.db <- lapply(base.dbname, function(db){
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(db,".db"))) # change this to proper var later

    RSQLite::dbExecute(conn, "CREATE TEMP TABLE formulas(formula text,
                                                charge TEXT)")
    RSQLite::dbWriteTable(conn, data.table::data.table(formula = formula,
                                                       charge = charge))
    res = RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname as name,
                                                     baseformula as formula,
                                                     identifier,
                                                     description,
                                                     structure
                                                     FROM base
                                                     JOIN formulas
                                                     ON base.baseformula = formulas.baseformula
                                                     AND base.basecharge = formulas.charge")

    res
  })

  merged.results <- data.table::rbindlist(table.per.db)

  return(merged.results)
}

searchRev <- function(structure, ext.dbname="extended", outfolder){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder,
                                                          paste0(ext.dbname,".db"))) # change this to proper var later
  result = RSQLite::dbSendStatement(conn,
                                    "SELECT DISTINCT fullmz, fullformula, finalcharge, adduct, isoprevalence
                                    FROM extended ext
                                    JOIN structures struct
                                    ON ext.struct_id = struct.struct_id
                                    WHERE struct.smiles = $structure")
  RSQLite::dbBind(result, list(structure = structure))
  res = RSQLite::dbFetch(result)
  RSQLite::dbClearResult(result)
  RSQLite::dbDisconnect(conn)
  # - - -
  res
}

searchCMMR <- function (cmm_url = "http://ceumass.eps.uspceu.es/mediator/api/v3/batch",
                        metabolites_type = "all-except-peptides", databases = "[\"all-except-mine\"]",
                        masses_mode = "mz", ion_mode = "positive", adducts = switch(ion_mode,
                                                                                    positive = "[\"M+H\", \"M+2H\", \"M+Na\", \"M+K\", \"M+NH4\", \"M+H-H2O\"]",
                                                                                    negative = "[\"M-H\", \"M+Cl\", \"M+FA-H\", \"M-H-H2O\"]"),
                        tolerance = 10, tolerance_mode = "ppm", unique_mz)
  {
  body <- cmmr::create_batch_body(metabolites_type, databases, masses_mode,
                            ion_mode, adducts, tolerance, tolerance_mode, unique_mz)
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
    json_file <- RJSONIO::fromJSON(httr::content(r, "text",
                                                 encoding = "UTF-8"))$results
    if (length(json_file) == 0) {
      cat("No compounds found in the database search.\n")
      return("No compounds found in the database search.")
    }
    pb <- progress::progress_bar$new(format = "  Parsing database search results [:bar] :percent in :elapsed",
                                     total = length(json_file) - 1, clear = FALSE, width = 100)
    df <- data.table::rbindlist(json_file,fill=T)
    return(df)
  }
  else {
    cat(paste0("Status: ", r$status_code, ", Fail to connect the API service!\n"))
    cat(paste0("Date: ", r$date, "\n"))
  }
}

searchMZonline <- function(mz=178.1219,
                           mode="positive",
                           adducts,
                           ppm=2,
                           which_db = "cmmr",
                           apikey = NULL){
  switch(which_db,
         cmmr={
           library(cmmr)
           results <- searchCMMR('http://ceumass.eps.uspceu.es/mediator/api/v3/batch',
                                 masses_mode = 'mz',
                                 ion_mode = 'positive',
                                 tolerance = ppm,
                                 tolerance_mode = 'ppm',
                                 unique_mz = mz)
            if(typeof(results) == "character"){
              data.table()
            }else{
              base.table = data.table::data.table(name = results$name,
                           baseformula = results$formula,
                           adduct = results$adduct,
                           `%iso` = c(100),
                           fullformula = c(NA),
                           finalcharge = c(NA),
                           dppm = results$error_ppm,
                           identifier = results$identifier,
                           description = c(""),
                           structure = results$inChiKey,
                           query_mz = results$EM)
              has.inchi <- which(!is.na(base.table$structure))
              if(length(has.inchi) >0){
                inchis <- base.table$structure[has.inchi]
                if(!is.null(apikey)){
                  inchi = pbapply::pbsapply(inchis, function(x) webchem::cs_convert(x, from = "inchikey", to = "inchi", apikey = apikey))
                  smiles = pbapply::pbsapply(inchi, function(x) webchem::cs_convert(x, from = "inchi", to = "smiles", apikey = apikey))
                  base.table$structure[has.inchi] <- smiles
                  base.table$structure[!has.inchi] <- c("")
                }
              }else{
                base.table$structure <- c("")
              }
              return(base.table)
            }
         })
}

