# === MOVE TO db_search.R ===

showAllBase <- function(outfolder, base.dbname){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(base.dbname,".db"))) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname as name, baseformula as formula, structure as structure, description as description, charge as charge FROM base")
  RSQLite::dbDisconnect(conn)
  # --- result ---
  result
}


searchMZ <- function(mzs, ionmodes, outfolder, base.dbname, ppm, ext.dbname="extended", append=F){

  # connect to extended DB
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(ext.dbname, ".db")))

  # create temp range tables for chosen mz (mimic mzvals, mzrangea
  mzvals = data.table::data.table(mzmed = as.numeric(mzs),
                                  foundinmode=ionmodes)
  mzranges = data.table::data.table(mzmin=sapply(as.numeric(mzs), function(mz) mz - (mz*(ppm*1e-6))),
                                    mzmax=sapply(as.numeric(mzs), function(mz) mz + (mz*(ppm*1e-6))))

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
  query <- gsubfn::fn$paste(strwrap(
    "CREATE TABLE unfiltered AS
    SELECT DISTINCT
    cpd.adduct as adduct,
    cpd.isoprevalence as isoprevalence,
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
    ON cpd.struct_id = struc.struct_id",width=10000, simplify=TRUE))
  RSQLite::dbExecute(conn, query)


  table.per.db <- lapply(base.dbname, function(db){
    dbpath = file.path(outfolder, paste0(db, ".db"))
    try({
      DBI::dbExecute(conn, gsubfn::fn$paste("DETACH base"))
    },silent=T)
    query = gsubfn::fn$paste("ATTACH '$dbpath' AS base")
    RSQLite::dbExecute(conn, query)
    results <- DBI::dbGetQuery(conn, strwrap("SELECT
                                             b.compoundname as name,
                                             b.baseformula,
                                             u.adduct,
                                             u.isoprevalence as perciso,
                                             u.dppm,
                                             b.identifier,
                                             b.description,
                                             u.structure,
                                             u.query_mz
                                             FROM unfiltered u
                                             JOIN base.base b
                                             ON u.structure = b.structure",width=10000, simplify=TRUE))
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

searchRev <- function(structure, ext.dbname="extended"){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder,
                                                          paste0(ext.dbname,".db"))) # change this to proper var later
  result = RSQLite::dbSendStatement(conn,
                                    "SELECT DISTINCT fullmz, adduct, isoprevalence
                                    FROM extended ext
                                    JOIN structures struct
                                    ON ext.struct_id = struct.struct_id
                                    WHERE struct.smiles = $structure")
  RSQLite::dbBind(result, list(structure = structure))
  res = RSQLite::dbFetch(result)
  RSQLite::dbClearResult(result)
  RSQLite::dbDisconnect(conn)
}
