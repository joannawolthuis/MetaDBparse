# === MOVE TO db_search.R ===

showAllBase <- function(outfolder, base.dbname){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(base.dbname,".db"))) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT * FROM base")
  # --- result ---
  result
}

searchMZ <- function(mzs, ionmodes, outfolder, base.dbname, ppm, append=F){

  # connect to extended DB
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, "extended.db"))

  # join base db
  base.db = file.path(outfolder, paste0(base.dbname, ".db"))
  query = gsubfn::fn$paste("ATTACH '$base.db' AS base")
  RSQLite::dbExecute(conn, query)

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

  if(!DBI::dbExistsTable(conn, "unfiltered") | !append)
  {
    RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS unfiltered')
    query <- gsubfn::fn$paste(strwrap(
      "CREATE TABLE unfiltered AS
          SELECT DISTINCT
          cpd.fullformula,
          cpd.adduct as adduct,
          cpd.isoprevalence as isoprevalence,
          struc.smiles as structure,
          mz.mzmed as query_mz,
          (ABS(mz.mzmed - cpd.fullmz) / mz.mzmed)/1e6 AS dppm
          FROM mzvals mz
          JOIN mzranges rng ON rng.ID = mz.ID
          JOIN extended cpd indexed by e_idx2
          ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
          AND mz.foundinmode = cpd.foundinmode
          JOIN structures struc
          ON cpd.struct_id = struc.struct_id",width=10000, simplify=TRUE))
    RSQLite::dbExecute(conn, query)
  }

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

  results$perciso <- round(results$perciso, 2)
  results$dppm <- signif(results$dppm, 2)

  colnames(results)[which(colnames(results) == "perciso")] <- "%iso"

  DBI::dbDisconnect(conn)

  results
}

searchFormula <- function(formula, outfolder, base.dbname){
  query =
    gsubfn::fn$paste(strwrap(
      "SELECT DISTINCT compoundname as name,
    baseformula as formula,
    identifier,
    description,
    structure
     FROM base
     WHERE baseformula = '$formula'"
      , width=10000, simplify=TRUE))

  conn <- RSQLite::dbConnect(RSQLite::SQLite(), file.path(outfolder, paste0(base.dbname,".db"))) # change this to proper var later
  res <- RSQLite::dbGetQuery(conn, query)

  return(res)
}
