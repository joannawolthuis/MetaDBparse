# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# define database

openBaseDB <- function(outfolder, dbname){
  if(!dir.exists(outfolder)) dir.create(outfolder)
  db <- file.path(normalizePath(outfolder), paste0(dbname))
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(conn, statement = "CREATE TABLE IF NOT EXISTS base(compoundname text,
                                                                    description text,
                                                                    baseformula text,
                                                                    identifier text,
                                                                    charge text,
                                                                    structure text)")
  #DBI::dbDisconnect(conn)
  conn
}

removeBaseDB <- function(outfolder, dbname){
  db <- file.path(normalizePath(outfolder), paste0(dbname))
  print(db)
  if(file.exists(db)) file.remove(db)
}

writeDB <- function(conn, table, tblname){
  DBI::dbWriteTable(conn, tblname, table, append=T)
}

makeExtDB <- function(outfolder){
  outfolder <- normalizePath(outfolder)
  data(isotopes, package = "enviPat")
  base.db <- file.path(outfolder, paste0(dbname, ".base.db"))
  full.db <<- file.path(outfolder, paste0("extended.db"))

  first.db = if(!file.exists(full.db)) T else F
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)

  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA foreign_keys = ON"))
}

removeExtDB <- function(outfolder){
  outfolder <- normalizePath(outfolder)
  db <- file.path(outfolder, paste0("extended.db"))
  if(file.exists(db)) file.remove(db)
}
