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

removeDB <- function(outfolder, dbname){
  db <- file.path(normalizePath(outfolder), paste0(dbname))
  if(file.exists(db)) file.remove(db)
}

writeDB <- function(conn, table, tblname){
  DBI::dbWriteTable(conn, tblname, table, append=T)
}

