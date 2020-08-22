#' @title Create/open and prepare SQLite database
#' @description Create/open and prepare SQLite database
#' @param outfolder Which folder are you building your databases in?
#' @param dbname What is the name of the database? (exclude .db)
#' @return Nothing, writes SQLITE database to outfolder
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[DBI]{dbExecute}}
#' @rdname openBaseDB
#' @export
#' @importFrom RSQLite dbConnect SQLite
#' @importFrom DBI dbExecute
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
  conn
}

#' @title Remove database
#' @description Removes database from disk completely.
#' @param outfolder Which folder are you building your databases in?
#' @param dbname What database do you want to remove? (exclude .db suffix)
#' @return Nothing, removes database from disk
#' @rdname removeDB
#' @export
removeDB <- function(outfolder, dbname){
  db <- file.path(normalizePath(outfolder), paste0(dbname))
  if(file.exists(db)) unlink(db)
}

#' @title Write table to SQLite database
#' @description Simple wrapper - take data table or data frame and append it to the given table in the given database.
#' @param conn Database connection (from DBI::dbConnect or similar)
#' @param table Data frame or data table
#' @param tblname SQLITE table name to append to
#' @return Nothing, writes table to SQLITE database
#' @seealso
#'  \code{\link[DBI]{dbWriteTable}}
#' @rdname writeDB
#' @export
#' @importFrom DBI dbWriteTable
writeDB <- function(conn, table, tblname){
  DBI::dbWriteTable(conn, tblname, table, append=TRUE)
}

#' @title Adjusted sdfStream version for databases that store their compounds in SDF format
#' @description Adjusted from existing function to extract the columns needed for MetaDBparse database format
#' @param input input SDF file
#' @param output output CSV file to write to
#' @param append Append to existing CSV file or start anew?, Default: FALSE
#' @param fct Function to apply to each object to get the wanted columns
#' @param Nlines Lines to read in one go?, Default: 10000
#' @param startline Start at line, Default: 1
#' @param restartNlines Restart after x lines, Default: 10000
#' @param silent Suppress warnings?, Default: FALSE
#' @param ... Other arguments
#' @return Nothing, writes csv version of given SDF files to disk
#' @rdname sdfStream.joanna
#' @export
sdfStream.joanna <- function(input, output, append = FALSE, fct, Nlines = 10000, startline = 1, restartNlines = 10000, silent = FALSE, ...) {
  stop <- FALSE
  f <- file(input, "r")
  n <- Nlines
  offset <- 0
  if (startline != 1) {
    fmap <- file(input, "r")
    shiftback <- 2
    chunkmap <- scan(fmap, skip = startline - shiftback, nlines = restartNlines, what = "a", blank.lines.skip = FALSE, quiet = TRUE, sep = "\n")
    startline <- startline + (which(grepl("^\\${4,4}", chunkmap, perl = TRUE))[1] + 1 - shiftback)
    if (is.na(startline)) {
      stop("Invalid value assigned to startline.")
    }
    dummy <- scan(f, skip = startline - 2, nlines = 1, what = "a", blank.lines.skip = FALSE, quiet = TRUE, sep = "\n")
    close(fmap)
    offset <- startline - 1
  }
  counter <- 0
  cmpid <- 1
  partial <- NULL
  while (!stop) {
    counter <- counter + 1
    chunk <- readLines(f, n = n)
    if (length(chunk) > 0) {
      if (length(partial) > 0) {
        chunk <- c(partial, chunk)
      }
      inner <- sum(grepl("^\\${4,4}", chunk, perl = TRUE)) < 2
      while (inner) {
        chunklength <- length(chunk)
        chunk <- c(chunk, readLines(f, n = n))
        if (chunklength == length(chunk)) {
          inner <- FALSE
        }
        else {
          inner <- sum(grepl("^\\${4,4}", chunk, perl = TRUE)) < 2
        }
      }
      y <- regexpr("^\\${4,4}", chunk, perl = TRUE)
      index <- which(y != -1)
      indexDF <- data.frame(start = c(1, index[-length(index)] + 1), end = index)
      complete <- chunk[1:index[length(index)]]
      if ((index[length(index)] + 1) <= length(chunk)) {
        partial <- chunk[(index[length(index)] + 1):length(chunk)]
      }
      else {
        partial <- NULL
      }
      index <- index + offset
      indexDF <- data.frame(SDFlineStart = c(offset + 1, index[-length(index)] + 1), SDFlineEnd = index)
      offset <- indexDF[length(indexDF[, 2]), 2]
      sdfset <- ChemmineR::read.SDFset(ChemmineR::read.SDFstr(complete), skipErrors = TRUE)
      if (length(indexDF[, 1]) == 1) {
        suppressWarnings(sdfset <- c(sdfset, sdfset))
        resultMA <- fct(sdfset, ...)
      }
      else {
        resultMA <- fct(sdfset, ...)
      }
      if (silent == FALSE) {
        print(rownames(resultMA))
      }
      if (counter == 1 & append != TRUE) {
        unlink(output)
        write.table(resultMA, output, quote = FALSE, col.names = NA, sep = "\t")
      }
      else {
        write.table(resultMA, output, quote = FALSE, append = TRUE, col.names = FALSE, sep = "\t")
      }
    }
    if (length(chunk) == 0) {
      stop <- TRUE
      close(f)
    }
  }
}

#' @title Is this item 'empty'?
#' @description Checks if given object is either NULL, NA, or just whitespace
#' @param item object to check
#' @return TRUE or FALSE
#' @examples
#'  is.empty(NA)
#' @rdname is.empty
#' @export
is.empty <- function(item) {
  if (is.null(item)) {
    return(TRUE)
  }
  else if (is.na(item)) {
    return(TRUE)
  }
  else if (gsub(item, pattern = " ", replacement = "") == "") {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

