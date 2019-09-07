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

sdfStream.joanna <- function (input, output, append = FALSE, fct, Nlines = 10000,
                              startline = 1, restartNlines = 10000, silent = FALSE, ...)
{
  require(ChemmineR)

  stop <- FALSE
  f <- file(input, "r")
  n <- Nlines
  offset <- 0
  if (startline != 1) {
    fmap <- file(input, "r")
    shiftback <- 2
    chunkmap <- scan(fmap, skip = startline - shiftback,
                     nlines = restartNlines, what = "a", blank.lines.skip = FALSE,
                     quiet = TRUE, sep = "\n")
    startline <- startline + (which(grepl("^\\${4,4}", chunkmap,
                                          perl = TRUE))[1] + 1 - shiftback)
    if (is.na(startline))
      stop("Invalid value assigned to startline.")
    dummy <- scan(f, skip = startline - 2, nlines = 1, what = "a",
                  blank.lines.skip = FALSE, quiet = TRUE, sep = "\n")
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
      inner <- sum(grepl("^\\${4,4}", chunk, perl = TRUE)) <
        2
      while (inner) {
        chunklength <- length(chunk)
        chunk <- c(chunk, readLines(f, n = n))
        if (chunklength == length(chunk)) {
          inner <- FALSE
        }
        else {
          inner <- sum(grepl("^\\${4,4}", chunk, perl = TRUE)) <
            2
        }
      }
      y <- regexpr("^\\${4,4}", chunk, perl = TRUE)
      index <- which(y != -1)
      indexDF <- data.frame(start = c(1, index[-length(index)] +
                                        1), end = index)
      complete <- chunk[1:index[length(index)]]
      if ((index[length(index)] + 1) <= length(chunk)) {
        partial <- chunk[(index[length(index)] + 1):length(chunk)]
      }
      else {
        partial <- NULL
      }
      index <- index + offset
      indexDF <- data.frame(SDFlineStart = c(offset +
                                               1, index[-length(index)] + 1), SDFlineEnd = index)
      offset <- indexDF[length(indexDF[, 2]), 2]
      sdfset <- read.SDFset(read.SDFstr(complete))
      #valid <- validSDF(sdfset)
      #sdfset <- sdfset[valid]
      #indexDForig <- indexDF
      #indexDF <- indexDF[valid, ]
      if (length(indexDF[, 1]) == 1) {
        suppressWarnings(sdfset <- c(sdfset, sdfset))
        resultMA <- fct(sdfset, ...)
        # resultMA <- cbind(as.data.frame(indexDF), as.data.frame(resultMA[1,
        #                                                                  , drop = FALSE]), row.names = row.names(resultMA)[1])
      }
      else {
        resultMA <- fct(sdfset, ...)
        # resultMA <- cbind(as.data.frame(indexDF), as.data.frame(resultMA),
        #                   row.names = row.names(resultMA))
      }
      #resultMA <- resultMA[names(valid), ]
      #if (any(is.na(resultMA))) {
      #  resultMA[, 1:2] <- indexDForig[, 1:2]
      #}
      #rownames(resultMA) <- paste("CMP", cmpid:(cmpid +
      #                                            length(resultMA[, 1]) - 1), sep = "")
      #cmpid <- cmpid + length(resultMA[, 1])
      if (silent == FALSE) {
        print(rownames(resultMA))
      }
      if (counter == 1 & append != TRUE) {
        unlink(output)
        write.table(resultMA, output, quote = FALSE,
                    col.names = NA, sep = "\t")
      }
      else {
        write.table(resultMA, output, quote = FALSE,
                    append = TRUE, col.names = FALSE, sep = "\t")
      }
    }
    if (length(chunk) == 0) {
      stop <- TRUE
      close(f)
    }
  }
}

#Overwrite RSQLite's sqliteWriteTable function with this one to be able to reject duplicates
#example:
#create table foo(x int, y int, z int, constraint id primary key(x,y));
#insert or ignore into foo(x,y,z) values(1,2,3);
#insert or ignore into foo(x,y,z) values(1,2,3);
#
#In R
#mysqliteWriteTable(con, "foo", data.frame(x=1,y=2,z=3), row.names=F, ignore=T)
mysqliteWriteTable = function (con, name, value, row.names = TRUE, overwrite = FALSE,
                               append = FALSE, field.types = NULL, ignore = FALSE, ...)
{
  if (overwrite && append)
    stop("overwrite and append cannot both be TRUE")
  if (length(RSQLite::dbListResults(con))) {
    new.con <- RSQLite::dbConnect(con)
    on.exit(RSQLite::dbDisconnect(new.con))
  }
  else {
    new.con <- con
  }
  foundTable <- RSQLite::dbExistsTable(con, name)
  new.table <- !foundTable
  createTable <- (new.table || foundTable && overwrite)
  removeTable <- (foundTable && overwrite)
  success <- RSQLite::dbBeginTransaction(con)
  if (!success) {
    warning("unable to begin transaction")
    return(FALSE)
  }
  if (foundTable && !removeTable && !append) {
    warning(paste("table", name, "exists in database: aborting dbWriteTable"))
    success <- FALSE
  }
  if (removeTable) {
    success <- tryCatch({
      if (RSQLite::dbRemoveTable(con, name)) {
        TRUE
      }
      else {
        warning(paste("table", name, "couldn't be overwritten"))
        FALSE
      }
    }, error = function(e) {
      warning(conditionMessage(e))
      FALSE
    })
  }
  if (!success) {
    RSQLite::dbRollback(con)
    return(FALSE)
  }
  if (row.names) {
    value <- cbind(row.names(value), value, stringsAsFactors = FALSE)
    names(value)[1] <- "row.names"
  }
  if (createTable) {
    sql <- RSQLite::dbBuildTableDefinition(new.con, name, value, field.types = field.types,
                                  row.names = FALSE)
    success <- tryCatch({
      RSQLite::dbGetQuery(new.con, sql)
      TRUE
    }, error = function(e) {
      warning(conditionMessage(e))
      FALSE
    })
    if (!success) {
      RSQLite::dbRollback(con)
      return(FALSE)
    }
  }
  valStr <- paste(rep("?", ncol(value)), collapse = ",")
  if (ignore) sql <- sprintf("insert or ignore into %s values (%s)", name, valStr)
  else sql <- sprintf("insert into %s values (%s)", name, valStr)
  success <- tryCatch({
    ret <- FALSE
    rs <- RSQLite::dbSendPreparedQuery(new.con, sql, bind.data = value)
    ret <- TRUE
  }, error = function(e) {
    warning(conditionMessage(e))
    ret <- FALSE
  }, finally = {
    if (exists("rs"))
      RSQLite::dbClearResult(rs)
    ret
  })
  if (!success)
    RSQLite::dbRollback(con)
  else {
    success <- RSQLite::dbCommit(con)
    if (!success) {
      warning(RSQLite::dbGetException(con)[["errorMsg"]])
      RSQLite::dbRollback(con)
    }
  }
  success
}

is.empty <- function(item){
  if(is.null(item)){
      return(TRUE)
    }else if(is.na(item)){
      return(TRUE)
    }else if(gsub(item, pattern = " ", replacement="") == ""){
        return(TRUE)
    }else{
      return(FALSE)
    }
  }

