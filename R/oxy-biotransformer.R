#' @title Download biotransformer jar file.
#' @description To predict enzymatic metabolite changes, you can enable Biotransformer processing in the base database creation. This downloads the necessary JAR file from the developer page.
#' @param outfolder In which folder do you want to save the .jar file?
#' @return File location
#' @seealso
#'  \code{\link[utils]{download.file}},\code{\link[utils]{unzip}}
#' @rdname downloadBT
#' @export
#' @examples
#'  \dontrun{jarloc = downloadBT(outfolder = tempdir())}
#' @importFrom utils download.file unzip
downloadBT <- function(outfolder) {
  file.url <- "https://bitbucket.org/djoumbou/biotransformerjar/get/f47aa4e3c0da.zip"
  zip.file <- file.path(outfolder, "biotransformer.zip")
  utils::download.file(file.url, zip.file, mode = "wb", cacheOK = TRUE, method = "auto")
  utils::unzip(zip.file, exdir = outfolder)
  newfolder <- file.path(outfolder, "djoumbou-biotransformerjar-f47aa4e3c0da")
  jarloc <- list.files(newfolder, pattern = "\\.jar", full.names = TRUE)
  system(gsubfn::fn$paste("chmod -x $jarloc"))
  return(jarloc)
}

#' @title Run Biotransformer on SMILES.
#' @description To predict enzymatic metabolite changes, you can enable Biotransformer processing in the base database creation. This runs the process on a given SMILES string.
#' @param smis SMILES string character vector
#' @param jarloc Location of Biotransformer .jar file.
#' @param opts -q command line options.
#' @param help Show help for opts only?
#' @param cl parallel::makeCluster object.
#' @return Data table with metabolites of given SMILES and which reaction it pertains.
#' @seealso
#'  \code{\link[data.table]{fread}},\code{\link[data.table]{rbindlist}}
#'  \code{\link[pbapply]{pblapply} }
#' @rdname doBT
#' @export
#' @examples
#' \dontrun{myFolder = tempdir()}
#' \dontrun{btLoc = downloadBT(outfolder = myFolder)}
#' \dontrun{doBT(smis = c("CCNC1=NC(=NC(=N1)Cl)NC(C)C"),
#'      jarloc = btLoc,
#'      opts = "cyp450:2; phaseII:1")}
#' @importFrom data.table fread rbindlist
#' @importFrom pbapply pblapply
doBT <- function(smis = c("CCNC1=NC(=NC(=N1)Cl)NC(C)C"),
                 jarloc,
                 opts = "cyp450:2; phaseII:1",
                 cl = 0,
                 help = FALSE){
  if (help == TRUE) {
    print("For opts, please specify the type of description: Type of Biotransformer - EC-based (ecbased), CYP450 (cyp450), Phase II (phaseII), Human gut microbial (hgut), human super transformer* (superbio, or allHuman), Environmental microbial (envimicro).")
  } else {
    oldDir <- getwd()
    if(getwd() != dirname(jarloc)){
      setwd(dirname(jarloc))
      on.exit({
        setwd(oldDir)
      })
    }
    btRows <- pbapply::pblapply(smis, cl = cl, function(smi) {
      mets <- data.table::data.table()
      try({
        smifile <- tempfile(fileext = ".csv")
        cmd <- gsubfn::fn$paste("java -jar $jarloc -ismi \"$smi\" -ocsv $smifile -k pred -q \"$opts\"")
        system(cmd, intern = TRUE)
        mets <- data.table::fread(smifile)
        mets$structure <- c(smi)
      })
      mets
    })
    data.table::rbindlist(btRows, fill = TRUE, use.names = TRUE)
  }
}
