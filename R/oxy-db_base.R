#' @title Get iatom containers from SMILES
#' @description FUNCTION_DESCRIPTION
#' @param smiles character vector of smiles
#' @param silent suppress warnings?, Default: TRUE
#' @param cl parallel::makeCluster object for multithreading, Default: 0
#' @return iatom containers for use in rcdk package
#' @examples
#' if(interactive()){
#'  smiles.to.iatom(c('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'))
#'  }
#' @seealso
#'  \code{\link[rJava]{jcall}}
#' @rdname smiles.to.iatom
#' @export
#' @importFrom rcdk parse.smiles do.aromaticity do.typing do.isotopes
#' @importFrom rJava .jcall
smiles.to.iatom <- function(smiles, silent=TRUE, cl=0){
  iatoms <- sapply(smiles, function(x, silent){
    mol=NULL
    try({
      mol = NULL
      try({
        mol = rcdk::parse.smiles(x)[[1]]
      })
      if(is.null(mol)){
        mol = rcdk::parse.smiles(x,kekulise = FALSE)[[1]]
      }
      rcdk::do.aromaticity(mol)
      rcdk::do.typing(mol)
      rcdk::do.isotopes(mol)
    }, silent=silent)
    mol
  },silent=silent)

  try({
    rJava::.jcall("java/lang/System","V","gc")
    gc()
  })

  return(iatoms)
}

#' @title Get SMILES from iatom container
#' @description This function takes an rcdk iatomcontainer and returns SMILES
#' @param iatoms list of iatomcontainers
#' @param smitype Which type of SMILES to export?, Default: 'Canonical'
#' @param silent Suppress warnings?, Default: TRUE
#' @return character vector of SMILES in the chosen format
#' @seealso
#'  \code{\link[rcdk]{get.smiles}},\code{\link[rcdk]{smiles.flavors}}
#'  \code{\link[rJava]{jcall}}
#' @rdname iatom.to.smiles
#' @export
#' @importFrom rcdk get.smiles smiles.flavors
#' @importFrom rJava .jcall
iatom.to.smiles <- function(iatoms, smitype="Canonical", silent=TRUE){

  if(!silent){
    cat("Valid SMILES output options include:\n\n")
    cat(c(" - - - \n" , "Absolute", "AtomAtomMap", "AtomicMass",
          "AtomicMassStrict", "Canonical", "Cx2dCoordinates", "Cx3dCoordinates",
          "CxAtomLabel", "CxAtomValue", "CxCoordinates", "CxFragmentGroup",
          "CxMulticenter", "CxPolymer", "CxRadical", "CxSmiles",
          "CxSmilesWithCoords", "Default", "Generic", "InChILabelling",
          "Isomeric", "Stereo", "StereoCisTrans", "StereoExTetrahedral",
          "StereoTetrahedral", "Unique", "UniversalSmiles", "UseAromaticSymbols\n", "- - - \n\n"))
    cat("Defaulting to 'Canonical'.")
  }

  new.smiles <- sapply(iatoms, function(mol, silent){
    smi = ""
    try({
      smi <- if(is.null(mol)) smi = "" else rcdk::get.smiles(mol, flavor = rcdk::smiles.flavors(smitype))
    }, silent=silent)
    smi
  },silent=silent)

  try({
    rJava::.jcall("java/lang/System","V","gc")
    gc()
  })

  return(new.smiles)
}

#' @title Get formal charge from iatomcontainer
#' @description This function takes iatomcontainer object and returns the formal charge.
#' @param iatoms list of rcdk iatomcontainers
#' @param silent suppress warnings?, Default: TRUE
#' @return Character vector of formal charges per iatomcontainer.
#' @seealso
#'  \code{\link[rcdk]{get.total.formal.charge}}
#'  \code{\link[rJava]{jcall}}
#' @rdname iatom.to.charge
#' @export
#' @importFrom rcdk get.total.formal.charge
#' @importFrom rJava .jcall
iatom.to.charge <- function(iatoms, silent=TRUE){

  new.charges <- sapply(iatoms, function(mol, silent){
    ch=0
    try({
      ch = rcdk::get.total.formal.charge(mol = mol)
    }, silent=silent)
    ch
  },silent=silent)

  try({
    rJava::.jcall("java/lang/System","V","gc")
    gc()
  })

  return(new.charges)
}

#' @title Get molecular formula from iatomcontainer
#' @description This function takes iatomcontainer object and returns the molecular formula.
#' @param iatoms list of rcdk iatomcontainers
#' @param silent suppress warnings?, Default: TRUE
#' @return Character vector of formulas per iatomcontainer.
#' @seealso
#'  \code{\link[rcdk]{get.mol2formula}}
#'  \code{\link[rJava]{jcall}}
#' @rdname iatom.to.formula
#' @export
#' @importFrom rcdk get.mol2formula
#' @importFrom rJava .jcall
iatom.to.formula <- function(iatoms, silent=TRUE){

  new.formulas <- sapply(iatoms, function(mol,silent){
    form = NULL
    try({
      form = rcdk::get.mol2formula(mol)@string
    }, silent=silent)
    form[[1]]
  },silent=silent)

  try({
    rJava::.jcall("java/lang/System","V","gc")
    gc()
  })

  return(new.formulas)
}

#' @title Uniformize database and remove invalid formulas/SMILES
#' @description This is a wrapper function to take a 'raw' input data table with compound information, uniformize the SMILES
#' @param db.formatted Data table with columns 'compoundname, structure, baseformula, charge, description'
#' @param cl parallel::makeCluster object for multithreading
#' @param silent Suppress warnings?
#' @param blocksize How many compounds to process per 'block'? Higher number means bigger memory spikes, but faster processing time.
#' @param smitype SMILES format, Default: 'Canonical'
#' @return Data table with SMILES in the correct format, and charge/formula re-generated from said SMILES if available.
#' @seealso
#'  \code{\link[parallel]{clusterApply}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[enviPat]{check_chemform}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname cleanDB
#' @export
#' @importFrom parallel clusterExport
#' @importFrom pbapply pblapply
#' @importFrom enviPat check_chemform
#' @importFrom data.table rbindlist
cleanDB <- function(db.formatted, cl, silent, blocksize, smitype='Canonical'){
  blocks = split(1:nrow(db.formatted), ceiling(seq_along(1:nrow(db.formatted))/blocksize))

  if(is.list(cl)){
    parallel::clusterExport(cl, varlist = c("db.formatted",
                                            "iatom.to.smiles",
                                            "smiles.to.iatom",
                                            "iatom.to.formula",
                                            "iatom.to.charge",
                                            "silent",
                                            "isotopes"),
                            envir = environment())
  }

  print("Uniformizing formulas/SMILES/charges and checking for mistakes...")

  db.fixed.rows <- pbapply::pblapply(blocks, cl=cl, function(block){

    db.form.block = db.formatted[block,]
    iats = smiles.to.iatom(db.form.block$structure,
                           silent=silent)

    valid.struct <- unlist(lapply(iats, function(x) !is.null(x)))

    new.smiles = iatom.to.smiles(iats[valid.struct], smitype = smitype, silent=silent)

    new.charge = iatom.to.charge(iats[valid.struct],silent=silent)

    new.formula = iatom.to.formula(iats[valid.struct],silent=silent)

    db.redone.struct <- db.form.block
    db.redone.struct$structure[valid.struct] <- new.smiles
    db.redone.struct$baseformula[valid.struct] <- new.formula
    db.redone.struct$charge[valid.struct] <- new.charge

    db.removed.invalid <- db.redone.struct
    formulas = as.character(db.redone.struct$baseformula)
    null.or.na <- which(is.null(formulas) | is.na(formulas) | formulas == "NULL")

    if(length(null.or.na) > 0){
      db.removed.invalid <- db.removed.invalid[-null.or.na,]
      valid.struct = valid.struct[-null.or.na]
    }

    print(head(db.removed.invalid))

    checked <- enviPat::check_chemform(isotopes,
                                       chemforms = as.character(db.removed.invalid$baseformula))

    db.removed.invalid$baseformula <- checked$new_formula

    invalid.struct <- !valid.struct

    # dummy structures
    db.removed.invalid$structure[which(invalid.struct)] <- paste0("[",
                                                                  db.removed.invalid$baseformula[invalid.struct],
                                                                  "]",
                                                                  db.removed.invalid$charge[invalid.struct])
    # remove these
    invalid.formula <- which(checked$warning)
    #no.structure.no.formula <- which(checked$warning & !valid.struct)
    if(length(invalid.formula) > 0){
      db.removed.invalid <- db.removed.invalid[-invalid.formula,]
    }

    return(db.removed.invalid)
  })

  return(data.table::rbindlist(db.fixed.rows))
}

#' @title Build the base database
#' @description This is a large wrapper function that calls upon all individual database parsers, cleans the resulting database and saves it in a SQLite database.
#' @param outfolder In which folder are you building your databases? Temp folders etc. will be put here.
#' @param dbname Which database do you want to build? Options: chebi,maconda,kegg,bloodexposome,dimedb,expoexplorer, foodb, drugbank, lipidmaps, massbank, metabolights, metacyc, phenolexplorer, respect, wikidata, wikipathways, t3db, vmh, hmdb, smpdb, lmdb, ymdb, ecmdb, bmdb, rmdb, stoff, nanpdb, mcdb, mvoc, pamdb
#' @param custom_csv_path PARAM_DESCRIPTION, Default: NULL
#' @param smitype Which SMILES format do you want?, Default: 'Canonical'
#' @param silent Suppress warnings?, Default: TRUE
#' @param cl parallel::makeCluster object for multithreading, Default: 0
#' @return Nothing, writes SQLite database to 'outfolder'.
#' @seealso
#'  \code{\link[data.table]{fread}},\code{\link[data.table]{as.data.table}}
#'  \code{\link[DBI]{dbDisconnect}}
#' @rdname buildBaseDB
#' @export
#' @importFrom data.table fread as.data.table data.table
#' @importFrom RSQLite dbExecute
#' @importFrom DBI dbDisconnect
buildBaseDB <- function(outfolder, dbname, custom_csv_path=NULL,
                        smitype = "Canonical", silent=TRUE, cl=0){

  removeDB(outfolder, paste0(dbname,".db"))
  conn <- openBaseDB(outfolder, paste0(dbname,".db"))
  if(is.null(custom_csv_path)){
    db.formatted.all <- switch(dbname,
                           chebi = build.CHEBI(outfolder),
                           maconda = build.MACONDA(outfolder, conn),
                           kegg = build.KEGG(outfolder),
                           bloodexposome = build.BLOODEXPOSOME(outfolder),
                           dimedb = build.DIMEDB(outfolder),
                           expoexplorer = build.EXPOSOMEEXPLORER(outfolder),
                           foodb = build.FOODB(outfolder),
                           drugbank = build.DRUGBANK(outfolder),
                           lipidmaps = build.LIPIDMAPS(outfolder),
                           massbank = build.MASSBANK(outfolder),
                           metabolights = build.METABOLIGHTS(outfolder),
                           metacyc = build.METACYC(outfolder),
                           phenolexplorer = build.PHENOLEXPLORER(outfolder),
                           respect = build.RESPECT(outfolder),
                           wikidata = build.WIKIDATA(outfolder),
                           #wikipathways = build.WIKIPATHWAYS(outfolder),
                           t3db = build.T3DB(outfolder),
                           vmh = build.VMH(outfolder),
                           hmdb = build.HMDB(outfolder),
                           smpdb = build.SMPDB(outfolder),
                           lmdb = build.LMDB(outfolder),
                           ymdb = build.YMDB(outfolder),
                           ecmdb = build.ECMDB(outfolder),
                           bmdb = build.BMDB(outfolder),
                           rmdb = build.RMDB(outfolder),
                           stoff = build.STOFF(outfolder),
                           nanpdb = build.NANPDB(outfolder),
                           mcdb = build.MCDB(outfolder),
                           mvoc = build.mVOC(outfolder),
                           pamdb = build.PAMDB(outfolder))
  }else{
    db.formatted.all <- list(db = data.table::fread(custom_csv_path, header=TRUE),
                             version = Sys.time())
  }

  if(dbname == "maconda") return(NA)

  print(db.formatted.all$version)
  print(head(db.formatted.all$db))

  db.formatted <- data.table::as.data.table(db.formatted.all$db)
  db.formatted <- data.frame(lapply(db.formatted, as.character), stringsAsFactors=FALSE)
  db.final <- data.table::as.data.table(cleanDB(db.formatted,
                                   cl = cl,
                                   silent = silent,
                                   blocksize=400))

  db.final <- db.final[, lapply(.SD, as.character)]

  # - - - - - - - - - - - - - - - - - -
  writeDB(conn, data.table::data.table(date = Sys.Date(),
                                       version = db.formatted.all$version),
          "metadata")
  writeDB(conn, table = db.final, "base")
  RSQLite::dbExecute(conn, "CREATE INDEX b_idx1 ON base(structure)")
  DBI::dbDisconnect(conn)
}

