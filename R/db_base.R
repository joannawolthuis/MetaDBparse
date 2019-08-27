smiles.to.iatom <- function(smiles, silent=T, cl=0){

  require(rcdk)

  iatoms <- sapply(smiles, function(x, silent){
    mol=NULL
    try({
      mol = rcdk::parse.smiles(x,kekulise = F)[[1]]
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

iatom.to.smiles <- function(iatoms, smitype="Canonical", silent=T){

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

  require(rcdk)

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

iatom.to.charge <- function(iatoms, silent=T){

  require(rcdk)

  new.charges <- sapply(iatoms, function(mol, silent){
    ch=0
    try({
      ch = rcdk::get.total.formal.charge(molecule = mol)
    }, silent=silent)
    ch
  },silent=silent)

  try({
    rJava::.jcall("java/lang/System","V","gc")
    gc()
  })

  return(new.charges)
}

iatom.to.formula <- function(iatoms, silent=T){

  require(rcdk)

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

# BIG BOI

buildBaseDB <- function(outfolder, dbname,
                        smitype = "Canonical", silent=T, cl=0){

  require(enviPat)
  data(isotopes)

  removeDB(outfolder, paste0(dbname,".db"))
  conn <- openBaseDB(outfolder, paste0(dbname,".db"))
  db.formatted <- switch(dbname,
                    chebi = build.CHEBI(outfolder),
                    maconda = build.MACONDA(outfolder),
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
                    wikipathways = build.WIKIPATHWAYS(outfolder),
                    t3db = build.T3DB(outfolder),
                    vmh = build.VMH(outfolder),
                    hmdb = build.HMDB(outfolder),
                    smpdb = build.SMPDB(outfolder),
                    supernatural = build.SUPERNATURAL(outfolder))


  db.formatted <- data.table::as.data.table(db.formatted)
  print(nrow(db.formatted))

  blocks = split(1:nrow(db.formatted), ceiling(seq_along(1:nrow(db.formatted))/1000))

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

    require(enviPat)

    new.smiles = iatom.to.smiles(iats[valid.struct], smitype = "Canonical",silent=silent)

    new.charge = iatom.to.charge(iats[valid.struct],silent=silent)

    new.formula = iatom.to.formula(iats[valid.struct],silent=silent)

    db.redone.struct <- db.form.block
    db.redone.struct$structure[valid.struct] <- new.smiles
    db.redone.struct$baseformula[valid.struct] <- new.formula
    db.redone.struct$charge[valid.struct] <- new.charge

    db.removed.invalid <- db.redone.struct
    formulas = as.character(db.redone.struct$baseformula)
    null.or.na <- which(is.null(formulas) | is.na(formulas) | formulas == "NULL")

    if(length(null.or.na)>0){
      db.removed.invalid <- db.removed.invalid[-null.or.na,]
      valid.struct = valid.struct[-null.or.na]
    }

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

    deuterated = which(grepl("D\\d*", x = db.removed.invalid$baseformula))

    if(length(deuterated)>0){
      nondeuterated = gsub("D(\\d)*", "H\\1", db.removed.invalid$baseformula[deuterated])
      matching = data.table::as.data.table(db.removed.invalid)[baseformula %in% nondeuterated,]
      if(nrow(matching)>0){
        print("in progress... merge descriptions and add a note for deuterated")
      }else{
        db.removed.invalid$baseformula[deuterated] <- gsub("D(\\d)*", "H\\1",
                                                           db.removed.invalid$baseformula[deuterated])
        db.removed.invalid$description[deuterated] <- paste0("THIS DESCRIPTION IS FOR A SPECIFIC ISOTOPE, LIKELY NOT THE 100 PEAK!",
                                                             db.removed.invalid$description[deuterated])
      }
    }
    return(db.removed.invalid)
  })

  db.final <- data.table::rbindlist(db.fixed.rows)
  # - - - - - - - - - - - - - - - - - -
  writeDB(conn, db.final, "base")
  DBI::dbDisconnect(conn)
}

