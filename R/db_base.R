smiles.to.iatom <- function(smiles){

  require(rcdk)

  iatoms <- sapply(smiles, function(x){
    mol=NULL
    try({
      mol = rcdk::parse.smiles(x,kekulise = F)[[1]]
      rcdk::do.aromaticity(mol)
      rcdk::do.typing(mol)
      rcdk::do.isotopes(mol)
    }, silent=T)
    mol
  })

  try({
    rJava::.jcall("java/lang/System","V","gc")
    gc()
  }, silent=T)

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

  new.smiles <- sapply(iatoms, function(mol){
    smi = ""
    try({
      smi <- if(is.null(mol)) smi = "" else rcdk::get.smiles(mol, flavor = rcdk::smiles.flavors(smitype))
    }, silent=T)
    smi
  })

  rJava::.jcall("java/lang/System","V","gc")
  gc()

  return(new.smiles)
}

iatom.to.charge <- function(iatoms){

  require(rcdk)

  new.charges <- sapply(iatoms, function(mol){
    ch = rcdk::get.total.formal.charge(molecule = mol)
  })

  rJava::.jcall("java/lang/System","V","gc")
  gc()

  return(new.charges)
}

iatom.to.formula <- function(iatoms){

  require(rcdk)

  new.formulas <- sapply(iatoms, function(mol){
    form = NULL
    try({
      form = rcdk::get.mol2formula(mol)@string
    }, silent=T)
    form
  })

  return(new.formulas)
}

# BIG BOI

buildBaseDB <- function(outfolder, dbname, smitype = "Canonical"){

  removeDB(outfolder, paste0(dbname,".db"))
  conn <- openBaseDB(outfolder, paste0(dbname,".db"))
  db.orig <- switch(dbname,
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

  iats = smiles.to.iatom(db.orig$structure)
  valid.struct <- unlist(lapply(iats, function(x) !is.null(x)))
  iats.valid <- iats[valid.struct]

  require(enviPat)

  new.smiles = iatom.to.smiles(iats.valid, smitype = "Canonical")
  new.charge = iatom.to.charge(iats.valid)
  new.formula = iatom.to.formula(iats.valid)

  db.redone.struct <- db.orig
  db.redone.struct$structure[valid.struct] <- new.smiles
  db.redone.struct$baseformula[valid.struct] <- new.formula
  db.redone.struct$charge[valid.struct] <- new.charge

  # - - - - - - - - -

  db.removed.invalid <- db.redone.struct
  checked <- enviPat::check_chemform(isotopes, chemforms = as.character(db.removed.invalid$baseformula))
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
  db.removed.invalid <- db.removed.invalid[-invalid.formula,]

  deuterated = which(grepl("D\\d*", x = db.removed.invalid$baseformula))
  nondeuterated = gsub("D(\\d)*", "H\\1", db.removed.invalid$baseformula[deuterated])
  matching = db.removed.invalid[baseformula %in% nondeuterated,]
  if(nrow(matching)>0){
    print("in progress... merge descriptions and add a note for deuterated")
  }else{
    db.removed.invalid$baseformula[deuterated] <- gsub("D(\\d)*", "H\\1", db.removed.invalid$baseformula[deuterated])
    db.removed.invalid$description[deuterated] <- paste0("THIS DESCRIPTION IS FOR A SPECIFIC ISOTOPE, LIKELY NOT THE 100 PEAK! ", db.removed.invalid$description[deuterated])
  }

  # - - - - - - - - - - - - - - - - - -
  writeDB(conn, db.removed.invalid, "base")
  DBI::dbDisconnect(conn)
}

