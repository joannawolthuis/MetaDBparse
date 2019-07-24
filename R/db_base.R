smiles.to.iatom <- function(smiles){

  require(rcdk)

  iatoms <- pbapply::pbsapply(smiles, function(x){
    mol=NULL
    try({
      mol = rcdk::parse.smiles(x,kekulise = F)[[1]]
      rcdk::do.aromaticity(mol)
      rcdk::do.typing(mol)
      rcdk::do.isotopes(mol)
    })
    mol
  })

  rJava::.jcall("java/lang/System","V","gc")
  gc()

  return(iatoms)
}

iatom.to.smiles <- function(iatoms, smitype="Canonical"){

  cat("Valid SMILES output options include:\n\n")
  cat(c(" - - - \n" , "Absolute", "AtomAtomMap", "AtomicMass",
        "AtomicMassStrict", "Canonical", "Cx2dCoordinates", "Cx3dCoordinates",
        "CxAtomLabel", "CxAtomValue", "CxCoordinates", "CxFragmentGroup",
        "CxMulticenter", "CxPolymer", "CxRadical", "CxSmiles",
        "CxSmilesWithCoords", "Default", "Generic", "InChILabelling",
        "Isomeric", "Stereo", "StereoCisTrans", "StereoExTetrahedral",
        "StereoTetrahedral", "Unique", "UniversalSmiles", "UseAromaticSymbols\n", "- - - \n\n"))
  cat("Defaulting to 'Canonical'.")

  require(rcdk)

  new.smiles <- pbapply::pbsapply(iatoms, function(mol){
    smi = ""
    try({
      smi <- if(is.null(mol)) smi = "" else rcdk::get.smiles(mol, flavor = rcdk::smiles.flavors(smitype))
    })
    smi
  })

  rJava::.jcall("java/lang/System","V","gc")
  gc()

  return(new.smiles)
}

iatom.to.charge <- function(iatoms){

  require(rcdk)

  new.charges <- pbapply::pbsapply(iatoms, function(mol){
    ch = rcdk::get.total.formal.charge(molecule = mol)
  })

  rJava::.jcall("java/lang/System","V","gc")
  gc()

  return(new.charges)
}

iatom.to.formula <- function(iatoms){

  require(rcdk)

  new.formulas <- pbapply::pbsapply(iatoms, function(mol){
    form = NULL
    try({
      form = rcdk::get.mol2formula(mol)@string
    })
    form
  })

  return(new.formulas)
}
