# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx16g")) # give java enough memory for smiles parsing

outfolder = "/Users/jwolthuis/MetaDBparse/dbtest/"

removeBaseDB(outfolder, "chebi.db")
conn <- openBaseDB(outfolder, "chebi.db")
db.orig <- build.CHEBI(outfolder)
db.orig <- build.MACONDA(outfolder)

iats = smiles.to.iatom(db.orig$structure)
valid.struct <- unlist(lapply(iats, function(x) !is.null(x)))
iats.valid <- iats[valid.struct]

# - - - - - - - -

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
checked <- enviPat::check_chemform(isotopes, as.character(db.removed.invalid$baseformula))
db.removed.invalid$baseformula <- checked$new_formula
# remove these
no.structure.no.formula <- which(checked$warning & !valid.struct)
db.removed.invalid <- db.removed.invalid[-no.structure.no.formula,]

writeDB(conn, db.removed.invalid, "base")
