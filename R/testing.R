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

# options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx16g")) # give java enough memory for smiles parsing
#
# outfolder = "/Users/jwolthuis/MetaDBparse/dbtest/"
#
outfolder <- "~/MetaboShiny/databases/"
buildBaseDB(outfolder, "chebi")
buildExtDB(outfolder,
           base.dbname = "chebi",
           cl = 0,
           blocksize = 600,
           mzrange = c(60,600),
           adduct_table = adducts,
           adduct_rules = adduct_rules)
#
# # === intermezzo ===
#
# adducts <- fread("~/Google Drive/MetaboShiny/backend/adducts/adduct_rule_table.csv", header = T) # V2 has di/trimers
# adduct_rules <- fread("~/Google Drive/MetaboShiny/backend/adducts/adduct_rule_smarts.csv", header = T) # V2 has di/trimers
#
# usethis::use_data(adducts, overwrite=T)
# usethis::use_data(adduct_rules, overwrite=T)
#
# # install package
# devtools::install()
