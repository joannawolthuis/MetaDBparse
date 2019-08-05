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
# outfolder <- "~/MetaboShiny/databases/"
# dbname="kegg"
#

# options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx16G")) # give java enough memory for smiles parsing
#
# dbs = c("lipidmaps", "kegg", "phenolexplorer", "chebi", "vmh", "dimedb",
#         "t3db", "expoexplorer", "wikidata", "massbank", "respect",
#         "hmdb", "metacyc", "drugbank", "foodb", "smpdb")
# for(db in dbs){
#   print(db)
#   try({
#     buildBaseDB(outfolder, db, silent=F)
#   })
# }
#
# #
# require(parallel)
# require(data.table)
# require(enviPat)
# data(isotopes)
#
# adduct_rules <- fread("~/Google Drive/MetaboShiny/backend/adducts/adduct_rule_smarts.csv")
# adducts <- fread("~/Google Drive/MetaboShiny/backend/adducts/adduct_rule_table.csv")
#
# usethis::use_data(adduct_rules, overwrite=T)
# usethis::use_data(adducts, overwrite=T)
#
# {
#   #file.remove(file.path(outfolder, "extended.db"))
#   try({
#     parallel::stopCluster(session_cl)
#   },silent=T)
#   session_cl <- parallel::makeCluster(3, outfile="")
#   parallel::clusterExport(session_cl, c("smiles.to.iatom",
#                                         "countAdductRuleMatches",
#                                         "checkAdductRule",
#                                         "doAdduct",
#                                         "iatom.to.smiles",
#                                         "adduct_rules",
#                                         "adducts",
#                                         "doIsotopes",
#                                         "isotopes"))
#   parallel::clusterEvalQ(cl = session_cl, expr = {
#     library(data.table)
#     library(enviPat)
#     library(pbapply)
#   })
#
#   buildExtDB(outfolder,
#              base.dbname = dbname,
#              cl = session_cl,
#              blocksize = 200,
#              mzrange = c(60,600),
#              adduct_table = adducts,
#              adduct_rules = adduct_rules)
#   }
