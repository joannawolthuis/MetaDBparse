# Hello, world!
#
#usethis::use_build_ignore(c("R/testing.R"))
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
#outfolder <- "~/MetaboShiny/databases/"
db="chebi"
#
cl = parallel::makeCluster(3)

# options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx16G")) # give java enough memory for smiles parsing
#
# dbs = c("lipidmaps", "kegg", "phenolexplorer", "chebi", "vmh", "dimedb",
#         "t3db", "expoexplorer", "wikidata", "massbank", "respect",
#         "hmdb", "metacyc", "drugbank", "foodb", "smpdb")
# for(db in dbs){
#   print(db)
#   try({
     buildBaseDB(outfolder = outfolder,
                 dbname = db,
                 cl=cl, silent=F)
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
#   try({
#     parallel::stopCluster(session_cl)
#   },silent=T)
#
#   session_cl <- parallel::makeCluster(3, outfile="")
#
#   parallel::clusterExport(session_cl, c("smiles.to.iatom",
#                                         "countAdductRuleMatches",
#                                         "checkAdductRule",
#                                         "doAdduct",
#                                         "iatom.to.smiles",
#                                         "smiles.to.iatom",
#                                         "iatom.to.formula",
#                                         "iatom.to.charge",
#                                         "adduct_rules",
#                                         "adducts",
#                                         "doIsotopes",
#                                         "isotopes"))
#
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
#              adduct_rules = adduct_rules,
#              silent = F,
#              ext.dbname = "extended")
# }
#
# dbs = gsub(basename(list.files(outfolder, pattern="\\.db$")), pattern = "\\.db", replacement="")
# dbs = dbs[dbs!="extended"]
#
# res = searchMZ(mzs = "110.071176455696",
#                ionmodes = "positive",
#                outfolder = outfolder,
#                base.dbname = dbs,
#                ppm = 5)
#
# "SELECT DISTINCT
# cpd.adduct as adduct,
# cpd.isoprevalence as isoprevalence,
# struc.smiles as structure,
# mz.mzmed as query_mz,
# (1e6*ABS(mz.mzmed - cpd.fullmz)/cpd.fullmz) AS dppm
# FROM mzvals mz
# JOIN mzranges rng ON rng.ID = mz.ID
# JOIN extended cpd INDEXED BY e_idx2
# ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
# JOIN adducts
# ON cpd.adduct = adducts.Name
# AND mz.foundinmode = adducts.Ion_Mode
# JOIN structures struc
# ON cpd.struct_id = struc.struct_id"
#
# dbs = list.files(outfolder, pattern = "\\.db$")
# dbs = gsub(dbs[dbs != "extended.db"], pattern="\\.db", replacement="")
#
# for(db in dbs){
#   print(db)
#   conn = openBaseDB(outfolder, db)
#   RSQLite::dbExecute(conn,"CREATE INDEX IF NOT EXISTS b_idx1 ON base(structure);")
#   RSQLite::dbDisconnect(conn)
# }
