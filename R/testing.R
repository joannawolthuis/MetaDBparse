# # # Hello, world!
# # #
# usethis::use_build_ignore(c("R/testing.R"))
# # # This is an example function named 'hello'
# # # which prints 'Hello, world!'.
# # #
# # # You can learn more about package authoring with RStudio at:
# # #
# # #   http://r-pkgs.had.co.nz/
# # #
# # # Some useful keyboard shortcuts for package authoring:
# # #
# # #   Install Package:           'Cmd + Shift + B'
# # #   Check Package:             'Cmd + Shift + E'
# # #   Test Package:             'Cmd + Shift + T'
# #
# # # options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx16g")) # give java enough memory for smiles parsing
# # #
# # # outfolder = "/Users/jwolthuis/MetaDBparse/dbtest/"
# # #
# # outfolder <- "~/MetaboShiny/databases/"
# #
# # #db="hmdb"
# #
# if(interactive()){
#
#   ncores = parallel::detectCores(all.tests = FALSE)
#
#   cl = parallel::makeCluster(ncores)
#
#   options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx16G")) # give java enough memory for smiles parsing
#
#   require(enviPat)
#   data(isotopes)
#   require(MetaDBparse)
#
#   adduct_rules <- data.table::fread("/Users/jwolthuis/MetaDBparse/inst/files/adduct_rules_deut.csv")
#   adducts <- data.table::fread("/Users/jwolthuis/MetaDBparse/inst/files/adducts_deut.csv")
#   adducts <- data.table::fread("/Users/jwolthuis/MetaboShiny/saves/admin/adducts.csv")
#   adducts <- adducts[Info != "deut"]
#
#   parallel::clusterExport(cl, c("smiles.to.iatom",
#                                 "countAdductRuleMatches",
#                                 "checkAdductRule",
#                                 "doAdduct",
#                                 "iatom.to.smiles",
#                                 "smiles.to.iatom",
#                                 "iatom.to.formula",
#                                 "iatom.to.charge",
#                                 "adduct_rules",
#                                 "adducts",
#                                 "doIsotopes",
#                                 "isotopes"))
#
#   parallel::clusterEvalQ(cl = cl, expr = {
#     library(data.table)
#     library(enviPat)
#     library(pbapply)
#   })
#
#   dbs = c("chebi",
#           "maconda",
#           "kegg",
#           "bloodexposome","dimedb",
#           "expoexplorer", "foodb",
#           "drugbank", "lipidmaps","massbank",
#           "metabolights",#<-ERRORS!!
#           "metacyc", "phenolexplorer",
#           "respect",
#           "wikidata",
#           "t3db", "vmh", "hmdb",
#           "smpdb", "lmdb", "ymdb",
#           "ecmdb", "bmdb", "rmdb",
#           "stoff", "nanpdb","mcdb",
#           "mvoc", "pamdb"
#   )
#
#   outfolder = normalizePath("~/MetaboShiny/databases")
#
#   registerDoParallel(cl, cores = detectCores() - 1)
#
#   library("parallel")
#   library("foreach")
#   library("doParallel")
#
#   successList=c()
#
#   for(db in dbs){
#     print(db)
#     success=F
#     try({
#       #=== INDEXING ===
#       conn <- openBaseDB(outfolder, paste0(db,".db"))
#       RSQLite::dbExecute(conn, "CREATE INDEX IF NOT EXISTS b_idx1 ON base(structure)")
#       RSQLite::dbDisconnect(conn)
#       #=== BUILD BASE ===
#       buildBaseDB(outfolder = outfolder,
#                   dbname = db,
#                   cl = cl, silent=F)
#       success=T
#       #=== BUILD EXTENDED ===
#       #buildExtDB(outfolder,
#       #           base.dbname = db,
#       #           cl = cl,
#       #           blocksize = 600,
#       #           mzrange = c(60,800),
#       #           adduct_table = adducts,
#       #           adduct_rules = adduct_rules,
#       #           silent = F,
#       #           ext.dbname = "extended_no_rules_2",
#       #           use.rules = F)
#     })
#     successList = c(successList, success)
#   }
#
#   # # # devtools::install()
#   # # #
#   # # # #
#   # # require(parallel)
#   # # require(data.table)
#   # # require(enviPat)
#   # # data(isotopes)
#   # #
#   # # #
#   # # usethis::use_data(adduct_rules, overwrite=T)
#   # # usethis::use_data(adducts, overwrite=T)
#   # # usethis::use_data(lmdb, overwrite=T)
#   # # usethis::use_data(adducts_no_deut, overwrite=T)
#   # #
#   # # #
#   # # # {
#   # # #   try({
#   # # #     parallel::stopCluster(session_cl)
#   # # #   },silent=T)
#   # # #
#   # #    session_cl <- parallel::makeCluster(3, outfile="")
#   # # #
#   # #   parallel::clusterExport(session_cl, c("smiles.to.iatom",
#   # #                                         "countAdductRuleMatches",
#   # #                                         "checkAdductRule",
#   # #                                         "doAdduct",
#   # #                                         "iatom.to.smiles",
#   # #                                         "smiles.to.iatom",
#   # #                                         "iatom.to.formula",
#   # #                                         "iatom.to.charge",
#   # #                                         "adduct_rules",
#   # #                                         "adducts",
#   # #                                         "doIsotopes",
#   # #                                         "isotopes"))
#   # # #
#   # #
#   # #     parallel::clusterEvalQ(cl = session_cl, expr = {
#   # #     library(data.table)
#   # #     library(enviPat)
#   # #     library(pbapply)
#   # #      library(ChemmineR)
#   # #   })
#   # # #
#   # #   dbname = "hmdb"
#   # #   buildExtDB(outfolder,
#   # #              base.dbname = dbname,
#   # #              cl = 0,#session_cl,
#   # #              blocksize = 200,
#   # #              mzrange = c(60,600),
#   # #              adduct_table = adducts,
#   # #              adduct_rules = adduct_rules,
#   # #              silent = F,
#   # #              ext.dbname = "extended",
#   # #              use.rules = T)
#   # #
#   # #   dbname = "chebi"
#   # #   buildExtDB(outfolder,
#   # #              base.dbname = dbname,
#   # #              cl = session_cl,
#   # #              blocksize = 200,
#   # #              mzrange = c(60,600),
#   # #              adduct_table = adducts,
#   # #              adduct_rules = adduct_rules,
#   # #              silent = F,
#   # #              ext.dbname = "extended_norules",
#   # #              use.rules = F)
#   # #
#   # #   dbname = "wikidata"
#   # #   buildExtDB(outfolder,
#   # #              base.dbname = dbname,
#   # #              cl = session_cl,
#   # #              blocksize = 200,
#   # #              mzrange = c(60,600),
#   # #              adduct_table = adducts,
#   # #              adduct_rules = adduct_rules,
#   # #              silent = F,
#   # #              ext.dbname = "extended_norules",
#   # #              use.rules = F)
#   # #
#   # # #
#   # # # dbs = gsub(basename(list.files(outfolder, pattern="\\.db$")), pattern = "\\.db", replacement="")
#   # # # dbs = dbs[dbs!="extended"]
#   # # #
#   # # # res = searchMZ(mzs = "110.071176455696",
#   # # #                ionmodes = "positive",
#   # # #                outfolder = outfolder,
#   # # #                base.dbname = dbs,
#   # # #                ppm = 5)
#   # # #
#   # # # "SELECT DISTINCT
#   # # # cpd.adduct as adduct,
#   # # # cpd.isoprevalence as isoprevalence,
#   # # # struc.smiles as structure,
#   # # # mz.mzmed as query_mz,
#   # # # (1e6*ABS(mz.mzmed - cpd.fullmz)/cpd.fullmz) AS dppm
#   # # # FROM mzvals mz
#   # # # JOIN mzranges rng ON rng.ID = mz.ID
#   # # # JOIN extended cpd INDEXED BY e_idx2
#   # # # ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
#   # # # JOIN adducts
#   # # # ON cpd.adduct = adducts.Name
#   # # # AND mz.foundinmode = adducts.Ion_Mode
#   # # # JOIN structures struc
#   # # # ON cpd.struct_id = struc.struct_id"
#   # # #
#   # # # dbs = list.files(outfolder, pattern = "\\.db$")
#   # # # dbs = gsub(dbs[dbs != "extended.db"], pattern="\\.db", replacement="")
#   # # #
#   # # # for(db in dbs){
#   # # #   print(db)
#   # # #   conn = openBaseDB(outfolder, db)
#   # # #   RSQLite::dbExecute(conn,"CREATE INDEX IF NOT EXISTS b_idx1 ON base(structure);")
#   # # #   RSQLite::dbDisconnect(conn)
#   # # # }
#   # #
#   #internal standards
#   # fn = "~/Documents/internal_umc_database_orig.txt"
#   # internal_standards = data.table::fread(fn)
#   # db.formatted <- data.table::data.table(
#   #    compoundname = as.character(internal_standards$CompoundName),
#   #    baseformula = as.character(internal_standards$Composition),
#   #    identifier = as.character(paste0("umcu", 1:nrow(internal_standards))),
#   #    structure = c(NA),
#   #    charge = c(0),
#   #    description = c("Internal database of UMCU")
#   # )
#   #
#   # db.formatted <- db.formatted[grepl(compoundname, pattern = "\\(IS\\)") & !grepl(compoundname, pattern = "Dimeric")]
#   # db.formatted$baseformula <- as.character(gsub(db.formatted$baseformula, pattern = "\\((\\d+)(\\w+)\\)", replacement="[\\1]\\2"))
#   # db.formatted$baseformula <- gsub("D", "[2]H", db.formatted$baseformula)
#   #
#   # data.table::fwrite(db.formatted, "/Users/jwolthuis/Desktop/internal_standards.csv")
#
#   # manually add smiles
#   # db.formatted <- data.table::fread("~/Documents/internal_standards_withstruct.csv")
#   #
#   # outfolder = "/Users/jwolthuis/Desktop"
#   # dbname="umc_internal"
#   # removeDB(outfolder, paste0(dbname,".db"))
#   # conn <- openBaseDB(outfolder, paste0(dbname,".db"))
#   #
#   # # NEEDS STRUCTURES
#   # db.final <- cleanDB(db.formatted,cl=0,silent=F,blocksize=50)
#   # db.final$structure <- as.character(db.final$structure)
#   # db.final$charge <- as.numeric(db.final$charge)
#   #
#   # writeDB(conn, db.final, "base")
#   # RSQLite::dbExecute(conn, "CREATE INDEX b_idx1 ON base(structure)")
#   # DBI::dbDisconnect(conn)
#   # #
#   # # # - - - WHICH ADDUCTS NEED DEUTERATED VARIANTS (ONLY DEDUCTION...) - - -
#   # #
#   # # dedH <- grepl(adducts$RemAt, pattern = "H\\d") | grepl(adducts$RemEx, pattern = "H\\d")
#   # # extra_adducts <- adducts[dedH,]
#   # #
#   # # adducts_deut <- data.table::fread("adducts_deut.csv")
#   # # # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   # #
#   # require(data.table)
#   # buildExtDB(outfolder,
#   #            base.dbname = dbname,
#   #            cl = 0,#session_cl,
#   #            blocksize = 20,
#   #            mzrange = c(60,600),
#   #            adduct_table = adducts,
#   #            adduct_rules = adduct_rules,
#   #            silent = F,
#   #            ext.dbname = "extended_is_only")
#   # #
#   # # # - - - - -
#   # #
#   # # formulas = getPredicted(mz = mz, ppm = 5)
#   # # results = searchFormulaWeb(formulas$baseformula,
#   # #                            search = c("chemspider"))
#   # # check.smi = TRUE
#   # # if(check.smi){
#   # #   mols = smiles.to.iatom(results_withsmi$structure)
#   # #   new.smi = iatom.to.smiles(mols)
#   # #   results_withsmi$structure <- new.smi
#   # # }
#   # # if(check.rules){
#   # #   rulematch = countAdductRuleMatches(mols, adduct_rules)
#   # #   structure.adducts.possible = checkAdductRule(rulematch,
#   # #                                                adduct_table)
#   # #   keep <- sapply(1:nrow(results_withsmi), function(i){
#   # #     adduct = results_withsmi[i, "adduct"][[1]]
#   # #     if(!is.na(adduct)){
#   # #       structure.adducts.possible[i, ..adduct][[1]]
#   # #     }
#   # #   })
#   # #   results_withsmi <- results_withsmi[keep,]
#   # # }
#   # # results_final <- rbind(results_nosmi,
#   # #                        results_withsmi)
#   # #
# }
