#' @title Weights of all isotopes of all elements used
#' @description Sourced from 'enviPat' package
"isotopes"

#' @title Remove structures where isotope generation failed.
#' @description Sometimes if a run crashes, or a structure is bugged somehow, it is still registered as 'done' in the extended database and cannot be redone. This function removes these structures. Warning: slow!
#' @param outfolder Which folder are your databases in?
#' @param ext.dbname Extended database name (without .db suffix), Default: 'extended'
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#' @rdname removeFailedStructures
#' @export
#' @importFrom RSQLite dbConnect SQLite dbExecute dbDisconnect
removeFailedStructures <- function(outfolder,
                                   ext.dbname = "extended"){
  outfolder <- normalizePath(outfolder)
  print(paste("Will remove failed isotope calculations, this may take a while. Removing structures from structure list that didn't make it into isotope calculation."))
  full.db <- file.path(outfolder, paste0(ext.dbname, ".db"))
  first.db <- !file.exists(full.db)
  if(first.db){
    print("Can't do this on a new db! :-(")
    return(NULL)
  }else{
    full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
    RSQLite::dbExecute(full.conn,
                       "CREATE INDEX IF NOT EXISTS ext_id ON extended(struct_id)")
    RSQLite::dbExecute(full.conn,
                       "DELETE FROM structures
                        WHERE struct_id IN(SELECT structures.struct_id
                                           FROM structures
                                           LEFT JOIN extended
                                           ON structures.struct_id=extended.struct_id
                                           WHERE extended.struct_id IS NULL)"
    )
    RSQLite::dbDisconnect(full.conn)
  }
}

#' @title Check which structures are OK according to given adduct rules
#' @description Calculate 'rules' for all compounds (requires iatom-ization)
#' @param iatoms Iatomcontainers with compounds
#' @param adduct_rules Adduct rule table (default is data(adduct_rules))
#' @return Table with all structures and if they pass the rules given for each adduct
#' @seealso
#'  \code{\link[rcdk]{matches}},\code{\link[rcdk]{get.total.formal.charge}}
#' @rdname countAdductRuleMatches
#' @export
#' @examples
#'  iatom = smiles.to.iatom(c('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'))
#'  data(adduct_rules)
#'  addScore <- countAdductRuleMatches(iatom, adduct_rules = adduct_rules)
#' @importFrom rcdk matches get.total.formal.charge
#' @importFrom data.table data.table
countAdductRuleMatches <- function(iatoms, adduct_rules) {
  smiles <- iatom.to.smiles(iatoms)
  res_cols <- lapply(1:nrow(adduct_rules), function(i) {
    curr <- adduct_rules$SHORT[i]
    query <- adduct_rules$SMARTS[i]
    all_matches <- lapply(iatoms, function(iatom) {
      if (curr != "Nch") {
        all_matches <- list(mapping = list(a = c(1:10)))
        try({
          rcdk::matches(query = query, target = iatom, TRUE)
        })
      }
      else {
        rcdk::get.total.formal.charge(mol = iatom)
      }
    })
    if (curr != "Nch") {
      vals <- as.numeric(sapply(all_matches, function(x) length(x[[1]]$mapping)))
    }
    else {
      vals <- as.numeric(unlist(all_matches))
    }
    add.col <- data.table::data.table(vals)
    colnames(add.col) <- curr
    add.col
  })
  res <- do.call("cbind", res_cols)
  cbind(structure = smiles, res)
}

#' @title Check for combined adduct rules
#' @description Sometimes multiple rules apply - this function checks if they all apply as noted in the rule table.
#' @param adduct_rule_scores Scores from countAdductRuleMatches.
#' @param adduct_table Adduct table
#' @return Table with TRUE/FALSE for each structure and adduct
#' @seealso
#'  \code{\link[data.table]{as.data.table}}
#' @rdname checkAdductRule
#' @export
#' @examples
#'  iatom = smiles.to.iatom(c('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'))
#'  data(adduct_rules)
#'  data(adducts)
#'  addScore <- countAdductRuleMatches(iatom, adduct_rules = adduct_rules)
#'  checkAdductRule(addScore, adduct_table = adducts)
#' @importFrom data.table as.data.table data.table
checkAdductRule <- function(adduct_rule_scores, adduct_table) {
  left_val <- NULL
  adduct.qualify.cols <- lapply(1:nrow(adduct_table), function(i) {
    row <- adduct_table[i, ]
    name <- row$Name
    ion_mode <- row$Ion_mode
    rules_raw <- row$Rule
    rules_split <- strsplit(rules_raw, "AND| AND ", )[[1]]
    qualified_per_rule <- data.table::as.data.table(sapply(rules_split, function(rule) {
      middle <- if (grepl("<", rule)) {
        "below"
      } else if (grepl(">", rule)) {
        "above"
      } else {
        "equals"
      }
      leftright <- strsplit(rule, ">|<|=")[[1]]
      left_val <- leftright[1]
      right <- as.numeric(leftright[2])
      left <- as.numeric(adduct_rule_scores[, left_val, with = FALSE][[1]])
      qualifies <- switch(middle, below = left < right, equals = left == right, above = left > right)
      data.table::as.data.table(qualifies)
    }))
    qualified.rule <- apply(qualified_per_rule, MARGIN = 1, FUN = all)
    add.row <- data.table::data.table(qualified.rule)
    colnames(add.row) <- row$Name
    add.row
  })
  qualified.per.adduct <- do.call("cbind", adduct.qualify.cols)
  data.table::data.table(structure = adduct_rule_scores$structure, qualified.per.adduct)
}

#' @title Generate adduct for given structure
#' @description Takes in formula, an adduct of interest, and returns adduct formulas and charges.
#' @param structure SMILES structure
#' @param formula Molecular formula
#' @param charge Initial charge
#' @param adduct_table Adduct table
#' @param query_adduct Adduct 'Name' of interest
#' @return Table with adducts of this compound
#' @seealso
#'  \code{\link[enviPat]{check_chemform}},\code{\link[enviPat]{mergeform}},\code{\link[enviPat]{check_ded}},\code{\link[enviPat]{subform}},\code{\link[enviPat]{multiform}}
#' @rdname doAdduct
#' @export
#' @examples
#'  data(adduct_rules)
#'  data(adducts)
#'  structure = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'
#'  doAdduct(structure = structure, formula="C6H12O6", charge=0,
#'  adduct_table=adducts, query_adduct="[M+H]1+")
#' @importFrom enviPat check_chemform mergeform check_ded subform multiform
#' @importFrom data.table data.table
doAdduct <- function(structure, formula, charge, adduct_table, query_adduct) {
  Name <- final.charge <- NULL
  row <- adduct_table[Name == query_adduct, ]
  name <- row$Name
  ion_mode <- row$Ion_mode
  checked <- enviPat::check_chemform(isotopes, chemforms = formula)
  phailed <- which(checked$warning)
  if (all(checked$warning)) {
    return(data.table::data.table())
  }
  unique_formulas <- unique(data.table::data.table(structure = structure, baseformula = checked$new_formula, charge = charge))
  if (length(phailed) > 0) {
    unique_formulas <- unique_formulas[-phailed, ]
  }
  adduct_before <- row$AddAt
  deduct_before <- row$RemAt
  adduct_before <- if (is.na(adduct_before)) {
    FALSE
  } else {
    adduct_before
  }
  deduct_before <- if (is.na(deduct_before)) {
    FALSE
  } else {
    deduct_before
  }
  if (!(adduct_before %in% c("", "FALSE", FALSE, NA))) {
    unique_formulas$adducted <- enviPat::mergeform(formula1 = unique_formulas$baseformula, formula2 = adduct_before)
  }
  else {
    unique_formulas$adducted <- unique_formulas$baseformula
  }
  if (!(deduct_before %in% c("", "FALSE", FALSE, NA))) {
    can.deduct <- which(!as.logical(enviPat::check_ded(formulas = unique_formulas$adducted, deduct = deduct_before)))
    if (length(can.deduct) == 0) {
      return(data.table::data.table())
    }
    deductibles <- unique_formulas$adducted[can.deduct]
    unique_formulas$adducted[can.deduct] <- enviPat::subform(deductibles, deduct_before)
    unique_formulas <- unique_formulas[can.deduct]
  }
  unique_formulas <- unique_formulas[which(unique_formulas$adducted != "NANA"), ]
  multiplier <- as.numeric(row$xM)
  if (multiplier > 1) {
    unique_formulas$adducted <- enviPat::multiform(unique_formulas$adducted, multiplier)
  }
  adduct_after <- row$AddEx
  deduct_after <- row$RemEx
  adduct_after <- if (is.na(adduct_after)) {
    FALSE
  } else {
    adduct_after
  }
  deduct_after <- if (is.na(deduct_after)) {
    FALSE
  } else {
    deduct_after
  }
  if (!(adduct_after %in% c("", "FALSE", FALSE, NA))) {
    unique_formulas$adducted <- enviPat::mergeform(formula1 = unique_formulas$adducted, formula2 = adduct_after)
  }
  if (!(deduct_after %in% c("", "FALSE", FALSE, NA))) {
    can.deduct <- which(!as.logical(enviPat::check_ded(formulas = unique_formulas$adducted, deduct = deduct_after)))
    if (length(can.deduct) == 0) {
      return(data.table::data.table())
    }
    unique_formulas$adducted[can.deduct] <- enviPat::subform(unique_formulas$adducted[can.deduct], deduct_after)
    unique_formulas <- unique_formulas[can.deduct]
  }
  unique_formulas$final <- unique_formulas$adducted
  if (nrow(unique_formulas) == 0) {
    return(data.table::data.table())
  }
  unique_formulas$final.charge <- c(as.numeric(unique_formulas$charge)) + c(as.numeric(row$Charge))
  unique_formulas <- unique_formulas[final.charge != 0]
  if (nrow(unique_formulas) == 0) {
    return(data.table::data.table())
  }
  unique_formulas
}

#' @title Generate isotopes for given formula
#' @description Takes in formula and returns isotope pattern m/z values.
#' @param formula Molecular formula
#' @param charge Final charge
#' @return Table with isotopes of this molecular formula
#' @examples
#'  doIsotopes(formula="C6H12O6", charge=0)
#' @seealso
#'  \code{\link[enviPat]{isopattern}}
#' @rdname doIsotopes
#' @export
#' @importFrom enviPat isopattern
#' @importFrom data.table data.table
doIsotopes <- function(formula, charge) {
  isotables <- enviPat::isopattern(isotopes, formula, threshold = 0.1, plotit = FALSE, charge = charge, verbose = FALSE, )
  isolist <- lapply(isotables, function(isotable) {
    if (isotable[[1]] == "error") {
      return(data.table())
    }
    iso.dt <- data.table::data.table(isotable, fill = TRUE)
    result <- iso.dt[, 1:2]
    names(result) <- c("fullmz", "isoprevalence")
    result
  })
  isolist.nonas <- isolist[!is.na(isolist)]
  isotable <- data.table::rbindlist(isolist.nonas)
  keep.isos <- names(isolist.nonas)
  charges <- charge[!is.na(isolist)]
  repeat.times <- c(unlist(lapply(isolist.nonas, FUN = function(list) nrow(list))))
  isotable$final <- rep(keep.isos, repeat.times)
  isotable$final.charge <- rep(charges, repeat.times)
  isotable
}

#' @title Build external database using a given base database
#' @description Wrapper function that takes a base database, an existing (or not) external database, and fills the extended database with adduct and isotope variants of the compounds in the base database.
#' @param outfolder Which folder are your databases in?
#' @param ext.dbname Extended database name (without .db suffix), Default: 'extended'
#' @param base.dbname Base database name (without .db suffix)
#' @param cl parallel::makeCluster object for multithreading, Default: 0
#' @param blocksize How many compounds to process simultanaously? Higher means more memory spikes but faster building, Default: 600
#' @param mzrange Range of m/zs to include in database, Default: c(60, 600)
#' @param adduct_table Adduct table, Default: adducts
#' @param adduct_rules Adduct rule table, Default: adduct_rules
#' @param silent Silence warnings?, Default: silent
#' @param use.rules Use adduct rules?, Default: TRUE
#' @seealso
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[gsubfn]{fn}}
#'  \code{\link[data.table]{as.data.table}},\code{\link[data.table]{rbindlist}},\code{\link[data.table]{fwrite}},\code{\link[data.table]{fread}}
#'  \code{\link[DBI]{dbWriteTable}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[enviPat]{check_chemform}}
#' @rdname buildExtDB
#' @export
#' @examples
#'  \dontrun{myFolder = tempdir()}
#'  \dontrun{buildBaseDB(outfolder = myFolder, "lmdb", test = TRUE)}
#'  \dontrun{file.remove(file.path(myFolder, "extended.db"))}
#'  \dontrun{data(adducts)}
#'  \dontrun{data(adduct_rules)}
#'  \dontrun{buildExtDB(outfolder = myFolder, base.dbname = "lmdb",
#'  silent=FALSE, adduct_table = adducts, adduct_rules = adduct_rules)}
#' @importFrom RSQLite dbConnect SQLite dbExecute dbExistsTable dbGetQuery dbDisconnect dbWriteTable dbReadTable
#' @importFrom gsubfn fn
#' @importFrom data.table as.data.table data.table rbindlist fwrite fread
#' @importFrom DBI dbWriteTable
#' @importFrom pbapply pblapply pbsapply
#' @importFrom enviPat check_chemform
buildExtDB <- function(outfolder, ext.dbname = "extended", base.dbname, cl = 0, blocksize = 600, mzrange = c(60, 600), adduct_table = adducts, adduct_rules = adduct_rules, silent = silent, use.rules = TRUE) {
  Name <- charge <- ..add <- NULL
  outfolder <- normalizePath(outfolder)
  print(paste("Will calculate adducts + isotopes for the", base.dbname, "database."))
  full.db <- file.path(outfolder, paste0(ext.dbname, ".db"))
  first.db <- !file.exists(full.db)
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  tempdir <- file.path(outfolder, paste0(ext.dbname, "_inprogress"))
  if (dir.exists(tempdir)) {
    tmpfiles.ext <- list.files(tempdir, full.names = TRUE, pattern = "_ext_")
    tmpfiles.struct <- list.files(tempdir, full.names = TRUE, pattern = "_str_")
    if (length(tmpfiles.ext) > 0) {
      prevDB <- gsub("_.*$|\\.csv", "", basename(tmpfiles.ext[1]))
      print(paste("! ! ! failed to finish", prevDB, "database ! ! !"))
      print("recovery mode, loading done structures into database... please stand by o wo ")
      pbapply::pbsapply(1:length(tmpfiles.ext), function(i) {
        try({
          extended_csv <- data.table::fread(tmpfiles.ext[i])
          RSQLite::dbWriteTable(full.conn, "extended", extended_csv, append = TRUE)
          if (length(tmpfiles.struct > 0)) {
            structures_csv <- data.table::fread(tmpfiles.struct[i])
            RSQLite::dbWriteTable(full.conn, "structures", structures_csv, append = TRUE)
          }
        })
      })
    }
    unlink(tempdir, recursive = TRUE)
  }
  dir.create(tempdir)
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA foreign_keys = ON"))
  RSQLite::dbExecute(full.conn, "CREATE TABLE IF NOT EXISTS structures(struct_id INT PRIMARY KEY,
                                                                       smiles TEXT,
                                                                       UNIQUE(struct_id, smiles))")
  if (RSQLite::dbExistsTable(full.conn, "extended")) {
    done_adducts <- RSQLite::dbGetQuery(full.conn, gsubfn::fn$paste("SELECT DISTINCT Name FROM adducts
                                                                    WHERE dbname == '$base.dbname'"))[, 1]
    new_adducts <- setdiff(adduct_table$Name, done_adducts)
  }
  else {
    new_adducts <- adducts$Name
  }
  RSQLite::dbExecute(full.conn, strwrap("CREATE TABLE IF NOT EXISTS extended(
                                        struct_id INT,
                                        fullformula text,
                                        finalcharge text,
                                        fullmz decimal(30,13),
                                        adduct text,
                                        isoprevalence float)", width = 10000, simplify = TRUE))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA auto_vacuum = 1;"))
  base.db <- normalizePath(file.path(outfolder, paste0(base.dbname, ".db")))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))
  if (first.db) {
    RSQLite::dbExecute(full.conn, "CREATE INDEX st_idx1 ON structures(smiles)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx1 on extended(struct_id)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx2 on extended(fullmz)")
    RSQLite::dbExecute(full.conn, "PRAGMA journal_mode=WAL;")
  }
  if (first.db) {
    to.do <- RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                             FROM tmp.base")
    to.do$struct_id <- c(NA)
    adduct_only <- FALSE
  }
  else {
    to.do <- RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                             FROM tmp.base LEFT JOIN structures str
                                             ON base.structure = str.smiles
                                             WHERE str.smiles IS NULL")
    if (length(new_adducts) > 0) {
      adduct_only <- if (nrow(to.do) == 0) {
        TRUE
      } else {
        FALSE
      }
      if (adduct_only) {
        print(paste0("New adducts: ", paste0(new_adducts, collapse = ",")))
        to.do <- RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula,
                                                 structure, charge
                                                 FROM tmp.base")
        to.do$struct_id <- c(NA)
        adduct_table <- data.table::as.data.table(adduct_table)[Name %in% new_adducts, ]
      }
    }
    else {
      adduct_only <- FALSE
    }
  }
  if (nrow(to.do) == 0) {
    print("all already done")
    RSQLite::dbDisconnect(full.conn)
    return(NULL)
  }
  done.structures <- if (first.db) {
    0
  } else {
    RSQLite::dbGetQuery(full.conn, "SELECT MAX(struct_id) FROM structures")[, 1]
  }
  to.do <- to.do[!grepl(pattern = " ", to.do$baseformula), ]
  start.id <- done.structures + 1
  mapper <- data.table::data.table(struct_id = seq(start.id, start.id + nrow(to.do) - 1, 1), smiles = to.do$structure, baseformula = to.do$baseformula, charge = as.numeric(to.do$charge))
  if (adduct_only) {
    print("Getting structure IDs...")
    mapper_now <- unique(mapper[, -"struct_id"])
    DBI::dbWriteTable(full.conn, "mapper", mapper_now, overwrite = TRUE)
    mapper_filled <- data.table::as.data.table(RSQLite::dbGetQuery(full.conn, "SELECT * FROM structures s JOIN mapper m ON s.smiles = m.smiles"))
    mapper <- mapper_filled[, c("struct_id", "smiles", "baseformula", "charge")]
  }
  mapper[is.na(charge)]$charge <- c(0)
  print(paste("Parsing", nrow(mapper), "new compounds..."))
  blocks <- split(mapper, ceiling(seq_along(1:nrow(mapper)) / blocksize))
  tmpfiles.ext <- sapply(1:length(blocks), function(i) file.path(tempdir, paste0(base.dbname, "_ext_", i, ".csv")))
  tmpfiles.struct <- sapply(1:length(blocks), function(i) file.path(tempdir, paste0(base.dbname, "_str_", i, ".csv")))
  RSQLite::dbDisconnect(full.conn)
  per.adduct.tables <- pbapply::pblapply(1:length(blocks), cl = cl, function(i, mzrange = mzrange, silent = silent, blocks = blocks, adduct_table = adduct_table, adduct_rules = adduct_rules, tmpfiles.ext = tmpfiles.ext, tmpfiles.struct = tmpfiles.struct, mapper = mapper, use.rules = use.rules) {
    print(i)
    block <- blocks[[i]]
    deut <- grep(block$baseformula, pattern = "D")
    if (length(deut) > 0) {
      block$baseformula[deut] <- gsub(block$baseformula[deut], pattern = "D", replacement = "[2]H")
    }
    if (use.rules) {
      iatoms <- smiles.to.iatom(block$smiles, silent = silent)
      bad.structures <- which(sapply(iatoms, is.null))
      good.structures <- which(!sapply(iatoms, is.null))
    }
    else {
      iatoms <- c()
      bad.structures <- 1:length(iatoms)
      good.structures <- c()
    }
    structure.adducts.possible <- data.table::data.table()
    if (length(good.structures) > 0 | use.rules) {
      parsable.iats <- iatoms[good.structures]
      structure.reactive.groups <- countAdductRuleMatches(parsable.iats, adduct_rules)
      structure.adducts.possible <- checkAdductRule(structure.reactive.groups, adduct_table)
    }
    per.adduct.tables <- lapply(adduct_table$Name, function(add) {
      has.structure.can.adduct <- data.table::data.table()
      if (use.rules) {
        if (nrow(structure.adducts.possible) > 0) {
          addCol <- unlist(structure.adducts.possible[, add, with = FALSE])
          has.structure.can.adduct <- block[good.structures[addCol]]
        }
        no.structure <- block[bad.structures]
        block.calc.adduct <- rbind(has.structure.can.adduct, no.structure)
      }
      else {
        block.calc.adduct <- block
      }
      if (nrow(block.calc.adduct) == 0) {
        return(data.table::data.table())
      }
      adducted <- doAdduct(structure = block.calc.adduct$smiles, formula = block.calc.adduct$baseformula, charge = block.calc.adduct$charge, adduct_table = adduct_table, query_adduct = add)
      if (nrow(adducted) == 0) {
        return(data.table::data.table())
      }
      checked <- enviPat::check_chemform(isotopes, adducted$final)
      checked$mz <- checked$monoisotopic_mass / abs(adducted$final.charge)
      keepers <- which(checked$mz %between% mzrange)
      if (length(keepers) == 0) {
        return(data.table::data.table())
      }
      else {
        isotable <- doIsotopes(formula = adducted[keepers, ]$final, charge = adducted[keepers, ]$final.charge)
        formula_plus_iso <- merge(adducted, isotable, by = c("final", "final.charge"), allow.cartesian = TRUE)
        adducted.plus.isotopes <- merge(adducted, formula_plus_iso, by = c("baseformula", "charge", "adducted", "final.charge", "final", "structure"), allow.cartesian = TRUE)
        meta.table <- data.table::data.table(fullmz = adducted.plus.isotopes$fullmz, fullformula = adducted.plus.isotopes$final, finalcharge = adducted.plus.isotopes$final.charge, adduct = c(add), isoprevalence = adducted.plus.isotopes$isoprevalence, structure = adducted.plus.isotopes$structure)
        ids <- mapper$struct_id[match(meta.table$structure, mapper$smiles)]
        meta.table$struct_id <- ids
        return(unique(meta.table[, -"structure"]))
      }
    })
    to.write <- data.table::rbindlist(per.adduct.tables)
    if (nrow(to.write) > 0) {
      structs <- unique(blocks[[i]][, c("struct_id", "smiles")])
      data.table::fwrite(to.write[, c("struct_id", "fullformula", "finalcharge", "fullmz", "adduct", "isoprevalence")], file = tmpfiles.ext[i], append = FALSE)
      if (!adduct_only) {
        data.table::fwrite(structs, file = tmpfiles.struct[i], append = FALSE)
      }
    }
  }, mzrange = mzrange, silent = silent, blocks = blocks, adduct_table = adduct_table, adduct_rules = adduct_rules, tmpfiles.ext = tmpfiles.ext, tmpfiles.struct = tmpfiles.struct, mapper = mapper, use.rules = use.rules)
  print("loading into database... please stand by - w- ")
  pbapply::pbsapply(1:length(blocks), function(i) {
    full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
    try({
      extended_csv <- data.table::fread(tmpfiles.ext[i])
      structures_csv <- data.table::fread(tmpfiles.struct[i])
      RSQLite::dbWriteTable(full.conn, "extended", extended_csv, append = TRUE)
      if (!adduct_only) {
        RSQLite::dbWriteTable(full.conn, "structures", structures_csv, append = TRUE)
      }
    })
    RSQLite::dbDisconnect(full.conn)
  })
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  adduct_table$dbname <- base.dbname
  if (!first.db) {
    existing_adducts <- RSQLite::dbReadTable(full.conn, "adducts")
    adduct_tbl_named <- unique(rbind(adduct_table, existing_adducts))
  }
  else {
    adduct_tbl_named <- adduct_table
  }
  RSQLite::dbWriteTable(full.conn, "adducts", as.data.frame(adduct_tbl_named), overwrite = TRUE)
  RSQLite::dbDisconnect(full.conn)
  unlink(c(tmpfiles.ext, tmpfiles.struct))
}
