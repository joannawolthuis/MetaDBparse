# calculate 'rules' for all compounds (requires iatom-ization)
countAdductRuleMatches <- function(iatoms, adduct_rules){
  smiles = iatom.to.smiles(iatoms)
  res_cols <- lapply(1:nrow(adduct_rules), function(i){
    curr = adduct_rules$SHORT[i]
    query <- adduct_rules$SMARTS[i]
    all_matches <- lapply(iatoms, function(iatom){
      if(curr != "Nch"){
        all_matches = list(mapping=list(a=c(1:10)))
        try({
          rcdk::matches(query = query, target = iatom, T)
        })
      }else{
        rcdk::get.total.formal.charge(iatom)
      }
    })
    if(curr != "Nch"){
      vals = as.numeric(sapply(all_matches, function(x) length(x[[1]]$mapping)))
    }else{
      vals = as.numeric(unlist(all_matches))
    }
    add.col = data.table::data.table(vals)
    colnames(add.col) <- curr
    add.col
  })
  res = do.call("cbind", res_cols)
  cbind(structure=smiles, res)
}

checkAdductRule <- function(adduct_rule_scores, adduct_table){
  adduct.qualify.rows = lapply(1:nrow(adduct_table), function(i){
    row <- adduct_table[i,]
    name <- row$Name
    ion_mode <- row$Ion_mode

    # - - CHECK IF ADDUCT CAN EVEN BE FORMED BASED ON STRUCTURE - - -
    rules_raw = row$Rule
    rules_split = strsplit(rules_raw, "AND| AND ",)[[1]]

    qualified_per_rule = data.table::as.data.table(sapply(rules_split, function(rule){
      middle = if(grepl("<", rule)) "below" else if(grepl(">", rule)) "above" else "equals"
      leftright = strsplit(rule, ">|<|=")[[1]]
      left_val = leftright[1]
      right = as.numeric(leftright[2])
      left = as.numeric(adduct_rule_scores[,..left_val][[1]])
      qualifies = switch(middle,
                         below = left < right,
                         equals = left == right,
                         above = left > right)
    }))
    qualified.rule = apply(qualified_per_rule, MARGIN = 1, all)
    add.row = data.table::data.table(qualified.rule)
    colnames(add.row) <- row$Name
    add.row
  })
  qualified.per.adduct <- do.call("cbind", adduct.qualify.rows)
  data.table::data.table(structure = adduct_rule_scores$structure,
                         qualified.per.adduct)
}

doAdduct <- function(structure, formula, charge, adduct_table, query_adduct){
  # query_adduct = "[M+H]1+"
  row = adduct_table[Name == query_adduct,]
  name <- row$Name
  ion_mode <- row$Ion_mode

  unique_formulas <- unique(data.table::data.table(
    structure = structure,
    baseformula = formula,
    charge = charge))

  adduct_before <- row$AddAt
  deduct_before <- row$RemAt
  adduct_before <- if(is.na(adduct_before)) FALSE else adduct_before
  deduct_before <- if(is.na(deduct_before)) FALSE else deduct_before

  if(!(adduct_before %in% c("", "FALSE", FALSE, NA))){
    unique_formulas$adducted <- enviPat::mergeform(formula1 = unique_formulas$baseformula,
                                                   formula2 = adduct_before)
  }else{
    unique_formulas$adducted <- unique_formulas$baseformula
  }

  # --- is deduction necessary? ---

  if(!(deduct_before %in% c("", "FALSE", FALSE, NA))){
    can.deduct <- which(!as.logical(enviPat::check_ded(formulas = unique_formulas$adducted,
                                                       deduct = deduct_before)))
    if(length(can.deduct) == 0) return(data.table::data.table())
    deductibles <- unique_formulas$adducted[can.deduct]
    unique_formulas$adducted[can.deduct] <- enviPat::subform(unique_formulas$adducted[can.deduct],
                                                             deduct_before)
    unique_formulas <- unique_formulas[can.deduct]
  }

  # --- multiplication ---
  multiplier <- as.numeric(row$xM)

  if(multiplier > 1){
    unique_formulas$adducted <- enviPat::multiform(unique_formulas$adducted,
                                                   multiplier)
  }

  # --- adduct after multiplication ---
  adduct_after <- row$AddEx
  deduct_after <- row$RemEx

  adduct_after <- if(is.na(adduct_after)) FALSE else adduct_after
  deduct_after <- if(is.na(deduct_after)) FALSE else deduct_after

  if(!(adduct_after %in% c("", "FALSE", FALSE, NA))){
    unique_formulas$adducted <- enviPat::mergeform(formula1 = unique_formulas$adducted,
                                                   formula2 = adduct_after)}

  # --- is deduction necessary? ---
  if(!(deduct_after  %in% c("", "FALSE", FALSE, NA))){
    can.deduct <- which(!as.logical(enviPat::check_ded(formulas = unique_formulas$adducted,
                                                       deduct = deduct_after)))

    if(length(can.deduct) == 0) return(data.table::data.table())
    unique_formulas$adducted[can.deduct] <- enviPat::subform(unique_formulas$adducted[can.deduct], deduct_after)
    unique_formulas <- unique_formulas[can.deduct]
  }

  unique_formulas$final <- unique_formulas$adducted

  # - - - fixing - -
  if(nrow(unique_formulas) == 0) return(data.table::data.table())
  unique_formulas$final.charge <- c(as.numeric(unique_formulas$charge)) + c(as.numeric(row$Charge))
  unique_formulas <- unique_formulas[final.charge != 0]
  if(nrow(unique_formulas) == 0) return(data.table::data.table())

  unique_formulas
}

doIsotopes <- function(formula, charge){
  #require(enviPat)
  #data(isotopes)
  # - - - - - -
  isotables <- enviPat::isopattern(
    isotopes,
    formula,
    threshold = 0.1,
    plotit = FALSE,
    charge = charge,
    verbose = FALSE,
  )

  isolist <- lapply(isotables, function(isotable){
    if(isotable[[1]] == "error"){
      return(data.table())
    }
    iso.dt <- data.table::data.table(isotable, fill=TRUE)
    result <- iso.dt[,1:2]
    names(result) <- c("fullmz", "isoprevalence")
    # --- return ---
    result
  })

  isolist.nonas <- isolist[!is.na(isolist)]
  isotable <- rbindlist(isolist.nonas)
  keep.isos <- names(isolist.nonas)

  charges <- charge[!is.na(isolist)]

  repeat.times <- c(unlist(lapply(isolist.nonas, FUN=function(list) nrow(list))))
  isotable$final <- rep(keep.isos, repeat.times)
  isotable$final.charge <- rep(charges, repeat.times)

  # - - -

  isotable

}

buildExtDB <- function(outfolder,
                       ext.dbname = "extended",
                       base.dbname,
                       cl = 0,
                       blocksize = 600,
                       mzrange = c(60,600),
                       adduct_table = adducts,
                       adduct_rules = adduct_rules,
                       silent = silent){
  # A
  outfolder <- normalizePath(outfolder)

  print(paste("Will calculate adducts + isotopes for the", base.dbname, "database."))
  full.db <- file.path(outfolder, paste0(ext.dbname, ".db"))
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA foreign_keys = ON"))
  RSQLite::dbExecute(full.conn, "CREATE TABLE IF NOT EXISTS structures(
                                 struct_id INT PRIMARY KEY,
                                 smiles TEXT,
                                 UNIQUE(struct_id, smiles))")
  if(RSQLite::dbExistsTable(full.conn, "extended")){
    new_adducts <- setdiff(adduct_table$Name,
                           RSQLite::dbGetQuery(full.conn,
                                               "SELECT DISTINCT adduct FROM extended")[,1])
  }else{
    new_adducts <- adduct_table$Name
  }

  RSQLite::dbExecute(full.conn, "DROP TABLE IF EXISTS adducts")
  RSQLite::dbWriteTable(conn, "adducts", adduct_table)


  RSQLite::dbExecute(full.conn, strwrap("CREATE TABLE IF NOT EXISTS extended(
                                         struct_id text,
                                         fullmz decimal(30,13),
                                         adduct text,
                                         isoprevalence float)", width=10000, simplify=TRUE))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA auto_vacuum = 1;"))

  # B
  first.db <- if(RSQLite::dbGetQuery(full.conn, "SELECT COUNT(*) FROM extended")[1,1] == 0) TRUE else FALSE
  base.db <- normalizePath(file.path(outfolder, paste0(base.dbname, ".db")))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))

  if(first.db){
    RSQLite::dbExecute(full.conn, "CREATE INDEX st_idx1 ON structures(smiles)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx1 on extended(struct_id)")
    RSQLite::dbExecute(full.conn, "CREATE INDEX e_idx2 on extended(fullmz)")
    RSQLite::dbExecute(full.conn, "PRAGMA journal_mode=WAL;")
    }
  if(first.db | length(new_adducts) > 0){
    to.do = RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                            FROM tmp.base")
    to.do$struct_id <- c(NA)
  }else{
    to.do = RSQLite::dbGetQuery(full.conn, "SELECT DISTINCT baseformula, structure, charge
                                            FROM tmp.base LEFT JOIN structures str
                                            ON base.structure = str.smiles
                                            WHERE str.smiles IS NULL")}
  if(nrow(to.do) == 0){
    print("all already done")
    RSQLite::dbDisconnect(full.conn)
    return(NULL)
  }

  done.structures = RSQLite::dbGetQuery(full.conn, "SELECT MAX(struct_id) FROM structures")[,1]
  start.id = done.structures + 1

  # new will be written
  mapper = data.table::data.table(struct_id = seq(start.id, start.id + nrow(to.do) - 1, 1),
                                  smiles = to.do$structure,
                                  baseformula = to.do$baseformula,
                                  charge = to.do$charge)

  print(paste("Parsing", nrow(mapper), "new compounds..."))

  blocks = split(mapper, ceiling(seq_along(1:nrow(mapper))/blocksize))

  #blocks = split(1:nrow(mapper), ceiling(seq_along(1:nrow(mapper))/blocksize))

  RSQLite::dbDisconnect(full.conn)

  # - - - - L O O P  T H I S - - - - -
  if(length(new_adducts) > 0){
    adduct_table <- data.table::as.data.table(adduct_table)[Name %in% new_adducts,]
  }

  # completely new compounds
  per.adduct.tables = pbapply::pblapply(1:length(blocks), cl=cl, function(i, mzrange=mzrange, silent=silent){

    full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
    RSQLite::dbExecute(full.conn, "PRAGMA busy_timeout=5000;")

    # for each block
    block = blocks[[i]]

    # convert to iatom
    iatoms <- smiles.to.iatom(block$smiles, silent=silent)

    # check which are invalid structures
    bad.structures = which(sapply(iatoms, is.null))
    good.structures = which(!sapply(iatoms, is.null))

    # = = = FOR THE GOOD STRUCTURES = = =
    parsable.iats = iatoms[good.structures]
    structure.reactive.groups = countAdductRuleMatches(parsable.iats,
                                                       adduct_rules)
    structure.adducts.possible = checkAdductRule(structure.reactive.groups,
                                                 adducts)

    # now loop through adducts
    per.adduct.tables = lapply(adduct_table$Name, function(add){# make cl session_cl later

      has.structure.can.adduct = block[good.structures[unlist(structure.adducts.possible[,..add])]]
      no.structure = block[bad.structures]

      block.calc.adduct = rbind(has.structure.can.adduct, no.structure)

      if(nrow(block.calc.adduct) == 0) return(data.table::data.table())

      adducted = doAdduct(structure = block.calc.adduct$smiles,
                          formula = block.calc.adduct$baseformula,
                          charge = block.calc.adduct$charge,
                          adduct_table = adduct_table,
                          query_adduct = add)

      if(nrow(adducted) == 0) return(data.table::data.table())

      # check mass
      checked <- enviPat::check_chemform(isotopes, adducted$final)
      checked$mz <- checked$monoisotopic_mass/abs(adducted$final.charge)

      keepers <- which(checked$mz %between% mzrange)

      if(length(keepers) == 0){
        return(data.table::data.table())
      }else{
        isotable <- doIsotopes(formula = adducted[keepers,]$final,
                               charge = adducted[keepers,]$final.charge)

        formula_plus_iso <- merge(adducted, isotable,
                                  by = c("final", "final.charge"),
                                  allow.cartesian = T)
        adducted.plus.isotopes <- merge(adducted, formula_plus_iso,
                                        by=c("baseformula", "charge"),
                                        allow.cartesian = T)

        # ==== MAKE FINAL TABLE ====

        meta.table <- data.table::data.table(fullmz = adducted.plus.isotopes$fullmz,
                                             adduct = c(add),
                                             isoprevalence = adducted.plus.isotopes$isoprevalence,
                                             structure = adducted.plus.isotopes$structure.x)

        # map SMILES to smile_id
        ids <- mapper$struct_id[match(meta.table$structure,
                                      mapper$smiles)]
        meta.table$struct_id <- ids

        # === RETURN ===

        return(unique(meta.table[,-"structure"]))
    }})

    to.write = data.table::rbindlist(per.adduct.tables)

    if(nrow(to.write)>0){
      repeat{
        rv <- try({
          RSQLite::dbAppendTable(full.conn,
                                 "extended",
                                 to.write)
          RSQLite::dbAppendTable(full.conn,
                                 "structures",
                                 unique(blocks[[i]][,c("struct_id", "smiles")]))
        },silent = silent)
        if(!is(rv, "try-error")){
          break
        }
      }
    }
    RSQLite::dbDisconnect(full.conn)
  },mzrange=mzrange, silent=silent)
}
