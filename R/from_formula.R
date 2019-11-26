getFormula <- function(mz, add_name, adducts, ppm, elements = c("C","H","N","O","P","S")
){
  def.ele = elements
  row <- adducts[Name == add_name]
  add.ele <- if(!is.na(row$AddEx)){
    ele = gsub(stringr::str_match_all(row$AddEx, pattern="[A-Z][a-z]?\\d*|\\(.*?\\)\\d+")[[1]][,1], pattern="\\d", replacement="")
  }else{
    c()
  }
  filter="."
  add.only.ele <- setdiff(add.ele,
                          def.ele)

  total.ele <- unique(c(def.ele,
                        add.ele))

  total.ele <- total.ele[total.ele != ""]

  # get which formulas are possible
  predicted = Rdisop::decomposeMass(as.numeric(mz),
                                    ppm = ppm,
                                    elements = Rdisop::initializeElements(names = total.ele),
                                    z = row$Charge
  )
  predicted
}

filterFormula = function(predicted, rules){
  require(enviPat)
  data(isotopes)

  corrected = enviPat::check_chemform(isotopes = isotopes, chemforms = predicted$formula)

  deconstructed = data.table::data.table(
    nrC = as.numeric(stringr::str_match(corrected$new_formula, pattern = "C(\\d*)")[,2]),
    nrH = as.numeric(stringr::str_match(corrected$new_formula, pattern = "H(\\d*)")[,2]),
    nrBr = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Br(\\d*)")[,2]),
    nrCl = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Cl(\\d*)")[,2]),
    nrF = as.numeric(stringr::str_match(corrected$new_formula, pattern = "F(\\d*)")[,2]),
    nrN = as.numeric(stringr::str_match(corrected$new_formula, pattern = "N(\\d*)")[,2]),
    nrO = as.numeric(stringr::str_match(corrected$new_formula, pattern = "O(\\d*)")[,2]),
    nrP = as.numeric(stringr::str_match(corrected$new_formula, pattern = "P(\\d*)")[,2]),
    nrS = as.numeric(stringr::str_match(corrected$new_formula, pattern = "S(\\d*)")[,2]),
    nrSi = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Si(\\d*)")[,2]),
    nrNa = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Na(\\d*)")[,2])
  )

  electron.per.atom <- data.table::data.table(
    nrC = 12,
    nrH = 1,
    nrBr = 79,
    nrCl = 35,
    nrF = 19,
    nrN = 14,
    nrO = 16,
    nrP = 31,
    nrS = 32,
    nrSi = 28,
    nrNa = 11
  )

  deconstructed[is.na(deconstructed)] <- 0
  deconstructed$nrAtoms <- rowSums(deconstructed)

  deconstructed$senior <- sapply(1:nrow(deconstructed), function(i){
    row = deconstructed[i,]
    with(row,{
      if((4*nrC+1*nrH+1*nrBr+1*nrCl+1*nrF+5*nrN+2*nrO+5*nrP+6*nrS+4*nrSi)>=(2*(nrAtoms-1))) TRUE else FALSE
    })
  })

  deconstructed$eminus <- sapply(1:nrow(deconstructed), function(i){
    row = deconstructed[i,]
    with(row,{
      # e- : =4*C4+D4+7*E4+7*F4+7*G4+5*H4+6*I4+5*J4+6*K4+4*L4
      4*nrC+nrH+7*nrBr+7*nrCl+7*nrF+5*nrN+6*nrO+5*nrP+6*nrS+4*nrSi
    })
  })

  deconstructed$lewis <- sapply(1:nrow(deconstructed), function(i){
    row = deconstructed[i,]
    with(row,{
      # 4*C4+ 1*D4 +1*E4 +1*F4 +1*G4 +3*H4 +2*I4 +3*J4 +2*K4 +4*L4
      lewis.sum = 4* nrC + 1* nrH +1* nrBr +1* nrCl +1* nrF +3* nrN +2* nrO +3* nrP +2* nrS +4* nrSi
      if(lewis.sum %% 2 == 0 & eminus >7) TRUE else FALSE
    })
  })

  ch_nops_chnops_rows <- lapply(1:nrow(deconstructed), function(i){
    row = deconstructed[i,]
    res= data.table::data.table(hc = F, chnops = F, nops = F)
    # T4, V4, X4, X4, Y4 -> H/C N/C	O/C	P/C	S/C
    # chnops
    #IF(AND(T4>=0.2,T4<=3,V4>=0,V4<=2,W4>=0,W4<=1.2,X4>=0,X4<=0.32,Y4>=0,Y4<=0.65),"YES","NO")
    # nops
    # IF(AND(V4>=0,V4<=4,W4>=0,W4<=3,X4>=0,X4<=2,Y4>=0,Y4<=3),"YES","NO"
    require(data.table)
    try({
      res = with(row,{
        HC = nrH/nrC
        NC = nrN/nrC
        OC = nrO/nrC
        PC = nrP/nrC
        SC = nrS/nrC
        hc = HC %between% c(0, 6)
        chnops = HC %between% c(0.2, 3) & NC %between% c(0, 2) & OC %between% c(0, 1.2) & PC %between% c(0, 0.32) & SC %between% c(0,65)
        nops = NC %between% c(0, 4) & OC %between% c(0,3) & PC %between% c(0,2) & SC %between% c(0,3)
        data.table::data.table(hc = hc, chnops = chnops, nops = nops)
      })
    })
    res
  })

  checks = data.table::rbindlist(ch_nops_chnops_rows)
  deconstructed <- cbind(deconstructed, checks)

  if(length(rules)>0){
    passes.checks <- rep(TRUE, nrow(deconstructed))
    for(rule in rules){
      passes.checks <- passes.checks & deconstructed[[rule]]
    }
    keep.candidates <- predicted$formula[passes.checks]
  }
  keep.candidates
}

revertAdduct <- function(formula, add_name, adduct_table=adducts){

  # C4H19N4S1Na1         C1H8N1.5S0.5    [2M+ACN+Na]1+  100

  row = adduct_table[Name == add_name]

  if(!is.na(row$AddEx) & row$AddEx != ""){
    if(!as.logical(enviPat::check_ded(formula, row$AddEx))){
      formula <- enviPat::subform(formula, row$AddEx)
    }else{
      return(NA)
    }
  }

  if(!is.na(row$RemEx) & row$RemEx != ""){
    formula <- enviPat::mergeform(formula, row$RemEx)
  }

  # remove multiplic
  if(row$xM > 1 & row$xM != ""){
    demultiplied = enviPat::multiform(formula, 1/row$xM)
    if(!grepl(demultiplied, pattern="\\.")){
      formula <- demultiplied
    }
  }

  # remove initial adduct
  if(!is.na(row$AddAt) & row$AddAt != ""){
    if(!as.logical(enviPat::check_ded(formula, row$AddAt))){
      formula <- enviPat::subform(formula, row$AddAt)
    }else{
      return(NA)
    }
  }
  if(!is.na(row$RemAt) & row$RemAt != ""){
    formula <- enviPat::mergeform(formula, row$RemAt)
  }
  formula
}

getPredicted <- function(mz,
                          ppm = 2,
                          mode = "positive",
                          rules = c("senior", "lewis", "hc", "chnops", "nops"),
                          elements = c("C","H","N","O","P","S"),
                          search= c("PubChem", "ChemSpider"),
                          detailed = T,
                          calc_adducts = adducts[Ion_mode == mode,]$Name){

  cat("
      _...._
    .`      `.    *             *
   / ***      \\            Predicting   *
  : **         :    *   molecular formulas...
  :            :
   \\          /
*   `-.,,,,.-'        *
     _(    )_             *
  * )        (                  *
   (          ) *
    `-......-`
")
  # find which mode we are in

  print("Considered adducts:")
  print(calc_adducts)

  per_adduct_results <- lapply(calc_adducts, function(add_name){

    row = adducts[Name == add_name]
    predicted = getFormula(mz = mz,
                           ppm = ppm,
                           add_name = add_name,
                           adducts = adducts,
                           elements = elements)
    if(is.null(predicted)) return(data.table::data.table())
    keep.candidates = filterFormula(predicted, rules = rules)
    if(length(keep.candidates) == 0) return(data.table::data.table())

    res = lapply(keep.candidates, function(formula, row){

        checked <- enviPat::check_chemform(isotopes, formula)
        new_formula <- checked[1,]$new_formula
        # check which adducts are possible
        theor_orig_formula = new_formula

        # remove last adduct
        theor_orig_formula = revertAdduct(theor_orig_formula, add_name)
        if(is.na(theor_orig_formula)){
          return(data.table())
        }else{
          data.table::data.table(query_mz = c(mz),
                                 name = new_formula,
                                 baseformula = theor_orig_formula,
                                 fullformula = new_formula,
                                 finalcharge = c(adducts[Name == add_name]$Charge),
                                 adduct = row$Name,
                                 `%iso` = 100,
                                 structure = NA,
                                 identifier = "???",
                                 description = "Predicted possible formula for this m/z value.",
                                 source = "magicball")
        }
      }, row = row)

      if(length(res) > 0){

        res_proc = MetaboShiny::flattenlist(res)

        tbl <- data.table::rbindlist(res_proc[!sapply(res_proc, is.null)])
    }

    uniques <- unique(tbl$baseformula)

    if(is.null(uniques)) return(data.table::data.table())

    unique(tbl)
  })
  total_tbl <- data.table::rbindlist(per_adduct_results[sapply(per_adduct_results, function(x)nrow(x)>0)], fill=T)
  total_tbl
}

searchFormulaWeb <- function(formulas,
                             search=c("PubChem", "ChemSpider"),
                             apikey="sp1pysTkYyC0wSETdkWjEeEK8eiXXFuG",
                             detailed=F){
  if(length(search)>0){
    i = 0
    count = length(formulas)
    if(count == 0) return(NULL)
    formulas <- unique(gsub(formulas, pattern = "\\[2\\]H", replacement = "D"))
    list_per_website <- lapply(search, function(db){
      switch(db,
             knapsack = {
               options(stringsAsFactors = F)
               rows = pbapply::pblapply(formulas, function(formula){
                 url = gsubfn::fn$paste("http://www.knapsackfamily.com/knapsack_jsp/result.jsp?sname=formula&word=$formula")
                 description = "No hits for this predicted formula."
                 rows = data.table::data.table(name = formula,
                                               baseformula = formula,
                                               structure = NA,
                                               description = description)
                 res = XML::readHTMLTable(url,header = T,)[[1]]
                 res.firstrow = colnames(res)
                 colnames(res) <- c("c_id", "cas_id", "metabolite", "formula", "mw", "organism")
                 res <- data.table::as.data.table(rbind(res.firstrow, res))
                 res.aggr = unique(res[, .(identifier = c_id, compoundname = metabolite, baseformula = formula,
                                    description = paste0("Found in organisms: ", paste0(organism, collapse = ", "))),by=c_id])
                 res.aggr = res.aggr[,-"c_id"]
                 uniques = unique(res.aggr$c_id)

                 if(detailed){
                   rows.detailed = pbapply::pblapply(uniques[!is.na(uniques)], function(id){
                     url = gsubfn::fn$paste("http://www.knapsackfamily.com/knapsack_jsp/information.jsp?word=$id")
                     tbl = XML::readHTMLTable(url,header = T,)[[1]]
                     flipped = t(tbl)
                     colnames(flipped) <- flipped[1,]
                     flipped = as.data.frame(flipped[-1,])
                     data.table::data.table(identifier = id,
                                            structure = flipped$SMILES[1])
                   })
                   res = merge(res.aggr, data.table::rbindlist(rows.detailed), by="identifier")
                 }else{
                   res = res.aggr
                 }
                 })
               data.table::rbindlist(rows)
             },
             pubchem = {
               print("Searching PubChem...")
               rows = pbapply::pblapply(formulas, function(formula){
                 url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/", formula, "/cids/JSON")
                 description = "No hits for this predicted formula."
                 rows = data.table::data.table(name = formula,
                                               baseformula = formula,
                                               structure = NA,
                                               description = description)

                 try({
                   pc_res <- jsonlite::read_json(url,simplifyVector = T)
                   ids <- pc_res$IdentifierList$CID
                   rows$description <- paste0("PubChem found these IDs: ",
                                                paste0(ids, collapse = ", "))
                   if(detailed){
                     rows = pubChemInfo(ids)
                   }
                 })
                 rows
               })
               data.table::rbindlist(rows, fill=T)
             },
             chemspider = {
               print("Searching ChemSpider...")

               cs_url = "https://api.rsc.org/compounds/v1/filter/formula/batch"
               formjson <- paste0('"', paste0(formulas, collapse='","'), '"')

               post_body = gsubfn::fn$paste('{
                                              "formulas": [
                                                  $formjson
                                              ],
                                              "orderBy": "",
                                              "orderDirection": ""
                                          }')
               r <- httr::POST(url = cs_url,
                               httr::add_headers('apikey' = apikey),
                               body = post_body,
                               httr::content_type("application/json"), encode="raw")
               qid <- RJSONIO::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
               done_url = gsubfn::fn$paste("https://api.rsc.org/compounds/v1/filter/formula/batch/$qid/status")
               res_url = gsubfn::fn$paste("https://api.rsc.org/compounds/v1/filter/formula/batch/$qid/results")
               done_r = httr::GET(done_url,httr::add_headers('apikey' = apikey))
               if(done_r$status_code == 200){
                 res_r = httr::GET(res_url,httr::add_headers('apikey' = apikey))
                 results = RJSONIO::fromJSON(httr::content(res_r, "text", encoding = "UTF-8"))
                 res_rows <- lapply(results[[1]], function(l){
                   if(length(l$results) == 0){
                     description = "No hits for this predicted formula."
                     data.table::data.table(name = l$formula,
                                            baseformula = l$formula,
                                            structure = NA,
                                            description = description)
                   }else{
                     pastedIds = paste0(l$results, collapse = ", ")
                     description <- paste0("ChemSpider found these IDs: ",
                                           pastedIds)
                     if(detailed){
                       res = chemspiderInfo(ids,
                                      apikey=apikey)
                       res$baseformula <- c(l$formula)
                       res
                     }else{
                       data.table::data.table(name = l$formula,
                                              baseformula = l$formula,
                                              structure = NA,
                                              `%iso` = c(100),
                                              description = description)
                     }
                   }
                 })
                data.table::rbindlist(res_rows,fill=T)
               }
             })
    })
    fin = data.table::rbindlist(list_per_website,fill=T)
    fin$identifier <- as.character(fin$identifier)
    fin$structure <- sapply(fin$structure, function(smi) if(is.na(smi)) "" else smi)
    fin
  }else{
    print("Please select at least one database to search in!")
    data.table::data.table()
  }
}

chemspiderInfo <- function(ids,
                           maxn=100,
                           apikey="sp1pysTkYyC0wSETdkWjEeEK8eiXXFuG"){
  split.ids = split(ids,
                    ceiling(seq_along(ids) / maxn))

  row.blocks <- lapply(split.ids, function(ids){
    pastedIds = paste0(ids, collapse = ", ")
    post_body = gsubfn::fn$paste('{
                                    "recordIds": [
                                        $pastedIds
                                    ],
                                    "fields": [
                                        "CommonName", "SMILES","PubMedCount","ReferenceCount","RSCCount"
                                    ]
                                }')
    cs_url = "https://api.rsc.org/compounds/v1/records/batch"
    r <- httr::POST(url = cs_url,
                    httr::add_headers('apikey' = apikey),
                    body = post_body)
    json = rawToChar(r$content)
    parsed = jsonlite::parse_json(json)[[1]]
    json_table = rbindlist(parsed)
    data.table::data.table(name = json_table$commonName,
                           structure = json_table$smiles,
                           identifier = json_table$id,
                           source = c("chemspider"),
                           `%iso` = c(100),
                           description = paste0("ChemSpider(", json_table$id, "). ",
                                                "Referenced ", json_table$referenceCount, " times. ",
                                                "Mentioned in PubMed ", json_table$pubMedCount, " times. ",
                                                "RSC count: ", json_table$rscCount, "."))
  })
  rbindlist(row.blocks)
}

pubChemInfo <- function(ids, maxn=30){
  # structural info
  split.ids = split(ids,
                     ceiling(seq_along(ids) / maxn))

  chunk.row.list <- lapply(split.ids, function(idgroup){

    url_struct = paste0("https://pubchem.ncbi.outnlm.nih.gov/rest/pug/compound/cid/",
                        paste0(idgroup, collapse=","),
                        "/property/MolecularFormula,CanonicalSMILES,Charge/JSON")

    struct_res <- jsonlite::fromJSON(url_struct,
                                     simplifyVector = T)

    # CHECK IF ORIGINAL CHARGE IS ZERO

    keep.ids <- which(struct_res$PropertyTable$Properties$Charge == 0)

    if(length(keep.ids)==0) return(NULL)

    idgroup = idgroup[keep.ids]

    rows <- struct_res$PropertyTable$Properties[keep.ids,]

    # descriptions
    url_desc = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                      paste0(idgroup, collapse=","),
                      "/description/JSON")

    try({

      desc_res <- jsonlite::fromJSON(url_desc,
                                     simplifyVector = T)

      descs <- as.data.table(desc_res$InformationList$Information)

      if("Description" %in% colnames(descs)){
        descs.adj <- descs[, list(name = Title[!is.na(Title)], Description = paste(Description[!is.na(Description)], collapse=" ")), by = CID]
      }else{
        descs.adj <- descs[, list(name = Title[!is.na(Title)], Description = c("No further description available")), by = CID]
      }

      descs.adj$Description <- gsub(descs.adj$Description,
                                    pattern = "</?a(|\\s+[^>]+)>",
                                    replacement = "",
                                    perl = T)

      if(any(descs.adj$Description == "")){
        descs.adj[Description == ""]$Description <- c("No further description available")
      }

      rows <- unique(merge(rows, descs.adj, by.x="CID", by.y="CID"))

    })

    if(is.null(rows)) return(NULL)
    if(nrow(rows) == 0) return(NULL)

    colnames(rows) <- c("identifier", "baseformula", "structure", "name","description")

    url_syn = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(idgroup, collapse=","),
                     "/synonyms/JSON")

    syn.adj = data.table()

    try({
      synonyms <- jsonlite::fromJSON(url_syn,
                                     simplifyVector = T)
      syn.adj = synonyms$InformationList$Information
    })

    if(nrow(syn.adj) > 0){
      rows.adj <- merge(rows, syn.adj, by.x="identifier", by.y="CID")
      rows.renamed <- lapply(1:nrow(rows.adj), function(i){
        row = rows.adj[i,]
        synonyms = row$Synonym[[1]]
        old.name <- row$name
        new.name <- synonyms[1]

        if(is.null(new.name)) new.name <- old.name

        desc.names <- synonyms[-1]

        row$name <- new.name

        row$description <- paste0(paste0("PubChem(", row$identifier, "). ",
                                         "Other names: ",
                                         paste0(if(length(desc.names) > 0) c(old.name, desc.names) else old.name, collapse="; "),
                                         ". "),
                                  row$description)
        row <- as.data.table(row)
        row[,-"Synonym"]

      })
      tbl.renamed <- data.table::rbindlist(rows.renamed, fill=T)
    }else{
      tbl.renamed <- rows
    }

    tbl.fin <- as.data.table(tbl.renamed)

    tbl.fin$source <- c("pubchem")

    # - - - return rows - - -

    result <- tbl.fin[,c("name", "baseformula", "structure", "description", "source")]
    result[,"%iso"] <- c(100)
    result

  })

  res <- chunk.row.list[sapply(chunk.row.list, function(x){
    returnme = F
    try({if(nrow(x)>0){
        returnme=T}})
    returnme
    })]
  data.table::rbindlist(res)
}

