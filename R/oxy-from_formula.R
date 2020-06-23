#' @title Find possible formulas for a given m/z
#' @description Using m/z and isotope distributions for each element, find possible molecular formulas within allowed error margin
#' @param mz M/z of interest
#' @param add_name Which adducts to consider
#' @param adducts Full adduct table (data(adducts) loads it into memory)
#' @param ppm Allowed error margin in parts per million
#' @param elements Considered elements in formula generation, Default: c("C", "H", "N", "O", "P", "S")
#' @return Table with found formulas, their adduct and isotope percentage.
#' @seealso
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[Rdisop]{initializeCHNOPS}}
#' @rdname getFormula
#' @export
#' @importFrom stringr str_match_all
#' @importFrom Rdisop decomposeMass initializeElements
#' @importFrom data.table data.table
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
  ppms = abs((predicted$exactmass - as.numeric(mz)) / (1e-6 * as.numeric(mz)))
  data.table::data.table(
    fullformula = predicted$formula,
    dppm = ppms
  )
}

#' @title Apply seven golden rules to a vector of formulas
#' @description Returns formulas that pass the user-given rules.
#' @param formulas Molecular formulas
#' @param rules Rules to apply. Default: c("senior", "lewis", "hc", "chnops", "nops")
#' @return Vector of formulas that pass rules
#' @seealso
#'  \code{\link[enviPat]{check_chemform}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[stringr]{str_match}}
#' @rdname filterFormula
#' @export
#' @importFrom enviPat check_chemform
#' @importFrom data.table data.table rbindlist
#' @importFrom stringr str_match
filterFormula = function(formulas, rules=c("senior", "lewis", "hc", "chnops", "nops")){
  data(isotopes, package = "enviPat", envir = environment())

  corrected = enviPat::check_chemform(isotopes = isotopes, chemforms = formulas)

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

  if(length(rules) > 0){
    passes.checks <- rep(TRUE, nrow(deconstructed))
    for(rule in rules){
      passes.checks <- passes.checks & deconstructed[[rule]]
    }
    keep.candidates <- formulas[passes.checks]
  }else{
    keep.candidates <- formulas
  }
  keep.candidates
}

#' @title Break formula apart into adduct and main formula
#' @description Used to use the formula creation function and consider adducts at the same time.
#' @param formula Formula of interest
#' @param add_name Adduct names to consider ('Name' column of adduct table)
#' @param adduct_table Full adduct table, Default: adducts
#' @return Table with formula, adduct, isotope
#' @seealso
#'  \code{\link[enviPat]{check_ded}},\code{\link[enviPat]{subform}},\code{\link[enviPat]{mergeform}},\code{\link[enviPat]{multiform}}
#' @rdname revertAdduct
#' @export
#' @importFrom enviPat check_ded subform mergeform multiform
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

#' @title Get predicted formulas and adducts from m/z value
#' @description Wrapper function to predict formulas and then consider adducts as well.
#' @param mz M/z of interest
#' @param ppm Error margin in parts per million, Default: 2
#' @param mode M/z found in positive or negative mode?, Default: 'positive'
#' @param rules Which golden rules to apply?, Default: c("senior", "lewis", "hc", "chnops", "nops")
#' @param elements Which elements to consider?, Default: c("C", "H", "N", "O", "P", "S")
#' @param search Check the found formulas on PubChem or ChemSpider?, Default: c("PubChem", "ChemSpider")
#' @param detailed Look up details like description etc. if hit found? Makes things slower!, Default: T
#' @param calc_adducts Which adducts to consider?, Default: adducts[Ion_mode == mode, ]$Name
#' @return Table of found matches and associated info
#' @seealso
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[enviPat]{check_chemform}}
#' @rdname getPredicted
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom enviPat check_chemform
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
    if(nrow(predicted) == 0) return(data.table::data.table())
    keep.candidates = filterFormula(predicted$fullformula, rules = rules)
    if(length(keep.candidates) == 0) return(data.table::data.table())

    iter.rows = predicted[fullformula %in% keep.candidates,]

    res = lapply(1:nrow(iter.rows), function(i, row){

      formula = iter.rows[i,]$fullformula
      checked <- enviPat::check_chemform(isotopes, formula)
      new_formula <- checked[1,]$new_formula
      # check which adducts are possible
      theor_orig_formula = new_formula

      # remove last adduct
      theor_orig_formula = revertAdduct(theor_orig_formula,
                                        add_name,
                                        adduct_table = adducts)
      if(is.na(theor_orig_formula)){
        return(data.table())
      }else{
        data.table::data.table(query_mz = c(mz),
                               name = theor_orig_formula,
                               baseformula = theor_orig_formula,
                               fullformula = new_formula,
                               basecharge = c(0),
                               finalcharge = c(adducts[Name == add_name]$Charge),
                               adduct = row$Name,
                               `%iso` = 100,
                               identifier = new_formula,
                               structure = paste0("[", new_formula, "]0"),
                               description = "Predicted possible formula for this m/z value.",
                               source = "magicball",
                               dppm = iter.rows[i,]$dppm)
      }
    }, row = row)

    if(length(res) > 0){

      flattenlist = function(x){
        morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
        out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE, use.names = T))
        if(sum(morelists)){
          Recall(out)
        }else{
          return(out)
        }
      }

      res_proc = flattenlist(res)

      tbl <- data.table::rbindlist(res_proc[!sapply(res_proc, is.null)])
    }

    uniques <- unique(tbl$baseformula)

    if(is.null(uniques)) return(data.table::data.table())

    unique(tbl)
  })
  total_tbl <- data.table::rbindlist(per_adduct_results[sapply(per_adduct_results, function(x)nrow(x) > 0)], fill=T)
  total_tbl
 }

#' @title Find web hits for a molecular formula
#' @description Takes molecular formula, and scours PubChem and/or ChemSpider for compounds matching that formula.
#' @param formulas Character vector of formulas to check
#' @param search Which databases to check?, Default: c("PubChem", "ChemSpider")
#' @param apikey API key for ChemSpider
#' @param detailed Find detailed results? Not just the compound name, but other associated info?, Default: F
#' @return Data table with match results.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hits = searchFormulaWeb(c("C6H12O6"), "PubChem", detailed = T)
#'  }
#' }
#' @seealso
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{as.data.table}},\code{\link[data.table]{rbindlist}}
#'  \code{\link[enviPat]{check_chemform}}
#'  \code{\link[gsubfn]{fn}}
#'  \code{\link[stringr]{str_extract}},\code{\link[stringr]{str_match}}
#'  \code{\link[XML]{readHTMLTable}}
#'  \code{\link[rlist]{list.clean}}
#'  \code{\link[jsonlite]{read_json}}
#'  \code{\link[httr]{POST}},\code{\link[httr]{add_headers}},\code{\link[httr]{content_type}},\code{\link[httr]{content}},\code{\link[httr]{GET}}
#'  \code{\link[RJSONIO]{fromJSON}}
#' @rdname searchFormulaWeb
#' @export
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table as.data.table rbindlist
#' @importFrom enviPat check_chemform
#' @importFrom gsubfn fn
#' @importFrom stringr str_extract str_match
#' @importFrom XML readHTMLTable
#' @importFrom rlist list.clean
#' @importFrom jsonlite read_json
#' @importFrom httr POST add_headers content_type content GET
#' @importFrom RJSONIO fromJSON
searchFormulaWeb <- function(formulas,
                             search = c("pubchem", "chemspider", "knapsack", "supernatural2"),
                             apikey = "sp1pysTkYyC0wSETdkWjEeEK8eiXXFuG",
                             detailed = T){

  if(length(search)>0){
    i = 0
    count = length(formulas)
    if(count == 0) return(NULL)
    formulas <- unique(gsub(formulas, pattern = "\\[2\\]H", replacement = "D"))
    list_per_website <- lapply(search, function(db){
      switch(db,
             supernatural2 = {
               print("Searching SUPER NATURAL II...")
               rows = pbapply::pblapply(formulas, function(formula){
                 rows = data.table::data.table(name = formula,
                                               baseformula = formula,
                                               structure = NA,
                                               description = "No hits for this predicted formula.")
                 # get molecular weight appoximately
                 molwt = enviPat::check_chemform(isotopes, formula)$monoisotopic_mass
                 # go .5da below and above
                 mzmin = molwt - 0.0001
                 mzmax = molwt + 0.0001
                 # generate url
                 theurl = gsubfn::fn$paste("http://bioinf-applied.charite.de/supernatural_new/index.php?site=compound_search&start=0&supplier=all&molwt1=$mzmin&molwt2=$mzmax&classification=all")
                 # fetch hits in range
                 html = readLines(theurl)
                 count = stringr::str_extract(html, "count_compounds=(\\d+)")
                 total = as.numeric(gsub(count[!is.na(count)], pattern = ".*=", replacement = ""))
                 if(length(total) > 0){
                   pages = total / 15
                   rows = lapply(0:ceiling(pages), function(i){
                     db.frag = data.table::data.table()
                     try({
                       theurl = gsub(theurl, pattern = "start=(\\d)", replacement = paste0("start=", i))
                       html = readLines(theurl)
                       # which formulas match?
                       formulas = stringr::str_extract(html, "(\\d+)")
                       identifiers = stringr::str_extract(html, "SN(\\d+)")
                       identifiers = unique(identifiers[!is.na(identifiers)])
                       smi_rows = grep(html, pattern = "SMILES") + 1
                       smiles = stringr::str_match(html[smi_rows], 'value="(.*)"\\/>')[,2]
                       identifiers = unique(identifiers[!is.na(identifiers)])
                       tables <- XML::readHTMLTable(theurl)
                       tables <- rlist::list.clean(tables, fun = is.null, recursive = FALSE)
                       rows = lapply(tables, function(tbl){
                         if("Name" %in% tbl$V1){
                           tbl = data.table::as.data.table(tbl)
                           data.table::data.table(
                             compoundname = gsub(tbl[V1 == "Name",V2], pattern = "[^\x01-\x7F]+.-", replacement = ""),
                             description = paste0("Toxicity class: ", tbl[V1 == "Tox-class",V2]),
                             baseformula = tbl[V1 == "Formula",V2],
                             identifier = "???",
                             charge = tbl[V1 == "Charge",V2],
                             structure = "",
                             source = "supernatural2"
                           )
                         }else{
                           data.table::data.table()
                         }
                       })
                       db.frag = data.table::rbindlist(rows)
                       db.frag$identifier = identifiers
                       db.frag$structure = smiles
                       keep = which(enviPat::check_chemform(isotopes, db.frag$baseformula)$new_formula == formula)
                       db.frag[keep,]
                     })
                     db.frag
                   })
                   rows = data.table::rbindlist(rows, fill=T)
                 }
                 rows
               })
               data.table::rbindlist(rows, fill=T)
             },
             knapsack = {
               print("Searching knapsack...")
               options(stringsAsFactors = F)
               rows = pbapply::pblapply(formulas, function(formula){
                 url = gsubfn::fn$paste("http://www.knapsackfamily.com/knapsack_core/result.php?sname=formula&word=$formula")
                 print(url)
                 description = "No hits for this predicted formula."
                 rows = data.table::data.table(name = formula,
                                               baseformula = formula,
                                               structure = NA,
                                               description = description,
                                               source = "magicball")
                 res = XML::readHTMLTable(url,header = T,)[[1]]
                 res.firstrow = colnames(res)
                 if(!is.null(res.firstrow)){
                   colnames(res) <- c("c_id", "cas_id", "metabolite",
                                      "formula", "mw", "organism")
                   res <- data.table::as.data.table(rbind(res.firstrow, res))
                   res.aggr = unique(res[, .(identifier = c_id, compoundname = metabolite, baseformula = formula, source = "knapsack",
                                             description = paste0("Found in organisms: ", paste0(organism, collapse = ", "))),by=c_id])
                   res.aggr = res.aggr[,-"c_id"]
                   uniques = unique(res.aggr$c_id)

                   if(detailed){
                     rows.detailed = pbapply::pblapply(uniques[!is.na(uniques)], function(id){
                       url = gsubfn::fn$paste("http://www.knapsackfamily.com/knapsack_core/information.php?word=$id")
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
                 }else{
                   rows
                 }
               })
               data.table::rbindlist(rows, fill=T)
             },
             pubchem = {
               print("Searching PubChem...")
               rows = pbapply::pblapply(formulas, function(formula){
                 url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/", formula, "/cids/JSON")
                 description = "No hits for this predicted formula."
                 rows = data.table::data.table(name = formula,
                                               baseformula = formula,
                                               structure = NA,
                                               description = description,
                                               source = "magicball")

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
               formblocks = split(formulas, ceiling(seq_along(formulas)/100))

               if(apikey == ""){
                 print("No ChemSpider API key supplied!")
                 return(data.table::data.table())
               }

               rows = pbapply::pblapply(formblocks, function(formgroup){
                 formjson <- paste0('"', paste0(formgroup,
                                                collapse='","'),
                                    '"')

                 post_body = gsubfn::fn$paste('{"formulas": [$formjson]}')
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
                                              description = description,
                                              source = "magicball")
                     }else{
                       pastedIds = paste0(l$results, collapse = ", ")
                       description <- paste0("ChemSpider found these IDs: ",
                                             pastedIds)
                       if(detailed){
                         print("...")
                         res = chemspiderInfo(l$results,
                                              apikey=apikey)
                         print("!")
                         res$baseformula <- c(l$formula)
                         res
                       }else{
                         data.table::data.table(name = l$formula,
                                                baseformula = l$formula,
                                                structure = NA,
                                                `%iso` = c(100),
                                                description = description,
                                                source = "chemspider")
                       }
                     }
                   })
                   data.table::rbindlist(res_rows, fill = T)
                 }
               })
               data.table::rbindlist(rows, fill = T)
             })
    })
    fin = data.table::rbindlist(list_per_website,fill=T)
    if(nrow(fin) > 0){
      fin$identifier <- as.character(fin$identifier)
      fin$structure <- sapply(fin$structure, function(smi) if(is.na(smi)) "" else smi)
    }
    unique(fin)
  }else{
    print("Please select at least one database to search in!")
    data.table::data.table()
  }
}

#' @title Find more info through ChemSpider.
#' @description Takes ChemSpider CIDs and finds name, SMILES, citations on pubmed/references.
#' @param ids Character vector of ChemSpider IDs.
#' @param maxn Max ids per batch (batch search is used), Default: 100
#' @param apikey ChemSpider API key
#' @return Data table with match results
#' @seealso
#'  \code{\link[gsubfn]{fn}}
#'  \code{\link[httr]{POST}},\code{\link[httr]{add_headers}}
#'  \code{\link[jsonlite]{read_json}}
#' @rdname chemspiderInfo
#' @export
#' @importFrom gsubfn fn
#' @importFrom httr POST add_headers
#' @importFrom jsonlite parse_json
#' @importFrom data.table data.table
chemspiderInfo <- function(ids,
                           maxn=100,
                           apikey){#="sp1pysTkYyC0wSETdkWjEeEK8eiXXFuG"){
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
    json_table = rbindlist(parsed, fill=T)
    data.table::data.table(name = json_table$commonName,
                           structure = json_table$smiles,
                           identifier = json_table$id,
                           `%iso` = c(100),
                           source = "chemspider",
                           description = paste0("ChemSpider(", json_table$id, "). ",
                                                "Referenced ", json_table$referenceCount, " times. ",
                                                "Mentioned in PubMed ", json_table$pubMedCount, " times. ",
                                                "RSC count: ", json_table$rscCount, "."))
  })
  rbindlist(row.blocks, fill=T)
}

#lmdb = lmdb[, lapply(.SD, function(x) textclean::replace_non_ascii(x, remove.nonconverted = T))]

#' @title Find additional info on a PubChem ID.
#' @description Takes PubChem ID and finds name, formula, smiles, charge
#' @param ids Vector of identifiers.
#' @param maxn Compounds searched per batch search, Default: 30
#' @return Table with additional info on PubChem IDs
#' @seealso
#'  \code{\link[data.table]{rbindlist}}
#' @rdname pubChemInfo
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
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
