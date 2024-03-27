#' Calculation of generalized additive mixed models (GAMMs)
#' @description Fits multiple generalized additive mixed models (GAMMs) to a longitudinal dataset. 
#' @param object an S4 object of class metime_analyser                      
#' @param which_data  a character defining the name of the dataset to be used.
#' @param verbose a logical on whether to print the calculation progress. Default set to FALSE.
#' @param cols_for_meta a character vector to define column names that are to be used for plotting purposes. Default set to NULL, therby not adding columns as metadata.
#' If you want automated facet wrapping option then set your new_columns as "facet_your_name"
#' @param name a character vector to define the index within the results. Should be equal to length of which_data. Default set to regression_gamm_1.
#' @param stratifications list to stratify data into a subset. Usage list(name=value). Default set to NULL, thereby not performing any type of stratification.
##' @param random a character vector defining which variables should be treated as random effects. Default set to "subject".
#' @param threshold a character of length 1 to define the type of threshold for significant interactions. 
#'      allowed inputs are "li", "FDR", "bonferroni" and "nominal"(cutoff p=0.05, set as Default)
##' @param interaction a character vector defining which interaction terms should be added to the model. Default set to NULL, with no interaction added.
#' @param num_cores numeric input to define the number of cores that you want to use for parallel computing. Default is set to NULL which is parallel::detectCores() -1.
#' @param k numeric input for setting the basis complexity for smoothing in gam. See more on [mgcv::bam]. Default is set to 10 as suggested by mgcv
#' @details The calculation function fits multiple generalized additive mixed models (GAMMs) on a longitudinal dataset. Here, one model fits one metabolite vs one trait. The degree of smoothness of a model term is estimated as part of the fitting. 
#' @return a S4 object of the class metime_analyzer with analysis results appended to the result section.
#' @export
setGeneric("calc_gamm", function(object, which_data, stratifications = NULL, cols_for_meta=NULL, 
  #random="subject", 
  threshold=c("none","nominal","li","fdr","bonferroni"), 
  #interaction = NULL,
  verbose=T,
  name="regression_gamm_1", 
  num_cores=NULL,
  k=10) standardGeneric("calc_gamm"))
setMethod("calc_gamm", "metime_analyser", function(object,
                                                   which_data,
                                                   stratifications = NULL,
                                                   cols_for_meta=NULL,
                                                   #random="subject",
                                                   threshold=c("none","nominal","li","fdr","bonferroni"),
                                                   #interaction=NULL,
                                                   verbose=T,
                                                   name="regression_gamm_1",
                                                   num_cores=NULL,
                                                   k=10) {
  # sanity checks ----
  ## check that covariates (cov) and type are set in the col_data
  if(!all(c("cov","type","interaction","random") %in% names(object@list_of_col_data[[which_data]]))) stop("calc_gamm() needs the columns cov and type to specify the model")
  
  while(grep(name, names(object@results)) %>% length() >= 1) {
    warning("name of the results was previously used, using a different name")
    index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    index <- c(0:9)[grep(index, 0:9)+1]
    name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  }
  
  # Data processing ----
  ## stratify data
  gamm_data <- get_stratified_data(which_data=which_data, 
                                 object=object, stratifications=stratifications)

  
    ## find list of metabolites
    my_met <- object@list_of_col_data[[which_data]] %>% 
      dplyr::filter(type=="met") %>% 
      dplyr::pull(id)
    
    ## find list of traits
    my_trait <- object@list_of_col_data[[which_data]] %>% 
      dplyr::filter(type=="trait") %>% 
      dplyr::pull(id)
    
    ## make list of formulas
    my_formula <- lapply(unique(my_trait),function(x) {
    object@list_of_col_data[[which_data]] %>% 
      dplyr::select(cov, type, id, random, interaction) %>% 
      dplyr::filter(id %in% my_met) %>% 
      dplyr::rename("met"="id") %>% 
      dplyr::mutate(trait=x,
                    cov = paste0(ifelse(is.na(cov), "", cov),
                                 object@list_of_col_data[[which_data]]$cov[which(object@list_of_col_data[[which_data]]$id==x)]),
                      random=object@list_of_col_data[[which_data]]$random[object@list_of_col_data[[which_data]]$id==x],
                      interaction=object@list_of_col_data[[which_data]]$interaction[object@list_of_col_data[[which_data]]$id==x])
      }) %>% 
    do.call(what=rbind.data.frame)

  # model calculation ----
  # changing mclapply to parLapply
  if(is.null(num_cores)) {
    cl <- parallel::makeCluster(spec = parallel::detectCores(all.tests = FALSE, logical = TRUE)-1, type="PSOCK")
  } else {
    cl <- parallel::makeCluster(spec = num_cores, type="PSOCK")
  }
  
  parallel::clusterExport(cl=cl, 
                          varlist=c("my_formula", "gamm_data", "verbose", "k"), # changed lmm_data to gamm_data
                          envir = environment())
  opb <- pbapply::pboptions(title="Running calc_gamm(): ", type="timer")
  on.exit(pbapply::pboptions(opb))
  results=pbapply::pblapply(cl=cl, 
                              1:nrow(my_formula),
                             #mc.cores=parallel::detectCores(all.tests = FALSE, logical = TRUE)-1,
                             #mc.preschedule = TRUE,
                             function(x) {
    if(verbose) cat(x, " , ") # report iteration
                               require(magrittr)
    # extract data 
    this_data <-  gamm_data$data %>% 
      dplyr::select(any_of(setdiff(names(gamm_data$data), c("subject","time")))) %>% 
      dplyr::mutate(id = rownames(.[])) %>%
      dplyr::left_join(by="id", gamm_data$row_data[,c("id","time","subject")]) %>% 
      dplyr::mutate(met=get(my_formula$met[x]), 
                    trait=get(my_formula$trait[x])) %>% 
      #dplyr::filter(!is.na(met), !is.na(trait)) %>% 
      #dplyr::drop_na(met) %>% dplyr::drop_na(trait) %>% 
      dplyr::mutate(time = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=time)),
                    subject = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=subject)))
    
    # extract formula information
      formula_met <- my_formula$met[x]

      formula_trait <- ifelse(is.na(my_formula$interaction[x]), 
                                                  "trait", 
                                                  paste0("trait", "*",
                                                     paste0(my_formula$interaction[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character() %>% .[!. %in% ""], collapse="*")))
      
      formula_cov <- my_formula$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character() %>% .[!. %in% c("", "NA", NA)]
       
      smoothing_cov <- sapply(formula_cov,function(y) ifelse(length(unique(this_data[[y]]))>3, paste0("s(",y, ", k=", k, ")"), y)) 
      
      formula_random <-  ifelse(!is.na(my_formula$random[x]), paste0("s(",my_formula$random[x], ", bs='re')"), NA)
      
      #formula_interaction <- ifelse(!is.null(interaction), paste0("s(",interaction, ")"), NA)
      
      
      # collapse formula
      formula = paste0("met ~ ",
                      ifelse(all(!is.na(formula_trait)), formula_trait, ""), 
                       ifelse(all(!is.na(smoothing_cov)), paste0("+",paste0(smoothing_cov,collapse = "+")), ""),
                       #ifelse(all(!is.na(formula_interaction)), paste0("+",paste0(formula_interaction,collapse = "+")), ""),
                       ifelse(all(!is.na(formula_random)), paste0("+",paste0(formula_random,collapse = "+")), "")
      )
    this_model = try(
      silent=T,
      mgcv::bam(
        formula=as.formula(formula), 
        data= this_data,
        method="REML"
      )
    )
    if(all(class(this_model)!="try-error")) {
      results_sum <- summary(this_model)
      out_this_model <- data.frame(
        stringsAsFactors = F,
        met = my_formula$met[x],
        trait = my_formula$trait[x],
        level = grep(x=names(results_sum$p.coeff), pattern="trait", value=T) %>% gsub(pattern="trait",replacement=my_formula$trait[x]),
        beta = results_sum$p.coeff[grep(x=names(results_sum$p.coeff), pattern="trait", value=T)] %>% as.numeric(),
        se = results_sum$se[grep(x=names(results_sum$se), pattern="trait", value=T)]%>% as.numeric(),
        pval = results_sum$p.pv[grep(x=names(results_sum$p.pv), pattern="trait", value=T)]%>% as.numeric(),
        tval = results_sum$p.t[grep(x=names(results_sum$p.t), pattern="trait", value=T)]%>% as.numeric(),
        formula=formula
      ) #%>% 
        #cbind(as.data.frame(t(results_sum$s.table[,1])))
    } else {
      out_this_model <- data.frame(
        stringsAsFactors = F,
        met = my_formula$met[x],
        trait = my_formula$trait[x],
        level=NA,
        beta= NA,
        se=NA,
        pval=NA,
        tval=NA,
        formula = formula
      )
    }
    return(out_this_model)
  })

  on.exit(parallel::stopCluster(cl))

  annotated_results <- plyr::rbind.fill(results)
  
  ## modify results to for pipeline
  annotated_results <- annotated_results %>% 
    dplyr::mutate(x=beta, 
                  xmin=beta-abs(se), 
                  xmax=beta+abs(se), 
                  y=met, 
                  color= "none",
                  type=ifelse(!is.na(level), level, trait))
  #rownames(annotated_results) <- annotated_results$y 
  

   # calculate the thresholds 
  thresh_bonferroni <- 0.05/length(my_met)
  eigenvals <- cor(object@list_of_data[[which_data]][,my_met], use="pairwise.complete.obs") %>%
    eigen()
  thresh_li <- 0.05/(sum(as.numeric(eigenvals$values >= 1) + (eigenvals$values - floor(eigenvals$values))))
  annotated_results <- annotated_results %>% dplyr::mutate(li_thresh=thresh_li, bonferroni_thresh=thresh_bonferroni)

  
  # for(i in intersect(c("none","nominal","li","fdr","bonferroni"),threshold)){
  #   if(i=="none"){
  #     annotated_results<-annotated_results %>% 
  #       dplyr::mutate(color="none") 
  #   }else if(i=="nominal"){
  #     annotated_results<-annotated_results %>% 
  #       dplyr::mutate(color=ifelse(pval<=0.05, "nominal",color))
  #   }else if(i == "li"){
  #     eigenvals <- cor(object@list_of_data[[which_data]][,my_met], use = "pairwise.complete.obs") %>%
  #       eigen()
  #     li_thresh <- 0.05/(sum(as.numeric(eigenvals$values >= 1) + (eigenvals$values - floor(eigenvals$values))))
  #     annotated_results<-annotated_results%>% 
  #       dplyr::mutate(color=ifelse(pval<=li_thresh, "li",color))
  #   }else if(i=="fdr"){
  #     annotated_results<-annotated_results %>% 
  #       dplyr::mutate(fdr=p.adjust(pval, method="BH")) %>% 
  #       dplyr::mutate(color=ifelse(fdr<=0.05, "fdr",color)) %>% 
  #       dplyr::select(-fdr)
  #   }else if(i=="bonferroni"){
  #     annotated_results=annotated_results %>% 
  #       dplyr::mutate(color=ifelse(pval<=0.05/length(my_met), "bonferroni",color))
  #   }
  # }
  # split into single results files ----
  out_results <- lapply(unique(annotated_results$type), function(y){
    annotated_results[which(annotated_results$type==y),] %>% 
      dplyr::mutate(qval = p.adjust(pval, method="BH")) %>% 
      dplyr::mutate(color = ifelse(pval <= 0.05, "nominally","none")) %>% 
      dplyr::mutate(color = ifelse(qval <= 0.05, "fdr",color)) %>% 
      dplyr::mutate(color = ifelse(pval <= thresh_li, "li",color)) %>% 
      dplyr::mutate(color = ifelse(pval <= thresh_bonferroni, "bonferroni",color)) %>% 
      `rownames<-`(.[,"met"])
  })
  
  # out_results <- lapply(seq_along(out_results), function(ind) {
  #         rownames(out_results[[ind]]) <- out_results[[ind]]$met
  #         return(out_results[[ind]]) 
  #   })
  names(out_results) <- unique(annotated_results$type)
  
  if(is.null(cols_for_meta)) {
    metadata <- NULL
  } else {
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta)
  }
     
  
  # continue here
  out <- get_make_results(object = object, 
                          data = out_results, 
                          metadata = metadata, 
                          calc_type = rep("regression",length(out_results)), 
                          calc_info = paste("Generalized additive mixed models for ", which_data,"_", names(out_results),sep=""),
                          name = name) %>%
    add_function_info(function_name = "calc_gamm", 
                      params = list(which_data = which_data, 
                                    verbose = verbose, 
                                    cols_for_meta = cols_for_meta, 
                                    #random=random,
                                    #interaction=interaction,
                                    num_cores=num_cores,
                                    threshold=threshold,
                                    name = name, 
                                    stratifications = stratifications))
  
  return(out)
  
})



