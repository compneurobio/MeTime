# Case Study: Disease Mechanism Discovery

## Case Study: Disease Mechanism Discovery

> One-line summary: end-to-end pipeline showing how MeTime supported the
> core results of the preprint: “A seven year old longitudinal study of
> the Alzheimer’s Disease metabolome”
> (<https://pubmed.ncbi.nlm.nih.gov/41409665/>)

### Preparing data

Download the MXPQuant 500 dataset of ADNI and follow the preprocessing
steps as described in the manuscript.

``` r


library(MeTime)
setwd("/path/for/analysis")

# load biocrates data here
my_analyser <- readRDS("/read/processed_data/from/imputation/script") # This data is not provided

which_data <- "biocrates_data" # change name based on the dataset you are analysing - this is an example for MXP Quant 500

my_analyser <- my_analyser %>%
  add_screening_vars(vars=c("PTEDUCAT", "PTGENDER", "APOEGrp")) %>% # this function only works for factor variables and not for continuous like age, bmi etc
  add_distribution_vars_to_rows(which_data=which_data, screening_vars=NULL, distribution_vars=c("DXGrp_long", "time_resolved"))   

path <- "write/path/for/analyser/objects/" # make sure it ends with "/"


# you can remove this chunk if you have already adjusted for this 
my_analyser <- my_analyser %>% 
          mod_filter(which_data=which_data, type="row_data", 
            !time_resolved %in% c("mix_CN_stable", "mix_MCI_stable", "mix_AD_stable", "MCI_CN_revert")) %>% # goal is to remove all samples with atypical diagnostic trajectories
          mod_mutate(which_data="phenotype_data", type="data",
             DXGrp_long=ifelse(time_resolved %in% c("CN_AD_convert", "MCI_AD_convert"), "AD_converter", time_resolved)) 


# filter for timepoints here 
# this depends on how time was encoded in the original dataset - both 0 and t0 are allowed so change it accordingly 
my_analyser_long <- my_analyser %>%
  mod_filter(which_data=which_data, type="row_data", time %in% c("0", "6", "12", "18", "24", "36", "48", "60", "72", "84"))
my_analyser_main_tp <- my_analyser %>%
  mod_filter(which_data=which_data, type="row_data", time %in% c("0", "12", "24"))



# Getting covariates info table
xlsx_path <- "biocrates_covariates.xlsx" # change this to current file path (this file is available in the same folder)
# This file is also shared to show how to approach more complicated mixed models
covariate_data <- lapply(excel_sheets(xlsx_path), function(x) read_excel(path = xlsx_path, sheet = x) %>% 
                           dplyr::select(any_of(c("id", "cov", "type", "interaction", "random")))) %>% 
  magrittr::set_names(excel_sheets(xlsx_path))


lmm <- covariate_data$lmm 
meta <- covariate_data$meta 
lm_model <- covariate_data$lm 
diag_lm <- covariate_data$diag_lm 
diag_meta <- covariate_data$diag_meta 
diag_lmm <- covariate_data$diag_lmm 

timepoints <- c("t0", "t6", "t12", "t18", "t24", "t36", "t48", "t60", "t72", "t84")


convert_analyser <- function(object) {
    new_object <- object %>%
        mod_mutate(which_data="phenotype_data", 
             regression_dxgrp = as.factor(DXGrp_longi)) %>%
        mod_merge_data(which_data=c("biocrates_data", "phenotype_data","medication_data", "nmr_data"), # change biocrates_data if the dataset is different 
                 filter_sample="biocrates_data", name="regression", append=TRUE,
                 cols_list=list(biocrates_data=names(object@list_of_data$biocrates_data), # same as above
                                phenotype_data=names(object@list_of_data$phenotype_data),
                                medication_data=names(object@list_of_data$medication_data),
                                nmr_data=c("TOTAL_C", "HDL_C"))) %>% 
        mod_mutate(which_data="regression", 
             "regression_dxgrp_long" = factor(DXGrp_long, levels=c("CN", "CN_MCI_convert","MCI_stable", "AD_converter", "AD")),
             "regression_scale(abeta)" =  CSF_Roche_Abeta42 %>% as.numeric() %>% scale(),
             "regression_entorhinal_thickness" = (LEntCtx_Thick + REntCtx_Thick) %>% as.numeric() %>% as.numeric(),
             "regression_entorhinal_volume" = (LEntCtx_vol + REntCtx_vol) %>% as.numeric() %>% as.numeric(),
             "regression_log2(ptau)" = CSF_Roche_PTau %>% as.numeric() %>% log2() %>% as.numeric(),
             "regression_log2(tau)" = CSF_Roche_Tau %>% as.numeric() %>% log2() %>% as.numeric(),
             "regression_scale(log2(CSF_Roche_Tau/CSF_Roche_Abeta42))" = as.numeric(CSF_Roche_Tau)/as.numeric(CSF_Roche_Abeta42) %>% log2() %>% as.numeric(),
             "regression_scale(sqrt(adas13))" = ADAS13 %>% sqrt() %>% as.numeric(),
             "regression_adni_mem" = ADNI_MEM,
             "regression_adni_ef" = ADNI_EF,
             "regression_adni_lan" = ADNI_LAN,
             "regression_hipp_vol" = HippVol,
             "regression_grey_matter_thick" = GlobalCtx_Thick,
             "regression_grey_matter_vol" = Global_vol,
             "regression_fdgpet" = (FDG_Cing_Mean + FDG_LAng_Mean + FDG_LTemp_Mean + FDG_RTemp_Mean + FDG_RAng_Mean) %>% as.numeric(),
             "regression_dxgrp_time" = as.numeric(ifelse(DXGrp %in% c(1,5), 0, ifelse(DXGrp %in% c(2,3), 1, 2)))) 

    new_object@list_of_data$regression <- new_object@list_of_data$regression %>% 
                                                dplyr::rename("Total_C"="TOTAL_C") # depends on how Total_C was encoded in the excel sheet 

    return(new_object)
}

# creating regression dataset now
my_analyser_long <- convert_analyser(my_analyser_long)
my_analyser_main_tp <- convert_analyser(my_analyser_main_tp)

phenotypes <- c("regression_scale(abeta)",
                "regression_entorhinal_thickness",
                "regression_entorhinal_volume",
                "regression_log2(ptau)",
                "regression_log2(tau)",
                "regression_hipp_vol",
                "regression_grey_matter_thick",
                "regression_grey_matter_vol",
                "regression_fdgpet",
                "regression_adni_mem",
                "regression_adni_ef",
                "regression_adni_lan")

others <- c("regression_dxgrp_long", "regression_dxgrp_time")


metabolites <- my_analyser_long@list_of_data[["biocrates_data"]] %>% colnames()
metabolites <- c(metabolites, "Total_C", "HDL_C") # check TOTAL_C encoding here again


# Updating regression data in by scaling them based on the type of data
for(i in c(phenotypes, metabolites)) {
    print(i)
    my_analyser_long@list_of_data[["regression"]][ ,i] <- scale(my_analyser_long@list_of_data[["regression"]][ ,i])
    my_analyser_main_tp@list_of_data[["regression"]][ ,i] <- scale(my_analyser_main_tp@list_of_data[["regression"]][ ,i])
}


    

filter_by_nonNAs <- function(object) {

  data <- object@list_of_data$regression %>% 
            dplyr::mutate(timepoint_col=rownames(.) %>% gsub(pattern="R[0-9]+_", replacement=""))
  # Identify the columns that start with "regression_"
  regression_cols <- grep("^regression_", names(data), value = TRUE)
  
  # Iterate over each regression column
  for (col in regression_cols) {
    # Iterate over each unique timepoint
    for (timepoint in unique(data$timepoint_col)) {
      # Subset the rows for the current timepoint
      timepoint_rows <- data$timepoint_col== timepoint
      
      # Count non-NA values in the current column at the current timepoint
      nonNA_count <- sum(!is.na(data[timepoint_rows, col]))
      
      # If fewer than 50 non-NAs, set all values to NA for this column at this timepoint
      if (nonNA_count < 50) {
        data[timepoint_rows, col] <- NA
      }
    }
  }
  
  data$timepoint_col <- NULL

  object@list_of_data$regression <- data
  return(object)
}


# filtering for timepoints with >50 trait values for regression
my_analyser_long <- filter_by_nonNAs(my_analyser_long)



# multi-timepoint meta-analysis for biocrates is shown here

std_meta <- my_analyser_main_tp %>% add_data(which_data="regression", type="col_data",x=meta, id="id") %>%
    calc_lmm(which_data="regression", 
            name="biocrates_meta", 
            stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
            num_cores=2)
saveRDS(std_meta, paste0(path,"biocrates_meta.rds")) # change name as per need

# Time-interaction analysis

std_lmm <- my_analyser_long %>% add_data(which_data="regression", type="col_data",x=lmm, id="id") %>%
    calc_lmm(which_data="regression", 
            name="biocrates_time_interaction", 
            stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
            num_cores=12)
saveRDS(std_lmm, paste0(path,"biocrates_time_interaction.rds")) # change name as per need
 
# Time-interaction analysis subset to main timepoints
std_lmm_main_tp <- my_analyser_main_tp %>% add_data(which_data="regression", type="col_data",x=lmm, id="id") %>%
    calc_lmm(which_data="regression", 
            name="biocrates_time_interaction_main_tp", 
            stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
            num_cores=12)
saveRDS(std_lmm_main_tp, paste0(path,"biocrates_time_interaction_main_tp.rds")) # change name as per need

# Cross-sectional single time point models

std_lm0 <- my_analyser_main_tp %>% add_data(which_data="regression", type="col_data",x=lm_model, id="id") %>%
    calc_lm(which_data="regression", 
             name="biocrates_lm0", 
             stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
             num_cores=2, timepoint="0")
saveRDS(std_lm0, paste0(path,"biocrates_lm0.rds")) # change name as per need

std_lm12 <- my_analyser_main_tp %>% add_data(which_data="regression", type="col_data",x=lm_model, id="id") %>%
    calc_lm(which_data="regression", 
             name="biocrates_lm12", 
             stratifications=NULL, 
             cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
             num_cores=2, timepoint="12")
saveRDS(std_lm12, paste0(path,"biocrates_lm12.rds")) # change name as per need


std_lm24 <- my_analyser_main_tp %>% add_data(which_data="regression", type="col_data",x=lm_model, id="id") %>%
    calc_lm(which_data="regression", 
             name="biocrates_lm24", 
             stratifications=NULL, 
             cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
             num_cores=2, timepoint="24")
saveRDS(std_lm24, paste0(path,"biocrates_lm24.rds")) # change name as per need


# Running the analysis for Diagnostic group comparisons
combinations <- list(c("CN", "MCI_stable"),
                    c("CN", "AD_converter"), 
                    c("CN", "AD"),
                    c("MCI_stable", "AD_converter"), 
                    c("AD", "AD_converter"), 
                    c("MCI_stable", "AD"))


# add on combinations

#combinations <- list()

for(i in seq_along(combinations)) {

    get_right_ones <- function(object, comb) {
        object <- object %>% mod_filter(which_data="regression", type="data", 
                regression_dxgrp_longi %in% comb)
        print(paste(comb, collapse="_v_"))
        print(object)
        return(object)
    }


    main_tp <- get_right_ones(object=my_analyser_main_tp, comb=combinations[[i]])
    long <- get_right_ones(object=my_analyser_long, comb=combinations[[i]])

    # multi-timepoint meta-analysis    

    std_meta <- main_tp %>% add_data(which_data="regression", type="col_data",x=diag_meta, id="id") %>%
        calc_lmm(which_data="regression", 
              name=paste(combinations[[i]], collapse="_v_"), 
              stratifications=NULL, 
              cols_for_meta=list(c(sub_pathway="Class")),
              num_cores=2)
    saveRDS(std_meta, paste0(path, "meta_diag_", paste(combinations[[i]], collapse="_v_"), ".rds"))

    # Time-interaction analysis

    std_lmm <- long %>% add_data(which_data="regression", type="col_data",x=diag_lmm, id="id") %>%
        calc_lmm(which_data="regression", 
             name=paste(combinations[[i]], collapse="_v_"), 
             stratifications=NULL, 
             cols_for_meta=list(c(sub_pathway="Class")),
             num_cores=2)
    saveRDS(std_lmm, paste0(path, "lmm_diag_", paste(combinations[[i]], collapse="_v_"), ".rds"))

    std_lmm_main_tp <- main_tp %>% add_data(which_data="regression", type="col_data",x=diag_lmm, id="id") %>%
        calc_lmm(which_data="regression", 
             name=paste(combinations[[i]], collapse="_v_"), 
             stratifications=NULL, 
             cols_for_meta=list(c(sub_pathway="Class")),
             num_cores=2)
    saveRDS(std_lmm_main_tp, paste0(path, "lmm_3tp_diag_", paste(combinations[[i]], collapse="_v_"), ".rds"))

    # linear models

    std_lm0 <- main_tp %>% add_data(which_data="regression", type="col_data",x=diag_lm, id="id") %>%
        calc_lm(which_data="regression", 
               name=paste(combinations[[i]], collapse="_v_"), 
               stratifications=NULL, 
              cols_for_meta=list(c(sub_pathway="Class")),
               num_cores=2, timepoint="0")
    saveRDS(std_lm0, paste0(path, "lm0_diag_", paste(combinations[[i]], collapse="_v_"), ".rds"))

    std_lm12 <- main_tp %>% add_data(which_data="regression", type="col_data",x=diag_lm, id="id") %>%
        calc_lm(which_data="regression", 
               name=paste(combinations[[i]], collapse="_v_"), 
               stratifications=NULL, 
               cols_for_meta=list(c(sub_pathway="Class")),
               num_cores=2, timepoint="12")
    saveRDS(std_lm12, paste0(path, "lm12_diag_", paste(combinations[[i]], collapse="_v_"), ".rds"))


    std_lm24 <- main_tp %>% add_data(which_data="regression", type="col_data",x=diag_lm, id="id") %>%
        calc_lm(which_data="regression", 
               name=paste(combinations[[i]], collapse="_v_"), 
               stratifications=NULL, 
               cols_for_meta=list(c(sub_pathway="Class")),
               num_cores=2, timepoint="24")
    saveRDS(std_lm24, paste0(path, "lm24_diag_", paste(combinations[[i]], collapse="_v_"), ".rds"))
    

}
```
