require(MeTime)

# Example regression workflows: meta-analysis, time-interaction, GAMM, and baseline LM.

# Loading the imputed analyser object.
load("adni_nmr_data")
# Load covariate sheets.
load("covariate_sheets")

meta <- covariate_sheets$meta
lmm <- covariate_sheets$lmm
lm_model <- covariate_sheets$lm
gamm <- covariate_sheets$lmm


which_data <- "nmr_data"
object <- adni_nmr_data
# In this example we show you different kinds of regressions you can perform in R.

adni_nmr_data <- adni_nmr_data %>%
  	add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long"), 
                                which_data=which_data) %>%
  	mod_filter(which_data=which_data, type="row_data", 
            !DXGrp_long %in% c("mix_CN_stable", "mix_MCI_stable", "mix_AD_stable", "MCI_CN_revert")) %>% # goal is to remove all samples with atypical diagnostic trajectories
    mod_mutate(which_data="phenotype_data", type="data",
             DXGrp_long=ifelse(DXGrp_long %in% c("CN_AD_convert", "MCI_AD_convert"), "AD_converter", DXGrp_long)) %>%
    mod_filter_timepoints(which_data=which_data, timepoints=c("t0", "t12", "t24")) %>%
    mod_trans_zscore(which_data=which_data) %>%
    mod_merge_data_and_row_data(which_data=which_data, 
        cols_list=list(data=NULL, row_data=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long")),
        name="nmr_data_merged") %>%
    mod_merge_data(which_data=c("phenotype_data","medication_data", "nmr_data_merged"), 
                 filter_sample="biocrates_data", name="regression", append=TRUE,
                 cols_list=list(nmr_data=names(object@list_of_data[["nmr_data_merged"]]), 
                                phenotype_data=c("CSF_Roche_Abeta42", "LEntCtx_vol", "LEntCtx_Thick", "REntCtx_vol", "REntCtx_Thick", 
                                				 "CSF_Roche_Tau", "CSF_Roche_PTau", "ADAS13", "ADNI_MEM", "ADNI_EF", "ADNI_LAN", "HippVol", 
                                				 "GlobalCtx_Thick", "Global_vol", "FDG_Cing_Mean", "FDG_LAng_Mean", "FDG_RAng_Mean", 
                                				 "FDG_LTemp_Mean", "FDG_RTemp_Mean", "DXGrp")
                                medication_data=names(object@list_of_data$medication_data),
                                nmr_data=c("TOTAL_C", "HDL_C"))) %>% 
    mod_mutate(which_data="regression", 
             "diagnostic_group" = factor(DXGrp_long, levels=c("CN", "CN_MCI_convert","MCI_stable", "AD_converter", "AD")),
             "abeta42" =  CSF_Roche_Abeta42 %>% as.numeric() %>% scale(),
             "entorhinal_thickness" = (LEntCtx_Thick + REntCtx_Thick) %>% as.numeric() %>% as.numeric() %>% scale(),
             "entorhinal_volume" = (LEntCtx_vol + REntCtx_vol) %>% as.numeric() %>% as.numeric() %>% scale(),
             "ptau" = CSF_Roche_PTau %>% as.numeric() %>% log2() %>% as.numeric() %>% scale(),
             "tau" = CSF_Roche_Tau %>% as.numeric() %>% log2() %>% as.numeric() %>% scale(),
             "tau_by_abeta42" = as.numeric(CSF_Roche_Tau)/as.numeric(CSF_Roche_Abeta42) %>% log2() %>% as.numeric() %>% scale(),
             "adas_cog13" = ADAS13 %>% sqrt() %>% as.numeric() %>% scale(),
             "memory_scare" = ADNI_MEM %>% scale(),
             "executive_function_score" = ADNI_EF %>% scale(),
             "language_sore" = ADNI_LAN %>% scale(),
             "hippocampal_volume" = HippVol %>% scale(),
             "grey_matter_thickness" = GlobalCtx_Thick %>% scale(),
             "grey_matter_volume" = Global_vol %>% scale(),
             "fdgpet" = (FDG_Cing_Mean + FDG_LAng_Mean + FDG_LTemp_Mean + FDG_RTemp_Mean + FDG_RAng_Mean) %>% as.numeric() %>% scale(),
             "diagnosis_at_visit" = as.numeric(ifelse(DXGrp %in% c(1,5), 0, ifelse(DXGrp %in% c(2,3), 1, 2))))

# Running regression analysis in different objects otherwise the objects would blow up in size.
# Multi-timepoint meta analysis.
std_meta <- adni_nmr_data %>% add_data(which_data="regression", type="col_data",x=meta, id="id") %>%
    calc_lmm(which_data="regression", 
            name="adni_nmr_meta", 
            stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Group")), # change column name as per need
            num_cores=12) %>%
    mod_generate_plots(type="regression", interactive=T)

# Time-interaction analysis.

std_lmm <- adni_nmr_data %>% add_data(which_data="regression", type="col_data",x=lmm, id="id") %>%
    calc_lmm(which_data="regression", 
            name="biocrates_time_interaction", 
            stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Group")), # change column name as per need
            num_cores=12) %>%
    mod_generate_plots(type="regression", interactive=T)

# GAMM analysis here - all continuous covariates are smoothed.

std_gamm <- adni_nmr_data %>% add_data(which_data="regression", type="col_data",x=lmm, id="id") %>%
    calc_gamm(which_data="regression", 
            name="biocrates_time_interaction", 
            stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Group")), # change column name as per need
            num_cores=12) %>%
    mod_generate_plots(type="regression", interactive=T)

# Linear models for cross-sectional analysis - example shown for baseline.

std_lm0 <- adni_nmr_data %>% add_data(which_data="regression", type="col_data",x=lm_model, id="id") %>%
    calc_lm(which_data="regression", 
             name="biocrates_lm0", 
             stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="Class")), # change column name as per need
             num_cores=2, timepoint="t0") %>%
    mod_generate_plots(type="regression", interactive=T)

# Run if you want to visualize reports.

std_meta %>% write_report(title="ADNI NMR meta-analysis results", file="meta_analysis.html")
std_lmm %>% write_report(title="ADNI NMR time-interaction results", file="time_interaction_analysis.html")
std_gamm %>% write_report(title="ADNI NMR GAMM results", file="GAMM_analysis.html")
std_lm0 %>% write_report(title="ADNI NMR baseline association results", file="baseline_analysis.html")
