require(MeTime)

#loading the imputed analyser object
load("adni_nmr_data")
which_data <- "nmr_data"
adni_nmr_data <- adni_nmr_data %>%
  add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long"), 
                                which_data=which_data) %>%
  mod_merge_data_and_row_data(which_data="nmr_data", 
        cols_list=list(data=NULL, row_data=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT")),
        name="ggm_data") %>%
  mod_mutate(which_data = "ggm_data", type="row_data", time=paste0("t", time)) %>% 
  mod_mutate(which_data = "ggm_data", type="data", 
                                          APOEGrp=as.numeric(APOEGrp), 
                                          PTGENDER=as.numeric(PTGENDER),
                                          PTEDUCAT=as.numeric(PTEDUCAT)) %>%
    calc_ggm_genenet(which_data = "ggm_data", threshold = "li", all=FALSE, 
                     cols_for_meta = list(biocrates_data=c(sub_pathway="Class")),
                     covariates = c("Age", "BMI", "APOEGrp", "PTGENDER", "PTEDUCAT", "TOTAL_C", "HDL_C"),
                     stratifications = list(time=c("t0", "t12", "t24")),
                     name="genenet_ggm_results") %>%
    mod_generate_plots(type="ggm_longitudinal")

adni_nmr_data %>% 
    write_report(title="ADNI NMR networks", file="networks.html")




