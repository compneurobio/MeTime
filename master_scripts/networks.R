require(MeTime)

#loading the imputed analyser object
load("adni_q500_data")
which_data <- "biocrates_data"
networks <- adni_q500_data %>%
  add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long"), 
                                which_data=which_data) %>%
  mod_merge_row_data_and_data(which_data=which_data, 
        cols_list=list(data=NULL, row_data=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", 'TOTAL_C', 'HDL_C')),
        name="ggm_data", append=T) %>%
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

networks %>% 
    write_report(title="ADNI NMR networks", file="networks.html")




