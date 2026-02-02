require(MeTime)

#loading the imputed analyser object

load("adni_q500_data")
# Dataset of interest
which_data <- "biocrates_data"

res_conservation_index <- adni_q500_data %>%
  add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_long", "PTGENDER", "Age", "BMI", "PTEDUCAT"), 
                                which_data=which_data) %>%
  mod_mutate(which_data = which_data, type="row_data", 
             APOEGrp=as.factor(APOEGrp),
             PTGENDER=as.factor(PTGENDER),
             PTEDUCAT=as.factor(PTEDUCAT),
             DXGrp_longi = ifelse(DXGrp_long == "CN_AD_convert", "AD_converter", DXGrp_long),
             DXGrp_longi = ifelse(DXGrp_long == "MCI_AD_convert", "AD_converter",DXGrp_long)) %>%
  mod_filter(which_data = which_data,
             type="row_data",
             !DXGrp_long %in% c("CN_MCI_convert","MCI_CN_revert","mix_AD_stable","mix_CN_stable","mix_MCI_stable")) %>%
  mod_trans_zscore(which_data=which_data) %>% 
  calc_conservation_metabotype(which_data=which_data, 
                               stratifications=list(time=c("t0", "t12", "t24")),
                               verbose=F, 
                               cols_for_meta=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_longi", "PTGENDER", "Age", "BMI"), 
                               name="Conservation_Metabotype_scaled") %>%
  mod_generate_plots(type="CI_metabotype", .interactive=T) %>%
  calc_conservation_metabolite(which_data=which_data, 
                               stratifications=list(time=c("t0", "t12", "t24")),
                               verbose=F, 
                               cols_for_meta = list(biocrates_data=c(id="id", sub_pathway="Class")), 
                               name="Conservation_Metabolite_scaled") %>%
  mod_generate_plots(type="CI_metabolite", .interactive=T)


res_conservation_index %>%
  write_report(file="conservation_index.html", title="ADNI Q500 Conservation index")
  



  