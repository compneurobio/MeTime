require(MeTime)

#loading the imputed analyser object

load("adni_nmr_data")
# Dataset of interest
which_data <- "nmr_data"


my_analyzer <- adni_nmr_data %>%
  add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long"), 
                                which_data=which_data) %>%
   mod_filter(which_data = which_data,
             type="row_data",
             !DXGrp_long %in% c("CN_MCI_convert","MCI_CN_revert","mix_AD_stable","mix_CN_stable","mix_MCI_stable")) %>%
  mod_mutate(which_data = which_data, type="row_data", 
             APOEGrp=as.factor(APOEGrp),
             PTGENDER=as.factor(PTGENDER),
             PTEDUCAT=as.factor(PTEDUCAT),
             DXGrp_long = ifelse(DXGrp_long %in% c("CN_AD_convert", "MCI_AD_convert"), "AD_converter", DXGrp_long),
             DXGrp_long = as.factor(DXGrp_long, levels=c("CN_stable", "CN_MCI_convert","MCI_stable", "AD_converter", "AD_stable"))) %>%
  mod_trans_zscore(which_data=which_data) %>%
  calc_dimensionality_reduction_samples(which_data = which_data, type="PCA", 
                                        cols_for_meta=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_long", "Age", "BMI"),
                                        name="PCA_samples") %>%
  mod_generate_plots(type="PCA") %>%
  calc_dimensionality_reduction_metabs(which_data = which_data,
                                       cols_for_meta = list(biocrates_data=c(id="id", sub_pathway="Class")),
                                       type="PCA",
                                       name="PCA_metabs") %>%
  mod_generate_plots(type="PCA") %>%
  calc_dimensionality_reduction_samples(which_data = which_data, type="UMAP", 
                                        cols_for_meta=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_long", "Age", "BMI"),
                                        name="UMAP_samples") %>%
  mod_generate_plots(type="UMAP") %>%
  calc_dimensionality_reduction_metabs(which_data = which_data,
                                       cols_for_meta = list(biocrates_data=c(id="id", sub_pathway="Class")),
                                       type="UMAP",
                                       name="UMAP_metabs") %>%
  mod_generate_plots(type="UMAP") %>%
  calc_dimensionality_reduction_samples(which_data = which_data, type="tSNE", 
                                        cols_for_meta=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_long", "Age", "BMI"),
                                        name="tSNE_samples") %>%
  mod_generate_plots(type="tSNE") %>%
  calc_dimensionality_reduction_metabs(which_data = which_data,
                                       cols_for_meta = list(biocrates_data=c(id="id", sub_pathway="Class")),
                                       type="tSNE",
                                       name="tSNE_metabs") %>%
  mod_generate_plots(type="tSNE")

my_analyzer %>%
  write_report(file="dimensionality_reduction.html", title="ADNI NMR Dimensionality Reduction")

