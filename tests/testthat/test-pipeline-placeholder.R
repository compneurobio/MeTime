test_that("placeholder: PCA pipeline runs on example data", {
  skip("PLACEHOLDER: enable after validating runtime profile in CI")

  data("humet_object")
  out <- humet_object %>%
    mod_trans_zscore(which_data = "humet_subset_data") %>%
    calc_dimensionality_reduction_samples(
      which_data = "humet_subset_data",
      type = "PCA",
      cols_for_meta = c("Factor.Challenge.Value.", "Factor.Challenge.Value.Day."),
      name = "PCA_samples"
    )

  expect_s4_class(out, "metime_analyser")
})
