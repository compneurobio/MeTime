test_that("MeTime loads bundled example object", {
  data("humet_object")
  expect_s4_class(humet_object, "metime_analyser")
  expect_true("humet_subset_data" %in% names(humet_object@list_of_data))
})
