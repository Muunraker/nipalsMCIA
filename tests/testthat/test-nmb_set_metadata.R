test_that("Incorrect metadata dimensions caught.", {
  bad_metadata <- data.frame(metadata_NCI60[1:20, ])
  expect_error(nmb_set_metadata(mcia_results_no_meta, bad_metadata))
})

test_that("Empty metadata caught.", {
  bad_metadata <- data.frame()
  expect_error(nmb_set_metadata(mcia_results_no_meta, bad_metadata))
})

test_that("Row name mismatch warning works.", {
  bad_metadata <- data.frame(metadata_NCI60[1:21, ])
  expect_warning(nmb_set_metadata(mcia_results_no_meta, bad_metadata))
})

test_that("Metadata embedded in object successfully.", {
  test_nmb <- nmb_set_metadata(mcia_results_no_meta, metadata_NCI60)
  expect_s3_class(test_nmb@metadata, "data.frame")
  expect_false(any(dim(test_nmb@metadata) == 0))
})