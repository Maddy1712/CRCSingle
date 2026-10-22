test_that("single sample equals batch", {
  x <- crc_example(); b <- predict_subtype(x); s <- predict_subtype(x[, 3, drop = FALSE])
  expect_equal(s$subtype, b$subtype[3])
})
test_that("output structure", {
  x <- crc_example(); b <- predict_subtype(x)
  expect_true(all(c("sample", "subtype", "margin", "call") %in% names(b)))
  expect_equal(nrow(b), ncol(x))
  expect_true(all(b$subtype %in% c("S1", "S2", "S3")))
})
test_that("lite and cohort modes run", {
  x <- crc_example()
  expect_equal(nrow(predict_subtype(x, model = "lite")), ncol(x))
  expect_equal(nrow(predict_subtype(x, cohort_adjust = TRUE)), ncol(x))
})
test_that("bad input fails clearly", {
  x <- crc_example()
  expect_error(predict_subtype(x[1:100, ]))
  expect_error(predict_subtype(x[, 1:2], cohort_adjust = TRUE))
})
