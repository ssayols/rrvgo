go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")

test_that("calculateSimMatrix returns a squared matrix", {
  expect_is(simMatrix, "matrix")
  expect_equal(nrow(simMatrix), ncol(simMatrix))
})
