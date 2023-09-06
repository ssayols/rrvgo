go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")

test_that("reduceSimMatrix returns a data.frame of the same length as the simMatrix", {
  expect_is(reducedTerms, "data.frame")
  expect_equal(nrow(reducedTerms), nrow(simMatrix))
})

test_that("reducedSimMatrix draws all parents from list", {
  expect_true(all(reducedTerms$parent %in% reducedTerms$go))
  expect_true(all(reducedTerms$parentTerm %in% reducedTerms$term))
})

test_that("reducedSimMatrix column types", {
  expect_is(reducedTerms$cluster, "integer")
  expect_is(reducedTerms$termUniqueness, "numeric")
  expect_is(reducedTerms$termDispensability, "numeric")
  expect_is(reducedTerms$size, "numeric")
})

test_that("Term uniquness and despisablility are anti-correlated",{
  expect_true(cor(reducedTerms$termUniqueness, reducedTerms$termDispensability) < 0)
})
