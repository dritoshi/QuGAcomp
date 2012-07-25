test_that("check metrics", {
	
  # demo data
  first  <- c(0,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,0)
  second <- c(0,0,0,0,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1)
  my.df <- data.frame(first=first, second=second)
  my.table <- table(my.df)

  first.Rle  <- Rle(first)
  second.Rle <- Rle(second)

  q <- qugacomp(first, second)
  expect_equal(
    cor(first, second),
    phiCoef(q)
  )
  expect_equal(
    cor(first, second),
    pearsonCoef(q)
  )

}
)