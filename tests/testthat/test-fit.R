test_that("blblm works with dataset", {
  a = blblm(speed ~ dist, cars, B = 10)
  expect_lte(a$estimates$`3`[[3]]$sigma, 10)
  expect_gt(a$estimates$`3`[[3]]$sigma, 0)
  #learned the trick where you use multiple $ in 141B
  expect_length(a$estimates, 10)
  expect_named(a, c("estimates", "formula"), ignore.case = TRUE)
})


test_that("blblm works with files", {
  car1 <- write.csv(cars, "car1.csv", row.names = FALSE)
  car2 <- write.csv(cars, "car2.csv", row.names = FALSE)
  b = blblm(speed ~ dist, c("car1.csv", "car2.csv"), B = 10)
  expect_lte(b$estimates[[2]][[8]]$sigma, 8)
  expect_gt(b$estimates[[1]][[1]]$coef[1], 0)
  expect_length(b$estimates[[2]], 10)
  expect_named(b, c("estimates", "formula"), ignore.case = TRUE)
  file.remove(c("car1.csv", "car2.csv"))
})

test_that("blblm works with log reg", {
  iris1 = iris[1:100,]
  iris1$Species <- as.numeric(iris1$Species)-1 #the -1 takes all the numbers from 1 or 2 to 0 or 1
  cc = blblm(Species ~ Sepal.Length, iris1, B = 50, m = 2, model_choice = "glm")
  expect_lte(cc$estimates[[2]][[8]]$sigma, 1)
  #expect_gt(blbsigma(cc, 0)) NOT THE SAME CLASS
  expect_length(cc$estimates[[2]], 50)
  expect_named(cc, c("estimates", "formula"), ignore.case = TRUE)
  expect_equal(cc$formula, Species ~ Sepal.Length)
})


test_that("blblm works with log reg on 2 clusters", {
  iris1 = iris[1:100,]
  iris1$Species <- as.numeric(iris1$Species)-1 #the -1 takes all the numbers from 1 or 2 to 0 or 1
  d = blblm(Species ~ Sepal.Length, iris1, B = 50, m = 2, model_choice = "glm", parallelization = TRUE, spec = 2)
  expect_lte(d$estimates[[2]][[8]]$sigma, 1)
  expect_lt(d$estimates[[1]][[1]]$coef[1], 0)
  expect_length(d$estimates[[2]], 50)
  expect_named(d, c("estimates", "formula"), ignore.case = TRUE)
  expect_equal(d$formula, Species ~ Sepal.Length)
})