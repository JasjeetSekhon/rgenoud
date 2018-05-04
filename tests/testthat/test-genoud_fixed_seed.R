test_that("Tests the old version of genoud() where the seed are given", {
  # here we test a lot of things for genoud()
  set.seed(746612)
  #maximize the sin function
  sin1 <-
    genoud(
      sin,
      nvars = 1,
      max = TRUE,
      unif.seed = 912821,
      int.seed = 93058,
      print.level = 0
    )
  expect_equal(sin1$value, 1)
  expect_equal(sin1$par, 7.85398, tolerance = 1e-6)
  expect_equal(sin1$gradients, sin1$gradients)
  
  #minimize the sin function
  sin2 <-
    genoud(
      sin,
      nvars = 1,
      max = FALSE,
      unif.seed = 912821,
      int.seed = 93058,
      print.level = 0
    )
  expect_equal(sin2$value, -1)
  expect_equal(sin2$par, -1.570796, tolerance=1e-4)
  expect_equal(sin2$gradients, 5.944754e-11)
  
  #maximize a univariate normal mixture which looks like a claw
  claw <- function(xx) {
    x <- xx[1]
    y <- (0.46 * (dnorm(x, -1.0, 2.0 / 3.0) + dnorm(x, 1.0, 2.0 / 3.0)) +
            (1.0 / 300.0) * (dnorm(x, -0.5, .01) + dnorm(x, -1.0, .01) + 
                               dnorm(x, -1.5, .01)) +
        (7.0 / 300.0) * (dnorm(x, 0.5, .07) + dnorm(x, 1.0, .07) + 
                           dnorm(x, 1.5, .07)))
    return(y)
  }
  claw1   <-
    genoud(
      claw,
      nvars = 1,
      pop.size = 3000,
      max = TRUE,
      unif.seed = 912821,
      int.seed = 93058,
      print.level = 0
    )
  expect_equal(claw1$value, 0.4113123, tolerance=1e-6)
  expect_equal(claw1$par, 0.9995032, tolerance=1e-6)
  expect_equal(claw1$gradients, -1.277157e-06, tolerance=1e-6)

  # Maximize a bivariate normal mixture which looks like a claw.
  biclaw <- function(xx) {
    mNd2 <- function(x1, x2, mu1, mu2, sigma1, sigma2, rho)
    {
      z1 <- (x1 - mu1) / sigma1
      z2 <- (x2 - mu2) / sigma2
      w <- (1.0 / (2.0 * pi * sigma1 * sigma2 * sqrt(1 - rho * rho)))
      w <- w * exp(-0.5 * (z1 * z1 - 2 * rho * z1 * z2 + z2 * z2) / 
                     (1 - rho * rho))
      return(w)
    }
    x1 <- xx[1] + 1
    x2 <- xx[2] + 1
    
    y <- (0.5 * mNd2(x1, x2, 0.0, 0.0, 1.0, 1.0, 0.0) +
            0.1 * (
              mNd2(x1, x2, -1.0, -1.0, 0.1, 0.1, 0.0) +
                mNd2(x1, x2, -0.5, -0.5, 0.1, 0.1, 0.0) +
                mNd2(x1, x2, 0.0, 0.0, 0.1, 0.1, 0.0) +
                mNd2(x1, x2, 0.5, 0.5, 0.1, 0.1, 0.0) +
                mNd2(x1, x2, 1.0, 1.0, 0.1, 0.1, 0.0)
            ))
    
    return(y)
  }
  biclaw1 <-
    genoud(
      biclaw,
      default.domains = 20,
      nvars = 2,
      pop.size = 5000,
      max = TRUE,
      unif.seed = 912821,
      int.seed = 93058,
      print.level = 0
    )
  expect_equal(biclaw1$value, 1.671127, tolerance=1e-6)
  expect_equal(biclaw1$par, c(-1, -1), tolerance=1e-6)
  expect_equal(biclaw1$gradients, c(2.491051e-08, -3.427565e-08), 
               tolerance=1e-6)
  
})
