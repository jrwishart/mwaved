context("Blur detection tests")

n <- 64
m <- 3
shape <- seq(from = 0.25, to = 0.75, length = m)
scale <- rep(0.25, m)
BAwidth <- 1/sqrt(c(89, 343, 199))
Gdirect <- directBlur(n, m)
Gboxcar <- boxcarBlur(n, BAwidth)
Gsmooth <- gammaBlur(n, shape, scale)

# Check correct blur types are identified
test_that("Correct identification", {
  expect_output(detectBlur(Gdirect),"direct")
  expect_output(detectBlur(Gboxcar),"box.car")
  expect_output(detectBlur(Gsmooth),"smooth")
})

# Warnings are thrown for incorrect specifications
test_that("Throw warnings when input inconsistent", {
  expect_warning(blurCheck(Gboxcar, 'smooth'))
})
