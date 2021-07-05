library(testthat)

source("./R/no_plasticity/AdaptSS_cervus_K_no_plasticity.R")


# Loading island data -----------------------------------------------------

island_area <- read.table("./data/island_area.txt", header = TRUE)


# Model constant parameters -----------------------------------------------

time <- 10
K.isl <- matrix(600, nrow = time, ncol = 1)
input_wrigth <- 0.1

# Running simulation in parallel ------------------------------------------

results <- lapply(1:time, function(i){
  try(run_generation( 
    K.isl = K.isl,
    area = island_area[,2],
    Precol = island_area[,3],
    time = time, 
    input_wrigth = input_wrigth
  )
  )
}
)

test_that("run_generation output size", {
  
  expect_equal(length(results), time)
  
})

test_that("run_generation simulation size", {
  
  expect_equal(length(results[[1]]), 2)
  
})

test_that("run_generation output", {
  
  expect_equal(unique(sapply(results, function(i) class(i$output))), "numeric")
  
})




