library(parallel)
library(snow)


# Loading AdaptSS and run_generation functions -----------------------------

source("./R/no_reproduction_system/AdaptSS_cervus_no_reproduction_system.R")


# Loading island data -----------------------------------------------------

island_area <- read.table("./data/island_area.txt", header = TRUE)


# Model constant parameters -----------------------------------------------

time <- nrow(island_area)
K.isl <- matrix(600, nrow = time, ncol = 1)
input_wrigth <- 0.1

# Setting clusters --------------------------------------------------------

cl <- makeCluster(2, type = "SOCK")
clusterExport(cl = cl, ls())


# Running simulation in parallel ------------------------------------------

results <- parLapplyLB(cl, X = 1:1000, function(i){
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


# Closing clusters --------------------------------------------------------

stopCluster(cl)

# Detect potential errors 
cores <- detectCores() - 2
stopped_runs <- which(sapply(results, function(i) class(i) == "try-error") == TRUE)

if(length(stopped_runs) > cores){
  cl <- makeCluster(cores)} else{
    cl <- makeCluster(length(stopped_runs))  
  }

clusterExport(cl = cl, ls())

# Re-running simulations with errors --------------------------------------

while(length(stopped_runs > 0)){
  
  stopped_rerun <- parLapplyLB(cl, X = 1:length(stopped_runs), function(i){
    try(run_generation( 
      K.isl = K.isl,
      area = island_area[,2],
      Precol = island_area[,3],
      time = time, 
      input_wrigth = input_wrigth
    ))
  }
  )
  
  
  # Save error fixing -------------------------------------------------------
  
  for(i in 1:length(stopped_runs)){
    results[[stopped_runs[i]]] <- stopped_rerun[[i]]
  }
  
  stopped_runs <- which(sapply(results, function(i) class(i) == "try-error") == TRUE)    
}

stopCluster(cl)


# Export output -----------------------------------------------------------

#Simulation results
resultado <- do.call("rbind", lapply(results, function(i) i[[1]]))
#write.table(resultado, "./baseline/Cervus_results__k.txt")

# Mean trait evolution through time
mean_P<- do.call("cbind", lapply(results, function(i) i[[2]]))
#write.table(resultado, "./baseline/Cervus_meanP.txt")
