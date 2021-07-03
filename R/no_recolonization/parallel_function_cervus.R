# Libraries
library(parallel)
library(snow)

# Loading adapt function
source("./no_recolonization/AdaptSS_cervus_K_no_recolonization.R")
island_area <- read.table("./island_area.txt", header = TRUE)
# Model constant parameters

nsim <- 1 
Dmax <-1
Dmin <-0.01
coef.k <-as.vector(lm(c(0.01,1) ~ c(log(50),log(1000)))$coefficients)
time <- nrow(island_area)
K.isl <- matrix(2000, nrow = time, ncol = 1)
input_wrigth <- 0.1

#Run parallel
#source("https://raw.githubusercontent.com/csdambros/R-functions/master/passSSHid.R")
#passSSHid()

#cl <- makeMPIcluster(14)
cl <- makeCluster(10, type = "SOCK")
clusterExport(cl = cl, ls())

teste <- parLapplyLB(cl, X = 1:1000, function(i){
  try(run_generation( 
    Dmax = Dmax, 
    Dmin = Dmin, 
    coef.k = coef.k, 
    K.isl = K.isl,
    area = island_area[,2],
    Precol = rep(0, nrow(island_area)),
    time = time, 
    input_wrigth = input_wrigth
    ))
}
)

stopCluster(cl)
# Detect errors 
cores <- detectCores() - 2
stopped_runs <- which(sapply(teste, function(i) class(i) == "try-error") == TRUE)
  if(length(stopped_runs) > cores){
    cl <- makeCluster(cores)} else{
    cl <- makeCluster(length(stopped_runs))  
    }

clusterExport(cl = cl, ls())
# Re-run error replicates
while(length(stopped_runs > 0)){

  stopped_rerun <- parLapplyLB(cl, X = 1:length(stopped_runs), function(i){
    try(run_generation( 
      Dmax = Dmax, 
      Dmin = Dmin, 
      coef.k = coef.k, 
      K.isl = K.isl,
      area = island_area[,2],
      Precol = island_area[,3],
      time = time, 
      input_wrigth = input_wrigth
    ))
  }
  )

# save error replicates

for(i in 1:length(stopped_runs)){
  teste[[stopped_runs[i]]] <- stopped_rerun[[i]]
}

stopped_runs <- which(sapply(teste, function(i) class(i) == "try-error") == TRUE)    
}
stopCluster(cl)

# Format results as a data.frame and export it
resultado <- do.call("rbind", lapply(teste, function(i) i))
write.table(resultado, "./no_recolonization/Cervus_results_no_recolonization.txt")

