library(testthat)
source("./R/no_reproduction_system/AdaptSS_cervus_no_reproduction_system.R")

h2m <-  0.6521528
vp0 <- 122.3638
w2 <- 8631.647
optO <-  195.0168
optN <- 62.9603
B <- 5
KN <- 600
coef_ID <- 1.44
Fw <- 0.1
NR <- 14
DR <- 0.19
vmu <- 0.02304844
bp <- 0.26
cv <- runif (1,0.04,0.06)
sdp <- cv*optO # sd equals the coefficient of variation of 5%
vp <- sdp*sdp # 6.25
va <- h2m*vp
Gv <- rnorm(200,optO,sqrt(va)) 
fs <-seq(1:length(Gv))
ib <-rep(0,length(Gv))
G <-cbind(Gv,fs,ib)

res <- AdaptSS(G,h2m,vp0,w2,optO,optN,B,KN,coef_ID,Fw,NR,DR,vmu,bp)


test_that("AdaptSS result size", {
  
  expect_equal(length(res), 2)
  
})

test_that("Result class", {
  
  expect_equal(dim(res[[1]])[2], 3 )
  
})




