AdaptSS <-function(G,h2m,vp0,w2,optO,optN,B,KN,coef_ID,Fw,NR,DR,vmu,bp){
  
  #G <-vector of genetic values
  #h2 <-heritability
  #vg0 <-inicial genetic variance
  #w2 >-width of adaptive surface
  #optO <-inicial peak
  #optN <-new adaptive peak
  #B <-fecundity (maximum population growth rate/replacemente rate per female)
  #KN <-carrying capacity
  #Fw <-inbreeding in the former generation
  #NR <-Number of immigrants  
  #DR <-probability of demographic rescue
  #vmu <-mutational variance (to be multiplied by vG; proportion)
  
  GS <-G[sample(nrow(G)),] 
  G <-GS[,1]
  FS <-GS[,2]
  IB <-GS[,3]+Fw
  
  
  #genetic and environmental values of the trait
  va <-var(G)
  vp <-va/h2m
  ve <- vp - va
  E <-rnorm(length(G),0,sqrt(ve))
  PG <- G + E
  
  #Plasticity
  react <-summary(lm(c(mean(G),optN)~c(0,1)))$coefficients[,1]
  max.slp <-react[2]
  #bp <-0 #50% de plasticidade
  react.b <-bp*max.slp
  #P <- (G + rnorm(length(G),react.b,(abs(react.b)*0.05)*1)) + E
  P <- (G + react.b*1) + E #Schmidt & Guillaume 2017 #slope fixed
  
  rh2 <-var(G)/(var(G)+var(E)) #this is a constant H2 model, with this realized
  NpK <-length(P)
  
  #calculating reduction in fitness due to inbreeding depression
  rW <- 0 + (IB * coef_ID) #increase in mortality due to ID for different inbreeding levels W = W * (1 - rW)
  rW <-ifelse(rW > 1,1,rW)
  rW <-ifelse(rW < 0,0.0001,rW) 
  
  #fitness per individual and survival 
  DPi <- abs(P - optN)
  wPi <- exp(-((DPi^2)/(2*w2))) 
  wPi <-wPi * (1-rW) #reducing fitness due to ID
  meanFit <-mean(wPi)
  grad <-summary(lm(wPi~P))$coefficients[2,1]
  
  # reduction in mean body size due to inbreeding (see Huismann and others...)
  rbs <-numeric()
  rbs <-ifelse(IB > 0.5,0.956,1)
  G <-G*rbs
  
  
  # survival depends on fitness in a continuous manner...with smaller than optN also truncated
  psurv <-runif(length(wPi),0,1)
  surv <-ifelse(wPi > psurv,1,0)
  #surv <-ifelse(P < (optN - (1.96*sqrt(vp0))),0,surv) #trava para bichos pequenos...
  #wPi <- ifelse(wPi > 0.5,wPi,0)
  
  P <-P[which(surv==1)]
  G <-G[which(surv==1)]
  wPi <-wPi[which(surv==1)] 
  FS <-FS[which(surv==1)] 
  #meanFit <-mean(wPi)
  N <-length(P)
  meanFit <-mean(wPi)
  
  #density dependent fitness in number of offspring
  NFG1 <-wPi * (B)^(1-(N/KN)) #Chevin & Lande 2011 Evolution #Density dependent
  lambda = mean(wPi)*(B)^(1-(N/KN))
  qw <-quantile(wPi,p=(floor(lambda) - lambda + 1))
  NFG1 <-ifelse(wPi < qw,floor(NFG1),ceiling(NFG1))
  
  
  #selecting first and last individuals as couples and getting mean fitness
  mdP <-numeric()
  mdG <-numeric()
  nF1 <-numeric()
  INB <-numeric()
  
  for(i in 1:floor(N/2)){
    mdP[i] <-(P[i]+P[(N/2)+i])/2
    mdG[i] <-(G[i]+G[(N/2)+i])/2
    #nF1[i] <- round(2*((NFG1[i]+NFG1[(N/2)+i])/2))
    nF1[i] <- (NFG1[i]+NFG1[(N/2)+i])
    INB[i] <-ifelse((FS[i] == FS[(N/2)+i]),1,0)
    
  }
  
  
  #INB[2] <-1 #SO PARA TESTAR LOGO EM t0
  
  nF1 <-ifelse(nF1 < 2,2,nF1) #?????
  nF1 <- na.omit(nF1)
  
  #Reproduction
  
  #New individuals in Population in t+1
  
  nF <-numeric()
  G0 <-1
  #Fw is the inbreending from previous generation
  #Fw <-0
  
  for(i in 1:length(nF1)){
    
    GF <-numeric()
    
    for(j in 1:nF1[i]){	
      GF_mean <- mdG[i]
      Gmean <- sum(G)/length(G)
      var_G <- sum((G-Gmean)^2/(length(G)-1))
      GF_var <- sqrt((var_G*(1-Fw))/2)
      GF[j] <- rnorm(1,GF_mean, GF_var) 
    }
    G0 <-c(G0,GF)
  }
  G1 <- G0[-1] #vector of G values in the offspring
  
  
  #Adding mutational variance (mutation kernel)
  vm <- vmu*var(G1) #Kemper et al. 2012
  M <-rnorm(length(G1),0,sqrt(vm))
  G1 <-G1 + M
  
  #New family structure
  c0 <-numeric()
  i0 <-numeric()
  c0[1] <-50
  i0[1] <-50
  
  for(i in 1:length(nF1)){
    fam <-rep(i,nF1[i]) 
    indF <-rep(INB[i],nF1[i])
    c0 <- c(c0,fam) 
    i0 <-c(i0,indF)
  }
  FS1 <-c0[-1]
  IB1 <-i0[-1] + Fw
  
  #inbreeding na popula??o parental (n?o baseado na possibilidade dos filhotes)
  csF <-ifelse(INB==0,1,INB*nF1)
  Fw0 <-Fw
  quad <-numeric()
  for(j in 1:length(nF1)){
    quad[j] <- ((csF[j]^2 - csF[j])/2) * 0.5
  }
  Fw <- (sum(quad)/ (((length(csF)^2)-length(csF))/2))+(Fw0*(1-(vmu^2)))
  
  Fw <-ifelse(Fw > 1,1-vmu,Fw)
  Fw <-ifelse(Fw < 0,0,Fw)
  
  #Selection gradients and differentials (for analysis)
  Resp <-mean(G1)-mean(G)
  Sdif <-Resp/rh2
  #beta.mean <-mean(G)*(Sdif/var(P))
  beta.mean <-grad*mean(G)
  
  
  
  N1 <-length(G1)
  
  #Final addition for F1 and vectors
  #Demographic rescue / recolonization of continent
  
  #DR <- 0.05 #probability of demographic rescue
  #NR <-1  #number of imigrants
  sd0 <-sqrt(vp0*h2m)
  
  if(runif(1,0,1) <= DR){
    G1rc <-c(G1,rnorm(NR,optO,sd0)) #final vector
    FS1rc <-c(FS1,rep(max(FS1)+1:NR))
    IB1rc <-c(IB1,rep(0,NR))
    G1F <-cbind(G1rc,FS1rc,IB1rc)
    #Fw <-(Fw*(length(G1)))/(length(G1)+NR) #the dispersal decreases inbreeding
    mig <- 1
    
  }else{
    G1F <-cbind(G1,FS1,IB1)
    mig <- 0
  }
  
  #G1 <- na.omit(G1)
  
  
  #OUTPUT
  
  #OUTPUT
  out <- list(G = G1F, parameters = c(inbreed = Fw,
                                      meanFit = meanFit,
                                      mrh2 = rh2,
                                      mig01 = mig,
                                      Sel.dif = Sdif,
                                      beta.mean = beta.mean,
                                      RL = N1/NpK,
                                      meanP = mean(P))
  )
  
  return(out)
}




run_generation <- function(Dmax, 
                           Dmin, 
                           coef.k, 
                           K.isl,
                           area,
                           Precol,
                           time, 
                           input_wrigth,
                           fixed = FALSE,
                           plot = FALSE
                           ){
  
  
  # Output vector
  out <-rep(NA, 18)
  
  # Sample parameters
  
  oldP <- runif(1,180, 220)
  newP <- runif(1, 33, 39)
  h2m <- runif(1,0.6,0.85)
  cv <- runif (1,0.04,0.06)
  sdp <- cv*oldP #sd equals the coefficient of variation of 5%
  vp <- sdp*sdp #6.25
  va <- h2m*vp
  vpi <- vp
  vm <- runif(1,0.02,0.04)
  w2 <- runif(1, (100*(vp*h2m)), (125*(vp*h2m)))  #weak to moderate Burger & Lynch 
  bp.max <-runif(1,0.1,0.5)
  coef.bp <-summary(lm(c(0,bp.max)~c(0,0.9)))$coefficients[,1] #ratio N/K
  
  
  #Sampling from a skewed mortality distribution due to ID
  mortID <- 0.72 #define a value for mortality...73 in Bittles & Neel
  coef_ID <-lm(c(0,mortID)~c(0,0.5))$coefficients[2]
  
  
  #fecund <-round(runif(1,6,16))
  Ni <-round(runif(1,250,500))
  Nrecol <-round(runif(1,1,20))
  #Precol <-runif(1,0.1,0.25) # migration rate
  
  Gv <-rnorm(Ni,oldP,sqrt(va)) 
  fs <-seq(1:length(Gv))
  ib <-rep(0,length(Gv))
  G <-cbind(Gv,fs,ib)
  
  #sets of outputs and vars
  
  Np <-numeric(time)
  va <-numeric(time)
  mean.bs <-numeric(time)
  mean.w <-numeric(time)
  varT <-numeric(time)
  adap <-rep(NA, time)
  sgrad1 <-numeric(time)
  meanP <- numeric(time)
  coef.peak <-summary(lm(c(oldP,newP)~c(0,0.9)))$coefficients[,1]
  
  #begin time series
  temp <- 0
  for(t in 1:time){
    # escrever num arquivo externo os passos da corrida (logfile)  
    #setTxtProgressBar(pb2, t)
    
    if(t == 1){
      input_wrigth <- input_wrigth
      Np[t] <-length(G[,1])
    }
    
    if(length(G[,1]) < 10){
      break
    }
    
    #Shifting peak
    Isl.Peak <-as.numeric(c(1,(Np[t]/mean(K.isl[1:t])))%*%coef.peak) # modifiquei K de 1:t
    ######################################################################
    Isl.Peak <-ifelse(Isl.Peak < newP,newP,Isl.Peak) # tem que tornar geral
    ######################################################################
    Isl.Peak <-ifelse(Isl.Peak > oldP,oldP,Isl.Peak)
    
    bp <- bp.max
   # bp <-as.numeric(c(1,(Np[t]/mean(K.isl[1:t])))%*%coef.bp) # modifiquei K de 1:t
    
    #Fecundity by generation, demographic stochasticity
    fecund <-rpois(1, 7)
    
    res <- try(AdaptSS(G,h2m,vpi,w2,oldP,Isl.Peak,fecund,round(K.isl[t]),coef_ID,input_wrigth,Nrecol,Precol[t],vm, bp))
    if(class(res) == "try-error"){
      break
    }
    Np[t+1] <-length(res$G[,1])
    input_wrigth <-res$parameters[1]
    mean.w[t] <- res$parameters[2]
    sgrad1[t] <-res$parameters[6]
    meanP[t] <- res$parameters[8]
    
    # k by body size #################################################
    
    if(t != time){K.isl[t + 1] <- (20.85366 - (0.079 * meanP[t])) * area[t]}
    
    ##################################################################
    G <-res$G
    #adap[t] <-ifelse(mean(G[,1]) <= newP +(1*sqrt(vpi)), 1, 0)
    adap[t] <- ifelse(mean(G[,1]) <= 40, 1, 0)
    
    #if(Np[t+1] <= 10){
    #  break
    #}
    
   temp <- temp + 1
    if(adap[t] == 1){
      break
    }
    if(plot == TRUE){
    plot(1:t, meanP[1:t], type = "l")
    }
    
    
    
  } #end time series
  
  
  
  #Output
  
  out[1:18] <- c(h2m, 
                 cv,
                 vm,
                 oldP,
                 mean(K.isl),
                 Ni,
                 Nrecol,
                 mean(Precol),
                 w2,
                 t, 
                 mean(sgrad1[1:temp]),
                 mean(G[,1]),
                 var(G[,1]),
                 length(G[,1]),
                 Isl.Peak,
                 meanP[temp],
                 bp.max,
                 bp)
  names(out) <-c("h2","cv","vm","Ancestral","meanK","Ni","Nrecol","Precol","w2","time_adap","selgrad","meanG","varG","N_end","F_Peak", "meanP", "bp.max", "bp")
  return(list(output = out, P = meanP))
}

library(compiler)

AdaptSS <- cmpfun(AdaptSS)
run_generation <- cmpfun(run_generation)
