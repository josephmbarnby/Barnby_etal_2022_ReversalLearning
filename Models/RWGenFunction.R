  tau   <-  1.3 # temperature
  lrc   <-  0.4 # learning rate
  reset <-  0 # reset (shrinks the final values at 30 (block1) toward neutral)
  sal   <-  0.005 # salience (Pierce-hall modification)
  mem   <-  0.999 # memory

  test_run <- list()
  for (i in 1:1000){
  test_run_dat<- gen_ppt_RW_win(plot = 0) # run function
  test_run[[i]] <- test_run_dat[[1]]
  }

  test_run_plot <- do.call(rbind, test_run)

  ggplot(test_run_plot %>% as.data.frame() %>% filter(Trial > 0))+
    stat_summary(geom = 'line', aes(Trial, Q1), color = 'red')+
    stat_summary(geom = 'line', aes(Trial, Q2), color = 'black')+
    stat_summary(geom = 'line', aes(Trial, Q3), color = 'blue')


gen_ppt_RW_win <- function(plot = 1){

    t_n = 61
  q <- matrix(
    rep(
      rep(0, 61),
      4),
    nrow = 61,
    ncol = 4) #empty matrix for loop
  colnames(q) <- c("Q1", "Q2", "Q3", "Trial")
  q[1,] <- c(2.5, 2.5, 2.5, 0) #set priors based on average of 10 and -5
  q[,4] <- 0:(t_n-1)

  #Loglikelihood initialise
  l <- rep(NA, 60) # initialise log likelihood#Simulated behavioural data initialised parameters
    simD <- l
    simP <- l
    simQ <- q
    simR <- l;     # will hold simulated outcome
    simS <- matrix(
    NA,
    nrow = 61,
    ncol = 4) #empty matrix for salience loop
    colnames(simS) <- c("Q1", "Q2", "Q3", "Trial")
    simL <- rep(NA, 61)
    simS[1, 1:4] <- c(rep(1, 3),0)
    simS[2:61,4] <- 1:60
    retp1 <- matrix(
      c(
        0.8,0.5,0.2,0.2,0.5,0.8
        ),
      3,2)   # prob of win per action in block 1
    retp2 <- matrix(
      c(
        0.5,0.2,0.8,0.5,0.8,0.2
        ),
      3,2)   # prob of win per action in block 1
    peSim <- rep(NA, 60) # initialise peSim

  for(t in 2:61) {
          if (t == 32) {
            simQ[t-1, 1:3]  <- simQ[t-1, 1:3] + (reset * (2.5 - simQ[t-1, 1:3])) #reset priors at block 2
          }

          Pi                <- pGibbs(simQ[t-1,1:3],tau); #(exp(simQ[t-1, 1:3]/tau))/sum(exp(simQ[t-1,1:3]/tau))
          #Pi                <- (zet/3) + ((1-zet) * simPmotiv)
          simA              <- sample(c(1, 2, 3), 1, prob = Pi)  # simulated action (to avoid unreferencing)
          simD[t]           <- simA                              # ... store it as decision
          simP[t]           <- Pi[simA]                          # to compare w. experimental in due

          if (t < 32){
             simR[t] <- sample(c(10,-5),1,prob=retp1[simA,])  # simulated outcome
          } else {  # simple but a bit esoteric reversing, relies on symmetry of two blocks
             simR[t] <- sample(c(10,-5),1,prob=retp2[simA,])  # simulated outcome
          }

          peSim[t] <- simR[t] - simQ[t-1,simA]

          simS[t, simA]          <- (sal * abs(peSim[t])) + ((1-sal) * simS[t-1, simA]) #Salience modifier of card
          simS[t, c(-simA, -4)]  <- simS[t-1, c(-simA, -4)] #update other rows with t-1 salience
          simL[t]                <- lrc * simS[t,simA] #learning rate based on lrc and salience modifier

          if(simL[t] > 1) {simL[t] <- 0.99999}
          if(simL[t] < 0) {simL[t] <- 0.00001}

          simQ[t,1:3]              <- simQ[t-1, 1:3] #set priors for Q using t-1
          simQ[t,simA]             <- simQ[t,simA] + simL[t] * peSim[t] #update Q at trial t
          simQ[t, c(-simA, -4)]    <- 2.5 - (mem * (2.5 - simQ[t-1, c(-simA, -4)])) #decay non chosen option by value of mem

  }

  data <- list()
  data[[1]] <- simQ
  if(plot){
  data[[2]] <- ggplot(simQ %>% as.data.frame() %>% filter(Trial > 0))+
    geom_line(aes(Trial, Q1), color = 'red')+
    geom_line(aes(Trial, Q2), color = 'black')+
    geom_line(aes(Trial, Q3), color = 'blue')+
    geom_point(aes(1:60, simD[2:61]))+
    geom_point(aes(1:60, simR[2:61]),
               color = ifelse(simR[2:61] == 10, 'red', 'grey'))
  }
    return(data)
}
