# Rescorla-Wagner and related models for reversal learning.
# Joe Barnby and Michael Moutoussis, Fall 2020 on.

# Simple RW model ---------------------------------------------------------

# some utilities useful for RW and related modelling:
#try(source("gen_ut.R"))
#try(source("/home/michael/Insync/michael.moutoussis@googlemail.com/Google Drive - Shared with me/ComputationalModels/ProbabilisticReasoning/gen_ut.R"))

reversal_RW <- function(par, data, detail = T) {

  # argument checks :
  par <- as.numeric(par);
  if (length(par) <3){ par[3] <- 0; }; # if no reset param provided, set it to 0
  if (length(par) <4){ par[4] <- 1e-6; }; # if no Pierce-Hall salience parameter, set it to almost 0
  if (length(par) <5){ par[5] <- 1-1e-6; }; # if no memory parameter provided, set it to almost 1.
  if (length(par) <6){ par[6] <- 1e-6; }; # if no lapse param (zet) provided, set it to almost 0.
  if (length(par) <7){ par[7] <- 0; }; # if no lapse param (dlrc) provided, set it to almost 0.
  # unchosen options.

  #Trial Number
  t_n <- 61

  #Parameters
  tau   <-  par[1] # temperature
  lrc1  <-  par[2] # learning rate
  reset <-  par[3] # reset (shrinks the final values at 30 (block1) toward neutral)
  sal   <-  par[4] # salience (Pierce-hall modification)
  mem   <-  par[5] # memory
  zet   <-  par[6] # lapse parameter
  dlrc  <-  par[7] # change in learning between blocks

  #Salience matrix to load salience parameter calculations
  salience <- matrix(
    NA,
    nrow = 61,
    ncol = 4) #empty matrix for salience loop
  colnames(salience) <- c("Q1", "Q2", "Q3", "Trial")

  salience[1, 1:4] <- c(rep(1, 3),0)
  salience[2:61,4] <- 1:60

  #Q matrix to store Q values
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
  l <- rep(NA, 60) # initialise log likelihood

  #Pe initialise
  pe <- rep(NA, 60) # initialise pe

  #PeSim initialise
  peSim <- rep(NA, 60) # initialise peSim

  #win switching vector
  ws <- rep(NA, 61)

  #lose stay vector
  ls <- rep(NA, 61)

  #reward vector
  r    <- rep(NA, 61)
  r[1] <- 0 #initialise reward

  #action vector
  action <- rep(NA, 61)
  action[1]<- 0 #initialise action

  #Learning vector
  learning <- rep(NA, 61)

  if (detail == T){   # make space for the gory detail variables
    #Simulated behavioural data initialised parameters
    simD <- l
    simP <- l
    simQ <- q
    simR <- l;     # will hold simulated outcome
    simS <- salience
    simL <- learning; simL[1] <- 0;
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
  }

  # Main loop over trials, for learning,
  # reversal, and likelihood estimates
  for (t in 2:t_n) {

      #Block to reset priors at block 2 (trial 31:60)
      if (t == 32) {
        q[t-1, 1:3] <- q[t-1, 1:3] + (reset * (2.5 - q[t-1, 1:3])) #reset priors at block 2
      }

      lrc = lrc1

      if (t >= 32) {
        lrc = lrc1 + dlrc
      }

      ## Learning block with RW equation ##

      #Prediction error calculation
      action[t]   <- data[t-1, 2] #store which action was taken
      r[t]        <- data[t-1, 3] #store the reward achieved
      pe[t]       <- r[t] - q[t-1, action[t]] #calculate the pe based on prior expectation

      #Salience and learning estimates
      salience[t, action[t]]          <- (sal * abs(pe[t])) + ((1-sal) * salience[t-1,action[t]]) #Salience modifier of card
      salience[t, c(-action[t], -4)]  <- salience[t-1, c(-action[t], -4)] #update other rows with t-1 salience
      learning[t]                     <- lrc * salience[t,action[t]] #learning rate based on lrc and salience modifier

      if(learning[t] > 1) {learning[t] <- 0.99999}
      if(learning[t] < 0) {learning[t] <- 0.00001}

      #Q updating
      q[t,1:3]                        <- q[t-1, 1:3] #update trial t of Q with priors
      q[t,action[t]]                  <- q[t,action[t]] + (learning[t] * pe[t]) #update Q with learning plus salience
      q[t, c(-action[t], -4)]         <- 2.5 - (mem * (2.5 - q[t-1, c(-action[t], -4)])) #decay non chosen option by value of mem

      #likelihood equation
      Pmotiv          <- pGibbs(q[t-1,1:3],tau,action[t]); # (exp(q[t-1, action]/tau))/sum(exp(q[t-1,1:3]/tau)) #softmax function
      Pr              <- (zet/3) + ((1-zet) * Pmotiv)
      l[t-1]          <- log(Pr)

      #calculate win-switch behaviour

      if (action[t] != action[t-1] & r[t-1] == 10){
        ws[t] <- 1
      } else {
        ws[t] <- 0
      }

      #calculate lose-stay behaviour
      if (action[t] == action[t-1] & r[t-1] == -5){
        ls[t] <- 1
      } else {
        ls[t] <- 0
      }

      #Generate simulated data in parralell
      if (detail == T) {
        if (t == 32) {
          simQ[t-1, 1:3]  <- simQ[t-1, 1:3] + (reset * (2.5 - simQ[t-1, 1:3])) #reset priors at block 2
        }

        simPmotiv         <- pGibbs(simQ[t-1,1:3],tau); #(exp(simQ[t-1, 1:3]/tau))/sum(exp(simQ[t-1,1:3]/tau))
        Pi                <- (zet/3) + ((1-zet) * simPmotiv)
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

    }

  #Save outcomes of learning loops into an output
  if (detail == F) {

      return(sum(l))

    } else {

      output <- matrix(
        NA,
        nrow = 61,
        ncol = 3 + 2 + 3 + 2 + 2 + 6
      )

      colnames(output) <- c("Trial", "Action", "Outcome",
                            "PE", "Learning",
                            "Q1", "Q2", "Q3",
                            "ll", "PEsim",
                            "ws", "ls",
                            "simQ1","simQ2","simQ3","simD", "simR", "simP")

      output[2:61, 1:3] <- as.matrix(data[,1:3])
      output[, 4]       <- pe
      output[2:61, 5]   <- learning[2:61]
      output[, 6:8]     <- q[,1:3] #Q learning
      output[2:61, 9]   <- l       #loglikelihood
      output[, 10]  <- peSim
      output[, 11]      <- ws      #win switch
      output[, 12]      <- ls      #lose-stay
      output[, 13:18]   <- cbind(simQ[,1:3],simD, simR, simP) #simulated outcomes

      return(output)

    }

} # end of function

#  WRAPPERS TO FACILITATE USE OF REGULARIZATION ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.a Wrapper for Inference Model A Log-likelihood ~~~~~~~~~~~~~~~~~~~

# We'll try to use beta distribution with scaled x values, so that instead of
# bet. 0 and 1, x obtains values between lo and hi, which is the interval of
# psychobiological interest. Can have lo > hi, in which case they are mirrored.
# dbetasc <- function(x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE)
# e.g. for uncertaint param, may use:
# x <- 0.02*(0:1000);  plot(x, dbetasc(x,1.2,3.6,0,25),t='l')
# for a pretty uniform probability or learning rate prior:
# x <- 0.002*(0:1000)-0.5;  y <-  dbetasc(x,1.05,1.05,0,1); plot(x,y,t='l')
# for a normal looking prior with mean 10 and sd about 108:
# x <- (0:1000)-500;  y <-  dbetasc(x,10.41,10,-500,500); plot(x,y,t='l')
# So each col. of ParM has the form ashape, bshape, lowlim, hilim.,
# This example has weakly informative priors for all:
#           tau   lrc  reset  sal    mem    zeta
# ashape   0.5    0.1  0.75   0.05   0.95   0.0001
# bshape   0.1   0.05  0.1    0.01   0.05   0.0001
# lowlim    0       0    0      0     0     0
# hilim     1       1    1      1     1     1

wrapper_RW <- function(par, data, scbeta0=-1,check=0){
  # return the NEGATIVE- or MINUS- log-posterior given either the
  # prior provided (or a default one, if scbeta0 < 0)
  ### par: parameter vector
  ### scbeta0 = scaled beta distribution. If < 0 it acts as a
  ###           flag to use 'default weak regularizing priors'.
  ###           REM the density dbetasc has arguments (x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE)
  ### check: flag to do somewhat time-consuming checking.

  parM <- as.numeric(par); # in case it's inputed in another format
  parN <- length(parM)
  maxp <- c(20,rep(1,(parN-1)))  # upper boundaries for params
  minp <- rep(0,parN)            # lower boundaries for params

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    #  For:     tau lrc reset sal mem zet dlrc
    scale1 <- c(1.1,1.5,1.1,1.25,1.1,1.05, 1.5);
    scale2 <- c(10 ,2.5,1.1,1.25,1.1,2.95, 2.5);
    mp     <- rep(0,parN);
    Mp     <- c(20,rep(1,parN-1))
    scbeta0 <- t(matrix(c(scale1[1:parN],scale2[1:parN],mp,Mp),nrow = parN, ncol = 4))
    if(check){
      colnames(scbeta0) <- c('tau','lrc','reset','sal','mem','zeta', 'dlrc')
      rownames(scbeta0) <- c('scale1','scale2','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  if (sum(par < minp) || sum(par > maxp)) {  # i.e. if at least one param out of range
    mRWPrior <- Inf
  } else {
    mRWPrior <- 0;
    if (length(scbeta0)>1){  # legit prior must have 3*parN elements or so!
      # dbetasc follows vectorization of dbeta, so can do all the elements of parM
      # at once, and als get the log-prior straight away:
      mRWPrior <-  - sum(dbetasc( parM,
                                scbeta0[1,1:parN],scbeta0[2,1:parN],
                                scbeta0[3,1:parN],scbeta0[4,1:parN], log=TRUE));
    }
  }

  if (mRWPrior == Inf){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculate the likelihood - it will be nonsense anyway.
    return(Inf);
  } else {
    return(mRWPrior - reversal_RW(par,data, detail = F))
  }


} # end of wrapper_RW

wrapper_WSLS <- function(par, data, scbeta0=-1,check=0){
  # return the NEGATIVE- or MINUS- log-posterior given either the
  # prior provided (or a default one, if scbeta0 < 0)
  ### par: parameter vector
  ### scbeta0 = scaled beta distribution. If < 0 it acts as a
  ###           flag to use 'default weak regularizing priors'.
  ###           REM the density dbetasc has arguments (x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE)
  ### check: flag to do somewhat time-consuming checking.

  parM <- as.numeric(par); # in case it's inputed in another format
  parM <- c(0.001, 0.999, parM) # add fixed parameters
  parN <- length(parM)-2
  maxp <- rep(1,parN)  # upper boundaries for params  MAY NEED ADJUSTMENT
  minp <- rep(0,parN)            # lower boundaries for params  MAY NEED ADJUSTMENT

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    #  For:     tau lrc reset sal
    scale1 <- c(1.1,1.25);
    scale2 <- c(1.1,1.25);
    mp     <- rep(0,2);
    Mp     <- rep(1,2);
    scbeta0 <- t(matrix(c(scale1,scale2,mp,Mp),2,4))
    if(check){
      colnames(scbeta0) <- c('reset','sal')
      rownames(scbeta0) <- c('scale1','scale2','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  if (sum(par < minp) || sum(par > maxp)) {  # i.e. if at least one param out of range
    mRWPrior <- Inf
  } else {
    mRWPrior <- 0;
    if (length(scbeta0)>1){  # legit prior must have 3*parN elements or so!
      # dbetasc follows vectorization of dbeta, so can do all the elements of parM
      # at once, and als get the log-prior straight away:
      mRWPrior <-  - sum(dbetasc( parM[3:4],
                                  scbeta0[1,1:parN],scbeta0[2,1:parN],
                                  scbeta0[3,1:parN],scbeta0[4,1:parN], log=TRUE));
    }
  }

  if (mRWPrior == Inf){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculate the likelihood - it will be nonsense anyway.
    return(Inf);
  } else {
    return(mRWPrior - reversal_RW(parM,data, detail = F))
  }


} # end of wrapper_RW

wrapper_PH <- function(par, data, scbeta0=-1,check=0){
  # return the NEGATIVE- or MINUS- log-posterior given either the
  # prior provided (or a default one, if scbeta0 < 0)
  ### par: parameter vector
  ### scbeta0 = scaled beta distribution. If < 0 it acts as a
  ###           flag to use 'default weak regularizing priors'.
  ###           REM the density dbetasc has arguments (x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE)
  ### check: flag to do somewhat time-consuming checking.

  parM <- as.numeric(par); # in case it's inputed in another format
  parM <- c(parM[1], parM[2], 1e-10, parM[3], 1-1e-6, 1e-6, 1e-10) # add fixed parameters
  parN <- length(parM)-4
  maxp <- c(20, 1, 1)  # upper boundaries for params  MAY NEED ADJUSTMENT
  minp <- rep(0,parN)            # lower boundaries for params  MAY NEED ADJUSTMENT

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    #  For:     tau lrc reset sal mem zet
    scale1 <- c(1.1,1.25,1.05);
    scale2 <- c(1.1,1.25,2.95);
    mp     <- rep(0,3);
    Mp     <- c(20, 1, 1);
    scbeta0 <- t(matrix(c(scale1,scale2,mp,Mp),3,4))
    if(check){
      colnames(scbeta0) <- c('tau', 'lrc', 'sal')
      rownames(scbeta0) <- c('scale1','scale2','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  if (sum(par < minp) || sum(par > maxp)) {  # i.e. if at least one param out of range
    mRWPrior <- Inf
  } else {
    mRWPrior <- 0;
    if (length(scbeta0)>1){  # legit prior must have 3*parN elements or so!
      # dbetasc follows vectorization of dbeta, so can do all the elements of parM
      # at once, and als get the log-prior straight away:
      mRWPrior <-  - sum(dbetasc( parM[3:6],
                                  scbeta0[1,1:parN],scbeta0[2,1:parN],
                                  scbeta0[3,1:parN],scbeta0[4,1:parN], log=TRUE));
    }
  }

  if (mRWPrior == Inf){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculate the likelihood - it will be nonsense anyway.
    return(Inf);
  } else {
    return(mRWPrior - reversal_RW(parM,data, detail = F))
  }


} # end of wrapper_RW

# Block 1 model ####

reversal_RW_block1 <- function(par, data, detail = T) {

  # argument checks :
  par <- as.numeric(par);
  if (length(par) <3){ par[3] <- 1e-6; }; # if no Pierce-Hall salience parameter, set it to almost 0
  if (length(par) <4){ par[4] <- 1-1e-6; }; # if no memory parameter provided, set it to almost 1.
  if (length(par) <5){ par[5] <- 1e-6; }; # if no lapse param (zet) provided, set it to almost 0.
  # unchosen options.

  #Trial Number
  t_n <- 31

  #Parameters
  tau   <-  par[1] # temperature
  lrc   <-  par[2] # learning rate
  sal   <-  par[3] # salience (Pierce-hall modification)
  mem   <-  par[4] # memory
  zet   <-  par[5] # lapse parameter

  #Salience matrix to load salience parameter calculations
  salience <- matrix(
    NA,
    nrow = 31,
    ncol = 4) #empty matrix for salience loop
  colnames(salience) <- c("Q1", "Q2", "Q3", "Trial")

  salience[1, 1:4] <- c(rep(1, 3),0)
  salience[2:31,4] <- 1:30

  #Q matrix to store Q values
  q <- matrix(
    rep(
      rep(0, 31),
      4),
    nrow = 31,
    ncol = 4) #empty matrix for loop
  colnames(q) <- c("Q1", "Q2", "Q3", "Trial")

  q[1,] <- c(2.5, 2.5, 2.5, 0) #set priors based on average of 10 and -5
  q[,4] <- 0:(t_n-1)

  #Loglikelihood initialise
  l <- rep(NA, 30) # initialise log likelihood

  #Pe initialise
  pe <- rep(NA, 30) # initialise pe

  #PeSim initialise
  peSim <- rep(NA, 30) # initialise peSim

  #win switching vector
  ws <- rep(NA, 31)

  #lose stay vector
  ls <- rep(NA, 31)

  #reward vector
  r    <- rep(NA, 31)
  r[1] <- 0 #initialise reward

  #action vector
  action <- rep(NA, 31)
  action[1]<- 0 #initialise action

  #Learning vector
  learning <- rep(NA, 31)

  if (detail == T){   # make space for the gory detail variables
    #Simulated behavioural data initialised parameters
    simD <- l
    simP <- l
    simQ <- q
    simR <- l;     # will hold simulated outcome
    simS <- salience
    simL <- learning
    retp <-c(
        0.8,0.5,0.2
      )   # prob of win per action in block 1
  }

  # Main loop over trials, for learning,
  # reversal, and likelihood estimates
  for (t in 2:t_n) {

    ## Learning block with RW equation ##

    #Prediction error calculation
    action[t]   <- data[t-1, 2] #store which action was taken
    r[t]        <- data[t-1, 3] #store the reward achieved
    pe[t]       <- r[t] - q[t-1, action[t]] #calculate the pe based on prior expectation

    #Salience and learning estimates
    salience[t, action[t]]          <- (sal * abs(pe[t])) + ((1-sal) * salience[t-1,action[t]]) #Salience modifier of card
    salience[t, c(-action[t], -4)]  <- salience[t-1, c(-action[t], -4)] #update other rows with t-1 salience
    learning[t]                     <- lrc * salience[t,action[t]] #learning rate based on lrc and salience modifier

    if(learning[t] > 1) {
      learning[t] <- 0.99999
    }
    if(learning[t] < 0) {
      learning[t] <- 0.00001
    }

    #Q updating
    q[t,1:3]                        <- q[t-1, 1:3] #update trial t of Q with priors
    q[t,action[t]]                  <- q[t,action[t]] + (learning[t] * pe[t]) #update Q with learning plus salience
    q[t, c(-action[t], -4)]         <- 2.5 - (mem * (2.5 - q[t-1, c(-action[t], -4)])) #decay non chosen option by value of mem

    #likelihood equation
    Pmotiv          <- pGibbs(q[t-1,1:3],tau,action[t]); # (exp(q[t-1, action]/tau))/sum(exp(q[t-1,1:3]/tau)) #softmax function
    Pr              <- (zet/3) + ((1-zet) * Pmotiv)
    l[t-1]          <- log(Pr)

    #calculate win-switch behaviour

    if (action[t] != action[t-1] & r[t-1] == 10){
      ws[t] <- 1
    } else {
      ws[t] <- 0
    }

    #calculate lose-stay behaviour
    if (action[t] == action[t-1] & r[t-1] == -5){
      ls[t] <- 1
    } else {
      ls[t] <- 0
    }

    #Generate simulated data in parralell
    if (detail == T) {

      simPmotiv         <- pGibbs(simQ[t-1,1:3],tau); #(exp(simQ[t-1, 1:3]/tau))/sum(exp(simQ[t-1,1:3]/tau))
      Pi                <- (zet/3) + ((1-zet) * simPmotiv)
      simA              <- sample(c(1, 2, 3), 1, prob = Pi)  # simulated action (to avoid unreferencing)
      simD[t]           <- simA                              # ... store it as decision
      simP[t]           <- Pi[simA]                          # to compare w. experimental in due


      simR[t] <- sample(c(10,-5),1,prob=c(retp[simA], 1-retp[simA]))  # simulated outcome


      peSim[t] <- simR[t] - simQ[t-1,simA]

      simS[t, simA]          <- (sal * abs(peSim[t])) + ((1-sal) * simS[t-1, simA]) #Salience modifier of card
      simS[t, c(-simA, -4)]  <- simS[t-1, c(-simA, -4)] #update other rows with t-1 salience
      simL[t]                <- lrc * simS[t,simA] #learning rate based on lrc and salience modifier

      if(simL[t] > 1) {
        simL[t] <- 0.99999
      }
      if(simL[t] < 0) {
        simL[t] <- 0.00001
      }

      simQ[t,1:3]              <- simQ[t-1, 1:3] #set priors for Q using t-1
      simQ[t,simA]             <- simQ[t,simA] + simL[t] * peSim[t] #update Q at trial t
      simQ[t, c(-simA, -4)]    <- 2.5 - (mem * (2.5 - simQ[t-1, c(-simA, -4)])) #decay non chosen option by value of mem

    }

  }

  #Save outcomes of learning loops into an output
  if (detail == F) {

    return(sum(l))

  } else {

    output <- matrix(
      NA,
      nrow = 31,
      ncol = 3 + 2 + 3 + 2 + 2 + 6
    )

    colnames(output) <- c("Trial", "Action", "Outcome",
                          "PE", "Learning",
                          "Q1", "Q2", "Q3",
                          "ll", "PEsim",
                          "ws", "ls",
                          "simQ1","simQ2","simQ3","simD", "simR", "simP")

    output[2:31, 1:3] <- as.matrix(data[,1:3])
    output[, 4]       <- pe
    output[2:31, 5]   <- learning[2:31]
    output[, 6:8]     <- q[,1:3] #Q learning
    output[2:31, 9]   <- l       #loglikelihood
    output[, 10]  <- peSim
    output[, 11]      <- ws      #win switch
    output[, 12]      <- ls      #lose-stay
    output[, 13:18]   <- cbind(simQ[,1:3],simD, simR, simP) #simulated outcomes

    return(output)

  }

} # end of function

wrapper_RW_30 <- function(par, data, scbeta0=-1,check=0){
  # return the NEGATIVE- or MINUS- log-posterior given either the
  # prior provided (or a default one, if scbeta0 < 0)
  ### par: parameter vector
  ### scbeta0 = scaled beta distribution. If < 0 it acts as a
  ###           flag to use 'default weak regularizing priors'.
  ###           REM the density dbetasc has arguments (x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE)
  ### check: flag to do somewhat time-consuming checking.

  parM <- as.numeric(par); # in case it's inputed in another format
  parN <- length(parM)
  maxp <- c(20,rep(1,(parN-1)))  # upper boundaries for params  MAY NEED ADJUSTMENT
  minp <- rep(0,parN)            # lower boundaries for params  MAY NEED ADJUSTMENT

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    #  For:     tau lrc reset sal mem zet dlrc
    scale1 <- c(1.1,1.5,1.25,1.1,1.05);
    scale2 <- c(10 ,2.5,1.25,1.1,2.95);
    mp     <- rep(0,5);
    Mp     <- c(20,rep(1,4))
    scbeta0 <- t(matrix(c(scale1,scale2,mp,Mp),5,4))
    if(check){
      colnames(scbeta0) <- c('tau','lrc','sal','mem','zeta')
      rownames(scbeta0) <- c('scale1','scale2','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  if (sum(par < minp) || sum(par > maxp)) {  # i.e. if at least one param out of range
    mRWPrior <- Inf
  } else {
    mRWPrior <- 0;
    if (length(scbeta0)>1){  # legit prior must have 3*parN elements or so!
      # dbetasc follows vectorization of dbeta, so can do all the elements of parM
      # at once, and als get the log-prior straight away:
      mRWPrior <-  - sum(dbetasc( parM,
                                  scbeta0[1,1:parN],scbeta0[2,1:parN],
                                  scbeta0[3,1:parN],scbeta0[4,1:parN], log=TRUE));
    }
  }

  if (mRWPrior == Inf){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculate the likelihood - it will be nonsense anyway.
    return(Inf);
  } else {
    return(mRWPrior - reversal_RW_block1(par,data, detail = F))
  }


} # end of wrapper_RW

#Utility
#
pGibbs <- function(q,T=1.0,ind=NA) {
# Gibbs softmax probabilities for vectors q and temperatures T.
# Expect q to be a matrix with columns the different actions, all
# subject to the same T.
# If ind is NA, return the whole pGibbs matrix;
# otherwise for cols. in ind.

# old comment - for future devel:
# Should use the result to create a number of decisions - in 1-D case
# policy <- matrix( pGibbs(c(-1,0,1)),1)
# choices <- as.vector( rMultinom( policy, 100));

# test with
# T <- c(0.1,0.3,1,2); q <- t(matrix(rep(c(0.1,0.5,0.6),4),3)); ind <- c(1,1,2,3);

if (is.vector(q)){  q <- matrix(q,1); }; # convert to 1-row matrix
nro = dim(q)[1]
nco = dim(q)[2]
if (!(is.vector(T))){ stop('T must be a vector in pGibbs'); };
if (!(length(T)==nro)){ stop('q and T incompatible');  };

T <- matrix(rep(T,nco),length(T)) ; # convert to matrix for ease
                                    # of vect. op.
# Express values relative to the biggest one to avoid
# exponentiation overflows. Underflows don't matter so much.
qmax   <- apply(q,1,max);
qmax   <- matrix(rep(qmax,nco),nro,nco);
unNorm <- exp((q-qmax)/T);
# Normalizing factor:
denoms <- rowSums(unNorm);
denoms <- matrix(rep(denoms,dim(q)[2]),dim(q)[1]); # again into matrix
# final probability:
pGibbs <- unNorm / denoms;

if (is.na(ind[1])) {
    return( pGibbs );
}
else {
    if (!(length(ind) == dim(q)[1])){stop('q and ind incompatible');}
    else{
      ind <- matrix(rep(ind,dim(q)[2]),length(ind));
      # Next row a bit ridiculous - can't I just select the round(ind) elements??
      return(  rowSums(  pGibbs*(col(ind) == round(ind)))   );
    }
}

}
