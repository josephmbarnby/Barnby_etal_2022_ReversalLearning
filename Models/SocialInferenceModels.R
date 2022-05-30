#######  MODELS ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(R.matlab)
library(tidyverse)
library(tidyquant)

### Dictator_Reversal_Functions.R ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.d.Parametric dependence of the attribute-policy map, with only
#     two-level policies like the actual experiment.
#     separate-rigidity (double eta) Inference Model A Log-likelihood ~~
#     Only 2 possible returns, .5 and 0, as per pt instructions
#     May 2021 on

# Define likelihood function for inference model. This gives
# the probability density of the responses given the participant's
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 1 partners x 20 trials each, in order of presentation, down the rows.

#Newer model with variable policy map
infHISIll_20d <- function(parms,d,details=0,plots=0,sim_only=1,tn = 10,phase = 2) {

  d <- as.matrix(d)

  tn    = tn;  # 10 trials
  phase = phase;  # one partner, two phases

  # Parameters to help describe personality and action. These are
  # not free parameters to be fitted, i.e. don't vary across pts.
  Nb = 9;     # number of bins for each personality dim.
  Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
  Na = 2;     # actual number of actions in Joe's expt.

  # i.e. proportions that may be received by the Other.
  # In general will be set to Nb.

  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(parms[3], parms[4],Nb); PHI0 = noisyBino(parms[1],parms[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0);

  if (plots){
    opar=par()
    anames = c('selfish','fair')
  }

  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = parms[5];     # convenience copy of policy variability (inverse precision)
  w0  = parms[6];     # This and two next will help map HI, SI onto policy (zero or half return)
  whi = parms[7];
  wsi = parms[8];

  if ((whi < 0) || (wsi < 0)){
    print('Params c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI)) :')
    print(parms)
    error("whi,wsi must all be non-negative")
  }

  if (length(parms) == 8){
  etaHI = 1
  etaSI = 1
  } else if (length(parms) == 9){
  eta = parms[9];
  } else if (length(parms) == 10){
  etaHI = parms[9];
  etaSI = parms[10]
  }

  err =0.02/(Nb*Nb) # arbitrary lapse-rate-like
  #}
  # param; Note scaling by the number of attribution states considered.

  # Set up the map between 'attributes' and actions :
  pi  = array(NA,c(Nb,Nb,Na));

  offs = (Nb+1)/2
  for (SI in 1:Nb){
    for (HI in 1:Nb){
      pi[SI,HI,1]  = invlogit(w0  + (wsi*(SI-offs)) + (whi*(HI-offs)))
      # the offs is to center at (usually) 5
      pi[SI,HI,2]  = 1 - pi[SI,HI,1]
    }
  }

  ### Run the inference

  # Here the posterior of one trial will form the prior for the next trial,
  # staring from the population prior beliefs PSIHI0.
  #

  sll = 0;
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)

  if (details ){
    llout = list();

    llout[[1]]=0;

    # [[2]] is the parameter vector
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','w0', 'wHI', 'wSI', 'etaHI','etaSI', 'err')
    llout[[2]] = c(parms[1:10],err);  names(llout[[2]]) <- hd;

    # [[3]] is the detailed evolution, evo:
    llout[[3]] = matrix(NA,tn*phase+1,10);
    llout[[3]][,1] = c(0,1:(tn*phase))
    colnames(llout[[3]]) = c('trial','ret','HI','SI','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*phase),2:4] <- as.matrix(d)

    # [[4]] will be a big 9 x 9 x 21 array with detailed policies
    # Hypothetical (attribn. reporting) policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*phase)))
    llout[[4]][,,1] <- pri0;
    llout[[5]] <- list();
    llout[[6]] <- list();
    names(llout) <- c('sll','par', 'evo','policy','attr2pol', 'plots')
    }

    post <- pri0; # this is the belief that will be updated with each trial

    # rows to be processed and convenience copies of other's actions,
    # malice and greed attributions:
    ro = 1:(phase*tn); # rows of data matrix
    as = d[ro,1];  aind = round((Na-1)*as+1)
    hi = d[ro,2];  hind = round((Nb-1)*hi+1)
    si = d[ro,3];  sind = round((Nb-1)*si+1)

    for (t in 1:(tn*phase)){  # loop

      #if(plots){
      #  heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
      #          margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
      #          main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
      #}


      if (t == (tn+1) & length(parms)==9){

        post = (pri0 * (1-eta)) + (post * eta);

      } else if (t == (tn+1) & length(parms)==10) {

        PSIPost <- rowSums(post);
        PHIPost <- colSums(post);
        post    <- ((PSI0 * (1-etaSI)) + (etaSI * PSIPost)) %*% t((PHI0 * (1-etaHI)) + (etaHI * PHIPost))        ;

      }

      pri = post;              # new prior is last posterior
      # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
      # over the character of the partner. This post is their post. beliefs
      post = pi[,,aind[t]] * pri
      post = post / sum(as.vector(post))  # Bayes

      # Now the probability of the response, incl. the lapse-rate-like err:
      pol = post^(1/upi);
      pol = pol/sum(as.vector(pol)); #renormalise

      pol = (pol+err)/(1+err*length(pol))
      lik = pol[sind[t],hind[t]];
      sll = sll + log(lik);         # accumulate sum log lik

      if (details ){

        llout$evo[(t+1),'lik'] <- lik
        llout$evo[(t+1),'ll'] <- log(lik)
        # find mode of pol
        c = max.col(pol);  # this finds the max. col. pos. for each row
        m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
        r=which.max(m);  # ... so was to find the row with the mode of the distro.
        llout$evo[(t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
        # Now sample randomly
        hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
        sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the phase from the corresponding conditional.
        # debug: barplot(colSums(pol));       barplot(pol[,hisim]);
        llout$evo[(t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
        llout$policy[,,(t+1)] <- pol

      }
    }

    if(plots){
      library(tidyquant)
      library(ggplot2)
      library(patchwork)
      initial_policy <- pi[,,1] %>%
        as.data.frame() %>%
        mutate(SI = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85)) %>%
        pivot_longer(1:9, 'HI', values_to = 'Probability') %>%
        ggplot(aes(HI, SI, fill = Probability)) +
        geom_tile() +
        scale_fill_gradient(low = '#98C1D9', high = '#C81D25')+
        scale_x_discrete(labels = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85))+
        theme_tq()
      plotobject <- llout$evo %>%
        as.data.frame() %>%
        na.omit()
      if(sim_only == 1){
      trial_wise <- plotobject  %>%
        pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
        mutate(ret = ifelse(ret == 0.5, 1, 0)) %>%
        ggplot(aes(trial, values, color = Attribution))+
          geom_line()+
          geom_point(aes(trial, ret), color = ifelse(ret == 1, 'black', 'dark red'))+
          scale_color_brewer(palette = 'Set1')+
          labs(x = 'Trial', y = 'Simulated Observation')+
          theme_tq()
      }else if (sim_only == 0){
      trial_wise <- plotobject  %>%
        pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
        pivot_longer(HI:SI, 'AttributionREAL', values_to = 'valuesREAL') %>%
        mutate(ret = ifelse(ret == 0.5, 1, 0)) %>%
        ggplot(aes(trial, values, color = Attribution))+
          geom_line(linetype = 1)+
          geom_line(aes(trial, valuesREAL, color = AttributionREAL), linetype = 2)+
          geom_point(aes(trial, ret), color = ifelse(ret == 1, 'black', 'dark red'))+
          scale_color_brewer(palette = 'Set1')+
          labs(x = 'Trial', y = 'Simulated Observation')+
          theme_tq()
      }
      fitness <- plotobject %>%
        ggplot(aes(ll))+
        geom_density()+
        geom_histogram(alpha = 0.5, binwidth = 0.2)+
        geom_vline(xintercept = -4.17)+
        labs(x = 'Trial-wise LL', y = 'Density')+
        theme_tq()
      plotter <- (trial_wise | (initial_policy/fitness)) + plot_annotation(tag_levels = 'A')
      llout[[6]] <- plotter
    }

    #if (plots) {
    #  heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
    #          margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
    #          main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
    #  par(opar)  # restore state of graphics where we found it.
    #}

  if (details ){llout$sll <- sll} else {llout <- sll}
  return(llout)

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.c separate-rigidity (double eta) Inference Model A Log-likelihood ~~
#     Only 2 possible returns, .5 and 0, as per pt instructions
#     May 2021 on

# Define likelihood function for inference model. This gives
# the probability density of the responses given the participant's
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 1 partners x 20 trials each, in order of presentation, down the rows.

#Older model (Barnby et al., 2020)
infHISIll_20b <- function(parms,d,details=0,plots=0) {

  d <- as.matrix(d)
  # rem parms is of the form c(pHI0,uHI0,pSI0,uSI0,upi,etaHI,etaSI)
  tn = 10;  # 10 trials
  phase = 2;  # one partner, two phases

  # Parameters to help describe personality and action. These are
  # not free parameters to be fitted, i.e. don't vary across pts.
  Nb = 9;     # number of bins for each personality dim.
  Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
  Na = 2;     # actual number of actions in Joe's expt.

  # i.e. proportions that may be received by the Other.
  # In general will be set to Nb.

  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(parms[3], parms[4],Nb); PHI0 = noisyBino(parms[1],parms[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0);

  if (plots){
    # Action names (was:  anames = c('selfish','v.tight','tight','fair','kind','v.kind','altruistic') )
    anames = c('selfish','fair')
  }

  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = parms[5];     # convenience copy of policy variability (inverse precision)

  if (length(parms) == 6){
    eta = parms[6];
  }else if (length(parms) == 7){
    etaHI = parms[6];
    etaSI = parms[7]
  }

  if (length(parms) > 7){piu=parms[8]} else {piu = 1} # arbitrary policy-inverse-uncertainty
  # param; could set to e.g. 1/upi to reference it to self, roughly.
  if (length(parms) > 8){err=parms[9]} else {err =0.02/(Nb*Nb)} # arbitrary lapse-rate-like
  # param; Note scaling by the number of attribution states considered.

  # Set up the map between 'attributes' and actions :
  pif = array(NA,c(Nb,Nb,Nf));    # This will hold all possible policies for full range of actions
  pi  = array(NA,c(Nb,Nb,Na));    # Possible policies for actual range of actions
  fp =  c(1,1,1,0,0,0,0,0,0);    # auxiliary vector to further bin down pi, here to a two-bin vector
                                 # Joe suggestions c(1,1,1,1,1,1,1,0,0),  c(1,1,0,0,0,0,0,0,0);
                                 # Original fp = c(1,1,1,0,0,0,0,0,0);
  f2p = t(matrix(c(fp,1-fp),Nb,Na));
  # To plot average policy / returns :
  #if (plots) {
    piav = array(NA,c(Nb,Nb)); # This will hold average policies
    pifav = piav;
  #}
  # pinit, pstep, bu, mu fine-tune the map ... see below.
  pinit = 0.1; pstep= (1-2*pinit)/(Nb-1);
  bu = 2.5; mu= -2*pstep;     # more handwavey constants to modulate u. Should have
  # bu+mu*Nb and bu+mu both valid values for u in the noisyBino. Negative values of
  # mu over-weigh high values of SI and HI, so that SI=1, HI=9 is skewed towards small
  # returns, while mu=0 would be symmetric around the middle.

  for (SI in 1:Nb){
    for (HI in 1:Nb) {
      x = noisyBino(pinit+(SI-1)*pstep, bu+mu*SI,Nf) *
        noisyBino(pinit+(HI-1)*pstep, bu+mu*HI,Nf);
      pif[SI,HI,] = fliplr(x^(1/piu) / sum(x^(1/piu)))
      pi[SI,HI,]  = as.vector(f2p %*% pif[SI,HI,]);      # further bin down!
      if (plots){
        piav[SI,HI] = sum(pi[SI,HI,]*c(0,0.5))
        pifav[SI,HI] = sum(pif[SI,HI,]*(1:Nb))
      }
    }
  }

  if (plots){           # Display the average policy / returns as a heatmap
    heatmap( piav,
             Rowv=NA, Colv=NA, col = topo.colors(512),
             scale="none", margins=c(5,8),asp=1,
             labRow=((0:(Nb-1))+0.5)*0.1,
             labCol=((0:(Nb-1))+0.5)*0.1,
             main = paste('\n attributes vs. mean policy'),
             xlab='HI',ylab='SI\n')
  }


  ### Run the inference

  # Here the posterior of one trial will form the prior for the next trial,
  # staring from the population prior beliefs PSIHI0.
  #

  sll = 0;
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)

  if (details ){
    llout = list();

    llout[[1]]=0;

    # [[2]] is the parameter vector
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','etaHI','etaSI','piu','err')
    llout[[2]] = c(parms[1:7],piu,err);  names(llout[[2]]) <- hd;

    # [[3]] is the detailed evolution, evao:
    llout[[3]] = matrix(NA,tn*phase+1,10);
    llout[[3]][,1] = c(0,1:(tn*phase))
    colnames(llout[[3]]) = c('trial','ret','HI','SI','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*phase),2:4] <- as.matrix(d)

    # [[4]] will be a big 9 x 9 x 21 array with detailed policies
    # Hypothetical (attribn. reporting) policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*phase)))
    llout[[4]][,,1] <- pri0;

    # [[5]] will be the HISI -> return policy (usef for the Other)
    llout[[5]] <- list();
    llout[[5]][[1]] <- pi; llout[[5]][[2]] <- pif;  llout[[5]][[3]] <- piav;  llout[[5]][[4]] <- pifav;
    names(llout[[5]]) <- c('piBinary','piFull','piBinAv','piFullAv')

    names(llout) <- c('sll','par', 'evo','policy','attr2pol')
  }

    post <- pri0; # this is the belief that will be updated with each trial

    # rows to be processed and convenience copies of other's actions,
    # malice and greed attributions:
    ro = 1:(phase*tn); # rows of data matrix
    as = d[ro,1];  aind = round((Na-1)*as+1)
    hi = d[ro,2];  hind = round((Nb-1)*hi+1)
    si = d[ro,3];  sind = round((Nb-1)*si+1)

    for (t in 1:(tn*phase)){  # loop

      if(plots){
        heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
                margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
                main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
      }

      if (t == 11 & length(parms) == 6){

        post <-  (pri0 * (1-eta)) + (post * eta)       ;

      } else if (t == 11 & length(parms) == 7) {

        PHIPost <- colSums(post);
        PSIPost <- rowSums(post);
        # Now do the exact equivalent of PSIHI0 = PSI0 %*% t(PHI0) above :
        post <-  ((PSI0 * (1-etaSI)) + (etaSI * PSIPost)) %*% t((PHI0 * (1-etaHI)) + (etaHI * PHIPost))        ;

      }

      pri = post;              # new prior is last posterior
      # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
      # over the character of the partner. This post is their post. beliefs
      post = pi[,,aind[t]] * pri; post = post / sum(as.vector(post))  # Bayes

      # Now the probability of the response, incl. the lapse-rate-like err:
      pol = post^(1/upi);  pol = pol/sum(as.vector(pol));
      pol = (pol+err)/(1+err*length(pol))
      lik = pol[sind[t],hind[t]];
      sll = sll + log(lik);         # accumulate sum log lik

      if (details ){

        llout$evo[t+1,'lik'] <- lik
        llout$evo[t+1,'ll'] <- log(lik)
        # find mode of pol
        c = max.col(pol);  # this finds the max. col. pos. for each row
        m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
        r=which.max(m);  # ... so was to find the row with the mode of the distro.
        llout$evo[(t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
        # Now sample randomly
        hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
        sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the phase from the corresponding conditional.
        # debug: barplot(colSums(pol));       barplot(pol[,hisim]);
        llout$evo[(t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
        llout$policy[,,(t+1)] <- pol

      }
    }

    if (plots) {
      heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
              margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
              main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
    }

  if (details ){llout$sll <- sll} else {llout <- sll}
  return(llout)

}

# 2.a Associative Model Log-likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     Single eta value.
#     Learning rate based on HI and SI separately
#     Only 2 possible returns, .5 and 0, as per pt instructions
#     May 2021 on

# Define likelihood function for associativ model. This gives
# the probability density of the responses given the participant's
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 1 partners x 20 trials each, in order of presentation, down the rows.

Assoc_model <- function(par, data, detail = F){


  data = as.data.frame(data)

  wHI0     <- par[1]   # harm intent report intercept
  wSI0     <- par[2]
  wHI      <- par[3]    #    "   "      "    slope w.r.t. HI expected social value
  wSI      <- par[4]
  alphaHI  <- par[5]    # lr for harm intent
  alphaSI  <- par[6]    # lr for selfishness intent
  SD       <- par[7]

  tn = 21

  ESVhi <- rep(NA, 21)
  ESVsi <- rep(NA, 21)
  ESVhi[1] <- par[8] # initial values of ESV for hi
  ESVsi[1] <- par[9] # initial values of ESV for si

  if(length(par) == 10){
  eta   <- par[10]
  }else if (length(par) == 11){
  etaHI   <- par[10]   # rigidity measure
  etaSI   <- par[11]
  }else{
  eta   <- 1
  }# rigidity measure

  SPEhi <- rep(NA, 21)
  SPEsi <- rep(NA, 21)
  HIsim <- rep(NA, 21)
  SIsim <- rep(NA, 21)
  sumll <- 0
  Bin   <- seq(0, 1, 0.1)
  Nb    <- length(Bin) - 2

  if(detail ){
    mat <- matrix(NA, tn, 13)
    colnames(mat) <- c('Trial', 'ret', 'SPEHI', 'SPESI', 'ESVHI', 'ESVSI', 'HI', 'SI', 'llHI', 'llSI', 'sumll', 'HIsim', 'SIsim')

    simD <- data; simD[,c(2,3)] <- NA;   # simulated data has same offers
  }

  for (t in 2:tn){

      ret   <- data[t-1,1]

      if (t == 12 & length(par) < 11){ # 'special' readjustment of values upon reversal at trial 11 i.e. t==12 :

        ESVhi[t-1] <- (eta*ESVhi[t-1] + (1-eta) * ESVhi[1])
        ESVsi[t-1] <- (eta*ESVsi[t-1] + (1-eta) * ESVsi[1])

      } else if (t == 12 & length(par) == 11){ # 'special' readjustment of values upon reversal at trial 11 i.e. t==12 :

        ESVhi[t-1] <- (etaHI*ESVhi[t-1] + (1-etaHI) * ESVhi[1])
        ESVsi[t-1] <- (etaSI*ESVsi[t-1] + (1-etaSI) * ESVsi[1])

      }

      SPEhi[t] <- ret - ESVhi[t-1]
      SPEsi[t] <- ret - ESVsi[t-1]

      # Associative update of expected social values for harm intent and selfishness intent.
      # In the case of reversal models, this should in due course include salience and resetting terms ...
      ESVhi[t] <- ESVhi[t-1] + alphaHI * SPEhi[t]
      ESVsi[t] <- ESVsi[t-1] + alphaSI * SPEsi[t]

      # Simplest linear map between ESV and expectation of reported value:
      expHI    <- wHI0 - (wHI * ESVhi[t])
      expSI    <- wSI0 - (wSI * ESVsi[t])

      if (expHI > 1){
        expHI = 0.999999
      }else if (expHI < 0) {
        expHI = 0.000001}

      if (expSI > 1){
        expSI = 0.999999
      }else if (expSI < 0) {
        expSI = 0.000001}

      # probability of actually observed data at this point:
      #prHI = dnorm( data[t-1,2], expHI , SD)
      #prSI = dnorm( data[t-1,3], expSI , SD)

      refHI = pnorm(c(Bin[2:10], 1), data[t-1,2], SD) - pnorm(c(0, Bin[2:10]), data[t-1,2], SD)
      refSI = pnorm(c(Bin[2:10], 1), data[t-1,3], SD) - pnorm(c(0, Bin[2:10]), data[t-1,3], SD)

      hind = round((Nb-1)*(expHI)+1);
      sind = round((Nb-1)*(expSI)+1);

      prHI = refHI[hind]
      prSI = refSI[sind]

      llHI = log(prHI)
      llSI = log(prSI)

      # if details asked for, generate some simulated data too:
      if (detail ){
        simD[t-1,2] <- rnorm(1, expHI , SD)
        simD[t-1,3] <- rnorm(1, expSI , SD)
        # clip off ridiculous values:
        if  (simD[t-1,2] < 0){ simD[t-1,2] <- 0 }
        if  (simD[t-1,2] > 1){ simD[t-1,2] <- 1 }
        if  (simD[t-1,3] < 0){ simD[t-1,3] <- 0 }
        if  (simD[t-1,3] > 1){ simD[t-1,3] <- 1 }
      }

      sumll = sumll + llHI + llSI

      if (detail ){
        mat[t,1]  <- t-1
        mat[t,2]  <- ret
        mat[t,3]  <- SPEhi[t]
        mat[t,4]  <- SPEsi[t]
        mat[t,5]  <- ESVhi[t]
        mat[t,6]  <- ESVsi[t]
        mat[t,7]  <- data[t-1, 2]
        mat[t,8]  <- data[t-1, 3]
        mat[t,9]  <- prHI
        mat[t,10] <- prSI
        mat[t,11] <- sumll
        mat[t,12] <- simD[t-1, 2]
        mat[t,13] <- simD[t-1, 3]

    }

  }   # main loop end

    if(detail ){
      res <- list();
      res[[1]] <- mat
      res[[2]] <- c(colSums(mat[2:tn,c(9:10)]), mat[tn,c(11)] )
      res[[3]] <- simD
      res[[4]] <- par
      names(res) <- c('evo','LLs','simD','par')
      return(res)
      } else {
    return(sumll)
    }
}

# 2.a Associative Model Log-likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     Single eta value.
#     Learning rate based on HI and SI separately
#     Only 2 possible returns, .5 and 0, as per pt instructions
#     May 2021 on

# Define likelihood function for associativ model. This gives
# the probability density of the responses given the participant's
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 1 partners x 20 trials each, in order of presentation, down the rows.

Assoc_model_singleESV <- function(par, data, detail = F){

  data = as.data.frame(data)

  wHI0     <- par[1]   # harm intent report intercept
  wSI0     <- par[2]
  wHI      <- par[3]    #    "   "      "    slope w.r.t. HI expected social value
  wSI      <- par[4]
  alpha    <- par[5]    # lr for harm intent
  SD       <- par[6]

  tn = 21

  ESV    <- rep(NA, 21)
  ESV[1] <- par[7] # initial values of ESV for hi

  if(length(par) == 8){
    eta   <- par[8]
  }else{
    eta   <- 1
  }# rigidity measure

  SPE   <- rep(NA, 21)
  HIsim <- rep(NA, 21)
  SIsim <- rep(NA, 21)
  sumll <- 0
  Bin   <- seq(0, 1, 0.1)
  Nb    <- length(Bin) - 2

  if(detail ){
    mat <- matrix(NA, tn, 11)
    colnames(mat) <- c('Trial', 'ret', 'SPE', 'ESV', 'HI', 'SI', 'llHI', 'llSI', 'sumll', 'HIsim', 'SIsim')

    simD <- data; simD[,c(2,3)] <- NA;   # simulated data has same offers
  }

  for (t in 2:tn){

    ret   <- data[t-1,1]

    if (t == 12){ # 'special' readjustment of values upon reversal at trial 11 i.e. t==12 :

      ESV[t-1] <- (eta*ESV[t-1]) + ((1-eta) * ESV[1])

    }

    SPE[t] <- ret - ESV[t-1]

    # Associative update of expected social values for harm intent and selfishness intent.
    # In the case of reversal models, this should in due course include salience and resetting terms ...
    ESV[t]   <- ESV[t-1] + alpha * SPE[t]

    # Simplest linear map between ESV and expectation of reported value:
    expHI    <- wHI0 - (wHI * ESV[t])
    expSI    <- wSI0 - (wSI * ESV[t])

    if (expHI > 0.999999){
      expHI = 0.999999
    }else if (expHI < 0.000001) {
      expHI = 0.000001}

    if (expSI > 0.999999){
      expSI = 0.999999
    }else if (expSI < 0.000001) {
      expSI = 0.000001}

    # probability of actually observed data at this point:
    #prHI = dnorm( data[t-1,2], expHI , SD)
    #prSI = dnorm( data[t-1,3], expSI , SD)

    refHI = pnorm(c(Bin[2:10], 1), data[t-1,2], SD) - pnorm(c(0, Bin[2:10]), data[t-1,2], SD)
    refSI = pnorm(c(Bin[2:10], 1), data[t-1,3], SD) - pnorm(c(0, Bin[2:10]), data[t-1,3], SD)

    hind = round((Nb*expHI)+1);
    sind = round((Nb*expSI)+1);

    prHI = refHI[hind]
    prSI = refSI[sind]

    llHI = log(prHI)
    llSI = log(prSI)

    # if details asked for, generate some simulated data too:
    if (detail ){
      simD[t-1,2] <- rnorm(1, expHI , SD)
      simD[t-1,3] <- rnorm(1, expSI , SD)
      # clip off ridiculous values:
      if  (simD[t-1,2] < 0){ simD[t-1,2] <- 0 }
      if  (simD[t-1,2] > 1){ simD[t-1,2] <- 1 }
      if  (simD[t-1,3] < 0){ simD[t-1,3] <- 0 }
      if  (simD[t-1,3] > 1){ simD[t-1,3] <- 1 }
    }

    sumll = sumll + llHI + llSI

    if (detail ){
      mat[t,1]  <- t-1
      mat[t,2]  <- ret
      mat[t,3]  <- SPE[t]
      mat[t,4]  <- ESV[t]
      mat[t,5]  <- data[t-1, 2]
      mat[t,6]  <- data[t-1, 3]
      mat[t,7]  <- prHI
      mat[t,8]  <- prSI
      mat[t,9]  <- sumll
      mat[t,10] <- simD[t-1, 2]
      mat[t,11] <- simD[t-1, 3]

    }

  }   # main loop end

  if(detail ){
    res <- list();
    res[[1]] <- mat
    res[[2]] <- c(colSums(mat[2:tn,c(7:8)]), mat[tn,c(9)] )
    res[[3]] <- simD
    res[[4]] <- par
    names(res) <- c('evo','LLs','simD','par')
    return(res)
  } else {
    return(sumll)
  }
}


# 2.c Associative Model Log-likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     Single Learnign rate.
#     Only 2 possible returns, .5 and 0, as per pt instructions
#     May 2021 on

# Define likelihood function for associativ model. This gives
# the probability density of the responses given the participant's
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 1 partners x 20 trials each, in order of presentation, down the rows.

Assoc_model_singleAlpha <- function(par, data, detail = F){

  data = as.data.frame(data)

  wHI0     <- par[1]   # harm intent report intercept
  wSI0     <- par[2]
  wHI      <- par[3]    #    "   "      "    slope w.r.t. HI expected social value
  wSI      <- par[4]
  alpha    <- par[5]    # lr for unfair outcome
  SD       <- par[6]

  tn = 21

  ESVhi <- rep(NA, 21)
  ESVsi <- rep(NA, 21)
  ESVhi[1] <- par[7] # initial values of ESV for hi
  ESVsi[1] <- par[8] # initial values of ESV for si

  eta   <- par[9]   # rigidity measure

  SPEhi <- rep(NA, 21)
  SPEsi <- rep(NA, 21)
  HIsim <- rep(NA, 21)
  SIsim <- rep(NA, 21)
  sumll <- 0
  Bin   <- seq(0, 1, 0.1)
  Nb    <- length(Bin) - 2

  if(detail ){
    mat <- matrix(NA, tn, 13)
    colnames(mat) <- c('Trial', 'ret', 'SPEHI', 'SPESI', 'ESVHI', 'ESVSI', 'HI', 'SI', 'llHI', 'llSI', 'sumll', 'HIsim', 'SIsim')

    simD <- data; simD[,c(2,3)] <- NA;   # simulated data has same offers
  }

  for (t in 2:tn){

    ret   <- data[t-1,1]

    if (t == 12){ # 'special' readjustment of values upon reversal at trial 11 i.e. t==12 :

      ESVhi[t-1] <- (eta*ESVhi[t-1] + (1-eta) * ESVhi[1])
      ESVsi[t-1] <- (eta*ESVsi[t-1] + (1-eta) * ESVsi[1])

    }

    SPEhi[t] <- ret - ESVhi[t-1]
    SPEsi[t] <- ret - ESVsi[t-1]

    # Associative update of expected social values for harm intent and selfishness intent.
      ESVhi[t] <- ESVhi[t-1] + alpha * SPEhi[t]
      ESVsi[t] <- ESVsi[t-1] + alpha * SPEsi[t]

    # Simplest linear map between ESV and expectation of reported value:
    expHI    <- wHI0 - (wHI * ESVhi[t])
    expSI    <- wSI0 - (wSI * ESVsi[t])

    if (expHI > 1){
      expHI = 0.999999
    }else if (expHI < 0) {
      expHI = 0.000001}

    if (expSI > 1){
      expSI = 0.999999
    }else if (expSI < 0) {
      expSI = 0.000001}

    # probability of actually observed data at this point:
    #prHI = dnorm( data[t-1,2], expHI , SD)
    #prSI = dnorm( data[t-1,3], expSI , SD)

    refHI = pnorm(c(Bin[2:10], 1), data[t-1,2], SD) - pnorm(c(0, Bin[2:10]), data[t-1,2], SD)
    refSI = pnorm(c(Bin[2:10], 1), data[t-1,3], SD) - pnorm(c(0, Bin[2:10]), data[t-1,3], SD)

    hind = round((Nb)*(expHI)+1);
    sind = round((Nb)*(expSI)+1);

    prHI = refHI[hind]
    prSI = refSI[sind]

    llHI = log(prHI)
    llSI = log(prSI)

    # if details asked for, generate some simulated data too:
    if (detail ){
      simD[t-1,2] <- rnorm(1, expHI , SD)
      simD[t-1,3] <- rnorm(1, expSI , SD)
      # clip off ridiculous values:
      if  (simD[t-1,2] < 0){ simD[t-1,2] <- 0 }
      if  (simD[t-1,2] > 1){ simD[t-1,2] <- 1 }
      if  (simD[t-1,3] < 0){ simD[t-1,3] <- 0 }
      if  (simD[t-1,3] > 1){ simD[t-1,3] <- 1 }
    }

    sumll = sumll + llHI + llSI


    if (detail ){
      mat[t,1]  <- t-1
      mat[t,2]  <- ret
      mat[t,3]  <- SPEhi[t]
      mat[t,4]  <- SPEsi[t]
      mat[t,5]  <- ESVhi[t]
      mat[t,6]  <- ESVsi[t]
      mat[t,7]  <- data[t-1, 2]
      mat[t,8]  <- data[t-1, 3]
      mat[t,9]  <- prHI
      mat[t,10] <- prSI
      mat[t,11] <- sumll
      mat[t,12] <- simD[t-1, 2]
      mat[t,13] <- simD[t-1, 3]

    }

  }   # main loop end

  if(detail ){
    res <- list();
    res[[1]] <- mat
    res[[2]] <- c(colSums(mat[2:tn,c(9:10)]), mat[tn,c(11)] )
    res[[3]] <- simD
    res[[4]] <- par
    names(res) <- c('evo','LLs','simD','par')
    return(res)
  } else {
    return(sumll)
  }
}

Assoc_model_salience <- function(par, data, detail = F){


  data = as.data.frame(data)

  wHI0     <- par[1]   # harm intent report intercept
  wSI0     <- par[2]
  wHI      <- par[3]    #    "   "      "    slope w.r.t. HI expected social value
  wSI      <- par[4]
  alphaHI  <- par[5]    # lr for harm intent
  alphaSI  <- par[6]    # lr for selfishness intent
  SD       <- par[7]

  tn = 21

  ESVhi <- rep(NA, 21)
  ESVsi <- rep(NA, 21)
  ESVhi[1] <- par[8] # initial values of ESV for hi
  ESVsi[1] <- par[9] # initial values of ESV for si

  salHI <- par[10]
  salSI <- par[11]

  if(length(par) < 12){
  eta   <- 1
  } else if (length(par) == 12){
  eta   <- par[12]   # rigidity measure
  }

  SPEhi <- rep(NA, 21)
  SPEsi <- rep(NA, 21)
  HIsim <- rep(NA, 21)
  SIsim <- rep(NA, 21)

  salience <- matrix(
    NA,
    nrow = 21,
    ncol = 3) #empty matrix for salience loop
  colnames(salience) <- c("HI", 'SI', "Trial")

  salience[1, 1:3] <- c(rep(1, 2),0)
  salience[2:21,3] <- 1:20

  sumll <- 0
  Bin   <- seq(0, 1, 0.1)
  Nb = length(Bin) - 2

  if(detail ){
    mat <- matrix(NA, tn, 13)
    colnames(mat) <- c('Trial', 'ret', 'SPEHI', 'SPESI', 'ESVHI', 'ESVSI', 'HI', 'SI', 'llHI', 'llSI', 'sumll', 'HIsim', 'SIsim')

    simD <- data; simD[,c(2,3)] <- NA;   # simulated data has same offers
  }

  for (t in 2:tn){

    ret   <- data[t-1,1]

    if (t == 12){ # 'special' readjustment of values upon reversal at trial 11 i.e. t==12 :

      ESVhi[t-1] <- (eta*ESVhi[t-1] + ((1-eta) * ESVhi[1]))
      ESVsi[t-1] <- (eta*ESVsi[t-1] + ((1-eta) * ESVsi[1]))

    }

    SPEhi[t] <- ret - ESVhi[t-1]
    SPEsi[t] <- ret - ESVsi[t-1]

    salience[t, 1]          <- (salHI * abs(SPEhi[t])) + ((1-salHI) * salience[t-1,1]) #Salience modifier of return
    learningHI              <- alphaHI * salience[t,1] #learning rate based on lrc and salience modifier
    salience[t, 2]          <- (salSI * abs(SPEhi[t])) + ((1-salSI) * salience[t-1,2]) #Salience modifier of return
    learningSI              <- alphaSI * salience[t,2] #learning rate based on lrc and salience modifier

    if(learningHI > 1) {learningHI <- 0.99999}
    if(learningHI < 0) {learningHI <- 0.00001}
    if(learningSI > 1) {learningSI <- 0.99999}
    if(learningSI < 0) {learningSI <- 0.00001}

    # Associative update of expected social values for harm intent and selfishness intent.
    # In the case of reversal models, this should in due course include salience and resetting terms ...
    ESVhi[t] <- ESVhi[t-1] + learningHI * SPEhi[t]
    ESVsi[t] <- ESVsi[t-1] + learningSI * SPEsi[t]

    # Simplest linear map between ESV and expectation of reported value:
    expHI    <- wHI0 - (wHI * ESVhi[t])
    expSI    <- wSI0 - (wSI * ESVsi[t])

    if (expHI > 1){
    expHI = 0.999999
    }else if (expHI < 0) {
    expHI = 0.000001}

    if (expSI > 1){
      expSI = 0.999999
    }else if (expSI < 0) {
      expSI = 0.000001}

    # probability of actually observed data at this point:
    #prHI = dnorm( data[t-1,2], expHI , SD)
    #prSI = dnorm( data[t-1,3], expSI , SD)

    refHI = pnorm(c(Bin[2:10], 1), data[t-1,2], SD) - pnorm(c(0, Bin[2:10]), data[t-1,2], SD)
    refSI = pnorm(c(Bin[2:10], 1), data[t-1,3], SD) - pnorm(c(0, Bin[2:10]), data[t-1,3], SD)

    hind = round((Nb)*(expHI)+1);
    sind = round((Nb)*(expSI)+1);

    prHI = refHI[hind]
    prSI = refSI[sind]

    llHI = log(prHI)
    llSI = log(prSI)

    # if details asked for, generate some simulated data too:
    if (detail ){
      simD[t-1,2] <- rnorm(1, expHI , SD)
      simD[t-1,3] <- rnorm(1, expSI , SD)
      # clip off ridiculous values:
      if  (simD[t-1,2] < 0){ simD[t-1,2] <- 0 }
      if  (simD[t-1,2] > 1){ simD[t-1,2] <- 1 }
      if  (simD[t-1,3] < 0){ simD[t-1,3] <- 0 }
      if  (simD[t-1,3] > 1){ simD[t-1,3] <- 1 }
    }

    sumll = sumll + llHI + llSI


    if (detail ){
      mat[t,1]  <- t-1
      mat[t,2]  <- ret
      mat[t,3]  <- SPEhi[t]
      mat[t,4]  <- SPEsi[t]
      mat[t,5]  <- ESVhi[t]
      mat[t,6]  <- ESVsi[t]
      mat[t,7]  <- data[t-1, 2]
      mat[t,8]  <- data[t-1, 3]
      mat[t,9]  <- prHI
      mat[t,10] <- prSI
      mat[t,11] <- sumll
      mat[t,12] <- simD[t-1, 2]
      mat[t,13] <- simD[t-1, 3]

    }

  }   # main loop end

  if(detail ){
    res <- list();
    res[[1]] <- mat
    res[[2]] <- c(colSums(mat[2:tn,c(9:10)]), mat[tn,c(11)] )
    res[[3]] <- simD
    res[[4]] <- par
    names(res) <- c('evo','LLs','simD','par')
    return(res)
  } else {
    return(sumll)
  }
}

#######  WRAPPERS TO FACILITATE USE OF REGULARIZATION ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.a Warpper for Inference Model A Log-likelihood ~~~~~~~~~~~~~~~~~~~

# REM               SEb, acc0max, acc0min,  eta,  gam,  wexp,  wrpe, Tresp,  sig
#      parMat <- c( 0.8,  0.8,    0.6,      0.1,  0.8,   0.3,   0.4,  0.2,   0.1)
#      if fixAcc (fixed acceptance proportions) are given, then reset eta to zero,
#      acc0min and max to NA and use fixAcc (for all ptN) in hapPol[1,expInd,ptN].
#      usu. leave fixAcc at default c(0.85,0.70,0.30,0.15)

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
# This example has weakly informative priors for all except piu, which is
# essentially fixed to the value of 2, and err, which is infomative:
#          pHI0   uHI0 pSI0  uSI0   upi   eta    piu    err
# ashape   1.05   1.2  1.1   1.2    1.2   1.1   100   1.05
# bshape   1.05   3.6  1.1   3.6    3.6   1.1   100    4
# lowlim    0      0    0      0     0     0      0     0
# hilim     1     25    1     25    25     1      4     1

msLPhisi1a_20 <- function(ParM, datAr, scbeta0=NA,details=0){

  parM <- as.vector(ParM); # in case it's inputed in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,
                        1.2, 3.6, 0, 25,
                        1.05,1.05,0, 1,
                        1.2, 3.6, 0, 25,
                        1.2, 3.6, 0, 25,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.2, 3.6, 0, 25,
                        1.1,1.1,0, 1  ), 4, 9)
    if(details){
      colnames(scbeta0) <- c('pHI0','uHI0','pSI0','uSI0','upi','etaHI', 'etaSI', 'piu','err')
      rownames(scbeta0) <- c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have 24 elements or so!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - infHISIll_20b(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
       res[[1]] <- Inf
    } else {
       res[[1]] <- mSLPrior - infHISIll_20b(ParM,datAr);
    }
    names(res) <- c('sumL','scbeta0','par','dat')
    return(res)
  }


} # end of msLPhisi1a

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.a. Wrapper for function wth attribution --> policy params

# params c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI) ...
msLPhisi_20d <- function(ParM, datAr, scbeta0=NA,details=0){

  parM <- as.vector(ParM) # in case it's inputted in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,     # 1  for pHI0
                        1.2, 3.6, 0, 25,    # 2  for uHI0
                        1.05,1.05,0, 1,     # 3  for pSI0
                        1.2, 3.6, 0, 25,    # 4  for uSI0
                        1.2, 3.6, 0, 25,    # 5  for upi
                        3.6, 3.6, -25, 25,  # 6  for w0
                        1.2, 3.6, 0, 25,    # 7  for whi
                        1.2, 3.6, 0, 25,    # 8  for wlo
                        1.05,1.05,0, 1,     # 9  for eta/etaHI
                        1.05,1.05,0, 1      # 10 for etaSI
                        ),
                        nrow=4, ncol=10)
    if(parn==10){
      scbeta0 = scbeta0
    }else if (parn == 9) {
      scbeta0 = scbeta0[,1:9]
    } else {
      scbeta0 = scbeta0[,1:8]
    }

    if(details){
      colnames(scbeta0) <-   c('pHI0','uHI0','pSI0','uSI0','upi','w0',   'whi',  'wsi','etaHI', 'etaSI')
      rownames(scbeta0) <-   c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have say 24 elements or more!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - infHISIll_20d(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf; res[[5]] <- NA;
    } else {
      res[[5]] <- infHISIll_20d(ParM,datAr);
      res[[1]] <- mSLPrior - res[[5]];
    }
    names(res) <- c('sLP','scbeta0','par','dat','sLL')
    return(res)
  }


} # end of msLPhisi_20d

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 Wrapper for function wth associative model

AssocWrapper <- function(ParM, datAr, scbeta0=NA,details=0){

  parM <- as.vector(ParM); # in case it's inputed in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1), 4, 11)

    if(details){
      colnames(scbeta0) <- c('wHI0','wSI0','wHI','wSI','alphaHI','alphaSI','SD', 'ESVhi','ESVsi', 'etaHI', 'etaSI')
      rownames(scbeta0) <- c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have 24 elements or so!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - Assoc_model(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf
    } else {
      res[[1]] <- mSLPrior - Assoc_model(ParM,datAr);
    }
    names(res) <- c('sumL','scbeta0','par','dat')
    return(res)
  }


} # end of msLPhisi1a

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 Wrapper for function wth associative model

AssocWrapper_singleESV <- function(ParM, datAr, scbeta0=NA,details=0){

  parM <- as.vector(ParM); # in case it's inputed in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1), 4, 8)

    if(details){
      colnames(scbeta0) <- c('wHI0','wSI0','wHI','wSI','alpha', 'SD', 'ESV', 'eta')
      rownames(scbeta0) <- c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have 24 elements or so!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - Assoc_model_singleESV(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf
    } else {
      res[[1]] <- mSLPrior - Assoc_model_singleESV(ParM,datAr);
    }
    names(res) <- c('sumL','scbeta0','par','dat')
    return(res)
  }


} # end of msLPhisi1a

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 Wrapper for function wth associative model

AssocWrapper_singleAlpha <- function(ParM, datAr, scbeta0=NA,details=0){

  parM <- as.vector(ParM); # in case it's inputed in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1), 4, 9)

    if(details){
      colnames(scbeta0) <- c('wHI0','wSI0','wHI','wSI','alpha', 'SD', 'ESVhi','ESVsi', 'eta')
      rownames(scbeta0) <- c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have 24 elements or so!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - Assoc_model_singleAlpha(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf
    } else {
      res[[1]] <- mSLPrior - Assoc_model_singleAlpha(ParM,datAr);
    }
    names(res) <- c('sumL','scbeta0','par','dat')
    return(res)
  }


} # end of msLPhisi1a


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 Wrapper for function wth associative model

AssocWrapper_salience <- function(ParM, datAr, scbeta0=-1,details=0){

  parM <- as.vector(ParM); # in case it's inputed in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1,
                        1.05,1.05,0, 1), 4, 12)

    if(parn==12){
      scbeta0 = scbeta0
    }else if (parn == 11) {
      scbeta0 = scbeta0[,1:11]
    }

    if(details){
      colnames(scbeta0) <- c('wHI0','wSI0','wHI','wSI','alphaHI', 'alphaSI','SD', 'ESVhi','ESVsi', 'SalHI', 'SalSI', 'eta')
      rownames(scbeta0) <- c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have 24 elements or so!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - Assoc_model_salience(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf
    } else {
      res[[1]] <- mSLPrior - Assoc_model_salience(ParM,datAr);
    }
    names(res) <- c('sumL','scbeta0','par','dat')
    return(res)
  }


} # end of msLPhisi1a


#######  Functions to Fit Data ###############

## FUNCTION TO RUN SANN ##
runSANN <- function(par, fn, n2do, ncores, filsav){
  source('/Volumes/GoogleDrive/My Drive/Dropbox/PhD/dictator_paranoia/MoutoussisRScript.R')
  source('/Volumes/GoogleDrive/My Drive/Dropbox/PhD/MOBS/MOBS2/ComputationalModels/SerialDictatorReversal/Dictator_Reversal_Functions.R')
  registerDoParallel(cores = ncores)
  x <- foreach(pt=1:n2do, .combine = rbind) %dopar% {
    fitAttempt <- NA; # clear the decks
    dat <- as.matrix(alld[[pt]][[1]][,3:5])
    fitAttempt   <- optim(par = par,
                          fn = fn,
                          gr = NULL,
                          datAr = dat,
                          scbeta0 = -1,
                          method  = 'SANN')
    L <- NA;         try(L <- -fn(fitAttempt$par,dat,scbeta0 = NA) );
    data.frame( ID        =  unlist(alld[[pt]][[2]][1]),
                one       =  fitAttempt$par[1],
                two       =  fitAttempt$par[2],
                three     =  fitAttempt$par[3],
                four      =  fitAttempt$par[4],
                five      =  fitAttempt$par[5],
                six       =  fitAttempt$par[6],
                seven     =  fitAttempt$par[7],
                eight     =  fitAttempt$par[8],
                nine      =  fitAttempt$par[9],
                ten       =  fitAttempt$par[10],
                eleven    =  fitAttempt$par[11],
                twelve    =  fitAttempt$par[12],
                ll        =  L)
  }
  write.csv(x, paste("/Volumes/GoogleDrive/My Drive/Dropbox/PhD/MOBS/MOBS2/ComputationalModels/SerialDictatorReversal/",filsav,".csv", sep = ""))
}

runOPTIM <- function(par, fn, n2do, sann2hisi, ncores, filsav, wrapperType){
  source('/Volumes/GoogleDrive/My Drive/Dropbox/PhD/dictator_paranoia/MoutoussisRScript.R')
  source('/Volumes/GoogleDrive/My Drive/Dropbox/PhD/MOBS/MOBS2/ComputationalModels/SerialDictatorReversal/Dictator_Reversal_Functions.R')
  registerDoParallel(cores = ncores)
  x <- foreach(pt=1:n2do, .combine = rbind) %dopar% {
    tryP <- as.numeric(sann2hisi[pt,2:length(sann2hisi)] );
    if(sum(is.na(tryP))){ tryP <- par }
    fitAttempt <- NA; # clear the decks and attempt the new optim fit
    prm        <- tryP * NA;
    L          <- NA;
    dat <- as.matrix(alld[[pt]][[1]][,3:5])
    fitAttempt   <- try(optim(par = tryP,
                          fn = fn,
                          datAr = dat,
                          scbeta0 = wrapperType))
    if(inherits(fitAttempt, "try-error"))
    {
      #error handling code
    fitAttempt   <- optim(par = par,
                          fn = fn,
                          datAr = dat,
                          scbeta0 = wrapperType)
    }
    prm        <- fitAttempt$par;
    if(sum(is.na(prm))){               # if fit failed, revert to best up to now:
      prm <- as.numeric( sann2hisi[pt,2:length(sann2hisi)] );
    }
    L <- -fn(prm,dat,NA);  # plain log-likelihood

    data.frame( ID        =  unlist(alld[[pt]][[2]][1]),
                one       =  prm[1],
                two       =  prm[2],
                three     =  prm[3],
                four      =  prm[4],
                five      =  prm[5],
                six       =  prm[6],
                seven     =  prm[7],
                eight     =  prm[8],
                nine      =  prm[9],
                ten       =  prm[10],
                eleven    =  prm[11],
                twelve    =  prm[12],
                ll        =  L,
                Persec    =  unlist(alld[[pt]][[3]][1,'Persec']),
                ICARTot   =  unlist(alld[[pt]][[3]][1,'ICARTot']),
                Order     =  unlist(alld[[pt]][[3]][1,'Order']))
  }
  write.csv(x, paste("/Volumes/GoogleDrive/My Drive/Dropbox/PhD/MOBS/MOBS2/ComputationalModels/SerialDictatorReversal/",filsav,".csv", sep = ""))
}


#######  Functions to Generate Data ###############

generateDataBayes <- function(fn, ncores, n2do, modelName){

  registerDoParallel(cores = ncores)
  simulated <- foreach (i = 1:n2do, .combine = rbind) %dopar% {

    #Generate data over raw parameter space
    parms           <- rnorm(10, mean = 0, sd = 2)
    #Transform data into sensible parameter space
    for (n in 1:length(parms)){
      if (n %in% c(1, 3, 9, 10)){
        parms[n] <- 1/(1+exp(-parms[n]))
      } else if (n %in% c(2, 4, 5, 7, 8)) {
        parms[n] <- exp(parms[n])
        if(parms[n] < 0.05){parms[n] <- 0.05}
        if(parms[n] > 25){parms[n] <- 25}
      } else {
        parms[n] <- parms[n]
      }
    }
    parms
    if (modelName == '1eta'){parms <- parms[1:9]} else if (modelName == '0eta'){parms <- parms[1:8]}

    test                <- fn(as.numeric(parms),
                              as.matrix(alld[[i]][[1]][,3:5]),
                              detail=T)
    sumll               <- try(fn(as.numeric(parms),
                                  as.matrix(alld[[i]][[1]][,3:5])))
    data.frame(
      ID       = alld[[i]][[1]][,'ID'],
      Trial    = as.numeric(1:20),
      ret      = as.numeric(test$evo[2:21,'ret']),
      HIsim    = as.numeric(test$evo[2:21,'HIsim']),
      SIsim    = as.numeric(test$evo[2:21,'SIsim']),
      ll       = as.numeric(sumll),
      one      =  parms[1],
      two      =  parms[2],
      three    =  parms[3],
      four     =  parms[4],
      five     =  parms[5],
      six      =  parms[6],
      seven    =  parms[7],
      eight    =  parms[8],
      nine     =  parms[9],
      ten      =  parms[10],
      model    =  modelName,
      Order    =  alld[[i]][[3]][,'Order'],
      id       = i
    )

  }
  return(simulated)
}

generateDataAssoc <- function(fn, ncores, n2do, modelName){

  registerDoParallel(cores = ncores)
  simulated <- foreach (i = 1:n2do, .combine = rbind) %dopar% {

    #Generate data over raw parameter space
    parms           <- rnorm(11, mean = 0, sd = 2)
    #Transform data into sensible parameter space
    for (n in 1:length(parms)){
        parms[n] <- 1/(1+exp(-parms[n]))
    }

    test                <- fn(as.numeric(parms),
                              as.matrix(alld[[i]][[1]][,3:5]),
                              detail=T)
    sumll               <- try(fn(as.numeric(parms),
                                  as.matrix(alld[[i]][[1]][,3:5])))
    data.frame(
      ID       = alld[[i]][[1]][,'ID'],
      Trial    = as.numeric(1:20),
      ret      = as.numeric(test$evo[2:21,'ret']),
      HIsim    = as.numeric(test$evo[2:21,'HIsim']),
      SIsim    = as.numeric(test$evo[2:21,'SIsim']),
      ll       = as.numeric(sumll),
      one      =  parms[1],
      two      =  parms[2],
      three    =  parms[3],
      four     =  parms[4],
      five     =  parms[5],
      six      =  parms[6],
      seven    =  parms[7],
      eight    =  parms[8],
      nine     =  parms[9],
      ten      =  parms[10],
      eleven   =  parms[11],
      model    =  modelName,
      Order    =  alld[[i]][[3]][,'Order'],
      id       = i
    )

  }
  return(simulated)
}

simulateData <- function(fn, filename, ncores, n2do, parms, modelName, plot = 0){
  hisi <- read.csv(paste("/Users/josephbarnby/My Drive/Dropbox/PhD/MOBS/MOBS2/ComputationalModels/SerialDictatorReversal/", filename,'.csv', sep = ""))
  hisi_nona <- hisi %>% dplyr::select(-X) %>% arrange(ID)
  hisi_nona[,2:(parms+1)] <- sapply(hisi_nona[,2:(parms+1)], as.numeric)
  colnames(hisi_nona)[1] <- "ID";
  #once we have the list for each individual, we can run the model over each individually
  registerDoParallel(cores = ncores)
  simulated <- foreach (i = 1:n2do, .combine = rbind) %dopar% {

    test                <- fn(as.numeric(hisi_nona[i,2:(parms+1)]),
                              as.matrix(alld[[i]][[1]][,3:5]),
                              detail=T)
    sumll               <- try(fn(as.numeric(hisi_nona[i,2:(parms+1)]),
                                  as.matrix(alld[[i]][[1]][,3:5])))
    data.frame(
      ID       = alld[[i]][[1]][,'ID'],
      Trial    = as.numeric(1:20),
      ret      = as.numeric(test$evo[2:21,'ret']),
      HIsim    = as.numeric(test$evo[2:21,'HIsim']),
      SIsim    = as.numeric(test$evo[2:21,'SIsim']),
      ll       = as.numeric(sumll),
      HI       = as.numeric(test$evo[2:21,'HI']),
      SI       = as.numeric(test$evo[2:21,'SI']),
      one      =  test$par[1],
      two      =  test$par[2],
      three    =  test$par[3],
      four     =  test$par[4],
      five     =  test$par[5],
      six      =  test$par[6],
      seven    =  test$par[7],
      eight    =  test$par[8],
      nine     =  test$par[9],
      ten      =  test$par[10],
      eleven   =  test$par[11],
      twelve   =  test$par[12],
      model    =  modelName,
      Order    =  alld[[i]][[3]][,'Order']
    )

  }

  if (plot == 1){

    return(ggplot(simulated %>%
                    pivot_longer(7:8, 'Attribution', values_to = 'Metric')) +
             stat_summary(aes(Trial, Metric, linetype = Attribution), color = 'black',  size = 1, geom = 'line') +
             stat_summary(aes(Trial, Metric, group = Attribution), fill = 'black',  size = 1, geom = 'ribbon', alpha = 0.5) +
             stat_summary(aes(Trial, HIsim), color = 'red', geom = 'line')+
             stat_summary(aes(Trial, SIsim), color = 'blue', geom = 'line')+
             scale_linetype_manual(name = 'True Attribution', values = c(1, 3))+
             labs(x = 'Trial', y = 'Attribution')+
             facet_wrap(~Order) +
             tidybayes::theme_tidybayes()+
             theme(axis.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   strip.text.x = element_text(size = 14),
                   strip.background.x = element_blank()))
  }
  return(simulated)
}

# Janitor --------------------------------------------------------------

sigmoid = function(x){1/(1+exp(-x))}

cleanBB <- function(x, y, i = 2){
  z <- x$cbm[,,1]$output[,,1]$parameters[[i]][[1]] %>%
    as.data.frame() %>%
    mutate(id = 1:length(x$cbm[,,1]$output[,,1]$parameters[[i]][[1]][,1])) %>%
    plyr::join(as.data.frame(y), by = 'id') %>%
    mutate(V1 = 1/(1+exp(-V1)),
           V2 = exp(V2),
           V3 = 1/(1+exp(-V3)),
           V4 = exp(V4),
           V5 = exp(V5),
           V7 = exp(V7),
           V8 = exp(V8),
           V9 = 1/(1+exp(-V9)))
  if(i == 3){z$V10 <- 1/(1+exp(-z$V10))}
  colnames(z)[1:9] <- c('pHI0', 'uHI0', 'pSI0', 'uSI0', 'uPi', 'w0', 'wHI', 'wSI', 'eta')
  return(z)
}

cleanAssoc <- function(x, y){
  z <- x$cbm[,,1]$output[,,1]$parameters[[2]][[1]] %>%
    as.data.frame() %>%
    mutate(id = 1:length(x$cbm[,,1]$output[,,1]$parameters[[2]][[1]][,1])) %>%
    plyr::join(as.data.frame(y), by = 'id') %>%
    mutate(across(1:8, .fns = sigmoid))
  colnames(z)[1:8] <- c('wHI0', 'wSI0', 'wHI', 'wSI', 'alpha', 'SD', 'ESV0', 'eta')
  return(z)
}

# AUXILLIARY --------------------------------------------------------------


infHISIll_20d0eta <- function(ptp,d,details=0,plots=0) {
  # rem ptp is of the form c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI))
  # rem the corresponding wrapper is msLPhisi_20d
  #
  # Test with:
  # toyDf  <- c(0.50, 0.00, 0.50, 0.00, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.00, 0.00, 0.50, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.40, 0.40, 0.40, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.65, 0.65, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80);  toyDf <- matrix(toyDf,20,3); colnames(toyDf) <- c('ret','hi','si')
  # params <-  c(0.4,2.00, 0.75, 2, 0.05, -1,  0.5,   1,   0.75 )
  # test <- infHISIll_20d(params,toyDf,details=1,plots=1)

  d <- as.matrix(d)

  tn = 10;  # 10 trials
  phase = 2;  # one partner, two phases

  # Parameters to help describe personality and action. These are
  # not free parameters to be fitted, i.e. don't vary across pts.
  Nb = 9;     # number of bins for each personality dim.
  Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
  Na = 2;     # actual number of actions in Joe's expt.

  # i.e. proportions that may be received by the Other.
  # In general will be set to Nb.

  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(ptp[3], ptp[4],Nb); PHI0 = noisyBino(ptp[1],ptp[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0);

  if (plots){
    opar=par()
    # Action names (was:  anames = c('selfish','v.tight','tight','fair','kind','v.kind','altruistic') )
    anames = c('selfish','fair')
  }

  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = ptp[5];     # convenience copy of policy variability (inverse precision)
  w0  = ptp[6];     # This and two next will help map HI, SI onto policy (zero or half return)
  whi = ptp[7];
  wsi = ptp[8];
  if ((whi < 0) | (wsi < 0)){
    print('Params c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi:')
    print(ptp)
    error("whi,wsi must all be non-negative")
  }

  # NOT in 1.d.: if (length(ptp) > 7){piu=ptp[8]} else {piu = 1} # arbitrary policy-inverse-uncertainty
  # param; could set to e.g. 1/upi to reference it to self, roughly.
  err =0.02/(Nb*Nb)
  # param; Note scaling by the number of attribution states considered.

  # Set up the map between 'attributes' and actions :
  # NOT in 1.d.: pif = array(NA,c(Nb,Nb,Nf));    # This will hold all possible policies for full range of actions
  pi  = array(NA,c(Nb,Nb,Na));    # Possible policies for actual range of actions

  offs = (Nb+1)/2
  for (SI in 1:Nb){
    for (HI in 1:Nb) {
      pi[SI,HI,1]  = invlogit(w0 + whi*(HI-offs) + wsi*(SI-offs))  # prob. of ufair offer goes up w. HI, SI
      # the offs is to center at (usually) 5
      pi[SI,HI,2]  = 1 - pi[HI,SI,1]
      # NOT in 1.d.: if (plots){
      # NOT in 1.d.: piav[HI,SI] = sum(pi[HI,SI,]*c(0,0.5))
      # NOT in 1.d.: pifav[HI,SI] = sum(pif[HI,SI,]*(1:Nb))
      # NOT in 1.d.: }
    }
  }

  if (plots){           # Display the pi(zeror return) returns as a heatmap
    heatmap( pi[,,1],
             Rowv=NA, Colv=NA, col = topo.colors(512),
             scale="none", margins=c(5,8),asp=1,
             labRow=((0:(Nb-1))+0.5)*0.1,
             labCol=((0:(Nb-1))+0.5)*0.1,
             main = paste('\n attributes vs. pi(zero return)'),
             xlab='HI',ylab='SI\n')
  }


  ### Run the inference

  # Here the posterior of one trial will form the prior for the next trial,
  # staring from the population prior beliefs PSIHI0.
  #

  sll = 0;
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)

  if (details ){
    llout = list();

    llout[[1]]=0;

    # [[2]] is the parameter vector
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','w0', 'wHI', 'wSI', 'piu','err')
    llout[[2]] = c(ptp[1:8],piu,err);  names(llout[[2]]) <- hd;

    # [[3]] is the detailed evolution, evao:
    llout[[3]] = matrix(NA,tn*phase+1,10);
    llout[[3]][,1] = c(0,1:(tn*phase))
    colnames(llout[[3]]) = c('trial','ret','HI','SI','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*phase),2:4] <- as.matrix(d)

    # [[4]] will be a big 9 x 9 x 21 array with detailed policies
    # Hypothetical (attribn. reporting) policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*phase)))
    llout[[4]][,,1] <- pri0;

    # [[5]] will be the HISI -> return policy (usef for the Other)
    llout[[5]] <- list();
    #llout[[5]][[1]] <- pi; llout[[5]][[2]] <- pif;  llout[[5]][[3]] <- piav;  llout[[5]][[4]] <- pifav;
    #names(llout[[5]]) <- c('piBinary','piFull','piBinAv','piFullAv')

    names(llout) <- c('sll','par', 'evo','policy','attr2pol')
  }

  #for (block in 1:2){    # counter of 2 different policy styles

  post <- pri0; # this is the belief that will be updated with each trial

  # rows to be processed and convenience copies of other's actions,
  # malice and greed attributions:
  ro = 1:(2*tn); # rows of data matrix
  as = d[ro,1];  aind = round((Na-1)*as+1)
  hi = d[ro,2];  hind = round((Nb-1)*hi+1)
  si = d[ro,3];  sind = round((Nb-1)*si+1)

  for (t in 1:(tn+tn)){  # loop
    if(plots){
      heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
              margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
              main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
    }

    pri = post;              # new prior is last posterior
    # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
    # over the character of the partner. This post is their post. beliefs
    post = pi[,,aind[t]] * pri; post = post / sum(as.vector(post))  # Bayes

    # Now the probability of the response, incl. the lapse-rate-like err:
    pol = post^(1/upi);  pol = pol/sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    lik = pol[sind[t],hind[t]];
    sll = sll + log(lik);         # accumulate sum log lik

    if (details ){

      llout$evo[(t+1),'lik'] <- lik
      llout$evo[(t+1),'ll'] <- log(lik)
      # find mode of pol
      c = max.col(pol);  # this finds the max. col. pos. for each row
      m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
      r=which.max(m);  # ... so was to find the row with the mode of the distro.
      llout$evo[(t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
      # Now sample randomly
      hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
      sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the phase from the corresponding conditional.
      # debug: barplot(colSums(pol));       barplot(pol[,hisim]);
      llout$evo[(t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
      llout$policy[,,(t+1)] <- pol

    }
  }

  if (plots) {
    heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
            margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
            main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
    par(opar)  # restore state of graphics where we found it.
  }
  #}

  if (details ){llout$sll <- sll} else {llout <- sll}
  return(llout)

}


infHISIll_20d2eta <- function(ptp,d,details=0,plots=0) {
  # rem ptp is of the form c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI))
  # rem the corresponding wrapper is msLPhisi_20d
  #
  # Test with:
  # toyDf  <- c(0.50, 0.00, 0.50, 0.00, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.00, 0.00, 0.50, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.40, 0.40, 0.40, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.65, 0.65, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80);  toyDf <- matrix(toyDf,20,3); colnames(toyDf) <- c('ret','hi','si')
  # params <-  c(0.4,2.00, 0.75, 2, 0.05, -1,  0.5,   1,   0.75 )
  # test <- infHISIll_20d(params,toyDf,details=1,plots=1)

  d <- as.matrix(d)

  tn = 10;  # 10 trials
  phase = 2;  # one partner, two phases

  # Parameters to help describe personality and action. These are
  # not free parameters to be fitted, i.e. don't vary across pts.
  Nb = 9;     # number of bins for each personality dim.
  Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
  Na = 2;     # actual number of actions in Joe's expt.

  # i.e. proportions that may be received by the Other.
  # In general will be set to Nb.

  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(ptp[3], ptp[4],Nb); PHI0 = noisyBino(ptp[1],ptp[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0);

  if (plots){
    opar=par()
    # Action names (was:  anames = c('selfish','v.tight','tight','fair','kind','v.kind','altruistic') )
    anames = c('selfish','fair')
  }

  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = ptp[5];     # convenience copy of policy variability (inverse precision)
  w0  = ptp[6];     # This and two next will help map HI, SI onto policy (zero or half return)
  whi = ptp[7];
  wsi = ptp[8];
  if ((whi < 0) || (wsi < 0)){
    print('Params c(pHI0,uHI0,pSI0,uSI0,upi,w0,whi,wsi,etaHI (or single one), (optionally etaSI)) :')
    print(ptp)
    error("whi,wsi must all be non-negative")
  }
  etaHI = ptp[9];   # Will allow for different rigidity along Harmful intent and SI - default is only one:
  etaSI = ptp[10];

  # NOT in 1.d.: if (length(ptp) > 7){piu=ptp[8]} else {piu = 1} # arbitrary policy-inverse-uncertainty
  # param; could set to e.g. 1/upi to reference it to self, roughly.
  err =0.02/(Nb*Nb)
  # param; Note scaling by the number of attribution states considered.

  # Set up the map between 'attributes' and actions :
  # NOT in 1.d.: pif = array(NA,c(Nb,Nb,Nf));    # This will hold all possible policies for full range of actions
  pi  = array(NA,c(Nb,Nb,Na));    # Possible policies for actual range of actions
  fp =  c(1,1,1,0,0,0,0,0,0);    # auxiliary vector to further bin down pi, here to a two-bin vector
  # Joe suggestions c(1,1,1,1,1,1,1,0,0),  c(1,1,0,0,0,0,0,0,0);
  # Original fp = c(1,1,1,0,0,0,0,0,0);
  f2p = t(matrix(c(fp,1-fp),Nb,Na));


  offs = (Nb+1)/2
  for (SI in 1:Nb){
    for (HI in 1:Nb) {
      pi[SI,HI,1]  = invlogit(w0 + wsi*(SI-offs) + whi*(HI-offs))  # prob. of ufair offer goes up w. HI, SI
      # the offs is to center at (usually) 5
      pi[SI,HI,2]  = 1 - pi[SI,HI,1]

    }
  }

  if (plots){           # Display the pi(zeror return) returns as a heatmap
    heatmap( pi[,,1],
             Rowv=NA, Colv=NA, col = topo.colors(512),
             scale="none", margins=c(5,8),asp=1,
             labRow=((0:(Nb-1))+0.5)*0.1,
             labCol=((0:(Nb-1))+0.5)*0.1,
             main = paste('\n attributes vs. pi(zero return)'),
             xlab='HI',ylab='SI\n')
  }


  ### Run the inference

  # Here the posterior of one trial will form the prior for the next trial,
  # staring from the population prior beliefs PSIHI0.
  #

  sll = 0;
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)

  if (details ){
    llout = list();

    llout[[1]]=0;

    # [[2]] is the parameter vector
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','w0', 'wHI', 'wSI', 'etaHI','etaSI','piu','err')
    llout[[2]] = c(ptp[1:10],piu,err);  names(llout[[2]]) <- hd;

    # [[3]] is the detailed evolution, evao:
    llout[[3]] = matrix(NA,tn*phase+1,10);
    llout[[3]][,1] = c(0,1:(tn*phase))
    colnames(llout[[3]]) = c('trial','ret','HI','SI','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*phase),2:4] <- as.matrix(d)

    # [[4]] will be a big 9 x 9 x 21 array with detailed policies
    # Hypothetical (attribn. reporting) policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*phase)))
    llout[[4]][,,1] <- pri0;

    # [[5]] will be the HISI -> return policy (usef for the Other)
    llout[[5]] <- list();
    #llout[[5]][[1]] <- pi; llout[[5]][[2]] <- pif;  llout[[5]][[3]] <- piav;  llout[[5]][[4]] <- pifav;
    #names(llout[[5]]) <- c('piBinary','piFull','piBinAv','piFullAv')

    names(llout) <- c('sll','par', 'evo','policy','attr2pol')
  }

  for (block in 1:2){    # counter of 2 different policy styles

    if (block > 1){

      PHIPost <- colSums(post);
      PSIPost <- rowSums(post);
      # Now do the exact equivalent of PSIHI0 = PSI0 %*% t(PHI0) above :
      pri0 <-  ((PSI0 * (1-etaSI)) + (etaSI * PSIPost)) %*% t((PHI0 * (1-etaHI)) + (etaHI * PHIPost))        ;
    }

    post <- pri0; # this is the belief that will be updated with each trial

    # rows to be processed and convenience copies of other's actions,
    # malice and greed attributions:
    ro = ((block-1)*tn+1):(block*tn); # rows of data matrix
    as = d[ro,1];  aind = round((Na-1)*as+1)
    hi = d[ro,2];  hind = round((Nb-1)*hi+1)
    si = d[ro,3];  sind = round((Nb-1)*si+1)

    for (t in 1:tn){  # loop
      if(plots){
        heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
                margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
                main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
      }

      pri = post;              # new prior is last posterior
      # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
      # over the character of the partner. This post is their post. beliefs
      post = pi[,,aind[t]] * pri; post = post / sum(as.vector(post))  # Bayes

      # Now the probability of the response, incl. the lapse-rate-like err:
      pol = post^(1/upi);  pol = pol/sum(as.vector(pol));
      pol = (pol+err)/(1+err*length(pol))
      lik = pol[sind[t],hind[t]];
      sll = sll + log(lik);         # accumulate sum log lik

      if (details ){

        llout$evo[((block-1)*tn+t+1),'lik'] <- lik
        llout$evo[((block-1)*tn+t+1),'ll'] <- log(lik)
        # find mode of pol
        c = max.col(pol);  # this finds the max. col. pos. for each row
        m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
        r=which.max(m);  # ... so was to find the row with the mode of the distro.
        llout$evo[(((block-1)*tn)+t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
        # Now sample randomly
        hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
        sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the phase from the corresponding conditional.
        # debug: barplot(colSums(pol));       barplot(pol[,hisim]);
        llout$evo[((block-1)*tn+t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
        llout$policy[,,(block-1)*tn+t+1] <- pol

      }
    }

    if (plots) {
      heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none",
              margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1),
              main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
      par(opar)  # restore state of graphics where we found it.
    }
  }

  if (details ){llout$sll <- sll} else {llout <- sll}
  return(llout)

}

msLPhisi_20d2eta<- function(ParM, datAr, scbeta0=NA,details=0){

  parM <- as.vector(ParM); # in case it's inputted in another format
  parn <- length(parM)

  if ((scbeta0[1] < 0) && !is.na(scbeta0)){
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,     # 1  for pHI0
                        1.2, 3.6, 0, 25,    # 2  for uHI0
                        1.05,1.05,0, 1,     # 3  for pSI0
                        1.2, 3.6, 0, 25,    # 4  for uSI0
                        1.2, 3.6, 0, 25,    # 5  for upi
                        3.6, 3.6, -25, 25,  # 6  for w0
                        1.2, 3.6, 0, 25,    # 7  for whi
                        1.2, 3.6, 0, 25,    # 8  for wlo
                        1.05,1.05,0, 1,     # 9  for etaHI
                        1.05,1.05,0, 1),     # 10 for etaSI),
                      nrow=4, ncol=10)
    if(details){
      colnames(scbeta0) <-   c('pHI0','uHI0','pSI0','uSI0','upi','w0',   'whi',  'wsi','etaHI','etaSI')
      # was:                      1      2      3      4     5     6       7       8      9      10      11     12
      # colnames(scbeta0) <- c('pHI0','uHI0','pSI0','uSI0','upi','etaHI','etaSI','piu','err')
      rownames(scbeta0) <- c('ashape','bshape','min','max')
    }
  }

  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have say 24 elements or more!
    mSLPrior <- mSLPrior - sum(dbetasc( parM,
                                        scbeta0[1,1:parn],scbeta0[2,1:parn],
                                        scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE));
  }

  if (!details){
    if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
      # do not attempt to calculate the likelihood - it will be nonsense anyway.
      return(Inf);
    } else {
      return(mSLPrior - infHISIll_20d2eta(ParM,datAr))
    }
  } else {
    res = list();
    res[[2]] <- scbeta0;
    res[[3]] <- ParM;        res[[4]] <- datAr;
    if (mSLPrior == Inf){
      res[[1]] <- Inf; res[[5]] <- NA;
    } else {
      res[[5]] <- infHISIll_20d2eta(ParM,datAr);
      res[[1]] <- mSLPrior - res[[5]];
    }
    names(res) <- c('sLP','scbeta0','par','dat','sLL')
    return(res)
  }


}


# Generate simulations ----------------------------------------------------

simulatedata_HISI <- function(x,
                              values,
                              samples,
                              plot = 1,
                              trials = 10,
                              partners = 2,
                              pHI0 = 0.5, uHI0 = 2, pSI0 = 0.5, uSI0 = 2, w0 = 0, wHI = 0.5, wSI = 0.5, upi = 2, eta = 0.5) {

  if(!x %in% c('pHI0', 'uHI0', 'pSI0', 'uSI0', 'upi', 'upi01', 'w0', 'wHI', 'wSI', 'eta')){
    stop("Incorrect input! Must be one of: 'pHI0', 'uHI0', 'pSI0', 'uSI0', 'upi', 'upi01', 'w0', 'wHI', 'wSI', 'eta'")
  }

  fair_d   = sample(c(0.5, 0), trials, prob = c(0.8, 0.2), replace = T)
  unfair_d = sample(c(0.5, 0), trials, prob = c(0.2, 0.8), replace = T)
  partial_d= sample(c(0.5, 0), trials, prob = c(0.5, 0.5), replace = T)

  d <- matrix(0, nrow = trials * partners, ncol = 3)

  test_upi    <- list()
  loop_k      <- list()
  loop_order  <- list()
  list_produce<- list()

  for (order in 1:partners){
    for (k in 1:values){
      for (i in 1:samples){

      if        (order == 1 & partners == 2) {d[,1] = c(unfair_d, fair_d)
      } else if (order == 2 & partners == 2) {d[,1] = c(fair_d, unfair_d)
      } else if (order == 1 & partners == 3) {d[,1] = c(unfair_d, fair_d, partial_d)
      } else if (order == 2 & partners == 3) {d[,1] = c(fair_d, partial_d, unfair_d)
      } else if (order == 3 & partners == 3) {d[,1] = c(partial_d, unfair_d, fair_d)
      }

      if(x == 'pHI0') {meanobject = infHISIll_20d(c(k/values, uHI0, pSI0    , uSI0, upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'uHI0') {meanobject = infHISIll_20d(c(pHI0    , k   , pSI0    , uSI0, upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'pSI0') {meanobject = infHISIll_20d(c(pHI0    , uHI0, k/values, uSI0, upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'uSI0') {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , k   , upi     , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'upi')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, k       , w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'upi01'){meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, k/values, w0, wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'w0')   {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , k , wHI     , wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'wHI')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, k/values, wSI     , eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'wSI')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, wHI     , k/values, eta)     , d, details = 1, tn = trials, phase = partners)}
      if(x == 'eta')  {meanobject = infHISIll_20d(c(pHI0    , uHI0, pSI0    , uSI0, upi     , w0, wHI     , wSI     , k/values), d, details = 1, tn = trials, phase = partners)}

      loop_k[[i]]     <- meanobject$evo %>%
                            as.data.frame() %>%
                            na.omit() %>%
                            mutate(loop = ifelse(x %in% c('pHI0', 'pSI0', 'wHI', 'wSI', 'upi01', 'eta'), k/values, k),
                                   id = i,
                                   Order = order)

      }

    loop_order[[k]] <- do.call(rbind,loop_k)
    cat(paste('\n Simulated ', k*samples, 'iterations out of ',values*samples,'| Order = ', order, '| Parameter = ', x))

    }

  test_upi[[order]] <- do.call(rbind, loop_order)

  }#end

  list_produce[[1]] <- do.call(rbind, test_upi)

    if(plot){
    list_produce[[2]] <- list_produce[[1]]  %>%
      pivot_longer(HIsim:SIsim, 'Attribution', values_to = 'values') %>%
      mutate(Attribution = ifelse(Attribution == 'HIsim', 'HI Simulated', 'SI Simulated'),
             Order = ifelse(Order == 1 & partners == 2, 'Unfair | Fair',
                            ifelse(Order == 2 & partners == 2, 'Fair | Unfair',
                                   ifelse(Order == 1 & partners == 3, 'Unfair | Fair | Partial',
                                          ifelse(Order == 2 & partners == 3, 'Fair | Partial | Unfair',
                                                 ifelse(Order == 3 & partners == 3, 'Partial | Unfair | Fair', NA))))),
             ret = ifelse(ret == 0.5, 1, 0))

     list_produce[[2]] <- ggplot(list_produce[[2]], aes(trial, values, color = loop, group = loop, linetype= Attribution))+
        stat_summary(geom = 'line')+
        scale_color_gradient(low = '#BFD7EA', high = '#0B3954', name = x)+
        facet_wrap(Attribution~Order, scales = 'free_y', nrow = 2, ncol = partners)+
        coord_cartesian(ylim = c(0, 1))+
        labs(x = 'Trial', y = 'Simulated Observation')+
        theme_tq()+
        theme(text = element_text(size = 14),
              strip.background.x = element_rect(fill = '#009DDC'),
              legend.key.width = unit(1, 'cm'))
    }

  return(list_produce)

}

