##########################################################################################################################################
### R script: Paternity certainty and parent-offspring conflict jointly explain restrictions on female premarital sex across societies ###
##########################################################################################################################################


### code for the analysis of confounding

setwd ("")

library (rethinking)

memory.limit (size=200000)

norms <- read.csv ("norms_final.csv", header=TRUE, row.names=1)
N_pop <- nrow (norms)

# standardize the variables

ext_aff_s <- (norms$ext_aff - mean (norms$ext_aff))/ (sd (norms$ext_aff)*2)
dir_care_s <- (norms$dir_care - mean (norms$dir_care))/ (sd (norms$dir_care)*2)
fem_cont_s <- (norms$fem_cont - mean (norms$fem_cont))/ (sd (norms$fem_cont)*2)
sex_rat_s <- (norms$sex_rat - mean (norms$sex_rat))/ (sd (norms$sex_rat)*2)
marr_arr_s <- (norms$marr_arr - mean (norms$marr_arr))/ (sd (norms$marr_arr)*2)
dec_mak_s <- (norms$dec_mak - mean (norms$dec_mak))/ (sd (norms$dec_mak)*2)

sex_norms <- as.integer (norms$sex_norms) # outcome remains unchanged, otherwise ulam cannot compute cutpoints

inher <- as.factor(norms$inher)
p_loc <- as.factor(norms$p_loc)
bride_pr <- as.factor (norms$bride_pr)
gift_ex <- as.factor (norms$gift_ex)
dowry <- as.factor (norms$dowry)

### create a data list

d_norm <- list (sex_norms=sex_norms,
                   ext_aff=ext_aff_s,
                   inher=inher,
                   dir_care=dir_care_s,
                   fem_cont=fem_cont_s,
                   sex_rat=sex_rat_s,
                   marr_arr=marr_arr_s,
                   dec_mak=dec_mak_s,
                   p_loc=p_loc,
                   bride_pr=bride_pr,
                   gift_ex=gift_ex,
                   dowry=dowry,
                   N_pop=N_pop,
                   society=1:N_pop,
                   Imat=diag(nrow(norms))
)

### model with extramarital sex excluded ###

m_norm_no_ext <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- bI*inher + bDC*dir_care + bFC*fem_cont + bSR*sex_rat + bMA*marr_arr + bDM*dec_mak + bPloc*p_loc + bBP*bride_pr + bGE*gift_ex + bD*dowry + k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- Imat*sigma,
    
    cutpoints ~ dnorm (0,1.5),
    
    c (bD,bGE,bBP,bPloc,bDM,bMA,bSR,bFC,bDC,bI) ~ dnorm (0,0.5),
    
    sigma ~ dexp (1)
    
  ), data=d_norm, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_no_ext, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bFC","bDC","bI","sigma"))

### get model summary

norm_no_ext <- precis (m_norm_no_ext, 2, omit=c("cutpoints","k","z"))
norm_no_ext <- data.frame (round (norm_no_ext, 3))
write.csv (data.frame (norm_no_ext), file="norm_no_ext_output.csv")

norm_post_no_ext <- extract.samples (m_norm_no_ext, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bFC","bDC","bI","sigma"))
write.csv (norm_post_no_ext, file="norm_post_no_ext.csv")


### model with all three marriage transactions predictors excluded ###

m_norm_no_trans <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- bEX*ext_aff + bI*inher + bDC*dir_care + bFC*fem_cont + bSR*sex_rat + bMA*marr_arr + bDM*dec_mak + bPloc*p_loc + k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- Imat*sigma,
    
    cutpoints ~ dnorm (0,1.5),
    
    c (bPloc,bDM,bMA,bSR,bFC,bDC,bI,bEX) ~ dnorm (0,0.5),
    
    sigma ~ dexp (1)
    
  ), data=d_norm, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_no_trans, pars=c("cutpoints","bPloc","bDM","bMA","bSR","bFC","bDC","bI","bEX","sigma"))

### get model summary

norm_no_trans <- precis (m_norm_no_trans, 2, omit=c("cutpoints","k","z"))
norm_no_trans <- data.frame (round (norm_no_trans, 3))
write.csv (data.frame (norm_no_trans), file="norm_no_trans_output.csv")

norm_post_no_trans <- extract.samples (m_norm_no_trans, pars=c("cutpoints","bPloc","bDM","bMA","bSR","bFC","bDC","bI","bEX","sigma"))
write.csv (norm_post_no_trans, file="norm_post_no_trans.csv")


### model with female contribution excluded ###

m_norm_no_cont <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- bEX*ext_aff + bI*inher + bDC*dir_care + bSR*sex_rat + bMA*marr_arr + bDM*dec_mak + bPloc*p_loc + bBP*bride_pr + bGE*gift_ex + bD*dowry + k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- Imat*sigma,
    
    cutpoints ~ dnorm (0,1.5),
    
    c (bD,bGE,bBP,bPloc,bDM,bMA,bSR,bDC,bI,bEX) ~ dnorm (0,0.5),
    
    sigma ~ dexp (1)
    
  ), data=d_norm, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_no_cont, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bDC","bI","bEX","sigma"))

### get model summary

norm_no_cont <- precis (m_norm_no_cont, 2, omit=c("cutpoints","k","z"))
norm_no_cont <- data.frame (round (norm_no_cont, 3))
write.csv (data.frame (norm_no_cont), file="norm_no_cont_output.csv")

norm_post_no_cont <- extract.samples (m_norm_no_cont, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bDC","bI","bEX","sigma"))
write.csv (norm_post_no_cont, file="norm_post_no_cont.csv")


### model with parental control excluded

m_norm_no_marr <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- bEX*ext_aff + bI*inher + bDC*dir_care + bFC*fem_cont + bSR*sex_rat + bDM*dec_mak + bPloc*p_loc + bBP*bride_pr + bGE*gift_ex + bD*dowry + k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- Imat*sigma,
    
    cutpoints ~ dnorm (0,1.5),
    
    c (bD,bGE,bBP,bPloc,bDM,bSR,bFC,bDC,bI,bEX) ~ dnorm (0,0.5),
    
    sigma ~ dexp (1)
    
  ), data=d_norm, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_no_marr, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bSR","bFC","bDC","bI","bEX","sigma"))

### get model summary

norm_no_marr <- precis (m_norm_no_marr, 2, omit=c("cutpoints","k","z"))
norm_no_marr <- data.frame (round (norm_no_marr, 3))
write.csv (data.frame (norm_no_marr), file="norm_no_marr_output.csv")

norm_post_no_marr <- extract.samples (m_norm_no_marr, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bSR","bFC","bDC","bI","bEX","sigma"))
write.csv (norm_post_no_marr, file="norm_post_no_marr.csv")


####################################
####################################