##########################################################################################################################################
### R script: Paternity certainty and parent-offspring conflict jointly explain restrictions on female premarital sex across societies ###
##########################################################################################################################################


###     

setwd ("F:/sex norms")

library (missForest)
library (rethinking)
library (ape)
library (geiger)


####################################
### imputation of missing values ###

norms <- read.csv ("sex_norms_data.csv", header=TRUE) # load the full non-imputed sample
str (norms)

norms2 <- norms[complete.cases(norms[, 4]), ] # remove societies with missing data for sex norms; it reduces to 145 societies
pop <- norms2$Population # store population column

norms2[c(18,27,29,30,32,33,35,36,43,44,45,46,49,53,54,55,119,143),12] <- 7 # change to indirect dowry according to Schlegel (1991)
norms2 <- norms2[,-c(1:3)] # remove population column and GPS coordinates

# convert from integers to factors, female contribution to numeric

cols <- c (1:6,8:12)
norms2[cols] <- lapply (norms2[cols], factor) # convert to factors
norms2$fem_cont <- as.numeric (norms2$fem_cont) # convert to numeric

norms_imp <- missForest (norms2) # impute
norms_imputed <- norms_imp$ximp

norms_imputed_pop <- cbind (pop, norms_imputed) # bind the population column with the imputed data
write.csv (norms_imputed_pop, file="norms_red_imputed.csv") # save the imputed data


######################
### data re-coding ###

norms <- read.csv ("norms_red_imputed.csv", header=TRUE) # load the imputed reduced sample
d_norms <- norms[,-1] # remove row index
str (d_norms)

# sex norms

d_norms$sex_norms <- as.factor (d_norms$sex_norms)
levels(d_norms$sex_norms) <- list ("1"="6", "2"="5", "3"="4", "4"="3", "5"="2", "6"="1") # invert the scale

# inheritance land

d_norms$inher_land <- as.factor (d_norms$inher_land)
levels(d_norms$inher_land) <- list ("1"=c("4","5","7"), "0"=c("1","2","3","6"))

# inheritance movable

d_norms$inher_mov <- as.factor (d_norms$inher_mov)
levels(d_norms$inher_mov) <- list ("1"=c("4","5","7"), "0"=c("1","2","3","6"))

# combine inheritance scores into a single predictor

d_norms$inher <- as.factor(paste(d_norms$inher_land, d_norms$inher_mov, sep=""))
levels(d_norms$inher) <- list ("1"="11", "0"=c("00","01","10"))

# direct care

d_norms$dir_care <- d_norms$dir_care_inf + d_norms$dir_care_child

# patrilocality

d_norms$p_loc <- as.factor (d_norms$trans_resid) 
levels (d_norms$p_loc) <- list ("1"="1", "0"=c("2","3","4"))

# bride-price & token bride-price; see Schlegel (1991)

d_norms$bride_pr <- as.factor (d_norms$marr_tr)
levels(d_norms$bride_pr) <- list ("1"=c("0","2"), "0"=c("1","3","4","5","6","7")) 

# gift exchange

d_norms$gift_ex <- as.factor (d_norms$marr_tr)
levels(d_norms$gift_ex) <- list ("1"="3", "0"=c("1","2","4","5","6","7"))

# dowry & indirect dowry

d_norms$dowry <- as.factor (d_norms$marr_tr)
levels(d_norms$dowry) <- list ("1"=c("6","7"), "0"=c("0","1","2","3","4","5"))

str (d_norms)
norms_data <- data.frame (d_norms)
write.csv (norms_data, "norms_red_final.csv")


#####################
### model fitting ###

### sex attitudes

memory.limit (size=200000) # allocate enough memory to R in order to avoid crashing

norms <- read.csv ("norms_red_final.csv", header=TRUE, row.names=2)
pop_i <- norms$X
N_pop <- nrow (norms)

# standardize the variables by centering and dividing by 2 SD; see Gelman (2008)

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

tree <- read.tree ("SPT.SCCS.tre") # load the phylogeny
name.check (tree, norms) # check the correspondence between tips of the phylogeny and societies in the data

# drop tips that are not in the data

tree <- drop.tip (tree, tip=c("Ajie","Badjau","Bambara","Banen","Bogo","Cayapa","Cayua","Comanche",
                                          "Cubeo_Tucano","Eyak","Hadza","Havasupai","Hidatsa","Huichol","Kenuzi_Nubians",
                                          "Khalka_Mongols","Kikuyu","Lengua","Lozi","Manchu","Massa_Masa","Mbau_Fijians",
                                          "Mbundu","Mende","New_Ireland","Papago","Pawnee","Pomo_Eastern","Popoluca",
                                          "Rhade","Romans","Saramacca","Saulteaux","Suku","Tobelorese","Toradja","Warrau",
                                          "Yokuts_Lake","Yukaghir","Yurok","Zuni"          
))

name.check(tree,norms) # OK

write.tree (tree,"SPT.SCCS_red.tre")
tree_red <- read.tree ("SPT.SCCS_red.tre") # load the phylogeny
tree_trimmed <- keep.tip (tree_red, pop_i) # combine tips of the phylogeny with societies in the sample
Dmat <- cophenetic (tree_trimmed) # compute distance matrix


### 'phylogeny' model; varying intercepts, no predictors (OU, de-centered) ### 

d1_norm <- list (sex_norms=sex_norms,
                 N_pop=N_pop,
                 society=1:N_pop
)
d1_norm$Dmat <- Dmat[pop_i,pop_i]/max(Dmat) # normalize the distances to scale from 0 to 1 and include it to the data list

m_norm_phylo <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- cov_GPL1 (Dmat, etasq, rhosq, 0.01),
    
    cutpoints ~ dnorm (0,1.5),
    
    etasq ~ dexp (1),
    rhosq ~ dexp (1)
    
  ), data=d1_norm, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_phylo, pars=c("cutpoints", "k")) # check if chains mixed well


### 'predictors' model (OU, de-centered) ###

d2_norm_g <- list (sex_norms=sex_norms,
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

m_norm_pred_g <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- bEX*ext_aff + bI*inher + bDC*dir_care + bFC*fem_cont + bSR*sex_rat + bMA*marr_arr + bDM*dec_mak + bPloc*p_loc + bBP*bride_pr + bGE*gift_ex + bD*dowry + k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- Imat*sigma,
    
    cutpoints ~ dnorm (0,1.5),
    
    c (bD,bGE,bBP,bPloc,bDM,bMA,bSR,bFC,bDC,bI,bEX) ~ dnorm (0,0.5),
    
    sigma ~ dexp (1)
    
  ), data=d2_norm_g, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_pred_g, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bFC","bDC","bI","bEX","sigma"))


### 'phylogeny + predictors' model (OU, de-centered) ###

d3_norm <- list (sex_norms=sex_norms,
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
                 Dmat=Dmat
)

m_norm_phylo_pred <- ulam(
  alist(
    sex_norms ~ dordlogit (phi, cutpoints),
    
    phi <- bEX*ext_aff + bI*inher + bDC*dir_care + bFC*fem_cont + bSR*sex_rat + bMA*marr_arr + bDM*dec_mak + bPloc*p_loc + bBP*bride_pr + bGE*gift_ex + bD*dowry + k[society],
    
    transpars> vector[N_pop]:k <<- L_SIGMA*z,
    vector [N_pop]:z ~ normal (0,1),
    
    transpars> matrix[N_pop,N_pop]:L_SIGMA <<- cholesky_decompose (SIGMA),
    transpars> matrix[N_pop,N_pop]:SIGMA <- cov_GPL1 (Dmat, etasq, rhosq, 0.01),
    
    cutpoints ~ dnorm (0,1.5),
    
    c (bD,bGE,bBP,bPloc,bDM,bMA,bSR,bFC,bDC,bI,bEX) ~ dnorm (0,0.5),
    
    etasq ~ dexp (1),
    rhosq ~ dexp (1)
    
  ), data=d3_norm, chains=4, cores=4, iter=2000, log_lik=TRUE)

traceplot (m_norm_phylo_pred, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bFC","bDC","bI","bEX","etasq","rhosq"))


### get posteriors ###

norm_model_comp <- compare (m_norm_phylo, m_norm_pred_g, m_norm_phylo_pred, func=WAIC) # model comparison
norm_model_comp <- data.frame (round (norm_model_comp, 3))
write.csv (norm_model_comp, file="norm_model_comp.csv")

norm_prec_g <- precis (m_norm_pred_g, 2, omit=c("cutpoints","k","z")) # 'predictors' model summary
norm_prec_g <- data.frame (round (norm_prec_g, 3))
write.csv (data.frame (norm_prec_g), file="norm_g_model_output.csv")

norm_prec_phylo <- precis (m_norm_phylo_pred, 2, omit=c("cutpoints","k","z")) # 'predictors + phylogeny' model summary
norm_prec_phylo <- data.frame (round (norm_prec_phylo, 3))
write.csv (data.frame (norm_prec_phylo), file="norm_phylo_model_output.csv")

norm_post_g <- extract.samples (m_norm_pred_g, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bFC","bDC","bI","bEX","k","bEX","sigma")) # 'predictors' model posterior
write.csv (norm_post_g, file="norm_post_g.csv") 

norm_post_phylo <- extract.samples (m_norm_phylo_pred, pars=c("cutpoints","bD","bGE","bBP","bPloc","bDM","bMA","bSR","bFC","bDC","bI","bEX","k","etasq","rhosq")) # 'predictors + phylogeny' model posterior
write.csv (norm_post_phylo, file="norm_post_phylo.csv") 


######################
### intervals & PP ###

### model without phylogeny

norm_post_g <- read.csv ("norm_post_g.csv", header=TRUE)

# compute 89% HPD intervals

norm_post_int <- c (round (HPDI (norm_post_g[,17], prob=0.89),2),
                    round (HPDI (norm_post_g[,16], prob=0.89),2),
                    round (HPDI (norm_post_g[,15], prob=0.89),2),
                    round (HPDI (norm_post_g[,14], prob=0.89),2), 
                    round (HPDI (norm_post_g[,13], prob=0.89),2),
                    round (HPDI (norm_post_g[,12], prob=0.89),2), 
                    round (HPDI (norm_post_g[,11], prob=0.89),2), 
                    round (HPDI (norm_post_g[,10], prob=0.89),2),
                    round (HPDI (norm_post_g[,9], prob=0.89),2),
                    round (HPDI (norm_post_g[,8], prob=0.89),2),
                    round (HPDI (norm_post_g[,7], prob=0.89),2)
                    
)

# compute the amount of posterior probability on the expected side of zero (PP)

norm_post_prob <- c (round ((sum (norm_post_g[,17] <0)/nrow(norm_post_g))*100,2),
                     round ((sum (norm_post_g[,16] >0)/nrow(norm_post_g))*100,2),
                     round ((sum (norm_post_g[,15] >0)/nrow(norm_post_g))*100,2),
                     round ((sum (norm_post_g[,14] <0)/nrow(norm_post_g))*100,2), 
                     round ((sum (norm_post_g[,13] >0)/nrow(norm_post_g))*100,2), 
                     round ((sum (norm_post_g[,12] >0)/nrow(norm_post_g))*100,2), 
                     round ((sum (norm_post_g[,11] <0)/nrow(norm_post_g))*100,2), 
                     round ((sum (norm_post_g[,10] >0)/nrow(norm_post_g))*100,2),
                     round ((sum (norm_post_g[,9] >0)/nrow(norm_post_g))*100,2),
                     round ((sum (norm_post_g[,8] >0)/nrow(norm_post_g))*100,2),
                     round ((sum (norm_post_g[,7] >0)/nrow(norm_post_g))*100,2)
                     
)

norm_post_int  <- matrix(norm_post_int, ncol=2, byrow=TRUE)
norm_post_prob <- matrix(norm_post_prob, ncol=1, byrow=TRUE)
norm_int_prob <- data.frame (cbind (norm_post_int,norm_post_prob))
colnames (norm_int_prob) <- c ("5.5%","94.5%","PP")
write.csv (norm_int_prob, file="norm_int_prob.csv")


### phylogenetic model

norm_post_phylo <- read.csv ("norm_post_phylo.csv", header=TRUE)

# compute 89% HPD intervals

norm_post_int <- c (round (HPDI (norm_post_phylo[,17], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,16], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,15], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,14], prob=0.89),2), 
                    round (HPDI (norm_post_phylo[,13], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,12], prob=0.89),2), 
                    round (HPDI (norm_post_phylo[,11], prob=0.89),2), 
                    round (HPDI (norm_post_phylo[,10], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,9], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,8], prob=0.89),2),
                    round (HPDI (norm_post_phylo[,7], prob=0.89),2)
                    
)

# compute the amount of posterior probability on the expected side of zero (PP)

norm_post_prob <- c (round ((sum (norm_post_phylo[,17] <0)/nrow(norm_post_phylo))*100,2),
                     round ((sum (norm_post_phylo[,16] >0)/nrow(norm_post_phylo))*100,2),
                     round ((sum (norm_post_phylo[,15] >0)/nrow(norm_post_phylo))*100,2),
                     round ((sum (norm_post_phylo[,14] <0)/nrow(norm_post_phylo))*100,2), 
                     round ((sum (norm_post_phylo[,13] >0)/nrow(norm_post_phylo))*100,2), 
                     round ((sum (norm_post_phylo[,12] >0)/nrow(norm_post_phylo))*100,2), 
                     round ((sum (norm_post_phylo[,11] <0)/nrow(norm_post_phylo))*100,2), 
                     round ((sum (norm_post_phylo[,10] >0)/nrow(norm_post_phylo))*100,2),
                     round ((sum (norm_post_phylo[,9] >0)/nrow(norm_post_phylo))*100,2),
                     round ((sum (norm_post_phylo[,8] >0)/nrow(norm_post_phylo))*100,2),
                     round ((sum (norm_post_phylo[,7] >0)/nrow(norm_post_phylo))*100,2)
                     
)

norm_post_int  <- matrix(norm_post_int, ncol=2, byrow=TRUE)
norm_post_prob <- matrix(norm_post_prob, ncol=1, byrow=TRUE)
norm_int_prob <- data.frame (cbind (norm_post_int,norm_post_prob))
colnames (norm_int_prob) <- c ("5.5%","94.5%","PP")
write.csv (norm_int_prob, file="norm_int_prob_phylo.csv")


####################################
####################################