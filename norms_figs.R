################################################################################################################################
### R script: Paternity uncertainty and parent-offspring conflict explain restrictions on female premarital sex across societies
################################################################################################################################


### code for reproducing the figures

setwd ("C:/R_folder/sex norms/norms resubmission/145")

library (tidyr)
library (grid)
library (ggplot2)
library (maps)
library (mapdata)
library (ggmap)
library (ggridges)
library (ggpubr)
library (grid)
library (gridExtra)
library (rethinking)


################
### Figure 1 ###

norms <- read.csv("norms_final.csv", header=TRUE)
norms <- subset(norms, sex_norms!="2" & sex_norms!="6") # remove societies with trial and early marriage (N=128)
oldvals <- c(1,3,4,5)
newvals <- as.integer (c(1,2,3,4))
norms$sex_norms <- newvals[match(norms$sex_norms, oldvals)]
norms$sex_norms <- as.factor(norms$sex_norms)

world <- map_data ("world")
labs_norm <- c ("Strongly sanctioned","Weakly sanctioned","Sanctions if pregnancy","No sanctions")
cols_norm <- c ("darkorange4","chocolate","goldenrod","yellow2")

map_norm <- ggplot() + geom_polygon(data=world, aes(x=long, y=lat, group=group), fill="grey85", color="grey85") + theme_classic() + labs (x="", y="") + coord_fixed(1.3) + geom_point(data=norms, aes(x=long, y=lat, colour=as.factor(norms[,4])), size=3) + ylim(-65,85) + theme (legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=16), axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), text=element_text(size=16), plot.margin=unit(c(0,1,0,0),"cm")) + scale_fill_manual (labels=labs_norm, values=cols_norm) + scale_colour_manual (labels=labs_norm, values=cols_norm)
map_norm

dev.off()


################
### Figure 2 ###

{
  
### phylogenetic covariance plot (Figure 2a)

norm_post <- read.csv ("norm_post_phylo.csv", header=TRUE)

layout.matrix <- matrix(c(1,0), nrow=1, ncol=2) # set the layout
layout(layout.matrix, widths=c(1,2))
par(mar=c(5.5,5.25,4.5,4.15))

plot (NULL, xlim=c(0,1), ylim=c(0,5), xlab="", ylab="", main="", xaxt="n", cex.axis=1.35)
axis (1, at=seq(from=0, to=1, length.out=3), labels=c("0.0","0.5","1.0"), cex.axis=1.35)
mtext("a)", side=3, line=0.25, cex=1.35, adj=0)
mtext ("covariance", 2, line=2.75, cex=1.5)
mtext ("phylogenetic distance", 1, line=2.75, cex=1.35)

x_seq <- seq (from=0, to=10, length.out=100)
pmcov <- sapply (x_seq, function(x) norm_post$etasq*exp(-norm_post$rhosq*x))
pmcov_mu <- apply (pmcov, 2, mean)

for (i in 1:100)
  curve (norm_post$etasq[i]*exp(-norm_post$rhosq[i]*x), add=TRUE, col=col.alpha("gray55",0.35))
lines (x_seq, pmcov_mu, lwd=3)


### density plots (Figure 2b)

norm_post_g <- read.csv ("norm_post_g.csv", header=TRUE)
norm_post_g <- norm_post_g[,5:15]
names(norm_post_g) <- c("dowry","gift exchange","bride-price","patrilocality","decision \nmakers","parental \ncontrol","sex ratio","female \ncontribution","direct care","inheritance","extramarital sex")

norm_mean <- apply(norm_post_g, 2, mean)

cols <- c("turquoise2","turquoise3","turquoise4","#663300","orange4","#993300","gold","yellow2","magenta4","magenta2","deeppink4")
names <- c("dowry","gift exchange","bride-price","patrilocality","decision \nmakers","parental \ncontrol","sex ratio","female \ncontribution","inheritance","direct care","extramarital sex")

norm_post_g <- gather(norm_post_g, key="var", value="est")
norm_post_g$var <- factor(norm_post_g$var, levels=names(norm_mean))
norm <- ggplot(norm_post_g, aes(x=est, y=var)) + labs(y="", x="posterior estimate") + geom_vline(xintercept=0, alpha=1, linetype="longdash") + geom_hline(yintercept=1:11, alpha=1, colour="gray85") + geom_density_ridges2(aes(fill=var, color=var), lwd=0.8, scale=2.5, rel_min_height=0.01, alpha=0.6) + scale_fill_manual(values=cols, labels=names) + scale_color_manual(values=cols, labels=names) + theme_classic(base_size=12) + theme(legend.position="none",axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), text=element_text(size=15),  plot.margin=unit(c(1.5,1,1.5,0.5),"cm")) + ggtitle ("b)") + scale_y_discrete (limits=names, labels=names, expand=c(0.05,0.0)) + coord_cartesian (clip="off") + annotate("text", x=2.6, y=1.3:11.3, colour=c("turquoise2","turquoise3","turquoise4","#663300","orange4","#993300","gold","yellow2","magenta2","magenta4","deeppink4"), size=6, fontface="bold", label=c("41.33","79.42","15.82","70.23","88.08","99.72","99.38","95.55","89.90","23.58","99.42"))

vp_norm <- viewport (height=unit(1,"npc"), width=unit(0.7,"npc"), just=c("left","top"), y=1, x=0.3)
print (norm, vp=vp_norm)

}

dev.off ()


################
### Figure 3 ###

norm_post_g <- read.csv ("norm_post_g.csv", header=TRUE)

seq_x <- seq (from=-2, to=2, length.out=30) # define a sequence for x axis
seq_f <- seq (from=0, to=1, length.out=30) # same for categorical predictors

### reverse the sign of the linear model inside the original 'pordlogit' function for more intuitive plots

pordlogit_rev <- function (x, phi, a, log=FALSE) 
{
  a <- c(as.numeric(a), Inf)
  if (length(phi)==1) {
    p <- logistic(a[x] + phi) # here
  }
  else {
    p <- matrix(NA, ncol=length(x), nrow=length(phi))
    for (i in 1:length(phi)) {
      p[i, ] <- logistic(a[x] + phi[i]) # and here
    }
  }
  if (log==TRUE) 
    p <- log(p)
  p
}


### posterior predictive plots

{
  
  par (mfrow=c(3,4), mar = c(5.15,5.55,3.5,4.15))
  
  ### extramarital sex
  
  phi_link <- function (ext_aff) norm_post_g$bEX*ext_aff # linear function
  phi <- sapply (seq_x, phi_link) # 
  
  plot (NULL, type="n", xlab="extramarital sex (std)", ylab="", xlim=c(-2,2), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=seq(-2,2,length.out = 5), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  mtext ("probability", 2, cex=1.25, line=2.75)
  
  ext_aff <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    ext_aff[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(ext_aff,2,PI), seq_x, col=col.alpha(cols[11],0.3))
  lines(x=seq_x, y=apply(ext_aff,2,mean), col=cols[11], lwd=2)
  
  ext_aff <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    ext_aff[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(ext_aff,2,PI), seq_x, col=col.alpha(cols[11],0.3))
  lines(x=seq_x, y=apply(ext_aff,2,mean), col=cols[11], lwd=2)
  
  ext_aff <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    ext_aff[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(ext_aff,2,PI), seq_x, col=col.alpha(cols[11],0.3))
  lines(x=seq_x, y=apply(ext_aff,2,mean), col=cols[11], lwd=2)
  
  
  ### direct care
  
  phi_link <- function (dir_care) norm_post_g$bDC*dir_care # linear function
  phi <- sapply (seq_x, phi_link) # 
  
  plot (NULL, type="n", xlab="direct care (std)", ylab="", xlim=c(-2,2), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=seq(-2,2,length.out = 5), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  dir_care <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    dir_care[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(dir_care,2,PI), seq_x, col=col.alpha(cols[9],0.3))
  lines(x=seq_x, y=apply(dir_care,2,mean), col=cols[9], lwd=2)
  
  dir_care <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    dir_care[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(dir_care,2,PI), seq_x, col=col.alpha(cols[9],0.3))
  lines(x=seq_x, y=apply(dir_care,2,mean), col=cols[9], lwd=2)
  
  dir_care <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    dir_care[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(dir_care,2,PI), seq_x, col=col.alpha(cols[9],0.3))
  lines(x=seq_x, y=apply(dir_care,2,mean), col=cols[9], lwd=2)
  
  
  ### inheritance
  
  phi_link <- function (inher) norm_post_g$bI*inher
  phi <- sapply (seq_f, phi_link)  
  
  plot (NULL, type="n", xlab="inheritance", ylab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=c(0.05,0.95), labels=c("absent", "present"), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  inher <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    inher[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(inher,2,PI), seq_f, col=col.alpha(cols[10],0.3))
  lines(x=seq_f, y=apply(inher,2,mean), col=cols[10], lwd=2)
  
  inher <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    inher[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(inher,2,PI), seq_f, col=col.alpha(cols[10],0.3))
  lines(x=seq_f, y=apply(inher,2,mean), col=cols[10], lwd=2)
  
  inher <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    inher[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(inher,2,PI), seq_f, col=col.alpha(cols[10],0.3))
  lines(x=seq_f, y=apply(inher,2,mean), col=cols[10], lwd=2)
  
  
  ### female contribution
  
  phi_link <- function (fem_cont) norm_post_g$bFC*fem_cont # linear function
  phi <- sapply (seq_x, phi_link) # 
  
  plot (NULL, type="n", xlab="female contribution (std)", ylab="", xlim=c(-2,2), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=seq(-2,2,length.out = 5), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  fem_cont <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    fem_cont[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(fem_cont,2,PI), seq_x, col=col.alpha(cols[8],0.3))
  lines(x=seq_x, y=apply(fem_cont,2,mean), col=cols[8], lwd=2)
  
  fem_cont <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    fem_cont[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(fem_cont,2,PI), seq_x, col=col.alpha(cols[8],0.3))
  lines(x=seq_x, y=apply(fem_cont,2,mean), col=cols[8], lwd=2)
  
  fem_cont <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    fem_cont[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(fem_cont,2,PI), seq_x, col=col.alpha(cols[8],0.3))
  lines(x=seq_x, y=apply(fem_cont,2,mean), col=cols[8], lwd=2)
  
  
  ### sex ratio
  
  phi_link <- function (sex_rat) norm_post_g$bSR*sex_rat # linear function
  phi <- sapply (seq_x, phi_link) # 
  
  plot (NULL, type="n", xlab="sex ratio (std)", ylab="", xlim=c(-2,2), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=seq(-2,2,length.out = 5), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  mtext ("probability", 2, cex=1.25, line=2.75)
  
  sex_rat <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    sex_rat[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(sex_rat,2,PI), seq_x, col=col.alpha(cols[7],0.3))
  lines(x=seq_x, y=apply(sex_rat,2,mean), col=cols[7], lwd=2)
  
  sex_rat <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    sex_rat[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(sex_rat,2,PI), seq_x, col=col.alpha(cols[7],0.3))
  lines(x=seq_x, y=apply(sex_rat,2,mean), col=cols[7], lwd=2)
  
  sex_rat <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    sex_rat[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(sex_rat,2,PI), seq_x, col=col.alpha(cols[7],0.3))
  lines(x=seq_x, y=apply(sex_rat,2,mean), col=cols[7], lwd=2)
  
  
  ### parental control
  
  phi_link <- function (marr_arr) norm_post_g$bMA*marr_arr # linear function
  phi <- sapply (seq_x, phi_link) # 
  
  plot (NULL, type="n", xlab="parental control (std)", ylab="", xlim=c(-2,2), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=seq(-2,2,length.out = 5), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  marr_arr <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    marr_arr[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(marr_arr,2,PI), seq_x, col=col.alpha(cols[6],0.3))
  lines(x=seq_x, y=apply(marr_arr,2,mean), col=cols[6], lwd=2)
  
  marr_arr <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    marr_arr[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(marr_arr,2,PI), seq_x, col=col.alpha(cols[6],0.3))
  lines(x=seq_x, y=apply(marr_arr,2,mean), col=cols[6], lwd=2)
  
  marr_arr <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    marr_arr[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(marr_arr,2,PI), seq_x, col=col.alpha(cols[6],0.3))
  lines(x=seq_x, y=apply(marr_arr,2,mean), col=cols[6], lwd=2)
  
  
  ### decision makers
  
  phi_link <- function (dec_mak) norm_post_g$bDM*dec_mak 
  phi <- sapply (seq_x, phi_link) # 
  
  plot (NULL, type="n", xlab="decision makers (std)", ylab="", xlim=c(-2,2), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=seq(-2,2,length.out = 5), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  dec_mak <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    dec_mak[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(dec_mak,2,PI), seq_x, col=col.alpha(cols[5],0.3))
  lines(x=seq_x, y=apply(dec_mak,2,mean), col=cols[5], lwd=2)
  
  dec_mak <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    dec_mak[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(dec_mak,2,PI), seq_x, col=col.alpha(cols[5],0.3))
  lines(x=seq_x, y=apply(dec_mak,2,mean), col=cols[5], lwd=2)
  
  dec_mak <- matrix(0, nrow=100, ncol=length(seq_x))
  for (i in 1:100) {
    dec_mak[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(dec_mak,2,PI), seq_x, col=col.alpha(cols[5],0.3))
  lines(x=seq_x, y=apply(dec_mak,2,mean), col=cols[5], lwd=2)
  
  
  ### patrilocality
  
  phi_link <- function (p_loc) norm_post_g$bPloc*p_loc
  phi <- sapply (seq_f, phi_link)  
  
  plot (NULL, type="n", xlab="patrilocality", ylab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=c(0.05,0.95), labels=c("absent", "present"), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  p_loc <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    p_loc[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(p_loc,2,PI), seq_f, col=col.alpha(cols[4],0.3))
  lines(x=seq_f, y=apply(p_loc,2,mean), col=cols[4], lwd=2)
  
  p_loc <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    p_loc[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(p_loc,2,PI), seq_f, col=col.alpha(cols[4],0.3))
  lines(x=seq_f, y=apply(p_loc,2,mean), col=cols[4], lwd=2)
  
  p_loc <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    p_loc[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(p_loc,2,PI), seq_f, col=col.alpha(cols[4],0.3))
  lines(x=seq_f, y=apply(p_loc,2,mean), col=cols[4], lwd=2)
  
  
  ### bride-price
  
  phi_link <- function (bride_pr) norm_post_g$bBP*bride_pr # linear function
  phi <- sapply (seq_f, phi_link) # 
  
  plot (NULL, type="n", xlab="bride-price", ylab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=c(0.05,0.95), labels=c("absent", "present"), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  mtext ("probability", 2, cex=1.25, line=2.75)
  
  bride_pr <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    bride_pr[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(bride_pr,2,PI), seq_f, col=col.alpha(cols[3],0.3))
  lines(x=seq_f, y=apply(bride_pr,2,mean), col=cols[3], lwd=2)
  
  bride_pr <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    bride_pr[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(bride_pr,2,PI), seq_f, col=col.alpha(cols[3],0.3))
  lines(x=seq_f, y=apply(bride_pr,2,mean), col=cols[3], lwd=2)
  
  bride_pr <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    bride_pr[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(bride_pr,2,PI), seq_f, col=col.alpha(cols[3],0.3))
  lines(x=seq_f, y=apply(bride_pr,2,mean), col=cols[3], lwd=2)
  
  
  ### gift exchange
  
  phi_link <- function (gift_ex) norm_post_g$bGE*gift_ex # linear function
  phi <- sapply (seq_f, phi_link) # 
  
  plot (NULL, type="n", xlab="gift exchange", ylab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=c(0.05,0.95), labels=c("absent", "present"), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  gift_ex <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    gift_ex[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(gift_ex,2,PI), seq_f, col=col.alpha(cols[2],0.3))
  lines(x=seq_f, y=apply(gift_ex,2,mean), col=cols[2], lwd=2)
  
  gift_ex <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    gift_ex[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(gift_ex,2,PI), seq_f, col=col.alpha(cols[2],0.3))
  lines(x=seq_f, y=apply(gift_ex,2,mean), col=cols[2], lwd=2)
  
  gift_ex <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    gift_ex[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(gift_ex,2,PI), seq_f, col=col.alpha(cols[2],0.3))
  lines(x=seq_f, y=apply(gift_ex,2,mean), col=cols[2], lwd=2)
  
  
  ### dowry
  
  phi_link <- function (dowry) norm_post_g$bD*dowry 
  phi <- sapply (seq_f, phi_link) 
  
  plot (NULL, type="n", xlab="dowry", ylab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n", cex.lab=1.75)
  axis(1, at=c(0.05,0.95), labels=c("absent", "present"), cex.axis=1.75)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.75)
  
  dowry <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    dowry[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,2])
  }
  
  shade(apply(dowry,2,PI), seq_f, col=col.alpha(cols[1],0.3))
  lines(x=seq_f, y=apply(dowry,2,mean), col=cols[1], lwd=2)
  
  dowry <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    dowry[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,3])
  }
  
  shade(apply(dowry,2,PI), seq_f, col=col.alpha(cols[1],0.3))
  lines(x=seq_f, y=apply(dowry,2,mean), col=cols[1], lwd=2)
  
  dowry <- matrix(0, nrow=100, ncol=length(seq_f))
  for (i in 1:100) {
    dowry[i,] <- pordlogit_rev(1, phi[i,], norm_post_g[i,4])
  }
  
  shade(apply(dowry,2,PI), seq_f, col=col.alpha(cols[1],0.3))
  lines(x=seq_f, y=apply(dowry,2,mean), col=cols[1], lwd=2)
  
}


dev.off()


################################################################################
################################################################################