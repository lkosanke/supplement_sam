

########################################
# Simulation Study 5
########################################

library(miceadds)
library(lavaan)

path <- "d:/_F/Projekte/SEM_2step"
pf1 <- file.path( path, "3factor")
pfc <- file.path( pf1, "Code")



stamp <- systime()[5]

name <- name0 <- "sam_3factor_cross_loadings"

# set.seed(998)

dfr <- NULL

lam1 <- .6; lam2 <- .65; lam3 <- .55
phi21 <- .6; phi31 <- .4; phi32 <- .2
phi21 <- .3; phi31 <- .3; phi32 <- .3
delta1 <- .4; delta2 <- .4; delta3 <- 0
phi21 <- 0.1; phi31 <- 0.3; phi32 <- 0.3

rho <- c(0, 0.1, .2, .3, .4, .5, .6)
delta <- c(0,.4)
des <- expand.grid( delta2=delta, delta3=delta, phi21=rho, phi31=rho, phi32=rho)
ND <- nrow(des)

for (dd in 1:ND){

    # dd <- 1
    phi21 <- des[dd,"phi21"]
    phi31 <- des[dd,"phi31"]
    phi32 <- des[dd,"phi32"]
    delta2 <- des[dd,"delta2"]
    delta3 <- des[dd,"delta3"]
        
    # dfr1
    
    
    
    LAM <- matrix(0, 9,3 )
    LAM[1:3,1] <- lam1
    LAM[4:6,2] <- lam2
    LAM[7:9,3] <- lam3
    LAM[1,2] <- delta1
    LAM[4,3] <- delta2
    LAM[7,1] <- delta3    
    PHI <- matrix(1,3,3)
    PHI[1,2] <- PHI[2,1] <- phi21
    PHI[1,3] <- PHI[3,1] <- phi31
    PHI[2,3] <- PHI[3,2] <- phi32
    
    TR <- LAM %*% PHI %*% t(LAM)
    PSI <- matrix(0, 9,9)
    diag(PSI) <- 1-diag(TR)
    
    S <- TR + PSI
    
    x0 <- rep(.5,3)
    
    
    source.all( pfc, "3factor")
    
    
    #**** SEM-ULS
    x <- c( rep(.5,9), .5, .5, .5)
    par <- names(x) <- c( paste0("lam" , 1:9), "phi21", "phi31", "phi32")
    true <- c( rep(lam1,3), rep(lam2,3), rep(lam3,3), phi21, phi31, phi32)
    upper <- .999 + 0*x    
    mod <- nlminb( start=x, objective=loss_fun_3factor, S=S, upper=upper)

    
    dfr1 <- data.frame( des_id=dd, lam1=lam1, lam2=lam2, lam3=lam3, phi21=phi21,
                phi31=phi31, phi32=phi32, delta1=delta1, delta2=delta2,
                delta3=delta3, par=par, true=true, est_sem=mod$par)
    
    
    #**** SAM
    ind <- 1:3; mod1 <- nlminb( start=x0, objective=loss_fun_1factor, S=S[ind,ind])
    ind <- 4:6; mod2 <- nlminb( start=x0, objective=loss_fun_1factor, S=S[ind,ind])
    ind <- 7:9; mod3 <- nlminb( start=x0, objective=loss_fun_1factor, S=S[ind,ind])
    
    x1 <- c( mod1$par, mod2$par, mod3$par)
    
    x <- rep(.5, 3)
    
    
    model <- "lsam_uls"
    mod <- nlminb( start=x, objective=loss_fun_3factor_lsam, S=S, x1=x1)
    dfr1$est_sam <- c(x1, mod$par)
    
    dfr <- rbind( dfr, dfr1 )
print(paste0(dd, " von ", ND));
}

dfr$bias_sem <- abs(dfr$est_sem - dfr$true)
dfr$bias_sam <- abs(dfr$est_sam - dfr$true)
# aggregate( dfr$bias_sem - dfr$bias_sam , list(dfr$par) , mean )

# define an indicator that SAM is more robust than SEM
bias_reduction <- .2
dfr$robustness <- "both no_bias"
dfr$robustness[ ( dfr$bias_sem > .05 ) & ( dfr$bias_sem > .05 ) ]  <- "both_bias"
dfr$robustness[ ( dfr$bias_sem > .05 ) & ( dfr$bias_sam / dfr$bias_sem < 1-bias_reduction ) ] <- "sam_better"
dfr$robustness[ ( dfr$bias_sam > .05 ) & ( dfr$bias_sem / dfr$bias_sam < 1-bias_reduction ) ] <- "sem_better"


dfr[ grep("lam", dfr$par), "robustness" ] <- NA

table(dfr$robustness)

save.data( dfr, paste0( name, "__RESULTS"), path=pf1, type="csv2", suffix=stamp, index=TRUE, systime=TRUE)


xtabs( ~ robustness + delta2 + delta3 , dfr )
c1 <- .25
xtabs( ~ robustness + delta2 + delta3 , dfr[ dfr$phi21 > c1 & dfr$phi31 > c1 & dfr$phi32 > c1 , ])
xtabs( ~ robustness + delta2 + delta3 , dfr[ dfr$phi21 < c1 & dfr$phi31 < c1 & dfr$phi32 < c1 , ])


# lavaan::lavParTable(mod1)

files_move(pfc)
files_move(pf1)

 ##  > xtabs( ~ robustness + delta2 + delta3 , dfr )
 ##  , , delta3 = 0
 ##  
 ##                delta2
 ##  robustness       0 0.4
 ##    both no_bias 488 149
 ##    both_bias    394 751
 ##    sam_better    97  69
 ##    sem_better    50  60
 ##  
 ##  , , delta3 = 0.4
 ##  
 ##                delta2
 ##  robustness       0 0.4
 ##    both no_bias 137   0
 ##    both_bias    698 978
 ##    sam_better    72  18
 ##    sem_better   122  33
 ##  
 ##  > c1 <- .25
 ##  > xtabs( ~ robustness + delta2 + delta3 , dfr[ dfr$phi21 > c1 & dfr$phi31 > c1 & dfr$phi32 > c1 , ])
 ##  , , delta3 = 0
 ##  
 ##                delta2
 ##  robustness       0 0.4
 ##    both no_bias  91  27
 ##    both_bias     76 140
 ##    sem_better    25  25
 ##  
 ##  , , delta3 = 0.4
 ##  
 ##                delta2
 ##  robustness       0 0.4
 ##    both no_bias  22   0
 ##    both_bias    123 186
 ##    sem_better    47   6
 ##  
 ##  > xtabs( ~ robustness + delta2 + delta3 , dfr[ dfr$phi21 < c1 & dfr$phi31 < c1 & dfr$phi32 < c1 , ])
 ##  , , delta3 = 0
 ##  
 ##                delta2
 ##  robustness      0 0.4
 ##    both no_bias 51  24
 ##    both_bias    18  46
 ##    sam_better   12  11
 ##  
 ##  , , delta3 = 0.4
 ##  
 ##                delta2
 ##  robustness      0 0.4
 ##    both no_bias 24   0
 ##    both_bias    47  77
 ##    sam_better   10   4
 ##  
