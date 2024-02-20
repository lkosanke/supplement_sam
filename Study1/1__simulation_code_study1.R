########################################
# Simulation Study 1
########################################

library(GPArotation)
library(CDM)
library(miceadds)
library(TAM)
library(sirt)
library(lavaan)


comp <- "1001"
# pf00 <- "d:/_F/Projekte/SEM_2step/"
pf00 <- "p:/Eigene_Projekte/SEM_2step/Analysen/"

pfR <- file.path( pf00, "R-Code")
files_move(pfR); source.all(pfR, "sam")


# simulation number
sim <- "212"

# path
pfo <- paste0( pf00 , "Simulation" , sim , "\\1__Durchfuehrung" )
pf14 <- paste0( pf00 , "Simulation" , sim )

files_move(pf14)



# source R functions
# source.all( file.path(pf00, "_R-Code" ), "link")

comp_index <- systime()[8]
comp_index <- paste0( "sim", sim, "_C" , comp , "_", comp_index )



# set.seed( round( as.numeric( Sys.time() ) %% 1E6 ) )

pf1 <- file.path( pfo , "Outputs" , comp_index )
dir.create(pf1)


# powers
pow <- c(2,1,.5, .25)

inf_val <- 1e100

des_list <- list( 
    N = c(100,250, 500, 1000, 2500, 100000), 
    delta_fac = c( .15, 0),
    rc = c(1,2)    
       )

design <- expand.grid(des_list)



# eliminations

design_elim <- function(design, elim)
{
    if (length(elim)>0){ design <- design[ - elim , ] }
    return(design)
}

# elim <- which( ( design$prop_dif > 0 ) & ( design$dif == 0 ) )
# design <- design_elim(design, elim)



rownames(design) <- NULL

ND <- nrow(design)





BMax <- 1000
bb <- 1

des_order <- 1:ND

# ind <- c(41,42)
# des_order <- c(ind, setdiff(1:ND, ind))

# des_order <- seq(ND, 1, -1)


write.csv2( design , file.path( pf1 ,"_DESIGN.csv" ) )





bb <- 1
dd <- 1
# dd <- 14; dd <- 3
dd <- 6
dd <- 18
# dd <- 24




for (bb in 1:BMax){

for (dd in des_order){

# dd <- 12
# dd <- 1

filename <- systime()[7]
filename <- paste0( comp_index, "__", filename , "_DES" , dd , "_REPL" , bb )
vars <- colnames(design)

for (vv in vars){
    Revalpr( paste0( vv , " <- design[dd,'", vv, "']") )
}
design_dd <- data.frame(desid=dd, design[dd,] )



#-----------------------------------------
#---- data generation --------------------


lam1 <- .55
lam2 <- .45
phi <- .60

LAM <- matrix(0, nrow=6, ncol=2)
LAM[1:3,1] <- lam1
LAM[4:6,2] <- lam2
PHI <- matrix(phi, nrow=2, ncol=2)
diag(PHI) <- 1

THETA <- diag( c( rep(1-lam1^2,3), rep(1-lam2^2,3) ))

delta <- round( (1-lam2^2)*delta_fac, 3 )
I0 <- min(3,rc)
for (ii in 1:I0){
    THETA[ii,ii+3] <- THETA[3+ii,ii] <- delta
}


if (rc>3){
for (ii in 1:( min(2,(rc-3)) ) ){
    THETA[ii,ii+4] <- THETA[4+ii,ii] <- delta
}
}
if (rc==6){
   ii <- 1; jj <- 4
       THETA[ii,jj] <- THETA[jj,ii] <- delta
}

# THETA[1,1] <- 1 ; THETA[1,2] <- THETA[2,1] <- .4

S0 <- S <- LAM %*% PHI %*% t(LAM) + THETA
rownames(S) <- colnames(S) <- c( paste0("X",1:3), paste0("Y",1:3) )

if (N<inf_val){
    dat <- MASS::mvrnorm(n=N, mu=rep(0,6), Sigma=S)
    S <- cov(dat) * (N-1) / N
}


# save example datasets
save.data( S , paste0("DATA_DES",dd ) , type="csv2" , path= pf1, row.names=TRUE)



dfr <- NULL


#*** estimates
lavmodel1 <- "
    FX =~ l1*X1 + l2*X2 + l3*X3
    X1 ~~ vX1*X1
    X2 ~~ vX2*X2    
    X3 ~~ vX3*X3
    FX ~~ 1*FX
    FY =~ l4*Y1 + l5*Y2 + l6*Y3
    FY ~~ 1*FY
    FX ~~ phi*FY    
    l1 > 0.01
    l2 > 0.01    
    l3 > 0.01
    l4 > 0.01
    l5 > 0.01
    l6 > 0.01       
    vX1 > 0.01
    vX2 > 0.01    
    vX3 > 0.01
    Y1 ~~ vY1*Y1
    Y2 ~~ vY2*Y2    
    Y3 ~~ vY3*Y3
    vY1 > 0.01
    vY2 > 0.01    
    vY3 > 0.01   
    phi < 0.99
    phi > -0.99
    "


pars <- scan.vec("phi")

    lavmodel <- lavmodel1 
    std.lv <- TRUE



model_name <- "SEM_ML"

define_dfr1 <- function()
{
    cp <- coef(mod1a)
    dfr1 <- data.frame( design_dd, model=model_name , par=pars, true=c(phi), delta=delta,
                        est=cp[pars] )
    return(dfr1)
}


mod1a <- lavaan::sem(lavmodel, sample.cov=S, sample.nobs=1e5, estimator="ML", std.lv=std.lv)
dfr1 <- define_dfr1()
dfr <- rbind(dfr, dfr1)

model_name <- "SEM_ULS"   
mod1a <- lavaan::sem(lavmodel, sample.cov=S, sample.nobs=1e5, estimator="ULS", std.lv=std.lv)
dfr1 <- define_dfr1()
dfr <- rbind(dfr, dfr1)


model_name <- "LSAM_ML"   
mod0a <- mod1a <- lavaan::sam(lavmodel, sample.cov=S, sample.nobs=1e5, sam.method="local", std.lv=std.lv,
                        local.options = list(M.method = "ML" ) )                       
dfr1 <- define_dfr1()
dfr <- rbind(dfr, dfr1)

 model_name <- "GSAM_ML"   
 mod1a <- lavaan::sam(lavmodel, sample.cov=S, sample.nobs=1e5, sam.method="global", std.lv=std.lv)
 dfr1 <- define_dfr1()
 dfr <- rbind(dfr, dfr1)





#*** multiple group estimation in sirt
x <- c( rep(lam1,3), rep(lam2,3), rep(1-lam1^2,3), rep(1-lam2^2,3), PHI[1,2] )
names(x) <- c( paste0("lam",1:6), paste0("theta",1:6), "phi")




fun <- function(x, p, eps=1e-4)
{
    LAM <- matrix( 0, nrow=6, ncol=2)
    LAM[1:3,1] <- x[1:3]
    LAM[4:6,2] <- x[4:6]
    THETA <- diag(rep(1,6))
    diag(THETA) <- x[7:12]
    PHI <- matrix(1,2,2)
    PHI[1,2] <- PHI[2,1] <- x["phi"]
    SIGMA0 <- LAM %*% PHI %*% t(LAM) + THETA
    resid <- S-SIGMA0
    sum( ( resid^2 + eps )^(p/2) )
}

lower <- 0.01+0*x
lower["phi"] <- -0.99
upper <- 1e5+0*x
upper["phi"] <- 0.99

x00 <- x


for (p in pow ){
    name3 <- paste0("ME_pow",p, "_P2")
    x0 <- x00
    for (eps in c(1e-1,1e-2,1e-3,1e-4,1e-5)){
        mod <- nlminb(start=x0, objective=fun, p=p,eps=eps, lower=lower, upper=upper)
        x0 <- mod$par
    }    
    x00 <- x0
    dfr2 <- dfr1
    dfr2$model <- name3
    dfr2$est <- mod$par[ c("phi")]
    dfr <- rbind(dfr, dfr2)
}


#@@@@@@ inclusion
#--- robust GSAM models

mm1 <- lavaan::sem( lav_mm1, sample.cov=S, sample.nobs=1e5, estimator="ML", std.lv=TRUE)
mm2 <- lavaan::sem( lav_mm2, sample.cov=S, sample.nobs=1e5, estimator="ML", std.lv=TRUE)
LAM_mm1 <- coef(mm1)[1:3]
LAM_mm2 <- coef(mm2)[1:3]

x00 <- c(.5)

LAM_mm <- matrix( 0 , 6, 2)
LAM_mm[1:3,1] <- LAM_mm1
LAM_mm[4:6,2] <- LAM_mm2


p <- 2

for (p in pow){
    name3 <- paste0("RGSAM_pow",p)
    x0 <- x00
    mod <- robust_gsam2(start=x00,  S=S, LAM_mm=LAM_mm, pow=p)
    dfr2 <- dfr1
    dfr2$model <- name3
    dfr2$est <- x00 <- mod$par
    dfr <- rbind(dfr, dfr2)
}



#@@@@@@ inclusion
#--- factor analysis models

TR <- matrix(0,6,6)
TR[1:3,4:6] <- 1
TR[4:6,1:3] <- 1

thresh <- .05
thresh <- 1
is_threshold <- TR

define_dfr2a <- function(){
    dfr2 <- dfr1
    dfr2$model <- model_name
    dfr2$est <- abs(m1$Phi[1,2])
    return(dfr2)
}

define_dfr2b <- function(){
    dfr2 <- dfr1
    dfr2$model <- paste0( model_name, "robme05")
    dfr2$est <- mod$phi
    return(dfr2)
}


thresh <- .03


loading_cutoff <- .15


for (thresh in c(1,.03) ){

res <- factanal_residual_threshold(covmat=S, n.obs=1e5, factors=2, is_threshold=TR, 
            thresh=thresh, a=2, maxiter=400, verbose=TRUE)                        
L <- res$LAMBDA

model_name <- paste0("efa_oblimin_thresh", thresh)
m1 <- oblimin(L=L)
dfr <- rbind( dfr, define_dfr2a() )

mod <- me_fa_ss(S=S, loadings=m1$loadings, p=0.5, loading_cutoff=loading_cutoff, eps=1e-4)
dfr <- rbind( dfr, define_dfr2b() )


model_name <- paste0("efa_geomin_thresh", thresh)
m1 <- geominQ(L=L)
dfr <- rbind( dfr, define_dfr2a() )
mod <- me_fa_ss(S=S, loadings=m1$loadings, p=0.5, loading_cutoff=loading_cutoff, eps=1e-4)
dfr <- rbind( dfr, define_dfr2b() )


model_name <- paste0("efa_quartimin_thresh", thresh)
m1 <- quartimin(L=L)
dfr <- rbind( dfr, define_dfr2a() )
mod <- me_fa_ss(S=S, loadings=m1$loadings, p=0.5, loading_cutoff=loading_cutoff, eps=1e-4)
dfr <- rbind( dfr, define_dfr2b() )

model_name <- paste0("efa_promax_thresh", thresh)
m1 <- fa_promax(L=L)
dfr <- rbind( dfr, define_dfr2a() )
mod <- me_fa_ss(S=S, loadings=m1$loadings, p=0.5, loading_cutoff=loading_cutoff, eps=1e-4)
dfr <- rbind( dfr, define_dfr2b() )


start <- c(.5,.5)
p <- 1
for (p in pow[-1]){
    model_name <- paste0("efa_lprot_p", p, "_thresh", thresh)
    m1 <- res <- lp_rotation2( start=start, L=L, pow=p)    
    start <- res$par
    dfr <- rbind( dfr, define_dfr2a() )

    mod <- me_fa_ss(S=S, loadings=m1$loadings, p=0.5, loading_cutoff=loading_cutoff, eps=1e-4)
    dfr <- rbind( dfr, define_dfr2b() )    
    
}
} # end thresh

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#--- output

save.data( dfr, filename , type="csv2" , path= pf1, suffix = "RESULTS")





}  # end dd
}  # end bb