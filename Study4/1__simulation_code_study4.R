
########################################
# Simulation Study 4
########################################

library(miceadds)
library(lavaan)


comp <- "1001"
# pf00 <- "d:/_F/Projekte/SEM_2step/"
pf00 <- "p:/Eigene_Projekte/SEM_2step/Analysen/"

pfR <- file.path( pf00, "R-Code")
files_move(pfR)
source.all(pfR, c("sam","rme") )

lavmodel <- readLines( file.path(pfR, "lavmodel_5factors.txt" ))
lavmodel <- paste0(paste0( lavmodel, collapse=" \n"), "\n")

# read true parameters
true_pars <- load.data( "_POPU", path=pfR, type="csv2")
true_pars <- true_pars[ grep("phi", true_pars[,1]) , ]

# read covariance matrices
setwd(pfR)
covmat_list <- list()
covmat_list[[1]] <- read.table("rl_simulation_true_model1__COVMAT.txt", header=TRUE)
covmat_list[[2]] <- read.table("rl_simulation_true_model2__COVMAT.txt", header=TRUE)
covmat_list[[3]] <- read.table("rl_simulation_true_model3__COVMAT.txt", header=TRUE)
NV <- ncol(covmat_list[[1]])


# simulation number
sim <- "501"

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
    dgm=c(1,2,3)
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
# dd <- 6
# dd <- 18
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


S <- covmat_list[[ dgm ]]




if (N<inf_val){
    dat <- MASS::mvrnorm(n=N, mu=rep(0,NV), Sigma=S)
    S <- cov(dat) * (N-1) / N
}


# save example datasets
save.data( S , paste0("DATA_DES",dd ) , type="csv2" , path= pf1, row.names=TRUE)



dfr <- NULL



# pars <- scan.vec("phi")
pars <- true_pars[,1]

#    lavmodel <- lavmodel1 
    std.lv <- TRUE



model_name <- "SEM_ML"

define_dfr1 <- function()
{
    cp <- coef(mod1a)
    dfr1 <- data.frame( design_dd, model=model_name , par=pars, true=true_pars[,2], 
                        est=cp[pars] )
    return(dfr1)
}


mod00 <- mod1a <- lavaan::sem(lavmodel, sample.cov=S, sample.nobs=1e5, estimator="ML", std.lv=std.lv)
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


model_name <- "RME"

# x <- c(x1, x2, x3)
x <- c11 <- coef(mod00)

lower <- .01 + 0*x
lower[ grep("phi", names(x)) ] <- -.99
upper <- 20 + 0*x
upper[ grep("phi", names(x)) ] <- .99


dfr1b <- dfr1
dfr1b$model <- model_name


start <- x

powers <- .5
for (p in powers){
    x0 <- start
    for (eps in c(1e-2, 1e-3)){    
        mod <- nlminb(start=x0, objective=fit_rme, S=S, p=p, eps=eps, lower=lower, upper=upper,
                    control=list(iter.max=1000000) )
       x0 <- mod$par  
       if (mod$convergence==0){
            c1 <- mod$par
       }       
    }
    start <- mod$par
}

dfr1b$est <- c1[pars]

dfr <- rbind(dfr, dfr1b)

#--- output

save.data( dfr, filename , type="csv2" , path= pf1, suffix = "RESULTS")





}  # end dd
}  # end bb
