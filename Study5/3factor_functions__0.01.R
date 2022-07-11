

comp_cov_1factor <- function(x)
{
    return( outer(x,x) )
}

loss_fun_1factor <- function(x, S)
{
    RES <- S - comp_cov_1factor(x=x)
    diag(RES) <- 0
    sum( RES^2 )
}


comp_cov_3factor <- function(x)
{
    LAM <- matrix(0, 9,3 )
    LAM[1:3,1] <- x[1:3]
    LAM[4:6,2] <- x[4:6]
    LAM[7:9,3] <- x[7:9]
    PHI <- matrix(1,3,3)
    PHI[1,2] <- PHI[2,1] <- x[10]
    PHI[1,3] <- PHI[3,1] <- x[11]
    PHI[2,3] <- PHI[3,2] <- x[12]
    COV <- LAM %*% PHI %*% t(LAM)
    return(COV)
}


loss_fun_3factor <- function(x, S)
{
    RES <- S - comp_cov_3factor(x=x)
    diag(RES) <- 0
    sum( RES^2 )
}

loss_fun_3factor_lsam <- function(x, x1, S)
{
    x <- c(x1, x)
    RES <- S - comp_cov_3factor(x=x)
    diag(RES) <- 0
    sum(RES^2)
}

