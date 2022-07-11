
#*** thresholding function
thresholding <- function(x, thresh, a=2)
{
	x0 <- x
	x <- abs(x)
	v2 <- a / (a-1) * ( x - thresh )
	y <- x*( x >= a*thresh ) + v2*(x > thresh )*(x < a*thresh)
	y <- sign(x0)*y
	return(y)
}

#*** factanal with residual thresholding
factanal_residual_threshold <- function(covmat, n.obs, factors, 
		is_threshold, rotation="varimax",
		thresh=1e9, a=2, maxiter=100, maxchange=1e-7, verbose=FALSE )
{	
	CM <- stats::cov2cor(covmat)
	
	#*** factor analysis without thresholding
	fa1 <- stats::factanal( covmat=CM, n.obs=n.obs, factors=factors)
	L0 <- L <- fa1$loadings
	resid <- CM - L %*% t(L)
	E <- thresholding( x=resid, thresh=thresh, a=a )*is_threshold

	#*** iterate
	iter <- 0
	E0 <- E
	print_digits_change <- round(abs(log10(maxchange))+3)
	iterate <- TRUE
	converged <- FALSE
	while(iterate){
		CM1 <- CM - E
		fa1 <- stats::factanal( covmat=CM1, n.obs=n.obs, factors=factors)
		L <- fa1$loadings
		resid <- CM - L %*% t(L)
		E <- thresholding( x=resid, thresh=thresh, a=a )*is_threshold
		change <- max( abs(E-E0) )
		iter <- iter + 1
		if (verbose){
			pp <- paste0("Iteration ", iter, " | maximum change = ", 
					round( change, print_digits_change) )
			cat(pp, "\n")
			utils::flush.console()
		}
		if ((change < maxchange)|(iter>maxiter)){
			iterate <- FALSE
			converged <- TRUE
		}
		E0 <- E
	}
	fa1 <- stats::factanal( covmat=CM1, n.obs=n.obs, factors=factors, rotation=rotation)
	
	#** unstandardized loadings
	L <- as.matrix(fa1$loadings) * sqrt(diag(covmat))

	#** covariance matrix
	R <- fa1$rotmat		
	CV <- t(R) %*% R
	PHI <- stats::cov2cor(CV)	
	h1 <- sqrt(diag(covmat))
	E <- E * outer(h1,h1)

	#--- output
	res <- fa1
	res$LAMBDA <- L
	res$Phi <- PHI
	res$E <- E	
	
	res$iter <- iter
	res$change <- change
	res$converged <- converged
	return(res)
}

#*** Lp rotation

define_rotmat2 <- function(x)
{
    a <- x[1]
    c <- x[2]
    A <- matrix( c( a, sqrt(1-a^2), sqrt(1-c^2), c), 2,2, byrow=TRUE)
    return(A)
}

rot_fun2 <- function(x, L, pow=1, eps=1e-4)
{
    A <- define_rotmat2(x=x)
    L1 <- L %*% MASS::ginv(A)
    sum( (L1^2+eps)^(pow/2) )
}



lp_rotation2 <- function(start, L, pow, eps=1e-4)
{
    #* optimize
    BB <- .995    
    one <- rep(1,2)
    res <- stats::nlminb( start, objective=rot_fun2, L=L, 
				pow=pow, eps=eps, lower=-BB*one, upper=BB*one)
    A <- define_rotmat2(x=res$par)
    Phi <- A %*% t(A)
    Lambda <- L %*% MASS::ginv(A)
    res$rotmat <- A
    res$Lambda <- Lambda
    res$Phi <- Phi
    return(res)
}

fa_promax <- function(L, m=4)
{
	res <- stats::promax(x=L, m=4)
	R <- res$rotmat
	res$Phi <- stats::cov2cor( t(R) %*% R )
	return(res)
}


lav_mm1 <- "
    FX =~ l1*X1 + l2*X2 + l3*X3
    FX ~~ 1*FX
    X1 ~~ v1*X1
    X2 ~~ v2*X2
    X3 ~~ v3*X3        
    l1 > 0.01
    l2 > 0.01
    l3 > 0.01
    v1 > 0.01
    v2 > 0.01
    v3 > 0.01
    "
lav_mm2 <- "
    FY =~ l1*Y1 + l2*Y2 + l3*Y3
    FY ~~ 1*FY
    Y1 ~~ v1*Y1
    Y2 ~~ v2*Y2
    Y3 ~~ v3*Y3        
    l1 > 0.01
    l2 > 0.01
    l3 > 0.01
    v1 > 0.01
    v2 > 0.01
    v3 > 0.01
    "
	

	
fun_cov2 <- function(x, pow, S, LAM_mm, eps=1e-4)
{    
	PHI <- matrix(1,2,2)
	PHI[1,2] <- PHI[2,1] <- x
	COV <- LAM_mm %*% PHI %*% t(LAM_mm)
	resid <- S[1:3,4:6] - COV[1:3,4:6]
	sum( ( resid^2 + eps )^(pow/2) )
}


robust_gsam2 <- function(start, pow, S, LAM_mm, eps=1e-4)
{
	res <- stats::nlminb(start=start, objective=fun_cov2, S=S, LAM_mm=LAM_mm, 
							pow=pow, eps=eps)
	return(res)
}



fun_cov2a <- function(x, pow, S, LAM_mm, eps=1e-4)
{    
	PHI <- matrix(1,2,2)
	PHI[1,1] <- x[2]
	PHI[2,2] <- x[3]
	PHI[1,2] <- PHI[2,1] <- x[1]*sqrt(x[2])*sqrt(x[3])
	COV <- LAM_mm %*% PHI %*% t(LAM_mm)
	resid <- S - COV
	diag(resid) <- 0
	sum( ( resid^2 + eps )^(pow/2) )
}


robust_gsam2a <- function(start, pow, S, LAM_mm, eps=1e-4)
{
	res <- stats::nlminb(start=start, objective=fun_cov2a, S=S, LAM_mm=LAM_mm, 
							pow=pow, eps=eps, lower=c(-Inf,.01,.01) )
	return(res)
}
