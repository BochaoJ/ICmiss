SimRegDat <- function(n = 100, p = 200, type = "indep", rate = 0.1)
{
  mu_x <- rep(0,p)
  sigma_x <- 2
  sigma_eps <- 1
  beta <- rep(0,p)
  beta_0 <- rep(1,n)
  beta[1:5] <- c(10,20,-15,-25,50)/10
  if(type=="indep"){
    cat("Generate independent covariates...\n")
    x <- rmvnorm(n,mu_x,diag(sigma_x^2,p))
    eps <- rnorm(n,0,sigma_eps)
    y <- beta_0 +x%*%beta+eps
    sec <- sample(1:length(x),size=rate*length(x))
    x[sec] <- NA
    result <- list(x=x,y=y,coef= beta)
  }else if(type=="dep"){
    cat("Generate dependent covariates...\n")
    C <- diag(1, p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (abs(j - i) == 1) 
          C[i, j] = 0.5 else if (abs(j - i) == 2) 
            C[i, j] = 0.25
      }
    }
    A <- diag(0, p)
    sigma <- solve(C)
    sigma2 <- sigma
    for (i in 1:p) {
      for (j in 1:p) {
        if (sigma[i, j] <= sigma[j, i]) {
          sigma2[i, j] <- sigma[j, i]
        } else {
          sigma2[i, j] <- sigma[i, j]
        }
        if (C[i, j] != 0 && i != j) 
          A[i, j] <- 1
      }
    }
    x <- rmvnorm(n, mu_x, sigma2)
    eps <- rnorm(n,0,sigma_eps)
    y <- beta_0 +x%*%beta+eps
    sec <- sample(1:length(x),size=rate*length(x))
    x[sec] <- NA
    result <- list(x=x,y=y,coef= beta)
  } else{
    stop("Only 'dep' or 'indep' types are provided !")
  }
  return(result)
}