library("WoodburyMatrix")
library("sandwich")
library("clubSandwich")
library("lme4")

# Import data and creat factors for clustering

load("ctdata.RData")
ctdata$fcluster = factor(ctdata$hdist)
ctdata$ftime = factor(ctdata$time)
ctdata$clustime = interaction(ctdata$fcluster,ctdata$ftime)

# fit models glmer object
rslt = lmer(ct ~ ftime + treat + (1 | fcluster) + (1 | clustime), data=ctdata)
rslt0 = glmer(ct ~ ftime + treat + (1 | fcluster), family=binomial, data=ctdata)
rslt1 = glmer(ct ~ ftime + treat + (1 | fcluster) + (1 | clustime), family=binomial, data=ctdata)
rslt2 = glmer(ct ~ ftime + treat + (treat | fcluster), family=binomial, data=ctdata)
rslt3 = glmer(ct ~ ftime + treat + (treat | fcluster) + (1 | clustime), family=binomial, data=ctdata)
rslt4 = glmer(ct ~ ftime + treat + (1 | fcluster) + (0 + treat | fcluster) + (1 | clustime), family=binomial, data=ctdata)


vcovCR.glmerModZY <- function(obj, cluster = NULL, type = "classic") {
  # put some check here to ensure the input is correct
  # If cluster not specified, will be set to attr(obj,"cluster")
  if (is.null(cluster)) {
    cluster <- attr(obj, "cluster")
  }
  
  link <- family(obj)$link
  
  beta <- matrix(fixef(obj), ncol = 1)
  
  np <- length(beta)
  
  gamma <- matrix(unlist(ranef(obj)), ncol = 1)
  
  nq <- length(gamma)
  
  X <- model.matrix(obj,type = "fixed")
  Z <- model.matrix(obj, type = "random")  # Z matrix for random effects
  Y <- obj@resp$y
  
  eta <- predict(obj, type = "link")
  ginv_eta <- predict(obj, type = "response")
  
  if (link == "identity"){
    delta = diag(nobs(obj)) # diag matrix
    deltainv = diag(1/diag(delta))# a more efficient way to get solve(delta)
  }
  
  else if (link == "logit"){
    delta <- diag(exp(eta)/(1+exp(eta))^2)
    deltainv <- diag(1/diag(delta)) 
  }
  
  else if (link == "log") {
    delta <- diag(exp(eta))
    deltainv <- diag(1/diag(delta))
  }
  
  P <- deltainv %*% (Y - ginv_eta) + eta
  e <- matrix(P - X %*% beta, ncol = 1)
  XtVX <- vcov(obj)  # model based variance
  #theta <- as.data.frame(VarCorr(obj)) # This is where the problem is, when you fit different models, this is different
  
  sigma2 <- sigma(obj)^2
  lambda <- getME(obj,"Lambda")
  R <- lambda %*% t(lambda) * sigma2
  
  G <- ngrps(obj)["fcluster"]
  sum <- matrix(0, np, np)
  
  for (g in 1:G){
    grp = ctdata$fcluster == g
    ng = sum(grp)
    
    # V = Z[grp,]%*%R%*%t(Z[grp,]) + deltainv[grp,grp]%*%Sigma%*%deltainv[grp,grp]
    # Vinv = solve(V)
    Sigma <- sigma2 * diag(ng)
    WB_A <- diag(1/diag(deltainv[grp,grp]%*%Sigma%*%deltainv[grp,grp]))

    WB_U <- Z[grp,]
    WB_C <- solve(R)
    WB_V <- t(Z[grp,])

    W <- WoodburyMatrix(A = WB_A, U = WB_U, B = WB_C, V = WB_V)
    Vinv <- solve(W)
    
    H = X[grp, ] %*% XtVX %*% t(X[grp, ]) %*% Vinv
    Q = t(X[grp, ]) %*% Vinv %*% X[grp, ] %*% XtVX
    
    # if loop, choose A, F and c based on robust variance form specified by the input 'type'
    F = diag(ng)
    A = diag(np)
    
    sum = sum + A %*% t(X[grp, ]) %*% Vinv %*% t(F) %*% e[grp, ] %*% t(e[grp, ]) %*% F %*% Vinv %*% X[grp, ] %*% A
  }
  
  c = 1
  robustVar <- c * XtVX %*% sum %*% XtVX
  diag(robustVar)  # Return diagonal elements representing the variances
}

#test the function

vcovCR.glmerModZY(rslt, type = "classic")
# vcovCR.glmerModZY(rslt0, type = "classic")
vcovCR.glmerModZY(rslt0,ctdata$fcluster, type = "classic") # manully input cluster info
vcovCR.glmerModZY(rslt1, type = "classic")
vcovCR.glmerModZY(rslt2, type = "classic")
vcovCR.glmerModZY(rslt3, type = "classic")
vcovCR.glmerModZY(rslt4, type = "classic")

# > vcovCR.glmerModZY(rslt, type = "classic")
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0001379472 0.0002188822 0.0002572821 0.0003607838 0.0003631747 0.0001417230 

#   > vcovCR.glmerModZY(rslt0,ctdata$fcluster, type = "classic") # manully input cluster info
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0001817311 0.0002142073 0.0002234002 0.0003933175 0.0004061382 0.0002053313 

# > vcovCR.glmerModZY(rslt1, type = "classic")
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0002074606 0.0002930088 0.0003169369 0.0005406408 0.0005345198 0.0002772918 

# > vcovCR.glmerModZY(rslt2, type = "classic")
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0003603288 0.0002999964 0.0002913427 0.0004321981 0.0004397180 0.0003357197 

# > vcovCR.glmerModZY(rslt3, type = "classic")
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0003608960 0.0003773237 0.0003381913 0.0005305571 0.0005124369 0.0003633539 

# > vcovCR.glmerModZY(rslt4, type = "classic")
# (Intercept)       ftime1       ftime2       ftime3       ftime4        treat 
# 0.0002081884 0.0003704127 0.0003326814 0.0005423342 0.0005245356 0.0002837417 
