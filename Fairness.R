eta2 <- function(x, gpe) {
  res <- varinter(x, gpe)/vartot(x)
  return(res)
}
varinter <- function(x, gpe) {
  moyennes <- tapply(x, gpe, mean)
  effectifs <- tapply(x, gpe, length)
  res <- (sum(effectifs * (moyennes - mean(x))^2))
  return(res)
}
vartot <- function(x) {
  res <- sum((x - mean(x))^2)
  return(res)}

library(FactoMineR)
library(cvTools)
library(factoextra)
library(pls)


delete.intercept <- function(mm) {
  ## Save the attributes prior to removing the intercept coloumn:
  saveattr <- attributes(mm)
  ## Find the intercept coloumn:
  intercept <- which(saveattr$assign == 0)
  ## Return if there was no intercept coloumn:
  if (!length(intercept)) return(mm)
  ## Remove the intercept coloumn:
  mm <- mm[,-intercept, drop=FALSE]
  ## Update the attributes with the new dimensions:
  saveattr$dim <- dim(mm)
  saveattr$dimnames <- dimnames(mm)
  ## Remove the assignment of the intercept from the attributes:
  saveattr$assign <- saveattr$assign[-intercept]
  ## Restore the (modified) attributes:
  attributes(mm) <- saveattr
  ## Return the model matrix:
  mm
}

p=60
sigma <- function(theta=0, lambda=c(1,1)) {
  cos.t <- cos(theta); sin.t <- sin(theta)
  a <- matrix(c(cos.t, sin.t, -sin.t, cos.t), ncol=2)
  t(a) %*% diag(lambda) %*% a
}
library(MASS)
n1 <- 20   # First group population
n2 <- 45   # Second group population
x <- rbind(mvrnorm(n1, c(-2,-1), sigma(0, c(1/2,1))),
           mvrnorm(n2, c(0,1), sigma(pi/3, c(1, 1/3))))

eps <- 0.25  # Error SD should be small compared to the SDs for the blobs
x <- cbind(x, matrix(rnorm(dim(x)[1]*(p-2), sd=eps), ncol=p-2))
rot <- qr.Q(qr(matrix(rnorm(p^2), p)))
y <- x %*% rot

vio <- c(rep(1,n1) ,rep(0,n2))
dat <- cbind(y,vio)

correlatedValue = function(x, r){
  r2 = r**2
  ve = 1 - r2
  SD = sqrt(ve)
  e  = rnorm(dim(x)[1], mean = 0, sd = SD)
  rvect <- rep(r,dim(x)[2])
  y  =  x %*% rvect + e
  return(y)
}

expl <- correlatedValue(y,0.9)

### fair pca 

groupes=cvFolds(n=dim(y)[1],K=5)
fibrosisScaled<- y
#penalty<-c(seq(0.000001,0.1,0.0001))
penalty<-c(seq(0.00000008,0.0000004,0.000001), seq(0.0000005,0.01,0.00001))
#penalty<-c(seq(0,0.0001,0.00001),seq(0.001,0.01,0.0001))
#dimension <- matrix(0,length(penalty),5)
RMSE <- c()
dimension<- matrix(0,ncol=5,nrow=length(penalty))
var_exp <- c()
for(i in 1:length(penalty)){
  rmse <- c()
  
  explained_var <- c()
  for(j in 1:5){
    test=groupes$subsets[groupes$which==j]
    train =groupes$subsets[groupes$which !=j]
    
    # Build X_train, y_train, X_test, y_test
    data <- cbind(expl, y)
    X_train <- data[train, -1]
    y_train <- data[train, "expl"]
    
    X_test <- data[test, -1]
    y_test <- data[test, "expl"]
    S_train  <- vio[train]
    
    res_pca_train <- PCA(X_train, ncp=dim(fibrosisScaled)[2],graph=F)
    
    result=c()
    for(l in 1:dim(res_pca_train$ind$coord)[2]){
      result[l]=eta2(res_pca_train$ind$coord[,l],S_train)
    }
    select <- which(result > penalty[i])
    
    if(length(select)!=dim(res_pca_train$ind$coord)[2]){
      eig.val <- get_eigenvalue(res_pca_train)
      explained_var[j] <- sum(eig.val[-select,2])
      
      train_proj <- res_pca_train$ind$coord[,-select]
      
      dimension[i,j] <- dim(res_pca_train$ind$coord)[2]-length(select)
      
      #dimension[i,j] <- dim( train_proj)[2]
      test_standard <- X_test
      # standardisation
      test_standard <- t(apply(X_test, MARGIN = 1, FUN = function(x) {(x-res_pca_train$call$centre )/res_pca_train$call$ecart.type} ))
      # projection (matrix multiplication between individuals and dimensions coordinates)
      test_proj <-test_standard %*%res_pca_train$svd$V[,-select]
      
    }else{
      eig.val <- get_eigenvalue(res_pca_train)
      index= which(result==min(result))
      explained_var[j] <- sum(eig.val[index,2])
      train_proj <- res_pca_train$ind$coord[,index]
      
      dimension[i,j] <- 1
      
      #dimension[i,j] <- dim( train_proj)[2]
      test_standard <- X_test
      # standardisation
      test_standard <- t(apply(X_test, MARGIN = 1, FUN = function(x) {(x-res_pca_train$call$centre )/res_pca_train$call$ecart.type} ))
      # projection (matrix multiplication between individuals and dimensions coordinates)
      test_proj <-test_standard %*%res_pca_train$svd$V[,index]
    }
    
    # convert to dataframe
    train_proj <- as.data.frame(train_proj)
    test_proj <- as.data.frame(test_proj)
    
    # add target column
    train_proj$expl <- expl[train]
    test_proj$expl <- expl[test]
    
    # rename variables to work with predict function
    colnames(test_proj) <- colnames(train_proj)
    
    # fit
    model_pcr <- lm(expl~., data=train_proj)
    
    # summary
    #summary(model_pcr)
    if(dim(test_proj)[2]<=2){
      projection <- test_proj[,-dim(test_proj)[2]]
      projection2 <- as.data.frame(projection)
      colnames(projection2)<- "train_proj"
      predict_pcr <- predict(model_pcr,projection2)
    }else{
      predict_pcr <- predict(model_pcr,test_proj[,-dim(test_proj)[2]])}
    
    minimal=min(test_proj$expl)
    maximal=max(test_proj$expl)
    squared_sums <- sum((test_proj$expl - predict_pcr)^2)
    mse <- squared_sums/length(test_proj$expl)
    rmsecalc <- sqrt(mse)
    rmse[j] <- rmsecalc/(maximal-minimal)
    
    
  }
  RMSE[i]=mean(rmse)
  
  var_exp[i]=mean(explained_var)
}
plot(penalty,RMSE,type="l",xlab="Threshold", ylab="NRMSE")
plot(penalty,var_exp,type="l",xlab="Threshold", ylab="Explained Variance")


####### fair PLS 


penalty<-c(seq(0.00000008,0.0000004,0.000001), seq(0.0000005,0.01,0.00001))

dimension <- matrix(0,length(penalty),5)
y <- data.frame(y)
fibrosisScaled<- y
groupes=cvFolds(n=dim(fibrosisScaled)[1],K=5)
RMSE <- c()
Var_ex <- c()
dimension <- matrix(0,ncol=5, nrow= length(penalty))
for(i in 1:length(penalty)){
  rmse <- c()
  exp_var <- c()
  for(j in 1:5){
    test=groupes$subsets[groupes$which==j]
    train =groupes$subsets[groupes$which !=j]
    
    # Build X_train, y_train, X_test, y_test
    data <- cbind(expl, fibrosisScaled)
    X_train <- data[train, -1]
    y_train <- data[train, "expl"]
    
    X_test <- data[test, -1]
    y_test <- data[test, "expl"]
    S_train  <- vio[train]
    data_train <- cbind(X_train, y_train)
    plsr_fit <- plsr(y_train ~ .,data = data_train, ncomp=50)
    
    
    result2=c()
    for(l in 1:dim(plsr_fit$score)[2]){
      result2[l]=eta2(plsr_fit$scores[,l],S_train)
    }
    select2 <- which(result2 > penalty[i])
    
    
    if(length(select2)!=50){
      exp_var[j] <- sum(explvar(plsr_fit)[-select2])
      # Predict
      components <- setdiff(c(1:50) ,select2)
      dimension[i,j] <- 50-length(select2)
    }else{
      index2= which(result2==min(result2))
      exp_var[j] <- sum(explvar(plsr_fit)[index2])
      # Predict
      components <- setdiff(c(1:50) ,index2)
      dimension[i,j]=1
      
    }
    
    
    pred <- predict(plsr_fit, comps=components,newdata=X_test)
    # RMSE
    minimal=min(y_test)
    maximal=max(y_test)
    
    squared_sums <- sum((y_test - pred)^2)
    mse <- squared_sums/length(y_test)
    rmsecalc <- sqrt(mse)
    rmse[j] <- rmsecalc/(maximal-minimal)
    
  }
  RMSE[i]=mean(rmse)
  Var_ex[i]=mean(exp_var)
}


plot(penalty,RMSE,type="l",xlab="Threshold", ylab="NRMSE")
plot(penalty,Var_ex,type="l",xlab="Threshold", ylab="Explained Variance")


## PLS penalized


penalty<-c(seq(0,5, 0.05))


fibrosisScaled <- y
groupes=cvFolds(n=dim(fibrosisScaled)[1],K=5)
RMSE <- c()
Var_ex <- c()
for(i in 1:length(penalty)){
  rmse <- c()
  exp_var <- c()
  for(j in 1:5){
    test=groupes$subsets[groupes$which==j]
    train =groupes$subsets[groupes$which !=j]
    
    # Build X_train, y_train, X_test, y_test
    data <- cbind(expl, fibrosisScaled)
    X_train <- data[train, -1]
    y_train <- data[train, "expl"]
    
    X_test <- data[test, -1]
    y_test <- data[test, "expl"]
    S_train  <- vio[train]
    y_new_train <- y_train - penalty[i] * S_train
    plsr_fit <- plsr(y_new_train ~ X_train, ncomp=1)
    
    exp_var[j] <- sum(explvar(plsr_fit))
    # Predict
    
    pred <- predict(plsr_fit,newdata=X_test)
    # RMSE
    minimal=min(y_test)
    maximal=max(y_test)
    
    squared_sums <- sum((y_test - pred)^2)
    mse <- squared_sums/length(y_test)
    rmsecalc <- sqrt(mse)
    rmse[j] <- rmsecalc/(maximal-minimal)
    
  }
  RMSE[i]=mean(rmse)
  Var_ex[i]=mean(exp_var)
}
plot(penalty,RMSE,type="l",xlab="Threshold", ylab="NRMSE")
plot(penalty,Var_ex,type="l",xlab="Threshold", ylab="Explained Variance")
