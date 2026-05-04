library(SomaDataIO)
library(ggfortify)
library(tidyverse)
library(factoextra)
library(pheatmap)
library(MASS)
library(psych)
library(truncnorm)
library(compositions)
library(splitTools)
library(pROC)
library(glmnet)
library(bnlearn)
library(xgboost)
library(mvtnorm)
library(ranger)
library(Boruta)
library(ropls)
library(GEOquery)
library(qpcR)
library(MASS)
library(modEvA)
library(glmtoolbox)
library(boot)
library(rpart)
library(DESeq2)
library(ciTools)

#Functions----------------------------------------------------------------------------------------

Simulate_OMICs_Data <- function(DF, DF_Response, transform, fast, Case_Rows, p, n) {
  
  Case <- DF[Case_Rows,]
  Control <- DF[-Case_Rows,]
  
  Case <- Case[,-which(colnames(DF) %in% DF_Response)]
  Control <- Control[,-which(colnames(DF) %in% DF_Response)]
  
  if(transform==T) {
    Case <- data.frame(apply(Case, 2, function(i) log2(i+0.01)))
    Control <- data.frame(apply(Control, 2, function(i) log2(i+0.01)))
  }
  
  #Generate----------------------------------------------------------------
  
  #Step 1 : Assume dataset has been powered for differential analysis

  Keep <- Get_Important_Features(DF, DF_Response, transform, p)
  
  case_n <- round(n*mean(DF[,colnames(DF) %in% DF_Response]),0)
  control_n <- round(n*(1-mean(DF[,colnames(DF) %in% DF_Response])),0)
  
  Case <- data.frame(apply(Case, 2, as.numeric))
  Control <- data.frame(apply(Control, 2, as.numeric))
  
  if(fast==F) {
    BNG_Case <- tabu(Case[,Keep]) 
  }
  if(fast==T) {
    BNG_Case <- rsmax2(Case[,Keep]) 
  }
  
  BNG_Case_fit <- bn.fit(BNG_Case, data=Case[,Keep])
  BNG_Case_Sim <- data.frame(rbn(BNG_Case_fit, case_n))
  
  if(fast==F) {
    BNG_Control <- tabu(Control[,Keep]) 
  }
  if(fast==T) {
    BNG_Control <- rsmax2(Control[,Keep]) 
  }
  BNG_Control_fit <- bn.fit(BNG_Control, data=Control[,Keep])
  BNG_Control_Sim <- data.frame(rbn(BNG_Control_fit, control_n))
  
  Data_Sim <- rbind(BNG_Case_Sim, BNG_Control_Sim)
  
  #Part 2 - Correlation Matrix of non-diff EX biomarkers
  
  if (transform==T) {
    Data_FeaturesOnly_TRANSFORM <- data.frame(apply(DF[,1:p], 2, function(i) scale(log2(i+0.01))))
  }
  
  else {
    Data_FeaturesOnly_TRANSFORM <- data.frame(apply(DF[,1:p], 2, function(i) scale(i)))
  }
  
  CorMat <- cov(Data_FeaturesOnly_TRANSFORM[,-which(colnames(Data_FeaturesOnly_TRANSFORM) %in% Keep)])
  
  CorMat_Smooth <- cor.smooth(CorMat)
  
  #Simulate log-normal data
  Data_Sim2 <- rmvnorm(n=n, mean=rep(0, ncol(Data_FeaturesOnly_TRANSFORM[,-which(colnames(Data_FeaturesOnly_TRANSFORM) %in% Keep)])),
                       sigma=CorMat_Smooth)
  
  Data_Sim2 <- data.frame(Data_Sim2)
  
  Data_Sim_scaled <- scale(Data_Sim)
  
  Data_Sim_Final <- cbind(Data_Sim_scaled, Data_Sim2)
  
  Data_Sim_Final$Group <- c(rep(1,case_n), rep(0,control_n))
  
  Data_Sim_Final
}



Simulate_OMICs_Data_SPLIT <- function(DF, DF_Response, Keep, transform, fast, Case_Rows, p, n) {
  
  Case <- DF[Case_Rows,]
  Control <- DF[-Case_Rows,]
  
  Case <- Case[,-which(colnames(DF) %in% DF_Response)]
  Control <- Control[,-which(colnames(DF) %in% DF_Response)]
  
  if(transform==T) {
    Case <- data.frame(apply(Case, 2, function(i) log2(i+0.01)))
    Control <- data.frame(apply(Control, 2, function(i) log2(i+0.01)))
  }
  
  #Generate----------------------------------------------------------------
  
  #Step 1 : Assume dataset has been powered for differential analysis
  
  case_n <- round(n*mean(DF[,colnames(DF) %in% DF_Response]),0)
  control_n <- round(n*(1-mean(DF[,colnames(DF) %in% DF_Response])),0)
  
  Case <- data.frame(apply(Case, 2, as.numeric))
  Control <- data.frame(apply(Control, 2, as.numeric))
  
  if(fast==F) {
    BNG_Case <- tabu(Case[,Keep]) 
  }
  if(fast==T) {
    BNG_Case <- rsmax2(Case[,Keep]) 
  }
  
  BNG_Case_fit <- bn.fit(BNG_Case, data=Case[,Keep])
  BNG_Case_Sim <- data.frame(rbn(BNG_Case_fit, case_n))
  
  if(fast==F) {
    BNG_Control <- tabu(Control[,Keep]) 
  }
  if(fast==T) {
    BNG_Control <- rsmax2(Control[,Keep]) 
  }
  BNG_Control_fit <- bn.fit(BNG_Control, data=Control[,Keep])
  BNG_Control_Sim <- data.frame(rbn(BNG_Control_fit, control_n))
  
  Data_Sim <- rbind(BNG_Case_Sim, BNG_Control_Sim)
  
  #Part 2 - Correlation Matrix of non-diff EX biomarkers
  
  if (transform==T) {
    
    GrandMeans <- apply(DF[,1:p], 2, function(i) mean(log2(i+0.01)))
    
    Data_FeaturesOnly_TRANSFORM <- data.frame(apply(DF[,1:p], 2, function(i) scale(log2(i+0.01))))
  }
  
  else {
    Data_FeaturesOnly_TRANSFORM <- data.frame(apply(DF[,1:p], 2, function(i) scale(i)))
  }
  
  #Correlation Matrix partition for computational speed
  
  Features <- p - length(Keep)
  
  n_partitions <- floor(Features / 1000) + 1
  remainder <- Features - ((n_partitions - 1)*1000)
  
  GrandMeans <- GrandMeans[-which(colnames(Data_FeaturesOnly_TRANSFORM) %in% Keep)]
  
  Data_FeaturesOnly_TRANSFORM <- Data_FeaturesOnly_TRANSFORM[,-which(colnames(Data_FeaturesOnly_TRANSFORM) %in% Keep)]
  
  
  Data_Sim2 <- NULL
  for (i in 1:n_partitions) {
    
    if (i == 1) {
    CorMat <- cov(Data_FeaturesOnly_TRANSFORM[,1:1000])
    Means <- GrandMeans[1:1000]
    }
    
    if (i == n_partitions) {
      CorMat <- cov(Data_FeaturesOnly_TRANSFORM[,(1000*(i-1) + 1):(1000*(i-1) + remainder)])
      Means <- GrandMeans[(1000*(i-1) + 1):(1000*(i-1) + remainder)]
    }
    
    if (i > 1 & i < n_partitions) {
    CorMat <- cov(Data_FeaturesOnly_TRANSFORM[,(1000*(i-1) + 1):(1000*(i))])
    Means <- GrandMeans[(1000*(i-1) + 1):(1000*(i))]
    }
    
    CorMat_Smooth <- cor.smooth(CorMat)
    
    #Simulate log-normal data
    Data_Sim2 <- cbind(Data_Sim2, rmvnorm(n=n, mean=Means,
                         sigma=CorMat_Smooth))
  
    
    }
  
  colnames(Data_Sim2) <- paste0("Y", 1:Features)
  
  Data_Sim2 <- data.frame(Data_Sim2)
  
  Data_Sim_Final <- cbind(Data_Sim, Data_Sim2)
  
  Data_Sim_Final$Group <- c(rep(1,case_n), rep(0,control_n))
  
  Data_Sim_Final
}



Simulate_OMICs_Data_BNG_Only <- function(DF, DF_Response, Keep, transform, fast, Case_Rows, p, n) {
  
  Case <- DF[Case_Rows,]
  Control <- DF[-Case_Rows,]
  
  Case <- Case[,-which(colnames(DF) %in% DF_Response)]
  Control <- Control[,-which(colnames(DF) %in% DF_Response)]
  
  if(transform==T) {
    Case <- data.frame(apply(Case, 2, function(i) log2(i+0.01)))
    Control <- data.frame(apply(Control, 2, function(i) log2(i+0.01)))
  }
  
  #Generate----------------------------------------------------------------
  
  #Step 1 : Assume dataset has been powered for differential analysis
  
  case_n <- round(n*mean(DF[,colnames(DF) %in% DF_Response]),0)
  control_n <- round(n*(1-mean(DF[,colnames(DF) %in% DF_Response])),0)
  
  Case <- data.frame(apply(Case, 2, as.numeric))
  Control <- data.frame(apply(Control, 2, as.numeric))
  
  if(fast==F) {
    BNG_Case <- tabu(Case[,Keep]) 
  }
  if(fast==T) {
    BNG_Case <- rsmax2(Case[,Keep]) 
  }
  
  BNG_Case_fit <- bn.fit(BNG_Case, data=Case[,Keep])
  BNG_Case_Sim <- data.frame(rbn(BNG_Case_fit, case_n))
  
  if(fast==F) {
    BNG_Control <- tabu(Control[,Keep]) 
  }
  if(fast==T) {
    BNG_Control <- rsmax2(Control[,Keep]) 
  }
  BNG_Control_fit <- bn.fit(BNG_Control, data=Control[,Keep])
  BNG_Control_Sim <- data.frame(rbn(BNG_Control_fit, control_n))
  
  Data_Sim <- rbind(BNG_Case_Sim, BNG_Control_Sim)
  
  Data_Sim$Group <- c(rep(1,case_n), rep(0,control_n))
  
  Data_Sim
}



plot_curve_full <- function(Curve, Type, FullAUC, Size) {
  
  if(Type=="PowerLaw") {
    
    test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = Curve, start = c(a=-0.5,b=-0.5, c=0.5), control = nls.control(maxiter = 1000))
    a <- summary(test)$coef[1,1]
    b <- summary(test)$coef[2,1]
    c <- summary(test)$coef[3,1]
    
    Model <- function(a, b, c, n) {a*(1/n^-b) + c}
    
    Values <- data.frame("Size"=1:Size, "Fitted"=Model(a=a, b=b, c=c, n=1:Size))
    
    Marker <- data.frame("N" = Get_N(Curve, FullAUC, Size), "AUC" = Values[Values==Get_N(Curve, FullAUC, Size), "Fitted"])
    
  }
  
  if(Type=="PowerLaw_Fixed") {
    c <- FullAUC
    
    test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = Curve, start = c(a=-0.5,b=-0.5), control = nls.control(maxiter = 1000))
    a <- summary(test)$coef[1,1]
    b <- summary(test)$coef[2,1]
    
    Model <- function(a, b, c, n) {a*(1/n^-b) + c}
    
    Values <- data.frame("Size"=1:Size, "Fitted"=Model(a=a, b=b, c=c, n=1:Size))
    
    Marker <- data.frame("N" = Get_N_Fixed(Curve, FullAUC, Size), "AUC" = Values[Values==Get_N_Fixed(Curve, FullAUC, Size), "Fitted"])
    
  }
  
  if(Type=="LOG") {
    
    test <- lm(data=Curve, Curve ~ log(n))
    
    B0 <- test$coefficients[1]
    
    B1 <- test$coefficients[2]
    
    Model <- function(B0, B1, n) {B0 + B1*log(n)}
    
    Values <- data.frame("Size"=1:Size, "Fitted"=Model(B0=B0, B1=B1, n=1:Size))
    
    Marker <- data.frame("N" = Get_N_LOG(Curve, FullAUC, Size), "AUC" = Values[Values==Get_N_LOG(Curve, FullAUC, Size), "Fitted"])
    
  }
  
  Data_forPlot_Raw <- data.frame("N" = Curve$n, "Raw" = Curve$Curve)
  
  Plot <- ggplot(data=Data_forPlot_Raw, aes(x=N, y=Raw)) + geom_point(color="blue", size=2.5) + ylim(c(min(Data_forPlot_Raw$Raw)-0.01, FullAUC+0.01)) + xlim(c(0,Size)) +
    geom_line(data=Values, aes(y=Fitted, x=Size), col="Blue", size=1) + 
    labs(x="Sample Size", y="AUC") + theme_bw() + theme(axis.title=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  
  
  
  Plot
  
  
}


Get_N <- function(DF, FullAUC, size) {
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5, c=0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  c <- summary(test)$coef[3,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:size)
  
  N <- which(Values - (FullAUC-0.02) >= 0)[1]
  
  N
  
}

Get_N_Th <- function(DF, FullAUC, size, Th) {
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5, c=0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  c <- summary(test)$coef[3,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:size)
  
  N <- which(Values - (FullAUC-Th) >= 0)[1]
  
  N
  
}

Get_N_Fixed <- function(DF, FullAUC, size) {
  c <- FullAUC
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:size)
  
  N <- which(Values - (FullAUC-0.02) >= 0)[1]
  
  N
  
}

Get_N_Fixed_Th <- function(DF, FullAUC, size, Th) {
  c <- FullAUC
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:size)
  
  N <- which(Values - (FullAUC-Th) >= 0)[1]
  
  N
  
}


Get_N_LOG <- function(DF, TargetAUC, Size) {
  
  test <- lm(data=DF, Curve ~ log(n))
  
  B0 <- test$coefficients[1]
  
  B1 <- test$coefficients[2]
  
  Model <- function(B0, B1, n) {B0 + B1*log(n)}
  
  Values <- data.frame("Size"=1:Size, "Fitted"=Model(B0=B0, B1=B1, n=1:Size))
  
  N <- which(Values$Fitted - (TargetAUC-0.02) >= 0)[1]
  
  N
  
}

Get_N_LOG_Th <- function(DF, TargetAUC, Size, Th) {
  
  test <- lm(data=DF, Curve ~ log(n))
  
  B0 <- test$coefficients[1]
  
  B1 <- test$coefficients[2]
  
  Model <- function(B0, B1, n) {B0 + B1*log(n)}
  
  Values <- data.frame("Size"=1:Size, "Fitted"=Model(B0=B0, B1=B1, n=1:Size))
  
  N <- which(Values$Fitted - (TargetAUC-Th) >= 0)[1]
  
  N
  
}

AUC_N <- function(DF, n) {
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5, c=0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  c <- summary(test)$coef[3,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:n)
  
  Values[n]
  
}

AUC_N_Fixed <- function(DF, FullAUC, n) {
  c <- FullAUC
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:n)
  
  Values[n]
  
}


AUC_N_LOG <- function(DF, n) {
  
  test <- lm(data=DF, Curve ~ log(n))
  
  B0 <- test$coefficients[1]
  
  B1 <- test$coefficients[2]
  
  Model <- function(B0, B1, n) {B0 + B1*log(n)}
  
  Values <- data.frame("Size"=1:n, "Fitted"=Model(B0=B0, B1=B1, n=1:n))
  
  Values[n,]$Fitted
  
}


#ML Functions -------------------------------------------------------------------------------------

Subsample_DF_random <- function(DF, DF_Response, n) {
  
  Reduced_Rows <- NULL
  
  for (i in 1:50) {
    Reduced_Rows <- append(Reduced_Rows, list(c(sample(rownames(DF), size=n, replace = F))
    ))
  }
  Reduced_Rows
}

XGBoost_Optimal <- function(DF, DF_Outcome, Test) {
  
  params <- list(booster="gbtree",
                 eta=0.3,
                 gamma=0,
                 max_depth=6,
                 lambda=0,
                 subsample=1)
    
    XGB_Example <- xgboost(data=as.matrix(DF[, -which(colnames(DF) %in% DF_Outcome)]),
                           label=DF[,colnames(DF) %in% DF_Outcome],
                           missing=NA, nrounds=25,
                           eval_metric="auc",verbose=F,
                           objective = "binary:logistic")
    
    Output <- predict(XGB_Example, newdata = as.matrix(Test[, -which(colnames(Test) %in% DF_Outcome)]),
                      type="response")
    
    AUC <- as.numeric(roc(as.factor(Test[,colnames(Test) %in% DF_Outcome]), Output, quiet=T)$auc)

    AUC
}

RF_Optimal_2 <- function(DF, DF_Outcome, Test) {
    
    RF_Example <- ranger(DF[,colnames(DF) %in% DF_Outcome] ~ ., data=DF[,-which(colnames(DF) %in% DF_Outcome)], probability = TRUE)
    
    Output <- predict(RF_Example, data = Test[, -which(colnames(Test) %in% DF_Outcome)],
                      type="response")$predictions[,2]
    
    AUC <- as.numeric(roc(as.factor(Test[,colnames(Test) %in% DF_Outcome]), Output, quiet=T)$auc)
    
    AUC
}

RF_Optimal <- function(DF, DF_Outcome) {
  
  AUC <- NULL
  
  Folds_inv <- create_folds(as.factor(DF[,colnames(DF) %in% DF_Outcome]), k=5, type="stratified", invert = T)
  
  for (j in 1:length(Folds_inv)) {
    HoldOut <- Folds_inv[[j]]
    
    RF_Example <- ranger(DF[-HoldOut,colnames(DF) %in% DF_Outcome] ~ ., data=DF[-HoldOut,-which(colnames(DF) %in% DF_Outcome)], probability = TRUE)
    
    Output <- predict(RF_Example, data = DF[HoldOut, -which(colnames(DF) %in% DF_Outcome)],
                      type="response")$predictions[,2]
    
    AUC <- append(AUC, as.numeric(roc(as.factor(DF[HoldOut,colnames(DF) %in% DF_Outcome]), Output, quiet=T)$auc))
  }
  c(mean(AUC), sd(AUC) / sqrt(5))
}

LR_Optimal <- function(DF, DF_Outcome) {
  
  AUC <- NULL
  
  Folds_inv <- create_folds(as.factor(DF[,colnames(DF) %in% DF_Outcome]), k=5, type="stratified", invert = T)
  
  for (j in 1:length(Folds_inv)) {
    HoldOut <- Folds_inv[[j]]
    
    LASSO_Example <- glmnet(x=DF[-HoldOut,-which(colnames(DF) %in% DF_Outcome)], 
                            y=DF[-HoldOut,colnames(DF) %in% DF_Outcome], 
                            family = "binomial", alpha = 1, lambda = 0)
    
    Output <- unname(predict(LASSO_Example, newx = as.matrix(DF[HoldOut, -which(colnames(DF) %in% DF_Outcome)]),
                             type="response")[,1])
    
    AUC <- append(AUC, as.numeric(roc(as.factor(DF[HoldOut,colnames(DF) %in% DF_Outcome]), Output, quiet=T)$auc))
  }
  c(mean(AUC), sd(AUC) / sqrt(5))
}

LR_Optimal_2 <- function(DF, DF_Outcome, Test) {

    LASSO_Example <- glmnet(x=DF[,-which(colnames(DF) %in% DF_Outcome)], 
                            y=DF[,colnames(DF) %in% DF_Outcome], 
                            family = "binomial", alpha = 1, lambda = 0)
    
    Output <- unname(predict(LASSO_Example, newx = as.matrix(Test[, -which(colnames(DF) %in% DF_Outcome)]),
                             type="response")[,1])
    
    AUC <- as.numeric(roc(as.factor(Test[,colnames(DF) %in% DF_Outcome]), Output, quiet=T)$auc)

    AUC
}




Filter_Pipeline <- function(DF, DF_Outcome, start, size, steps, n_predictors, True) {
  
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  p <- n_predictors
  
  for (n in samps) {
    
    PercentRecovered <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
    
    for (i in 1:25) {
      
      Reduced <- DF[Reduced_Rows[[i]],]
      
      #Step 1 - Filter based differential expression
      DiffEx <- Get_Important_Features(Reduced, DF_Response, transform=T, n_predictors)
      
      if (length(DiffEx)==0) {
        
        PercentRecovered <- append(PercentRecovered, 0)
      }
      
      else {
        
        PercentRecovered <- append(PercentRecovered, ifelse(mean(True %in% DiffEx)==1, 1, 0))
        
      }
    }
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve"=mean(PercentRecovered, na.rm=T))))
    
    
  }
  Curve
}

Filter_Pipeline_FAST <- function(DF, DF_Outcome, start, size, steps, n_predictors, True) {
  
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  p <- n_predictors
  
  for (n in samps) {
    
    PercentRecovered <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
    
    for (i in 1:1) {
      
      Reduced <- DF[Reduced_Rows[[i]],]
      
      #Step 1 - Filter based differential expression
      names <- colnames(Reduced)[1:p]
      pvalues <- p.adjust(sapply(1:p, function(k) {t.test(Reduced[,k] ~ Reduced[,which(colnames(Reduced) %in% DF_Outcome)])$p.value}), "fdr")
      
      DF_Diff_exp <- data.frame(names, pvalues)
      
      DiffEx <- DF_Diff_exp[DF_Diff_exp$pvalues<0.05, "names"]
      
      if (length(DiffEx)==0) {
        
        PercentRecovered <- append(PercentRecovered, 0)
      }
      
      else {
        
        PercentRecovered <- append(PercentRecovered, ifelse(mean(True %in% DiffEx)==1, 1, 0))
        
      }
    }
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve"=mean(PercentRecovered, na.rm=T))))
    
    
  }
  Curve
}


XGBoost_Pipeline_Curve <- function(DF, DF_Outcome, start, size, steps, True) {
  
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    PercentRecovered <- NULL
    
    FoldChange_Recovered <- NULL
    
    AUC <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=2*n)
    
    for (i in 1:50) {
      
      Reduced2 <- DF[Reduced_Rows[[i]][(1:n)],c(True, DF_Outcome)]
      
      params <- list(booster="gbtree",
                     eta=0.3,
                     gamma=0,
                     max_depth=6,
                     lambda=0,
                     subsample=1)
      
      Test <- DF[Reduced_Rows[[i]][((n+1):(2*n))],c(True, DF_Outcome)]
      
      XGB_Example <- xgboost(data=as.matrix(Reduced2[, -which(colnames(Reduced2) %in% DF_Outcome)]),
                             label=Reduced2[,colnames(Reduced2) %in% DF_Outcome],
                             missing=NA, nrounds=25,
                             eval_metric="auc",verbose=F,
                             objective = "binary:logistic")
      
      Output <- predict(XGB_Example, newdata = as.matrix(Test[, -which(colnames(Test) %in% DF_Outcome)]),
                        type="response")
      
      AUC <- append(AUC, as.numeric(roc(as.factor(Test[,colnames(Test) %in% DF_Outcome]), Output, quiet=T)$auc))
    }
    
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(AUC, na.rm=T))))
    
  }
  Curve
}

XGBoost_Pipeline_Curve_2 <- function(DF, DF_Outcome, start, size, steps, True, Test) {
  
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    PercentRecovered <- NULL
    
    AUC <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
    
    for (i in 1:50) {
      
      Reduced2 <- DF[Reduced_Rows[[i]][(1:n)],c(True, DF_Outcome)]
      
      params <- list(booster="gbtree",
                     eta=0.3,
                     gamma=0,
                     max_depth=6,
                     lambda=0,
                     subsample=1)
      
      Test <- Test
      
      XGB_Example <- xgboost(data=as.matrix(Reduced2[, -which(colnames(Reduced2) %in% DF_Outcome)]),
                             label=Reduced2[,colnames(Reduced2) %in% DF_Outcome],
                             missing=NA, nrounds=25,
                             eval_metric="auc",verbose=F,
                             objective = "binary:logistic")
      
      Output <- predict(XGB_Example, newdata = as.matrix(Test[, -which(colnames(Test) %in% DF_Outcome)]),
                        type="response")
      
      AUC <- append(AUC, as.numeric(roc(as.factor(Test[,colnames(Test) %in% DF_Outcome]), Output, quiet=T)$auc))
    }
    
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(AUC, na.rm=T))))
    
  }
  Curve
}

XGBoost_Curve_Evaluate <- function(DF, DF_Outcome, n) {
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
    Reduced <- DF[Reduced_Rows[[i]],]
    
    AUC <- NULL
    
    params <- list(booster="gbtree",
                   eta=0.3,
                   gamma=0,
                   max_depth=6,
                   lambda=0,
                   subsample=1)
    
    #if (n < 100) {
    
    #Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=2, type="stratified", invert = T)
    
    #}
    
    #else {
    
    Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=5, type="stratified", invert = T)
    
    
    #}
    
    for (j in 1:length(Folds_inv)) {
      HoldOut <- Folds_inv[[j]]
      
      XGB_Example <- xgboost(data=as.matrix(Reduced[-HoldOut, -which(colnames(Reduced) %in% DF_Outcome)]),
                             label=Reduced[-HoldOut,colnames(Reduced) %in% DF_Outcome],
                             missing=NA, nrounds=25,
                             eval_metric="auc",verbose=F,
                             objective = "binary:logistic")
      
      Output <- predict(XGB_Example, newdata = as.matrix(Reduced[HoldOut, -which(colnames(Reduced) %in% DF_Outcome)]),
                        type="response")
      
      AUC <- append(AUC, as.numeric(roc(as.factor(Reduced[HoldOut,colnames(Reduced) %in% DF_Outcome]), Output, quiet=T)$auc))
    }
    
    ComboAUC <- append(ComboAUC, mean(AUC))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}

RF_Curve_Evaluate <- function(DF, DF_Outcome, n) {
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:50) {
    
    Reduced <- DF[Reduced_Rows[[i]],]
    
    AUC <- NULL
    
    #if (n < 100) {
    
    #Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=2, type="stratified", invert = T)
    
    #}
    
    #else {
    
    Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=5, type="stratified", invert = T)
    
    
    #}
    
    for (j in 1:length(Folds_inv)) {
      HoldOut <- Folds_inv[[j]]
      
      RF_Example <- ranger(Reduced[-HoldOut,colnames(Reduced) %in% DF_Outcome] ~ ., data=Reduced[-HoldOut,-which(colnames(Reduced) %in% DF_Outcome)], probability = TRUE)
      
      Output <- predict(RF_Example, data = Reduced[HoldOut, -which(colnames(Reduced) %in% DF_Outcome)],
                        type="response")$predictions[,2]
      AUC <- append(AUC, as.numeric(roc(as.factor(Reduced[HoldOut,colnames(Reduced) %in% DF_Outcome]), Output, quiet=T)$auc))
    }
    
    ComboAUC <- append(ComboAUC, mean(AUC))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}


RF_Pipeline_Curve <- function(DF, DF_Outcome, start, size, steps, True) {
  
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    AUC <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=2*n)
    
    for (i in 1:50) {
      
      Reduced2 <- DF[Reduced_Rows[[i]][(1:n)],c(True, DF_Outcome)]
      
      Test <- DF[Reduced_Rows[[i]][((n+1):(2*n))],c(True, DF_Outcome)]
      
      RF_Example <- ranger(Reduced2[,colnames(Reduced2) %in% DF_Outcome] ~ ., data=Reduced2[,-which(colnames(Reduced2) %in% DF_Outcome)], probability = TRUE)
      
      Output <- predict(RF_Example, data = Test[, -which(colnames(Test) %in% DF_Outcome)],
                        type="response")$predictions[,2]
      
      AUC <- append(AUC, as.numeric(roc(as.factor(Test[,colnames(Test) %in% DF_Outcome]), Output, quiet=T)$auc))
      
    }
    
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(AUC, na.rm=T))))
    
  }
  Curve
}

RF_Pipeline_Curve_2 <- function(DF, DF_Outcome, start, size, steps, True, Test) {
  
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    PercentRecovered <- NULL
    
    AUC <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
    
    for (i in 1:50) {
      
      Reduced2 <- DF[Reduced_Rows[[i]][(1:n)],c(True, DF_Outcome)]
      
      Test <- Test
      
      RF_Example <- ranger(Reduced2[,colnames(Reduced2) %in% DF_Outcome] ~ ., data=Reduced2[,-which(colnames(Reduced2) %in% DF_Outcome)], probability = TRUE)
      
      Output <- predict(RF_Example, data = Test[, -which(colnames(Test) %in% DF_Outcome)],
                        type="response")$predictions[,2]
      
      AUC <- append(AUC, as.numeric(roc(as.factor(Test[,colnames(Test) %in% DF_Outcome]), Output, quiet=T)$auc))
    }
    
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(AUC, na.rm=T))))
    
  }
  Curve
}



NN_Optimal <- function(DF, DF_Outcome, h) {
  
  DF <- data.frame(cbind(scale(log2(DF[,-ncol(DF)]+0.01)), "Group"=DF[,which(colnames(DF) %in% DF_Outcome)]))
  
  DF$Group <- as.factor(DF$Group)
  
  DF <- as.h2o(DF)
  
  NN_Example <- h2o.deeplearning(y="Group",
                                 training_frame = DF, nfolds=5, standardize = F, fold_assignment="Stratified",
                                 keep_cross_validation_predictions=T,  reproducible = TRUE, hidden = c(h), epochs = 10,
                                 seed = 1)
  
  h2o.auc(NN_Example, valid=T)
}

NN_Optimal_2 <- function(DF, DF_Outcome, Test, h) {
  
  DF <- data.frame(cbind(scale((DF[,-ncol(DF)])), "Group"=DF[,which(colnames(DF) %in% DF_Outcome)]))
  
  DF$Group <- as.factor(DF$Group)
  
  DF <- as.h2o(DF)
  
  Test <- data.frame(cbind(scale((Test[,-ncol(Test)])), "Group"=Test[,which(colnames(Test) %in% DF_Outcome)]))
  
  Test$Group <- as.factor(Test$Group)
  
  Test <- as.h2o(Test)
  
  NN_Example <- h2o.deeplearning(y="Group",
                                 training_frame = DF, validation_frame = Test, standardize = F, 
                                 reproducible = F, hidden = c(h), epochs = 10, seed=1)
  
  h2o.auc(NN_Example, valid=T)
}

NN_Curve_Evaluate <- function(DF, DF_Outcome, n, True, h) {
  
  DF <- data.frame(cbind(scale(log2(DF[,-ncol(DF)]+0.01)), "Group"=DF[,which(colnames(DF) %in% DF_Outcome)]))
  
  DF$Group <- as.factor(DF$Group)
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n*2)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
    Reduced2 <- DF[Reduced_Rows[[i]][(1:n)],c(True, DF_Outcome)]
    
    Test <- Test <- DF[Reduced_Rows[[i]][((n+1):(2*n))],c(True, DF_Outcome)]
    
    #if (n < 100) {
    
    #Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=2, type="stratified", invert = T)
    
    #}
    
    #else {
    
    Reduced2 <- as.h2o(Reduced2)
    
    Test <- as.h2o(Test)
    
    NN_Example <- h2o.deeplearning(y="Group",
                                   training_frame = Reduced2, validation_frame = Test, standardize = F, 
                                   reproducible = FALSE, hidden = c(h), epochs = 10)
    
    ComboAUC <- append(ComboAUC, h2o.auc(NN_Example, valid = T))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}

NN_Curve_Evaluate_2 <- function(DF, DF_Outcome, n, True, Test, h) {
  
  DF <- data.frame(cbind(scale((DF[,-ncol(DF)])), "Group"=DF[,which(colnames(DF) %in% DF_Outcome)]))
  
  DF$Group <- as.factor(DF$Group)
  
  Test <- data.frame(cbind(scale((Test[,-ncol(Test)])), "Group"=Test[,which(colnames(Test) %in% DF_Outcome)]))
  
  Test$Group <- as.factor(Test$Group)
  
  Test <- as.h2o(Test)
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
    Reduced2 <- DF[Reduced_Rows[[i]],c(True, DF_Outcome)]
    
    Reduced2 <- as.h2o(Reduced2)
    
    NN_Example <- h2o.deeplearning(y="Group",
                                   training_frame = Reduced2, validation_frame = Test, standardize = F, 
                                   reproducible = FALSE, hidden = c(h), epochs = 10)
    
    ComboAUC <- append(ComboAUC, h2o.auc(NN_Example, valid = T))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}


#INPUT FOR THIS SHOULD BE PRE-PROCESSED ORIGINAL DATASET WITH RAW COUNTS
Feature_Curve <- function(DF, a, b, stepsize) {
  
  #Simulate DF from raw counts
  Sim <- Simulate_OMICs_Data_SPLIT(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                               p=ncol(DF)-1, n=5000)
  
  
  #Back-transform to neg-bin distribution
  Sim <- cbind((2^Sim[,-ncol(Sim)])-0.01, "Group"=Sim$Group)
  
  Sim[,1:(ncol(Sim)-1)] <- round(Sim[,1:(ncol(Sim)-1)], 0)

  
  #DESEQ Curves
  a <- start
  b <- size
  
  stepsize <- round((b - a) / steps)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > steps) {samps <- samps[1:steps]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    FeatureIndicator <- NULL
    
    Reduced_Rows <- Subsample_DF_random(DF=Sim, DF_Response="Group", n=n)
    
    for (i in 1:10) {
      
    Reduced <- Sim[Reduced_Rows[[i]],]
      
    ExpDesign <- data.frame(row.names=rownames(Reduced),
                              Group = Reduced[,"Group"])
      
    DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Reduced[,-which(colnames(Reduced) %in% "Group")])), 
                                             colData=ExpDesign, design=~Group))
      
      
    DiffEx <- results(DiffEx, alpha=0.05)
    DiffEx <- na.omit(DiffEx)

    FeatureIndicator <- ifelse((sum(Keep %in% rownames(DiffEx[DiffEx$padj<0.05,])) > length(Keep)*0.95) == T, 1, 0)
      
    }
    
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(FeatureIndicator, na.rm=T))))
    
  }
  Curve
  
  
}


AvgReadCount <- function(DS) {
  AvgCounts <- apply(DS[,-which(colnames(DS) %in% "Group")], 2, mean)
  AvgCounts
  
}


#Functions to get DS Characterisitcs ---------------------------------------------------------------

#DISPERSION

#START POINT HERE SHOULD ALWAYS BE RAW COUNTS
GetDispersion <- function(DF, UseWhole, Keep) {
  
  if(UseWhole==F) {
  
  DF <- DF[,c(Keep, "Group")]
    
  }
  
  ExpDesign <- data.frame(row.names=rownames(DF))
  
  DiffEx <- estimateSizeFactors(DESeqDataSetFromMatrix(countData = as.matrix(t(DF[,-which(colnames(DF) %in% "Group")])), 
                                         colData=ExpDesign, design=~1))
  
  
  DIS <- dispersions(estimateDispersions(DiffEx))
  
  DIS <- DIS[which(colnames(DF) %in% Keep)]
  DIS
}


#START POINT HERE SHOULD ALWAYS BE RAW COUNTS
GetDispersion_Sim <- function(DF, UseWhole, Keep) {
  
  Sim <- Simulate_OMICs_Data_BNG_Only(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                                      p=ncol(DF)-1, n=5000)
  
  #Simulate DF from raw counts
  if(UseWhole==T) {
  
    Sim <- Simulate_OMICs_Data_SPLIT(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                                        p=ncol(DF)-1, n=5000)
  
  }
  
  
  #Back-transform to neg-bin distribution
  Sim <- cbind((2^Sim[,-ncol(Sim)])-0.01, "Group"=Sim$Group)
  
  Sim[,1:(ncol(Sim)-1)] <- round(Sim[,1:(ncol(Sim)-1)], 0)
  
  ExpDesign <- data.frame(row.names=rownames(Sim))
  
  DiffEx <- estimateSizeFactors(DESeqDataSetFromMatrix(countData = as.matrix(t(Sim[,-which(colnames(Sim) %in% "Group")])), 
                                                       colData=ExpDesign, design=~1))
  
  
  DIS <- dispersions(estimateDispersions(DiffEx))
  
  DIS <- DIS[which(colnames(Sim) %in% Keep)]
  DIS
}



#CORRELATION

#START POINT HERE SHOULD ALWAYS BE UN-TRANSFORMED DESEQ2 NORMALIZED COUNTS
GetCorrelation <- function(DF, Keep) {
  
  DF <- DF[,c(Keep, "Group")]
  
  DF <- cbind(log2(DF[,-ncol(DF)]+0.01), "Group"=DF$Group) #Transform to normal distribution
  
  Cor <- abs(cor(DF[,1:(ncol(DF)-1)]))
  
  Cor[lower.tri(Cor)]
  
}


GetLFC <- function(DF, Keep) {
    
  DF <- DF[,c(Keep, "Group")]
  
  DF <- cbind(log2(DF[,-ncol(DF)]+0.01), "Group"=DF$Group)
  
  LFC_True <- sapply(1:(ncol(DF)-1), function(i) {abs(mean(DF[DF$Group==1,i]) - mean(DF[DF$Group==0,i]))})
  
  c("median"=median(LFC_True), "min"=min(LFC_True), "max"=max(LFC_True))
  
}


#Data----------------------------------------------------------------------------------------------

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Raw Count Data/Aim 2 Raw Count Data")

#NAFLD
NAFLD <- read.csv("GSE167523_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE167523")
NAFLD <- NAFLD[,-1]
NAFLD <- data.frame(t(NAFLD))

NAFLD$Group <- ifelse(Meta$GSE167523_series_matrix.txt.gz$`disease subtype:ch1`=="NASH", 1, 0)

Vars <- apply(NAFLD[,-ncol(NAFLD)], 2, var)
NAFLD <- NAFLD[,-which(colnames(NAFLD) %in% names(Vars[Vars==0]))]
NAFLD <- data.frame(apply(NAFLD, 2, as.numeric))



ExpDesign <- data.frame(row.names=rownames(NAFLD),
                        Group = NAFLD[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(NAFLD[,-which(colnames(NAFLD) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))

DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

NAFLD[,1:(ncol(NAFLD)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=NAFLD[,-which(colnames(NAFLD) %in% "Group")], y=as.factor(NAFLD[,which(colnames(NAFLD) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]


DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

DIS <- DIS[which(colnames(NAFLD) %in% Keep)]
median(DIS)

mean(NAFLD$Group)

NAFLD_Keep <- NAFLD[,c(Keep, "Group")]

NAFLD_Keep <- cbind(log2(NAFLD_Keep[,-ncol(NAFLD_Keep)]+0.01), "Group"=NAFLD_Keep$Group)

LFC_True <- sapply(1:(ncol(NAFLD_Keep)-1), function(i) {abs(mean(NAFLD_Keep[NAFLD_Keep$Group==1,i]) - mean(NAFLD_Keep[NAFLD_Keep$Group==0,i]))})



set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(NAFLD, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(NAFLD[which(NAFLD$Group==1),])), 
                                     p=ncol(NAFLD)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
LFC_Sim
plot(LFC_Sim, LFC_True)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

#DISPERSION on the simulated Data---------------
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)


Test_Round <- Test
Test_Round[,1:(ncol(Test_Round)-1)] <- round(Test_Round[,1:(ncol(Test_Round)-1)] , 0)

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test_Round[,-which(colnames(Test_Round) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DIS_Sim <- dispersions(DiffEx)
plot(DIS_Sim, DIS)
median(DIS_Sim)

#-----------------------------------------------


set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

NN_Optimal(Test, "Group", 3*52)

#Curves
set.seed(2024)
NAFLD_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)
set.seed(2024)
NAFLD_RF <- RF_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

plot_curve_full(NAFLD_XGB, "PowerLaw", 0.999, 250)
plot_curve_full(NAFLD_RF, "PowerLaw", 0.999, 250)

Get_N(NAFLD_XGB, 0.999, 250)
Get_N(NAFLD_RF, 0.999, 250)


#HCC
HCC <- read.csv("GSE114564_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE114564")
HCC <- HCC[,-1]
HCC <- data.frame(t(HCC))

HCC$Group <- Meta$GSE114564_series_matrix.txt.gz$`disease state:ch1`
HCC <- HCC[HCC$Group %in% c("Early HCC", "Advanced HCC"),]
HCC$Group <- ifelse(HCC$Group=="Advanced HCC", 1, 0)

Vars <- apply(HCC[,-ncol(HCC)], 2, var)
HCC <- HCC[,-which(colnames(HCC) %in% names(Vars[Vars==0]))]
HCC <- data.frame(apply(HCC, 2, as.numeric))



1-mean(HCC$Group)

ExpDesign <- data.frame(row.names=rownames(HCC),
                        Group = HCC[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(HCC[,-which(colnames(HCC) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

HCC[,1:(ncol(HCC)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=HCC[,-which(colnames(HCC) %in% "Group")], y=as.factor(HCC[,which(colnames(HCC) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

HCC_Keep <- HCC[,c(Keep, "Group")]

Vars_True <- apply(HCC_Keep[,-ncol(HCC_Keep)], 2, var)

HCC_Keep <- cbind(log2(HCC_Keep[,-ncol(HCC_Keep)]+0.01), "Group"=HCC_Keep$Group)

LFC_True <- sapply(1:(ncol(HCC_Keep)-1), function(i) {abs(mean(HCC_Keep[HCC_Keep$Group==1,i]) - mean(HCC_Keep[HCC_Keep$Group==0,i]))})



set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(HCC, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(HCC[which(HCC$Group==1),])), 
                                     p=ncol(HCC)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
LFC_Sim
cor(LFC_Sim, LFC_True)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "HCC_Sim.csv")



set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

NN_Optimal(Test, "Group", 9*3)

set.seed(2024)
HCC_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

set.seed(2024)
HCC_RF <- RF_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

plot_curve_full(HCC_XGB, "PowerLaw", 0.995, 250)
plot_curve_full(HCC_RF, "PowerLaw", 0.994, 250)

Get_N(HCC_XGB, 0.995, 250)
Get_N(HCC_RF, 0.994, 250)

#Kidney
  Kidney <- read.csv("GSE124685_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
  Meta <- getGEO("GSE124685")
  Kidney <- Kidney[,-1]
  Kidney <- data.frame(t(Kidney))
  
  Kidney$Group <- ifelse(Meta$GSE124685_series_matrix.txt.gz$`disease state:ch1`=="IPF", 1, 0)

Vars <- apply(Kidney[,-ncol(Kidney)], 2, var)
Kidney <- Kidney[,-which(colnames(Kidney) %in% names(Vars[Vars==0]))]
Kidney <- data.frame(apply(Kidney, 2, as.numeric))



ExpDesign <- data.frame(row.names=rownames(Kidney),
                        Group = Kidney[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Kidney[,-which(colnames(Kidney) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Kidney[,1:(ncol(Kidney)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Kidney[,-which(colnames(Kidney) %in% "Group")], y=as.factor(Kidney[,which(colnames(Kidney) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])
1-mean(Kidney$Group)

Kidney_Keep <- Kidney[,c(Keep, "Group")]

Vars_True <- apply(Kidney_Keep[,-ncol(Kidney_Keep)], 2, var)

Kidney_Keep <- cbind(log2(Kidney_Keep[,-ncol(Kidney_Keep)]+0.01), "Group"=Kidney_Keep$Group)

LFC_True <- sapply(1:(ncol(Kidney_Keep)-1), function(i) {abs(mean(Kidney_Keep[Kidney_Keep$Group==1,i]) - mean(Kidney_Keep[Kidney_Keep$Group==0,i]))})


set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Kidney, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Kidney[which(Kidney$Group==1),])), 
                                     p=ncol(Kidney)-1, n=5000)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "IPF_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
LFC_Sim
max(LFC_Sim)
plot(LFC_Sim, LFC_True)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
Kidney_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 25, 150, steps=10, True=Keep)

set.seed(2024)
Kidney_RF <- RF_Pipeline_Curve(Test, "Group", 10, 50, steps=10, True=Keep)

plot_curve_full(Kidney_XGB, "PowerLaw_Fixed", 1, 150)
Get_N_Fixed(Kidney_XGB, 1, 1000)

plot_curve_full(Kidney_RF, "LOG", 1, 50)
Get_N_LOG(Kidney_RF, 1, 150)




#IPF
IPF <- read.csv("GSE218048_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE218048")
IPF <- IPF[,-1]
IPF <- data.frame(t(IPF))

IPF$Group <- rep(0, nrow(IPF))
IPF[base::grep("c", Meta$GSE218048_series_matrix.txt.gz$title, "ifta"),"Group"] <- 1


Vars <- apply(IPF[,-ncol(IPF)], 2, var)
IPF <- IPF[,-which(colnames(IPF) %in% names(Vars[Vars==0]))]
IPF <- data.frame(apply(IPF, 2, as.numeric))



ExpDesign <- data.frame(row.names=rownames(IPF),
                        Group = as.factor(IPF[,"Group"]))

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(IPF[,-which(colnames(IPF) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

IPF[,1:(ncol(IPF)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=IPF[,-which(colnames(IPF) %in% "Group")], y=as.factor(IPF[,which(colnames(IPF) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])
mean(IPF$Group)

IPF_Keep <- IPF[,c(Keep, "Group")]

Vars_True <- apply(IPF_Keep[,-ncol(IPF_Keep)], 2, var)

IPF_Keep <- cbind(log2(IPF_Keep[,-ncol(IPF_Keep)]+0.01), "Group"=IPF_Keep$Group)

LFC_True <- sapply(1:(ncol(IPF_Keep)-1), function(i) {abs(mean(IPF_Keep[IPF_Keep$Group==1,i]) - mean(IPF_Keep[IPF_Keep$Group==0,i]))})


set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(IPF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(IPF[which(IPF$Group==1),])), 
                                     p=ncol(IPF)-1, n=5000)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "Kidney_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
LFC_Sim
median(LFC_Sim)
plot(LFC_Sim, LFC_True)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)
Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
IPF_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

set.seed(2024)
IPF_RF <- RF_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

plot_curve_full(IPF_XGB, "LOG", 0.848, 1000)
Get_N_LOG(IPF_XGB, 0.848, 1000)

plot_curve_full(IPF_RF, "LOG", 0.850, 1000)
Get_N_LOG(IPF_RF, 0.850, 1000)


#Tuberculosis
Tub <- read.csv("GSE89403_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE89403")
Tub <- Tub[,-1]
Tub <- data.frame(t(Tub))

Tub <- Tub[which(Meta$GSE89403_series_matrix.txt.gz$`time:ch1`=="DX"),]
Tub$Group <- Meta$GSE89403_series_matrix.txt.gz$`treatmentresult:ch1`[which(Meta$GSE89403_series_matrix.txt.gz$`time:ch1`=="DX")]

Tub$Group <- ifelse(Tub$Group=="NA", NA, ifelse(Tub$Group=="Not Cured", 0, 1))
Tub <- na.omit(Tub)

Vars <- apply(Tub[,-ncol(Tub)], 2, var)
Tub <- Tub[,-which(colnames(Tub) %in% names(Vars[Vars==0]))]
Tub <- data.frame(apply(Tub, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(Tub),
                        Group = Tub[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Tub[,-which(colnames(Tub) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Tub[,1:(ncol(Tub)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Tub[,-which(colnames(Tub) %in% "Group")], y=as.factor(Tub[,which(colnames(Tub) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DIS <- dispersions(DiffEx)

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

DIS <- DIS[which(colnames(Tub) %in% Keep)]
median(DIS)

Tub_Keep <- Tub[,c(Keep, "Group")]

Vars_True <- apply(Tub_Keep[,-ncol(Tub_Keep)], 2, var)

Tub_Keep <- cbind(log2(Tub_Keep[,-ncol(Tub_Keep)]+0.01), "Group"=Tub_Keep$Group)

LFC_True <- sapply(1:(ncol(Tub_Keep)-1), function(i) {abs(mean(Tub_Keep[Tub_Keep$Group==1,i]) - mean(Tub_Keep[Tub_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Tub, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Tub[which(Tub$Group==1),])), 
                                     p=ncol(Tub)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})

cor(LFC_Sim, LFC_True)
LFC_Sim


#DISPERSION on the simulated Data---------------
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "Tuberculosis_Sim.csv")

Test_Round <- Test
Test_Round[,1:(ncol(Test_Round)-1)] <- round(Test_Round[,1:(ncol(Test_Round)-1)] , 0)

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test_Round[,-which(colnames(Test_Round) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DIS_Sim <- dispersions(DiffEx)
cor(DIS_Sim, DIS)
max(DIS_Sim)

#-----------------------------------------------


set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2025)
Tub_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 75, 500, steps=10, True=Keep)

set.seed(2025)
Tub_RF <- RF_Pipeline_Curve(Test, "Group", 75, 500, steps=10, True=Keep)

plot_curve_full(Tub_XGB, "PowerLaw_Fixed", 0.995, 500)
plot_curve_full(Tub_RF, "PowerLaw", 0.992, 500)

Get_N_Fixed(Tub_XGB, 0.995, 500)
Get_N(Tub_RF, 0.992, 500)

#AVSC
AVSC <- read.csv("GSE218474_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE218474")
AVSC <- AVSC[,-1]
AVSC <- data.frame(t(AVSC))

Meta <- data.frame("ID"=Meta$GSE218474_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE218474_series_matrix.txt.gz$`avsc:ch1`=="Yes", 1, 0))

AVSC$ID <- rownames(AVSC)

AVSC <- left_join(AVSC, Meta, "ID")
AVSC <- AVSC[,-39377]


Vars <- apply(AVSC[,-ncol(AVSC)], 2, var)
AVSC <- AVSC[,-which(colnames(AVSC) %in% names(Vars[Vars==0]))]
AVSC <- data.frame(apply(AVSC, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(AVSC),
                        Group = AVSC[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(AVSC[,-which(colnames(AVSC) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

AVSC[,1:(ncol(AVSC)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=AVSC[,-which(colnames(AVSC) %in% "Group")], y=as.factor(AVSC[,which(colnames(AVSC) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])


AVSC_Keep <- AVSC[,c(Keep, "Group")]


Vars_True <- apply(AVSC_Keep[,-ncol(AVSC_Keep)], 2, var)

AVSC_Keep <- cbind(log2(AVSC_Keep[,-ncol(AVSC_Keep)]+0.01), "Group"=AVSC_Keep$Group)

LFC_True <- sapply(1:(ncol(AVSC_Keep)-1), function(i) {abs(mean(AVSC_Keep[AVSC_Keep$Group==1,i]) - mean(AVSC_Keep[AVSC_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(AVSC, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(AVSC[which(AVSC$Group==1),])), 
                                     p=ncol(AVSC)-1, n=5000)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "AVSC_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
plot(LFC_Sim, LFC_True)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

NN_Optimal(Test, "Group", 3*3)

set.seed(2024)
AVSC_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 1500, steps=10, True=Keep)

set.seed(2024)
AVSC_RF <- RF_Pipeline_Curve(Test, "Group", 50, 1500, steps=10, True=Keep)

plot_curve_full(AVSC_XGB, "PowerLaw", 0.852, 1500)

plot_curve_full(AVSC_RF, "PowerLaw", 0.849, 1500)

Get_N(AVSC_RF, 0.852, 1500)
Get_N(AVSC_XGB, 0.852, 1500)

#RA
RA <- read.csv("GSE117769_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE117769")
RA <- RA[,-1]
RA <- data.frame(t(RA))

Meta <- data.frame("ID"=Meta$GSE117769_series_matrix.txt.gz$geo_accession,
                   "Group"=Meta$GSE117769_series_matrix.txt.gz$`phenotype:ch1`)

RA$ID <- rownames(RA)

RA <- left_join(RA, Meta, "ID")
RA <- RA[,-39377]

RA <- RA[RA$Group %in% c("Healthy", "Rheumatoid arthritis"),]
RA$Group <- ifelse(RA$Group=="Healthy", 0, 1)

Vars <- apply(RA[,-ncol(RA)], 2, var)
RA <- RA[,-which(colnames(RA) %in% names(Vars[Vars==0]))]
RA <- data.frame(apply(RA, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(RA),
                        Group = RA[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(RA[,-which(colnames(RA) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

RA[,1:(ncol(RA)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=RA[,-which(colnames(RA) %in% "Group")], y=as.factor(RA[,which(colnames(RA) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

RA_Keep <- RA[,c(Keep, "Group")]

Vars_True <- apply(RA_Keep[,-ncol(RA_Keep)], 2, var)

RA_Keep <- cbind(log2(RA_Keep[,-ncol(RA_Keep)]+0.01), "Group"=RA_Keep$Group)

LFC_True <- sapply(1:(ncol(RA_Keep)-1), function(i) {abs(mean(RA_Keep[RA_Keep$Group==1,i]) - mean(RA_Keep[RA_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(RA, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(RA[which(RA$Group==1),])), 
                                     p=ncol(RA)-1, n=5000)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "RA_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
plot(LFC_Sim, LFC_True)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)


set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
RA_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 750, steps=10, True=Keep)

set.seed(2024)
RA_RF <- RF_Pipeline_Curve(Test, "Group", 50, 750, steps=10, True=Keep)

plot_curve_full(RA_XGB, "PowerLaw", 0.981, 750)
plot_curve_full(RA_RF, "PowerLaw", 0.977, 750)

Get_N(RA_XGB, 0.981, 750)
Get_N(RA_RF, 0.977, 750)


#MISC
MISC <- read.csv("GSE178491_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE178491")
MISC <- MISC[,-1]
MISC <- data.frame(t(MISC))

Meta <- data.frame("ID"=Meta$`GSE178491-GPL20301_series_matrix.txt.gz`$geo_accession,
                   "Group"=Meta$`GSE178491-GPL20301_series_matrix.txt.gz`$`disease:ch1`)

MISC$ID <- rownames(MISC)

MISC <- left_join(MISC, Meta, "ID")
MISC <- MISC[,-39377]
MISC <- na.omit(MISC)

MISC$Group <- ifelse(MISC$Group=="Febrile control", 0, 1)

Vars <- apply(MISC[,-ncol(MISC)], 2, var)
MISC <- MISC[,-which(colnames(MISC) %in% names(Vars[Vars==0]))]
MISC <- data.frame(apply(MISC, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(MISC),
                        Group = MISC[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(MISC[,-which(colnames(MISC) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

MISC[,1:(ncol(MISC)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=MISC[,-which(colnames(MISC) %in% "Group")], y=as.factor(MISC[,which(colnames(MISC) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

MISC_Keep <- MISC[,c(Keep, "Group")]

Vars_True <- apply(MISC_Keep[,-ncol(MISC_Keep)], 2, var)

MISC_Keep <- cbind(log2(MISC_Keep[,-ncol(MISC_Keep)]+0.01), "Group"=MISC_Keep$Group)

LFC_True <- sapply(1:(ncol(MISC_Keep)-1), function(i) {abs(mean(MISC_Keep[MISC_Keep$Group==1,i]) - mean(MISC_Keep[MISC_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(MISC, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(MISC[which(MISC$Group==1),])), 
                                     p=ncol(MISC)-1, n=5000)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "MISC_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})

plot(LFC_Sim, LFC_True, xlim=c(-0.5, 5), ylim=c(-0.5, 5), pch=19, xlab="Simulated Effect Sizes (LFC)", ylab="True Effect Sizes (LFC)") 
lines(x=c(0,4.25), y=c(0,4.25), col="blue", lwd=2)

cor(scale(MISC_Keep[,-42]))

pheatmap(cor(scale(MISC_Keep[,-42])), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
pheatmap(cor(scale(Test[,-42])), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)

max(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
MISC_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

set.seed(2024)
MISC_RF <- RF_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

plot_curve_full(MISC_XGB, "PowerLaw_Fixed", 0.999, 250)
plot_curve_full(MISC_RF, "PowerLaw_Fixed", 0.998, 250)

Get_N_Fixed(MISC_XGB, 0.999, 250)
Get_N_Fixed(MISC_RF, 0.998, 250)


#Test DS
Cancer <- read.csv("GSE183635_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE183635")


Meta <- data.frame("ID"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$geo_accession, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$geo_accession),
                   "Group"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$`patient group:ch1`, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$`patient group:ch1`))


Cancer <- Cancer[,-1]
Cancer <- data.frame(t(Cancer))

Cancer$ID <- rownames(Cancer)

Cancer <- left_join(Cancer, Meta, "ID")
Cancer <- Cancer[,-39377]

Cancer <- Cancer[Cancer$Group %in% c("Asymptomatic Controls", "Multiple Sclerosis"),]
Cancer$Group <- ifelse(Cancer$Group=="Asymptomatic Controls", 0, 1)

Vars <- apply(Cancer[,-ncol(Cancer)], 2, var)
Cancer <- Cancer[,-which(colnames(Cancer) %in% names(Vars[Vars==0]))]
MS <- data.frame(apply(Cancer, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(MS),
                        Group = MS[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(MS[,-which(colnames(MS) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

MS[,1:(ncol(MS)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Cancer[,-which(colnames(Cancer) %in% "Group")], y=as.factor(Cancer[,which(colnames(Cancer) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Cancer_Keep <- Cancer[,c(Keep, "Group")]

Vars_True <- apply(Cancer_Keep[,-ncol(Cancer_Keep)], 2, var)

Cancer_Keep <- cbind(log2(Cancer_Keep[,-ncol(Cancer_Keep)]+0.01), "Group"=Cancer_Keep$Group)

LFC_True <- sapply(1:(ncol(Cancer_Keep)-1), function(i) {abs(mean(Cancer_Keep[Cancer_Keep$Group==1,i]) - mean(Cancer_Keep[Cancer_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Cancer, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Cancer[which(Cancer$Group==1),])), 
                                     p=ncol(Cancer)-1, n=5000)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "MS_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
plot(LFC_Sim, LFC_True)
max(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
MS_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
MS_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(MS_XGB, "PowerLaw", 0.883, 2000)
Get_N(MS_XGB, 0.883, 2000)

plot_curve_full(MS_RF, "PowerLaw", 0.883, 2000)
Get_N(MS_RF, 0.883, 2000)

#CANCER 2
Cancer <- read.csv("GSE183635_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE183635")


Meta <- data.frame("ID"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$geo_accession, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$geo_accession),
                   "Group"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$`patient group:ch1`, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$`patient group:ch1`))


Cancer <- Cancer[,-1]
Cancer <- data.frame(t(Cancer))

Cancer$ID <- rownames(Cancer)

Cancer <- left_join(Cancer, Meta, "ID")
Cancer <- Cancer[,-39377]

Cancer <- Cancer[Cancer$Group %in% c("Pulmonary Hypertension", "Angina Pectoris"),]
Cancer$Group <- ifelse(Cancer$Group=="Pulmonary Hypertension", 0, 1)

Vars <- apply(Cancer[,-ncol(Cancer)], 2, var)
Cancer <- Cancer[,-which(colnames(Cancer) %in% names(Vars[Vars==0]))]
Hypertension <- data.frame(apply(Cancer, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(Hypertension),
                        Group = Hypertension[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Hypertension[,-which(colnames(Hypertension) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Hypertension[,1:(ncol(Hypertension)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Cancer[,-which(colnames(Cancer) %in% "Group")], y=as.factor(Cancer[,which(colnames(Cancer) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Cancer_Keep <- Cancer[,c(Keep, "Group")]

Vars_True <- apply(Cancer_Keep[,-ncol(Cancer_Keep)], 2, var)

Cancer_Keep <- cbind(log2(Cancer_Keep[,-ncol(Cancer_Keep)]+0.01), "Group"=Cancer_Keep$Group)

LFC_True <- sapply(1:(ncol(Cancer_Keep)-1), function(i) {abs(mean(Cancer_Keep[Cancer_Keep$Group==1,i]) - mean(Cancer_Keep[Cancer_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Cancer, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Cancer[which(Cancer$Group==1),])), 
                                     p=ncol(Cancer)-1, n=5000)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "Hypertension_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
max(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
HTN_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

set.seed(2024)
HTN_RF <- RF_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

plot_curve_full(HTN_XGB, "PowerLaw_Fixed", 0.999, 250)
Get_N_Fixed(HTN_XGB, 0.999, 250)

plot_curve_full(HTN_RF, "PowerLaw_Fixed", 0.999, 250)
Get_N_Fixed(HTN_RF, 0.999, 250)


#CANCER 3
Cancer <- read.csv("GSE183635_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE183635")


Meta <- data.frame("ID"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$geo_accession, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$geo_accession),
                   "Group"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$`patient group:ch1`, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$`patient group:ch1`))


Cancer <- Cancer[,-1]
Cancer <- data.frame(t(Cancer))

Cancer$ID <- rownames(Cancer)

Cancer <- left_join(Cancer, Meta, "ID")
Cancer <- Cancer[,-39377]

Cancer <- Cancer[Cancer$Group %in% c("Cholangiocarcinoma", "Colorectal Cancer"),]
Cancer$Group <- ifelse(Cancer$Group=="Colorectal Cancer", 0, 1)

Vars <- apply(Cancer[,-ncol(Cancer)], 2, var)
Cancer <- Cancer[,-which(colnames(Cancer) %in% names(Vars[Vars==0]))]
CCA <- data.frame(apply(Cancer, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(CCA),
                        Group = CCA[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(CCA[,-which(colnames(CCA) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

CCA[,1:(ncol(CCA)-1)] <- DDS


set.seed(2024)
B <- Boruta(x=Cancer[,-which(colnames(Cancer) %in% "Group")], y=as.factor(Cancer[,which(colnames(Cancer) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Cancer_Keep <- Cancer[,c(Keep, "Group")]

Vars_True <- apply(Cancer_Keep[,-ncol(Cancer_Keep)], 2, var)

Cancer_Keep <- cbind(log2(Cancer_Keep[,-ncol(Cancer_Keep)]+0.01), "Group"=Cancer_Keep$Group)

LFC_True <- sapply(1:(ncol(Cancer_Keep)-1), function(i) {abs(mean(Cancer_Keep[Cancer_Keep$Group==1,i]) - mean(Cancer_Keep[Cancer_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Cancer, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Cancer[which(Cancer$Group==1),])), 
                                     p=ncol(Cancer)-1, n=5000)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "CCA_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
plot(LFC_Sim, LFC_True)
max(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
CCA_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

set.seed(2024)
CCA_RF <- RF_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

plot_curve_full(CCA_XGB, "PowerLaw", 0.924, 1000)
Get_N(CCA_XGB, 0.924, 1000)

plot_curve_full(CCA_RF, "PowerLaw", 0.922, 1000)
Get_N(CCA_RF, 0.922, 1000)


#CANCER 4
Cancer <- read.csv("GSE183635_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE183635")


Meta <- data.frame("ID"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$geo_accession, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$geo_accession),
                   "Group"=c(Meta$`GSE183635-GPL16791_series_matrix.txt.gz`$`patient group:ch1`, Meta$`GSE183635-GPL20301_series_matrix.txt.gz`$`patient group:ch1`))


Cancer <- Cancer[,-1]
Cancer <- data.frame(t(Cancer))

Cancer$ID <- rownames(Cancer)

Cancer <- left_join(Cancer, Meta, "ID")
Cancer <- Cancer[,-39377]

Cancer <- Cancer[Cancer$Group %in% c("Glioma", "Epilepsy"),]
Cancer$Group <- ifelse(Cancer$Group=="Epilepsy", 0, 1)

Vars <- apply(Cancer[,-ncol(Cancer)], 2, var)
Cancer <- Cancer[,-which(colnames(Cancer) %in% names(Vars[Vars==0]))]
Glioma <- data.frame(apply(Cancer, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(Glioma),
                        Group = Glioma[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Glioma[,-which(colnames(Glioma) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Glioma[,1:(ncol(Glioma)-1)] <- DDS


set.seed(2024)
B <- Boruta(x=Cancer[,-which(colnames(Cancer) %in% "Group")], y=as.factor(Cancer[,which(colnames(Cancer) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Cancer_Keep <- Cancer[,c(Keep, "Group")]

Vars_True <- apply(Cancer_Keep[,-ncol(Cancer_Keep)], 2, var)

Cancer_Keep <- cbind(log2(Cancer_Keep[,-ncol(Cancer_Keep)]+0.01), "Group"=Cancer_Keep$Group)

LFC_True <- sapply(1:(ncol(Cancer_Keep)-1), function(i) {abs(mean(Cancer_Keep[Cancer_Keep$Group==1,i]) - mean(Cancer_Keep[Cancer_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Cancer, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Cancer[which(Cancer$Group==1),])), 
                                     p=ncol(Cancer)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
max(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "Glioma_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
Glioma_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

set.seed(2024)
Glioma_RF <- RF_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

plot_curve_full(Glioma_XGB, "PowerLaw", 0.964, 1000)
Get_N(Glioma_XGB, 0.964, 1000)

plot_curve_full(Glioma_RF, "PowerLaw", 0.963, 1000)
Get_N(Glioma_RF, 0.963, 1000)

#Lithium
Bipolar <- read.csv("GSE124326_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE124326")

Bipolar <- Bipolar[,-1]
Bipolar <- data.frame(t(Bipolar))

Meta <- data.frame("ID"=c(Meta$GSE124326_series_matrix.txt.gz$geo_accession),
                   "Group"=ifelse(Meta$GSE124326_series_matrix.txt.gz$`bipolar disorder diagnosis:ch1`=="Control", NA,
                                  ifelse(Meta$GSE124326_series_matrix.txt.gz$`lithium use (non-user=0, user = 1):ch1`==1, 1, 0)))

Bipolar$ID <- rownames(Bipolar)

Bipolar <- left_join(Bipolar, Meta, "ID")
Bipolar <- Bipolar[,-39377]

Bipolar <- na.omit(Bipolar)

Vars <- apply(Bipolar[,-ncol(Bipolar)], 2, var)
Bipolar <- Bipolar[,-which(colnames(Bipolar) %in% names(Vars[Vars==0]))]
Bipolar <- data.frame(apply(Bipolar, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(Bipolar),
                        Group = Bipolar[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Bipolar[,-which(colnames(Bipolar) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Bipolar[,1:(ncol(Bipolar)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Bipolar[,-which(colnames(Bipolar) %in% "Group")], y=as.factor(Bipolar[,which(colnames(Bipolar) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Bipolar_Keep <- Bipolar[,c(Keep, "Group")]

Vars_True <- apply(Bipolar_Keep[,-ncol(Bipolar_Keep)], 2, var)

Bipolar_Keep <- cbind(log2(Bipolar_Keep[,-ncol(Bipolar_Keep)]+0.01), "Group"=Bipolar_Keep$Group)

LFC_True <- sapply(1:(ncol(Bipolar_Keep)-1), function(i) {abs(mean(Bipolar_Keep[Bipolar_Keep$Group==1,i]) - mean(Bipolar_Keep[Bipolar_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Bipolar, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Bipolar[which(Bipolar$Group==1),])), 
                                     p=ncol(Bipolar)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
median(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "Bipolar_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
Bipolar_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
Bipolar_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(Bipolar_XGB, "PowerLaw", 0.933, 2000)
Get_N(Bipolar_XGB, 0.933, 2000)

plot_curve_full(Bipolar_RF, "PowerLaw", 0.930, 2000)
Get_N(Bipolar_RF, 0.930, 2000)


###
MDD <- read.csv("GSE80655_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE80655")

MDD <- MDD[,-1]
MDD <- data.frame(t(MDD))

Meta <- data.frame("ID"=c(Meta$GSE80655_series_matrix.txt.gz$geo_accession),
                   "Group"=ifelse(Meta$GSE80655_series_matrix.txt.gz$`clinical diagnosis:ch1` %in% c("Bipolar Disorder", "Schizophrenia"), NA,
                                  ifelse(Meta$GSE80655_series_matrix.txt.gz$`clinical diagnosis:ch1`=="Control", 0, 1)))

MDD$ID <- rownames(MDD)

MDD <- left_join(MDD, Meta, "ID")
MDD <- MDD[,-39377]

MDD <- na.omit(MDD)

Vars <- apply(MDD[,-ncol(MDD)], 2, var)
MDD <- MDD[,-which(colnames(MDD) %in% names(Vars[Vars==0]))]
MDD <- data.frame(apply(MDD, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(MDD),
                        Group = MDD[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(MDD[,-which(colnames(MDD) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

MDD[,1:(ncol(MDD)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=MDD[,-which(colnames(MDD) %in% "Group")], y=as.factor(MDD[,which(colnames(MDD) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

MDD_Keep <- MDD[,c(Keep, "Group")]

Vars_True <- apply(MDD_Keep[,-ncol(MDD_Keep)], 2, var)

MDD_Keep <- cbind(log2(MDD_Keep[,-ncol(MDD_Keep)]+0.01), "Group"=MDD_Keep$Group)

LFC_True <- sapply(1:(ncol(MDD_Keep)-1), function(i) {abs(mean(MDD_Keep[MDD_Keep$Group==1,i]) - mean(MDD_Keep[MDD_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(MDD, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(MDD[which(MDD$Group==1),])), 
                                     p=ncol(MDD)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
median(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "MDD_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2025)
MDD_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
MDD_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(MDD_XGB, "PowerLaw", 0.862, 2000)
Get_N_LOG(MDD_XGB, 0.862, 2000)

plot_curve_full(MDD_RF, "PowerLaw", 0.857, 2000)
Get_N(MDD_RF, 0.857, 2000)

###Crohn
Crohn <- read.csv("GSE137344_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE137344")

Crohn <- Crohn[,-1]
Crohn <- data.frame(t(Crohn))


Meta <- data.frame("ID"=c(Meta$GSE137344_series_matrix.txt.gz$geo_accession),
                   "Group"=ifelse(Meta$GSE137344_series_matrix.txt.gz$`diagnosis:ch1`=="Control",0,1))

Crohn$ID <- rownames(Crohn)

Crohn <- left_join(Crohn, Meta, "ID")
Crohn <- Crohn[,-39377]

Crohn <- na.omit(Crohn)

Vars <- apply(Crohn[,-ncol(Crohn)], 2, var)
Crohn <- Crohn[,-which(colnames(Crohn) %in% names(Vars[Vars==0]))]
Crohn <- data.frame(apply(Crohn, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(Crohn),
                        Group = Crohn[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Crohn[,-which(colnames(Crohn) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Crohn[,1:(ncol(Crohn)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Crohn[,-which(colnames(Crohn) %in% "Group")], y=as.factor(Crohn[,which(colnames(Crohn) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Crohn_Keep <- Crohn[,c(Keep, "Group")]

Vars_True <- apply(Crohn_Keep[,-ncol(Crohn_Keep)], 2, var)

Crohn_Keep <- cbind(log2(Crohn_Keep[,-ncol(Crohn_Keep)]+0.01), "Group"=Crohn_Keep$Group)

LFC_True <- sapply(1:(ncol(Crohn_Keep)-1), function(i) {abs(mean(Crohn_Keep[Crohn_Keep$Group==1,i]) - mean(Crohn_Keep[Crohn_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Crohn, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Crohn[which(Crohn$Group==1),])), 
                                     p=ncol(Crohn)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
max(LFC_Sim)

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "Crohn_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
Crohn_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
Crohn_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(Crohn_XGB, "PowerLaw", 0.903, 2000)
Get_N(Crohn_XGB, 0.903, 2000) #613

plot_curve_full(Crohn_RF, "PowerLaw", 0.896, 2000)
Get_N(Crohn_RF, 0.896, 2000) #191

1-mean(Crohn$Group)


#COVID
COVID <- read.csv("GSE184610_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE184610")

COVID  <- COVID [,-1]
COVID  <- data.frame(t(COVID))

Meta <- data.frame("ID"=c(Meta$GSE184610_series_matrix.txt.gz$geo_accession),
                   "Group"=ifelse(Meta$GSE184610_series_matrix.txt.gz$`disease state:ch1`=="Control",0,1))

COVID$ID <- rownames(COVID)

COVID <- left_join(COVID, Meta, "ID")
COVID <- COVID[,-39377]

COVID <- na.omit(COVID)

Vars <- apply(COVID[,-ncol(COVID)], 2, var)
COVID <- COVID[,-which(colnames(COVID) %in% names(Vars[Vars==0]))]
COVID <- data.frame(apply(COVID, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(COVID),
                        Group = COVID[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(COVID[,-which(colnames(COVID) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

COVID[,1:(ncol(COVID)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=COVID[,-which(colnames(COVID) %in% "Group")], y=as.factor(COVID[,which(colnames(COVID) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

COVID_Keep <- COVID[,c(Keep, "Group")]

Vars_True <- apply(COVID_Keep[,-ncol(COVID_Keep)], 2, var)

COVID_Keep <- cbind(log2(COVID_Keep[,-ncol(COVID_Keep)]+0.01), "Group"=COVID_Keep$Group)

LFC_True <- sapply(1:(ncol(COVID_Keep)-1), function(i) {abs(mean(COVID_Keep[COVID_Keep$Group==1,i]) - mean(COVID_Keep[COVID_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(COVID, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(COVID[which(COVID$Group==1),])), 
                                     p=ncol(COVID)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
median(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "COVID_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

set.seed(2024)
COVID_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
COVID_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(COVID_XGB, "PowerLaw", 0.936, 2000)
Get_N(COVID_XGB, 0.936, 2000) #613

plot_curve_full(COVID_RF, "PowerLaw", 0.936, 2000)
Get_N(COVID_RF, 0.936, 2000) #191

mean(COVID$Group)


#Therapy Response
ALS <- read.csv("GSE124439_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE124439")

ALS  <- ALS[,-1]
ALS  <- data.frame(t(ALS))

ALS$Group <- ifelse(Meta$GSE124439_series_matrix.txt.gz$`sample group:ch1` == "Non-Neurological Control", 0,
                    ifelse(Meta$GSE124439_series_matrix.txt.gz$`sample group:ch1`=="Other Neurological Disorders", NA, 1))

ALS <- na.omit(ALS)

Vars <- apply(ALS[,-ncol(ALS)], 2, var)
ALS <- ALS[,-which(colnames(ALS) %in% names(Vars[Vars==0]))]
ALS <- data.frame(apply(ALS, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(ALS),
                        Group = ALS[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(ALS[,-which(colnames(ALS) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

ALS[,1:(ncol(ALS)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=ALS[,-which(colnames(ALS) %in% "Group")], y=as.factor(ALS[,which(colnames(ALS) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

ALS_Keep <- ALS[,c(Keep, "Group")]

Vars_True <- apply(ALS_Keep[,-ncol(ALS_Keep)], 2, var)

ALS_Keep <- cbind(log2(ALS_Keep[,-ncol(ALS_Keep)]+0.01), "Group"=ALS_Keep$Group)

LFC_True <- sapply(1:(ncol(ALS_Keep)-1), function(i) {abs(mean(ALS_Keep[ALS_Keep$Group==1,i]) - mean(ALS_Keep[ALS_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(ALS, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(ALS[which(ALS$Group==1),])), 
                                     p=ncol(ALS)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
median(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "ALS_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

1-mean(ALS$Group)

set.seed(2)
ALS_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2)
ALS_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(ALS_XGB, "PowerLaw", 0.949, 2000)
Get_N(ALS_XGB, 0.949, 2000) #1048

plot_curve_full(ALS_RF, "PowerLaw", 0.943, 2000)
Get_N(ALS_RF, 0.943, 2000) #970

mean(ALS$Group)


##PrePost

PrePost <- read.csv("GSE229083_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE229083")

PrePost  <- PrePost[,-1]
PrePost  <- data.frame(t(PrePost))

Meta <- data.frame("ID"=c(Meta$GSE229083_series_matrix.txt.gz$geo_accession),
                   "Group"=ifelse(Meta$GSE229083_series_matrix.txt.gz$`classification group:ch1`=="PreCapillary",0,1))

PrePost$ID <- rownames(PrePost)

PrePost <- left_join(PrePost, Meta, "ID")
PrePost <- PrePost[,-39377]

PrePost <- na.omit(PrePost)

Vars <- apply(PrePost[,-ncol(PrePost)], 2, var)
PrePost <- PrePost[,-which(colnames(PrePost) %in% names(Vars[Vars==0]))]
PrePost <- data.frame(apply(PrePost, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(PrePost),
                        Group = PrePost[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(PrePost[,-which(colnames(PrePost) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

PrePost[,1:(ncol(PrePost)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=PrePost[,-which(colnames(PrePost) %in% "Group")], y=as.factor(PrePost[,which(colnames(PrePost) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

PrePost_Keep <- PrePost[,c(Keep, "Group")]

Vars_True <- apply(PrePost_Keep[,-ncol(PrePost_Keep)], 2, var)

PrePost_Keep <- cbind(log2(PrePost_Keep[,-ncol(PrePost_Keep)]+0.01), "Group"=PrePost_Keep$Group)

LFC_True <- sapply(1:(ncol(PrePost_Keep)-1), function(i) {abs(mean(PrePost_Keep[PrePost_Keep$Group==1,i]) - mean(PrePost_Keep[PrePost_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(PrePost, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(PrePost[which(PrePost$Group==1),])), 
                                     p=ncol(PrePost)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
median(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "PrePost_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

1-mean(PrePost$Group)

set.seed(2024)
PrePost_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 500, steps=10, True=Keep)

set.seed(2024)
PrePost_RF <- RF_Pipeline_Curve(Test, "Group", 50, 500, steps=10, True=Keep)

plot_curve_full(PrePost_XGB, "PowerLaw", 0.998, 500)
Get_N(PrePost_XGB, 0.998, 500) #154

plot_curve_full(PrePost_RF, "PowerLaw", 0.996, 500)
Get_N(PrePost_RF, 0.996, 500) #93

mean(PrePost$Group)


#HBV

HBV <- read.csv("GSE173897_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE173897")

HBV  <- HBV[,-1]
HBV  <- data.frame(t(HBV))

Meta$GSE173897_series_matrix.txt.gz$`viral load:ch1`

HBV$Group <- ifelse(Meta$GSE173897_series_matrix.txt.gz$`ethnicity:ch1`=="Hmong", 1, 0)

HBV <- na.omit(HBV)


Vars <- apply(HBV[,-ncol(HBV)], 2, var)
HBV <- HBV[,-which(colnames(HBV) %in% names(Vars[Vars==0]))]
HBV <- data.frame(apply(HBV, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(HBV),
                        Group = HBV[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(HBV[,-which(colnames(HBV) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

HBV[,1:(ncol(HBV)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=HBV[,-which(colnames(HBV) %in% "Group")], y=as.factor(HBV[,which(colnames(HBV) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

HBV_Keep <- HBV[,c(Keep, "Group")]

HBV_Keep <- cbind(log2(HBV_Keep[,-ncol(HBV_Keep)]+0.01), "Group"=HBV_Keep$Group)

LFC_True <- sapply(1:(ncol(HBV_Keep)-1), function(i) {abs(mean(HBV_Keep[HBV_Keep$Group==1,i]) - mean(HBV_Keep[HBV_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(HBV, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(HBV[which(HBV$Group==1),])), 
                                     p=ncol(HBV)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
max(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "HBV_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

1-mean(HBV$Group)

set.seed(2024)
HBV_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

set.seed(2024)
HBV_RF <- RF_Pipeline_Curve(Test, "Group", 50, 1000, steps=10, True=Keep)

plot_curve_full(HBV_XGB, "LOG", 0.883, 1000)
Get_N_LOG(HBV_XGB, 0.883, 1000) #529

plot_curve_full(HBV_RF, "LOG", 0.884, 1000)
Get_N_LOG(HBV_RF, 0.884, 1000) #106

mean(HBV$Group)

#NSCLC

NSCLC <- read.csv("GSE207586_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE207586")

NSCLC  <- NSCLC[,-1]
NSCLC  <- data.frame(t(NSCLC))

Meta <- data.frame("ID"=c(Meta$`GSE207586-GPL16791_series_matrix.txt.gz`$geo_accession, Meta$`GSE207586-GPL20301_series_matrix.txt.gz`$geo_accession),
                   "Group"=c(Meta$`GSE207586-GPL16791_series_matrix.txt.gz`$`patient group:ch1`, Meta$`GSE207586-GPL20301_series_matrix.txt.gz`$`patient group:ch1`))

Meta$Group <- ifelse(Meta$Group=="Control", 0, 1)

NSCLC$ID <- rownames(NSCLC)

NSCLC <- left_join(NSCLC, Meta, "ID")

NSCLC <- na.omit(NSCLC)

NSCLC <- NSCLC[,-c(39377)]



Vars <- apply(NSCLC[,-ncol(NSCLC)], 2, var)
NSCLC <- NSCLC[,-which(colnames(NSCLC) %in% names(Vars[Vars==0]))]
NSCLC <- data.frame(apply(NSCLC, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(NSCLC),
                        Group = NSCLC[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(NSCLC[,-which(colnames(NSCLC) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

NSCLC[,1:(ncol(NSCLC)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=NSCLC[,-which(colnames(NSCLC) %in% "Group")], y=as.factor(NSCLC[,which(colnames(NSCLC) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

NSCLC_Keep <- NSCLC[,c(Keep, "Group")]

Vars_True <- apply(NSCLC_Keep[,-ncol(NSCLC_Keep)], 2, var)

NSCLC_Keep <- cbind(log2(NSCLC_Keep[,-ncol(NSCLC_Keep)]+0.01), "Group"=NSCLC_Keep$Group)

LFC_True <- sapply(1:(ncol(NSCLC_Keep)-1), function(i) {abs(mean(NSCLC_Keep[NSCLC_Keep$Group==1,i]) - mean(NSCLC_Keep[NSCLC_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(NSCLC, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(NSCLC[which(NSCLC$Group==1),])), 
                                     p=ncol(NSCLC)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
max(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "NSCLC_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

1-mean(NSCLC$Group)

set.seed(2024)
NSCLC_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 500, steps=10, True=Keep)

set.seed(2024)
NSCLC_RF <- RF_Pipeline_Curve(Test, "Group", 50, 500, steps=10, True=Keep)

plot_curve_full(NSCLC_XGB, "PowerLaw", 0.995, 500)
Get_N(NSCLC_XGB, 0.995, 500) #235

plot_curve_full(NSCLC_RF, "PowerLaw", 0.993, 500)
Get_N(NSCLC_RF, 0.993, 500) #189

LR_Optimal(cbind(log2(Test[,-ncol(Test)]+0.01), "Group"=Test$Group), "Group")[1]


#GBM

GBM <- read.csv("GSE188812_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE188812")

GBM  <- GBM[,-1]
GBM  <- data.frame(t(GBM))

GRP <- NULL
for (i in 1:32) {
GRP <- append(GRP, str_split(Meta$GSE188812_series_matrix.txt.gz$`overall survival:ch1`, " ")[[i]][1])
}

GRP <- as.numeric(GRP)

GBM$Group <- ifelse(GRP>20, 1, 0)

Vars <- apply(GBM[,-ncol(GBM)], 2, var)
GBM <- GBM[,-which(colnames(GBM) %in% names(Vars[Vars==0]))]
GBM <- data.frame(apply(GBM, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(GBM),
                        Group = GBM[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(GBM[,-which(colnames(GBM) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

GBM[,1:(ncol(GBM)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=GBM[,-which(colnames(GBM) %in% "Group")], y=as.factor(GBM[,which(colnames(GBM) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

GBM_Keep <- GBM[,c(Keep, "Group")]

Vars_True <- apply(GBM_Keep[,-ncol(GBM_Keep)], 2, var)

GBM_Keep <- cbind(log2(GBM_Keep[,-ncol(GBM_Keep)]+0.01), "Group"=GBM_Keep$Group)

LFC_True <- sapply(1:(ncol(GBM_Keep)-1), function(i) {abs(mean(GBM_Keep[GBM_Keep$Group==1,i]) - mean(GBM_Keep[GBM_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(GBM, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(GBM[which(GBM$Group==1),])), 
                                     p=ncol(GBM)-1, n=5000)

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True)
max(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "GBM_Sim.csv")

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

mean(GBM$Group)

set.seed(2024)
GBM_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

set.seed(2024)
GBM_RF <- RF_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

plot_curve_full(GBM_XGB, "PowerLaw", 1, 250)
Get_N(GBM_XGB, 1, 500) #54

plot_curve_full(GBM_RF, "PowerLaw", 1, 250)
Get_N(GBM_RF, 1, 500) #50

LR_Optimal(cbind(log2(Test[,-ncol(Test)]+0.01), "Group"=Test$Group), "Group")[1]


#

TB <- read.csv("GSE103147_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE103147")

TB  <- TB[,-1]
TB  <- data.frame(t(TB))


Meta <- data.frame("ID"=Meta$GSE103147_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE103147_series_matrix.txt.gz$`qft:ch1`=="positive", 1, 
                                  ifelse(Meta$GSE103147_series_matrix.txt.gz$`qft:ch1`=="NA", NA, 0)),
                   "Time"=Meta$GSE103147_series_matrix.txt.gz$`timepoint:ch1`)

Meta <- Meta[Meta$Time=="0",]

Meta <- Meta[,-3]

TB$ID <- rownames(TB)

TB <- left_join(TB, Meta, "ID")

TB <- na.omit(TB)

TB <- TB[,-c(39377)]

Vars <- apply(TB[,-ncol(TB)], 2, var)
TB <- TB[,-which(colnames(TB) %in% names(Vars[Vars==0]))]
TB <- data.frame(apply(TB, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(TB),
                        Group = TB[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(TB[,-which(colnames(TB) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

TB[,1:(ncol(TB)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=TB[,-which(colnames(TB) %in% "Group")], y=as.factor(TB[,which(colnames(TB) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])


TB_Keep <- TB[,c(Keep, "Group")]

Vars_True <- apply(TB_Keep[,-ncol(TB_Keep)], 2, var)

TB_Keep <- cbind(log2(TB_Keep[,-ncol(TB_Keep)]+0.01), "Group"=TB_Keep$Group)

LFC_True <- sapply(1:(ncol(TB_Keep)-1), function(i) {abs(mean(TB_Keep[TB_Keep$Group==1,i]) - mean(TB_Keep[TB_Keep$Group==0,i]))})

LFC_True

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(TB, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(TB[which(TB$Group==1),])), 
                                     p=ncol(TB)-1, n=5000)



LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True) #0.994
min(LFC_Sim)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

set.seed(2024)
XGBoost_Optimal(Test, "Group") #0.893

set.seed(2024)
RF_Optimal(Test, "Group") #0.895

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "TB_Sim.csv")

set.seed(2025)
TB_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
TB_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(TB_XGB, "PowerLaw", 0.893, 2000)
Get_N(TB_XGB, 0.893, 2000) #904

plot_curve_full(TB_RF, "PowerLaw", 0.895, 2000)
Get_N(TB_RF, 0.895, 2000) #1005

LR_Optimal(cbind(log2(Test[,-ncol(Test)]+0.01), "Group"=Test$Group), "Group")[1] #0.841

1-mean(TB$Group) #17.5

#Ovarian
Ovarian <- read.csv("GSE230522_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE230522")



Ovarian  <- Ovarian[,-1]
Ovarian  <- data.frame(t(Ovarian))


Meta <- data.frame("ID"=Meta$GSE230522_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE230522_series_matrix.txt.gz$`asp_twogroup_recoded:ch1`=="current", 1, 
                                  ifelse(Meta$GSE230522_series_matrix.txt.gz$`asp_twogroup_recoded:ch1`=="NA", NA, 0)))


Ovarian$ID <- rownames(Ovarian)

Ovarian <- left_join(Ovarian, Meta, "ID")

Ovarian <- na.omit(Ovarian)

Ovarian <- Ovarian[,-c(39377)]

Vars <- apply(Ovarian[,-ncol(Ovarian)], 2, var)
Ovarian <- Ovarian[,-which(colnames(Ovarian) %in% names(Vars[Vars==0]))]
Ovarian <- data.frame(apply(Ovarian, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(Ovarian),
                        Group = Ovarian[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Ovarian[,-which(colnames(Ovarian) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Ovarian[,1:(ncol(Ovarian)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=Ovarian[,-which(colnames(Ovarian) %in% "Group")], y=as.factor(Ovarian[,which(colnames(Ovarian) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])


Ovarian_Keep <- Ovarian[,c(Keep, "Group")]

Vars_True <- apply(Ovarian_Keep[,-ncol(Ovarian_Keep)], 2, var)

Ovarian_Keep <- cbind(log2(Ovarian_Keep[,-ncol(Ovarian_Keep)]+0.01), "Group"=Ovarian_Keep$Group)

LFC_True <- sapply(1:(ncol(Ovarian_Keep)-1), function(i) {abs(mean(Ovarian_Keep[Ovarian_Keep$Group==1,i]) - mean(Ovarian_Keep[Ovarian_Keep$Group==0,i]))})

LFC_True

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(Ovarian, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(Ovarian[which(Ovarian$Group==1),])), 
                                     p=ncol(Ovarian)-1, n=5000)



LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True) #>0.999
median(LFC_Sim) #0.60 (0.34, 2.24)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

set.seed(2024)
XGBoost_Optimal(Test, "Group") #0.893

set.seed(2024)
RF_Optimal(Test, "Group") #0.895

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "Ovarian_Sim.csv")

set.seed(2025)
Ovarian_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

set.seed(2024)
Ovarian_RF <- RF_Pipeline_Curve(Test, "Group", 50, 2000, steps=10, True=Keep)

plot_curve_full(Ovarian_XGB, "PowerLaw", 0.845, 2000)
Get_N(Ovarian_XGB, 0.845, 2000) #713

plot_curve_full(Ovarian_RF, "PowerLaw", 0.846, 2000)
Get_N(Ovarian_RF, 0.846, 2000) #215

LR_Optimal(cbind(log2(Test[,-ncol(Test)]+0.01), "Group"=Test$Group), "Group")[1] #0.831

mean(Ovarian$Group) #39.9%


#EDS
EDS <- read.csv("GSE218012_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE218012")


EDS  <- EDS[,-1]
EDS  <- data.frame(t(EDS))


Meta <- data.frame("ID"=Meta$GSE218012_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE218012_series_matrix.txt.gz$`disease state:ch1`=="affected", 1, 0)
)


EDS$ID <- rownames(EDS)

EDS <- left_join(EDS, Meta, "ID")

EDS <- na.omit(EDS)

EDS <- EDS[,-c(39377)]

Vars <- apply(EDS[,-ncol(EDS)], 2, var)
EDS <- EDS[,-which(colnames(EDS) %in% names(Vars[Vars==0]))]
EDS <- data.frame(apply(EDS, 2, as.numeric))

ExpDesign <- data.frame(row.names=rownames(EDS),
                        Group = EDS[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(EDS[,-which(colnames(EDS) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

EDS[,1:(ncol(EDS)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=EDS[,-which(colnames(EDS) %in% "Group")], y=as.factor(EDS[,which(colnames(EDS) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])


EDS_Keep <- EDS[,c(Keep, "Group")]

Vars_True <- apply(EDS_Keep[,-ncol(EDS_Keep)], 2, var)

EDS_Keep <- cbind(log2(EDS_Keep[,-ncol(EDS_Keep)]+0.01), "Group"=EDS_Keep$Group)

LFC_True <- sapply(1:(ncol(EDS_Keep)-1), function(i) {abs(mean(EDS_Keep[EDS_Keep$Group==1,i]) - mean(EDS_Keep[EDS_Keep$Group==0,i]))})

LFC_True

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(EDS, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(EDS[which(EDS$Group==1),])), 
                                     p=ncol(EDS)-1, n=5000)



LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
cor(LFC_Sim, LFC_True) #>0.999
max(LFC_Sim) #0.86 (0.20, 5.16)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

set.seed(2024)
XGBoost_Optimal(Test, "Group") #1

set.seed(2024)
RF_Optimal(Test, "Group") #1

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
#write.csv(Test, "EDS_Sim.csv")

set.seed(2024)
EDS_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

set.seed(2024)
EDS_RF <- RF_Pipeline_Curve(Test, "Group", 50, 250, steps=10, True=Keep)

plot_curve_full(EDS_XGB, "PowerLaw", 1, 250)
Get_N(EDS_XGB, 1, 250) #132

plot_curve_full(EDS_RF, "PowerLaw", 1, 250)
Get_N(EDS_RF, 1, 250) #41

LR_Optimal(cbind(log2(Test[,-ncol(Test)]+0.01), "Group"=Test$Group), "Group")[1] #0.998

1-mean(EDS$Group) #50%

#TCGA----------------------------------------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(here)
setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")

LGG_Query <- GDCquery(project= "TCGA-LGG",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

GDCdownload(LGG_Query)

TCGA_LGG_Data <- GDCprepare(LGG_Query)

saveRDS(TCGA_LGG_Data, "TCGA_LGG_Data_SummarizedExperiment.rds")

TCGA_LGG_Data <- readRDS("TCGA_LGG_Data_SummarizedExperiment.rds")

TCGA_LGG_Mat <- assay(TCGA_LGG_Data)

#Number 2
GBM_Query <- GDCquery(project= "TCGA-GBM",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

GDCdownload(GBM_Query)

TCGA_GBM_Data <- GDCprepare(GBM_Query)

saveRDS(TCGA_GBM_Data, "TCGA_GBM_Data_SummarizedExperiment.rds")

TCGA_GBM_Data <- readRDS("TCGA_GBM_Data_SummarizedExperiment.rds")

TCGA_GBM_Mat <- assay(TCGA_GBM_Data)

library(org.Hs.eg.db)

TCGA_GBM_Genes <- rownames(TCGA_GBM_Mat) %>% tibble::enframe() %>% mutate(ENSEMBL= stringr::str_replace(value, "\\.[0-9]+", ""))

TCGA_GBM_gene_map <- clusterProfiler::bitr(TCGA_GBM_Genes$ENSEMBL,
                                            fromType = "ENSEMBL",
                                            toType = "SYMBOL",
                                            OrgDb = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

TCGA_GBM_gene_map <- TCGA_GBM_gene_map %>% left_join(TCGA_GBM_Genes)

head(TCGA_GBM_gene_map)

TCGA_LGG_Mat <- TCGA_LGG_Mat[TCGA_GBM_gene_map$value,]
row.names(TCGA_LGG_Mat) <- TCGA_GBM_gene_map$SYMBOL

TCGA_GBM_Mat <- TCGA_GBM_Mat[TCGA_GBM_gene_map$value,]
row.names(TCGA_GBM_Mat) <- TCGA_GBM_gene_map$SYMBOL

dim(TCGA_GBM_Mat)
dim(TCGA_LGG_Mat)

#Combined LUSC and LUAD
TCGA_Brain_Mat <- cbind(TCGA_LGG_Mat, TCGA_GBM_Mat)

TCGA_Brain_Meta <- data.frame("Group"=c(rep("LGG", ncol(TCGA_LGG_Mat)), 
                                       rep("GBM", ncol(TCGA_GBM_Mat))))

TCGA_BRAIN <- data.frame(t(TCGA_Brain_Mat))

TCGA_BRAIN$Group <- TCGA_Brain_Meta$Group

TCGA_BRAIN$Group <- ifelse(TCGA_BRAIN$Group=="GBM", 1, 0)

TCGA_BRAIN <- apply(TCGA_BRAIN, 2, as.numeric)
TCGA_BRAIN <- data.frame(TCGA_BRAIN)

dim(TCGA_BRAIN)


Vars <- apply(TCGA_BRAIN[,-ncol(TCGA_BRAIN)], 2, var)
TCGA_BRAIN <- TCGA_BRAIN[,-which(colnames(TCGA_BRAIN) %in% names(Vars[Vars==0]))]
TCGA_BRAIN <- data.frame(apply(TCGA_BRAIN, 2, as.numeric))

set.seed(2025)
Test <- TCGA_BRAIN[sample(1:nrow(TCGA_BRAIN), 125), ]
TCGA_BRAIN <- TCGA_BRAIN[-which(rownames(TCGA_BRAIN) %in% rownames(Test)),]

median(GetDispersion(TCGA_BRAIN, T, Keep=colnames(TCGA_BRAIN_Keep)[-ncol(TCGA_BRAIN_Keep)])) #1.49

#Predictors
ExpDesign <- data.frame(row.names=rownames(TCGA_BRAIN),
                        Group = TCGA_BRAIN[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(TCGA_BRAIN[,-which(colnames(TCGA_BRAIN) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

TCGA_BRAIN[,1:(ncol(TCGA_BRAIN)-1)] <- TCGA_BRAIN

#Test Transformation
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test[,-which(colnames(Test) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Test[,1:(ncol(Test)-1)] <- DDS


#Find Predictors in IBD Set + Simulate BNG Set
set.seed(2024)
B <- Boruta(x=TCGA_BRAIN[,-which(colnames(TCGA_BRAIN) %in% "Group")], y=as.factor(TCGA_BRAIN[,which(colnames(TCGA_BRAIN) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Keep


#Analysis
TCGA_BRAIN_Keep <- TCGA_BRAIN[,c(Keep, "Group")]

TCGA_BRAIN_Keep <- cbind(log2(TCGA_BRAIN_Keep[,-ncol(TCGA_BRAIN_Keep)]+0.01), "Group"=TCGA_BRAIN_Keep$Group)

LFC_True <- sapply(1:(ncol(TCGA_BRAIN_Keep)-1), function(i) {abs(mean(TCGA_BRAIN_Keep[TCGA_BRAIN_Keep$Group==1,i]) - mean(TCGA_BRAIN_Keep[TCGA_BRAIN_Keep$Group==0,i]))})

#TCGA_BRAIN_Keep <- TCGA_BRAIN_Keep[,c(which(LFC_True<=1.25), 103)]

max(LFC_True) #8.62
median(LFC_True) #3.39
min(LFC_True) #0.43
mean(TCGA_BRAIN_Keep$Group) #41%

colnames(TCGA_BRAIN_Keep)

LFC_True

Test_Keep <- Test[,c(colnames(TCGA_BRAIN_Keep)[-ncol(TCGA_BRAIN_Keep)], "Group")]

Test_Keep <- cbind(log2(Test_Keep[,-ncol(Test_Keep)]+0.01), "Group"=Test_Keep$Group)

XGBoost_Optimal(TCGA_BRAIN_Keep, "Group", Test_Keep) #0.973

set.seed(2025)
RF_Optimal_2(TCGA_BRAIN_Keep, "Group", Test_Keep) #0.978

set.seed(2025)
XGB_True_TCGA_BRAIN <- XGBoost_Pipeline_Curve_2(TCGA_BRAIN_Keep, 'Group', start=25, size=250, steps=10, True=colnames(TCGA_BRAIN_Keep)[-ncol(TCGA_BRAIN_Keep)], Test=Test_Keep)

set.seed(2025)
RF_True_TCGA_BRAIN <- RF_Pipeline_Curve_2(TCGA_BRAIN_Keep, 'Group', start=25, size=250, steps=10, True=colnames(TCGA_BRAIN_Keep)[-ncol(TCGA_BRAIN_Keep)], Test=Test_Keep)

LR_Optimal_2(TCGA_BRAIN_Keep, "Group", Test_Keep) 

plot_curve_full(XGB_True_TCGA_BRAIN, "PowerLaw", 1, 250)
Get_N(XGB_True_TCGA_BRAIN, 1, 250) #68
AUC_N(XGB_True_TCGA_BRAIN, 25) #94.2

plot_curve_full(RF_True_TCGA_BRAIN, "PowerLaw_Fixed", 1, 250)
Get_N_Fixed(RF_True_TCGA_BRAIN, 1, 250) #25
AUC_N_Fixed(RF_True_TCGA_BRAIN, 1, 25) #99.4

NN_Optimal_2(TCGA_BRAIN_Keep, "Group", Test_Keep, h=100) #1

#write.csv(TCGA_BRAIN_Keep, "TCGA_Brain.csv")
#write.csv(Test_Keep, "TCGA_Brain_Test.csv")

#write.csv(XGB_True_TCGA_BRAIN, "TCGA_Brain_XGB.csv")
#write.csv(RF_True_TCGA_BRAIN, "TCGA_Brain_RF.csv")

TCGA_BRAIN_NN <- NULL
for (i in round(seq(25, 500, by=475/10),0)) {
  TCGA_BRAIN_NN <- rbind(TCGA_BRAIN_NN, NN_Curve_Evaluate_2(TCGA_BRAIN_Keep, "Group", n=i, True=colnames(TCGA_BRAIN_Keep)[-ncol(TCGA_BRAIN_Keep)], Test=Test_Keep, h=100))
}
TCGA_BRAIN_NN

write.csv(TCGA_BRAIN_NN, "TCGA_BRAIN_NN.csv")

plot_curve_full(TCGA_BRAIN_NN, "PowerLaw", 1, 500)
Get_N(TCGA_BRAIN_NN, 1, 500) #25
AUC_N(TCGA_BRAIN_NN, 25) #99.0

median(AvgReadCount(cbind(2^(TCGA_LUNG_Keep[,-ncol(TCGA_LUNG_Keep)]-0.01), "Group"=TCGA_LUNG_Keep$Group))) #1128.37
median(AvgReadCount(cbind(2^(TCGA_BRAIN_Keep[,-ncol(TCGA_BRAIN_Keep)]-0.01), "Group"=TCGA_BRAIN_Keep$Group))) #35.49
median(AvgReadCount(cbind(2^(TCGA_BRCA_Keep[,-ncol(TCGA_BRCA_Keep)]-0.01), "Group"=TCGA_BRCA_Keep$Group))) #615.84

median(GetCorrelation(cbind(2^(TCGA_LUNG_Keep[,-ncol(TCGA_LUNG_Keep)]-0.01), "Group"=TCGA_LUNG_Keep$Group), Keep=colnames(TCGA_LUNG_Keep)[-ncol(TCGA_LUNG_Keep)])) #0.56
median(GetCorrelation(cbind(2^(TCGA_BRAIN_Keep[,-ncol(TCGA_BRAIN_Keep)]-0.01), "Group"=TCGA_BRAIN_Keep$Group), Keep=colnames(TCGA_BRAIN_Keep)[-ncol(TCGA_BRAIN_Keep)])) #0.60
median(GetCorrelation(cbind(2^(TCGA_BRCA_Keep[,-ncol(TCGA_BRCA_Keep)]-0.01), "Group"=TCGA_BRCA_Keep$Group), Keep=colnames(TCGA_BRCA_Keep)[-ncol(TCGA_BRCA_Keep)])) #0.25

#LUNG
max(LFC_True) #11.19
median(LFC_True) #4.42
min(LFC_True) #0.21




#LUNG ------------------------------------------------------------------------------------------------------------------------------
LUAD_Query <- GDCquery(project= "TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

GDCdownload(LUAD_Query)

TCGA_LUAD_Data <- GDCprepare(LUAD_Query)

saveRDS(TCGA_LUAD_Data, "TCGA_LUAD_Data_SummarizedExperiment.rds")

TCGA_LUAD_Data <- readRDS("TCGA_LUAD_Data_SummarizedExperiment.rds")

TCGA_LUAD_Mat <- assay(TCGA_LUAD_Data)

#Number 2
LUSC_Query <- GDCquery(project= "TCGA-LUSC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

GDCdownload(LUSC_Query)

TCGA_LUSC_Data <- GDCprepare(LUSC_Query)

saveRDS(TCGA_LUSC_Data, "TCGA_LUSC_Data_SummarizedExperiment.rds")

TCGA_LUSC_Data <- readRDS("TCGA_LUSC_Data_SummarizedExperiment.rds")

TCGA_LUSC_Mat <- assay(TCGA_LUSC_Data)

library(org.Hs.eg.db)

TCGA_LUAD_Genes <- rownames(TCGA_LUAD_Mat) %>% tibble::enframe() %>% mutate(ENSEMBL= stringr::str_replace(value, "\\.[0-9]+", ""))

TCGA_LUAD_gene_map <- clusterProfiler::bitr(TCGA_LUAD_Genes$ENSEMBL,
                                           fromType = "ENSEMBL",
                                           toType = "SYMBOL",
                                           OrgDb = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

TCGA_LUAD_gene_map <- TCGA_LUAD_gene_map %>% left_join(TCGA_LUAD_Genes)

head(TCGA_LUAD_gene_map)

TCGA_LUSC_Mat <- TCGA_LUSC_Mat[TCGA_LUAD_gene_map$value,]
row.names(TCGA_LUSC_Mat) <- TCGA_LUAD_gene_map$SYMBOL

TCGA_LUAD_Mat <- TCGA_LUAD_Mat[TCGA_LUAD_gene_map$value,]
row.names(TCGA_LUAD_Mat) <- TCGA_LUAD_gene_map$SYMBOL

dim(TCGA_LUSC_Mat)
dim(TCGA_LUAD_Mat)

#Combined LUSC and LUAD
TCGA_LUNG_Mat <- cbind(TCGA_LUAD_Mat, TCGA_LUSC_Mat)

TCGA_LUNG_Meta <- data.frame("Group"=c(rep("LUAD", ncol(TCGA_LUAD_Mat)), 
                                        rep("LUSC", ncol(TCGA_LUSC_Mat))))

TCGA_LUNG <- data.frame(t(TCGA_LUNG_Mat))

TCGA_LUNG$Group <- TCGA_LUNG_Meta$Group

TCGA_LUNG$Group <- ifelse(TCGA_LUNG$Group=="LUSC", 1, 0)

TCGA_LUNG <- apply(TCGA_LUNG, 2, as.numeric)
TCGA_LUNG <- data.frame(TCGA_LUNG)

dim(TCGA_LUNG)
mean(TCGA_LUNG$Group)

Vars <- apply(TCGA_LUNG[,-ncol(TCGA_LUNG)], 2, var)
TCGA_LUNG <- TCGA_LUNG[,-which(colnames(TCGA_LUNG) %in% names(Vars[Vars==0]))]
TCGA_LUNG <- data.frame(apply(TCGA_LUNG, 2, as.numeric))

median(GetDispersion(TCGA_LUNG, T, Keep=colnames(TCGA_LUNG_Keep)[-ncol(TCGA_LUNG_Keep)])) #2.86

set.seed(2025)
Test <- TCGA_LUNG[sample(1:nrow(TCGA_LUNG), 162), ]
TCGA_LUNG <- TCGA_LUNG[-which(rownames(TCGA_LUNG) %in% rownames(Test)),]


#Predictors
ExpDesign <- data.frame(row.names=rownames(TCGA_LUNG),
                        Group = TCGA_LUNG[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(TCGA_LUNG[,-which(colnames(TCGA_LUNG) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

TCGA_LUNG[,1:(ncol(TCGA_LUNG)-1)] <- TCGA_LUNG

#Test Transformation
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test[,-which(colnames(Test) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Test[,1:(ncol(Test)-1)] <- DDS


#Find Predictors in IBD Set + Simulate BNG Set
set.seed(2024)
B <- Boruta(x=TCGA_LUNG[,-which(colnames(TCGA_LUNG) %in% "Group")], y=as.factor(TCGA_LUNG[,which(colnames(TCGA_LUNG) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Keep


#Analysis
TCGA_LUNG_Keep <- TCGA_LUNG[,c(Keep, "Group")]

TCGA_LUNG_Keep <- cbind(log2(TCGA_LUNG_Keep[,-ncol(TCGA_LUNG_Keep)]+0.01), "Group"=TCGA_LUNG_Keep$Group)

LFC_True <- sapply(1:(ncol(TCGA_LUNG_Keep)-1), function(i) {abs(mean(TCGA_LUNG_Keep[TCGA_LUNG_Keep$Group==1,i]) - mean(TCGA_LUNG_Keep[TCGA_LUNG_Keep$Group==0,i]))})

#TCGA_LUNG_Keep <- TCGA_LUNG_Keep[,c(which(LFC_True<=1.5), 104)]

max(LFC_True) #11.19
median(LFC_True) #4.42
min(LFC_True) #0.21
mean(TCGA_LUNG_Keep$Group) #0.49

colnames(TCGA_LUNG_Keep)

#write.csv(TCGA_LUNG_Keep, "TCGA_Lung.csv")

LFC_True

Test_Keep <- Test[,c(colnames(TCGA_LUNG_Keep)[-ncol(TCGA_LUNG_Keep)], "Group")]

Test_Keep <- cbind(log2(Test_Keep[,-ncol(Test_Keep)]+0.01), "Group"=Test_Keep$Group)

write.csv(Test_Keep, "TCGA_Lung_Test.csv")

XGBoost_Optimal(TCGA_LUNG_Keep, "Group", Test_Keep) #0.985

set.seed(2025)
RF_Optimal_2(TCGA_LUNG_Keep, "Group", Test_Keep) #0.983

set.seed(2025)
XGB_True_TCGA_LUNG <- XGBoost_Pipeline_Curve_2(TCGA_LUNG_Keep, 'Group', start=25, size=300, steps=10, True=colnames(TCGA_LUNG_Keep)[-ncol(TCGA_LUNG_Keep)], Test=Test_Keep)

set.seed(2025)
RF_True_TCGA_LUNG <- RF_Pipeline_Curve_2(TCGA_LUNG_Keep, 'Group', start=25, size=300, steps=10, True=colnames(TCGA_LUNG_Keep)[-ncol(TCGA_LUNG_Keep)], Test=Test_Keep)

LR_Optimal_2(TCGA_LUNG_Keep, "Group", Test_Keep) #0.949

plot_curve_full(XGB_True_TCGA_LUNG, "PowerLaw", 0.985, 300)
Get_N(XGB_True_TCGA_LUNG, 0.985, 300) #86
AUC_N(XGB_True_TCGA_LUNG, 25) #89.0

plot_curve_full(RF_True_TCGA_LUNG, "PowerLaw_Fixed", 0.983, 300)
Get_N_Fixed(RF_True_TCGA_LUNG, 0.983, 700) #25
AUC_N_Fixed(RF_True_TCGA_LUNG, 0.983, 25) #96.5


#write.csv(XGB_True_TCGA_LUNG, "TCGA_Lung_XGB.csv")
#write.csv(RF_True_TCGA_LUNG, "TCGA_Lung_RF.csv")


NN_Optimal_2(TCGA_LUNG_Keep, "Group", Test_Keep, h=100) #0.981


TCGA_LUNG_NN <- NULL
for (i in round(seq(25, 500, by=475/10),0)) {
  TCGA_LUNG_NN <- rbind(TCGA_LUNG_NN, NN_Curve_Evaluate_2(TCGA_LUNG_Keep, "Group", n=i, True=colnames(TCGA_LUNG_Keep)[-ncol(TCGA_LUNG_Keep)], Test=Test_Keep, h=100))
}
TCGA_LUNG_NN

plot_curve_full(TCGA_LUNG_NN, "PowerLaw", 0.981, 500)
Get_N(TCGA_LUNG_NN, 0.981, 800) #69

write.csv(TCGA_LUNG_NN, "TCGA_LUNG_NN.csv")

AUC_N(TCGA_LUNG_NN, 25) #94.1



#TCGA BRCA--------------------------------------------------------------------------------------------------------------------------
BRCA_Query <- GDCquery(project= "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

GDCdownload(BRCA_Query)

TCGA_BRCA_Data <- GDCprepare(BRCA_Query)

saveRDS(TCGA_BRCA_Data, "TCGA_BRCA_Data_SummarizedExperiment.rds")

TCGA_BRCA_Data <- readRDS("TCGA_BRCA_Data_SummarizedExperiment.rds")

TCGA_BRCA_Mat <- assay(TCGA_BRCA_Data)


library(org.Hs.eg.db)

TCGA_BRCA_Genes <- rownames(TCGA_BRCA_Mat) %>% tibble::enframe() %>% mutate(ENSEMBL= stringr::str_replace(value, "\\.[0-9]+", ""))

TCGA_BRCA_gene_map <- clusterProfiler::bitr(TCGA_BRCA_Genes$ENSEMBL,
                                            fromType = "ENSEMBL",
                                            toType = "SYMBOL",
                                            OrgDb = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

TCGA_BRCA_gene_map <- TCGA_BRCA_gene_map %>% left_join(TCGA_BRCA_Genes)

head(TCGA_BRCA_gene_map)

TCGA_BRCA_Mat <- TCGA_BRCA_Mat[TCGA_BRCA_gene_map$value,]
row.names(TCGA_BRCA_Mat) <- TCGA_BRCA_gene_map$SYMBOL


#Get meta?

TCGA_BRCA <- data.frame(t(TCGA_BRCA_Mat))

TCGA_BRCA$Group <- TCGA_BRCA_Data$paper_BRCA_Pathology
#table(TCGA_BRCA_Data$`paper_CNV Clusters`)

TCGA_BRCA <- na.omit(TCGA_BRCA)

TCGA_BRCA$Group <- ifelse(TCGA_BRCA$Group %in% c("IDC"), 1, ifelse(TCGA_BRCA$Group=="NA", NA, 0))

TCGA_BRCA <- na.omit(TCGA_BRCA)
mean(TCGA_BRCA$Group)

TCGA_BRCA <- apply(TCGA_BRCA, 2, as.numeric)
TCGA_BRCA <- data.frame(TCGA_BRCA)

dim(TCGA_BRCA)


Vars <- apply(TCGA_BRCA[,-ncol(TCGA_BRCA)], 2, var)
TCGA_BRCA <- TCGA_BRCA[,-which(colnames(TCGA_BRCA) %in% names(Vars[Vars==0]))]
TCGA_BRCA <- data.frame(apply(TCGA_BRCA, 2, as.numeric))

median(GetDispersion(TCGA_BRCA, T, Keep=colnames(TCGA_BRCA_Keep)[-ncol(TCGA_BRCA_Keep)])) #0.84

set.seed(2025)
Test <- TCGA_BRCA[sample(1:nrow(TCGA_BRCA), 50), ]
TCGA_BRCA <- TCGA_BRCA[-which(rownames(TCGA_BRCA) %in% rownames(Test)),]


#Predictors
ExpDesign <- data.frame(row.names=rownames(TCGA_BRCA),
                        Group = TCGA_BRCA[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(TCGA_BRCA[,-which(colnames(TCGA_BRCA) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

TCGA_BRCA[,1:(ncol(TCGA_BRCA)-1)] <- TCGA_BRCA

#Test Transformation
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test[,-which(colnames(Test) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Test[,1:(ncol(Test)-1)] <- DDS


#Find Predictors in IBD Set + Simulate BNG Set
set.seed(2024)
B <- Boruta(x=TCGA_BRCA[,-which(colnames(TCGA_BRCA) %in% "Group")], y=as.factor(TCGA_BRCA[,which(colnames(TCGA_BRCA) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Keep


#Analysis
TCGA_BRCA_Keep <- TCGA_BRCA[,c(Keep, "Group")]

TCGA_BRCA_Keep <- cbind(log2(TCGA_BRCA_Keep[,-ncol(TCGA_BRCA_Keep)]+0.01), "Group"=TCGA_BRCA_Keep$Group)

LFC_True <- sapply(1:(ncol(TCGA_BRCA_Keep)-1), function(i) {abs(mean(TCGA_BRCA_Keep[TCGA_BRCA_Keep$Group==1,i]) - mean(TCGA_BRCA_Keep[TCGA_BRCA_Keep$Group==0,i]))})
max(LFC_True) #1.90
median(LFC_True) #0.83
min(LFC_True) #0.18
1-mean(TCGA_BRCA$Group) #39.6%

colnames(TCGA_BRCA_Keep)

LFC_True

Test_Keep <- Test[,c(Keep, "Group")]
Test_Keep <- Test[,c(colnames(TCGA_BRCA_Keep)[-ncol(TCGA_BRCA_Keep)], "Group")]

Test_Keep <- cbind(log2(Test_Keep[,-ncol(Test_Keep)]+0.01), "Group"=Test_Keep$Group)

XGBoost_Optimal(TCGA_BRCA_Keep, "Group", Test_Keep) #0.718

#write.csv(TCGA_BRCA_Keep, "TCGA_BRCA.csv")
#write.csv(Test_Keep, "TCGA_BRCA_Test.csv")

set.seed(2025)
RF_Optimal_2(TCGA_BRCA_Keep, "Group", Test_Keep) #0.75

set.seed(2025)
XGB_True_TCGA_BRCA <- XGBoost_Pipeline_Curve_2(TCGA_BRCA_Keep, 'Group', start=25, size=800, steps=10, True=Keep, Test=Test_Keep)

set.seed(2025)
RF_True_TCGA_BRCA <- RF_Pipeline_Curve_2(TCGA_BRCA_Keep, 'Group', start=25, size=800, steps=10, True=Keep, Test=Test_Keep)

LR_Optimal_2(TCGA_BRCA_Keep, "Group", Test_Keep) #0.711

library(h2o)
h2o.init()
NN_Optimal_2(TCGA_BRCA_Keep, "Group", Test_Keep, h=200) #0.724
mean(c(0.724, 0.716, 0.721, 0.726, 0.769, 0.744, 0.714, 0.716, 0.698, 0.711))

TCGA_BRCA_NN <- NULL
for (i in round(seq(25, 800, by=775/10),0)) {
  TCGA_BRCA_NN <- rbind(TCGA_BRCA_NN, NN_Curve_Evaluate_2(TCGA_BRCA_Keep, "Group", n=i, True=colnames(TCGA_BRCA_Keep)[-ncol(TCGA_BRCA_Keep)], Test=Test_Keep, h=100))
}
TCGA_BRCA_NN

plot_curve_full(TCGA_BRCA_NN, "LOG", 0.724, 800)
Get_N_LOG(TCGA_BRCA_NN, 0.724, 800) #298
AUC_N_LOG(TCGA_BRCA_NN, 25) #).673

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
#write.csv(TCGA_BRCA_NN, "TCGA_BRCA_NN.csv")


plot_curve_full(XGB_True_TCGA_BRCA, "PowerLaw_Fixed", 0.718, 800)
Get_N_Fixed(XGB_True_TCGA_BRCA, 0.718, 800) #607
AUC_N_Fixed(XGB_True_TCGA_BRCA, 0.718, 25) #0.642

plot_curve_full(RF_True_TCGA_BRCA, "LOG", 0.75, 800)
Get_N_LOG(RF_True_TCGA_BRCA, 0.75, 800) #419
AUC_N_Fixed(RF_True_TCGA_BRCA, 0.75, 25) #0.671

#write.csv(XGB_True_TCGA_BRCA, "TCGA_Breast_XGB.csv")
#write.csv(RF_True_TCGA_BRCA, "TCGA_Breast_RF.csv")


#Analysis-------------------------------------------------------------------------------------------------

#Read in from CSV saved
setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")


NAFLD_Sim <- read.csv("NAFLD_Sim.csv")[,-1]
HCC_Sim <- read.csv("HCC_Sim.csv")[,-1]
Kidney_Sim <- read.csv("Kidney_Sim.csv")[,-1]
IPF_Sim <- read.csv("IPF_Sim.csv")[,-1]
Tuberculosis_Sim <- read.csv("Tuberculosis_Sim.csv")[,-1]
RA_Sim <- read.csv("RA_Sim.csv")[,-1]
MS_Sim <- read.csv("MS_Sim.csv")[,-1]
Glioma_Sim <- read.csv("Glioma_Sim.csv")[,-1]
CCA_Sim <- read.csv("CCA_Sim.csv")[,-1]
Bipolar_Sim <- read.csv("Bipolar_Sim.csv")[,-1]
MDD_Sim <- read.csv("MDD_Sim.csv")[,-1]
MISC_Sim <- read.csv("MISC_Sim.csv")[,-1]
Hypertension_Sim <- read.csv("Hypertension_Sim.csv")[,-1]
AVSC_Sim <- read.csv("AVSC_Sim.csv")[,-1]
Crohn_Sim <- read.csv("Crohn_Sim.csv")[,-1] 
COVID_Sim <- read.csv("COVID_Sim.csv")[,-1]
ALS_Sim <- read.csv("ALS_Sim.csv")[,-1]
PrePost_Sim <- read.csv("PrePost_Sim.csv")[,-1]
Gastric_Sim <- read.csv("Gastric_Sim.csv")[,-1]
HBV_Sim <- read.csv("HBV_Sim.csv")[,-1]
NSCLC_Sim <- read.csv("NSCLC_Sim.csv")[,-1]
GBM_Sim <- read.csv("GBM_Sim.csv")[,-1]
TB_Sim <- read.csv("TB_Sim.csv")[,-1]
Ovarian_Sim <- read.csv("Ovarian_Sim.csv")[,-1]
EDS_Sim <- read.csv("EDS_Sim.csv")[,-1]
TCGA_Brain <- read.csv("TCGA_Brain.csv")[,-1]
TCGA_Breast <- read.csv("TCGA_BRCA.csv")[,-1]
TCGA_Lung <- read.csv("TCGA_Lung.csv")[,-1]
TCGA_Brain_Test <- read.csv("TCGA_Brain_Test.csv")[,-1]
TCGA_Breast_Test <- read.csv("TCGA_BRCA_Test.csv")[,-1]
TCGA_Lung_Test <- read.csv("TCGA_Lung_Test.csv")[,-1]

IBD <- read.csv("IBD_Keep.csv")[,-1]
IBD_Test <- read.csv("IBD_Test.csv")[,-1]
PDAC <- read.csv("PDAC.csv")[,-1]
PDAC_Test <- read.csv("PDAC_Test.csv")[,-1]



#FullAUC with SE
XGBoost_Optimal(NAFLD_Sim, "Group")[2]*2.25 #0.0001
XGBoost_Optimal(HCC_Sim, "Group")[2]*2.25 #0.0009
XGBoost_Optimal(Kidney_Sim, "Group")[2]*2.25 #0.002
XGBoost_Optimal(IPF_Sim, "Group")[2]*2.25 #0
XGBoost_Optimal(Tuberculosis_Sim, "Group")[2]*2.25 #0.001
XGBoost_Optimal(RA_Sim, "Group")[2]*2.25 #0.001
XGBoost_Optimal(MS_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(Glioma_Sim, "Group")[2]*2.25 #0.002
XGBoost_Optimal(CCA_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(Bipolar_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(MDD_Sim, "Group")[2]*2.25 #0.003
XGBoost_Optimal(MISC_Sim, "Group")[2]*2.25 #0.0003
XGBoost_Optimal(Hypertension_Sim, "Group")[2]*2.25 #0.0002
XGBoost_Optimal(AVSC_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(Crohn_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(COVID_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(ALS_Sim, "Group")[2]*2.25 #0.001
XGBoost_Optimal(PrePost_Sim, "Group")[2]*2.25 #0.0004
XGBoost_Optimal(Gastric_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(HBV_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(NSCLC_Sim, "Group")[2]*2.25 #0.0007
XGBoost_Optimal(GBM_Sim, "Group")[2]*2.25 #0.0001
XGBoost_Optimal(TB_Sim, "Group")[2]*2.25 #0.003
XGBoost_Optimal(Ovarian_Sim, "Group")[2]*2.25 #0.004
XGBoost_Optimal(EDS_Sim, "Group")[2]*2.25 #0.0001


RF_Optimal(NAFLD_Sim, "Group")[2]*2.25 #0.0002
RF_Optimal(HCC_Sim, "Group")[2]*2.25 #0.0004
RF_Optimal(Kidney_Sim, "Group")[2]*2.25 #0.003
RF_Optimal(IPF_Sim, "Group")[2]*2.25 #0
RF_Optimal(Tuberculosis_Sim, "Group")[2]*2.25 #0.003
RF_Optimal(RA_Sim, "Group")[2]*2.25 #0.002
RF_Optimal(MS_Sim, "Group")[2]*2.25 #0.005
RF_Optimal(Glioma_Sim, "Group")[2]*2.25 #0.002
RF_Optimal(CCA_Sim, "Group")[2]*2.25 #0.005
RF_Optimal(Bipolar_Sim, "Group")[2]*2.25 #0.004
RF_Optimal(MDD_Sim, "Group")[2]*2.25 #
RF_Optimal(MISC_Sim, "Group")[2]*2.25 #
RF_Optimal(Hypertension_Sim, "Group")[2]*2.25 #
RF_Optimal(AVSC_Sim, "Group")[2]*2.25 #
RF_Optimal(Crohn_Sim, "Group")[2]*2.25 #
RF_Optimal(COVID_Sim, "Group")[2]*2.25 #
RF_Optimal(ALS_Sim, "Group")[2]*2.25 #
RF_Optimal(PrePost_Sim, "Group")[2]*2.25 #
RF_Optimal(Gastric_Sim, "Group")[2]*2.25 #
RF_Optimal(HBV_Sim, "Group")[2]*2.25 #
RF_Optimal(NSCLC_Sim, "Group")[2]*2.25 #
RF_Optimal(GBM_Sim, "Group")[2]*2.25 #
RF_Optimal(TB_Sim, "Group")[2]*2.25 #
RF_Optimal(Ovarian_Sim, "Group")[2]*2.25 #
RF_Optimal(EDS_Sim, "Group")[2]*2.25 #

#Example for NN (if needed)
NN_SE_Curve_Evaluate_2(TCGA_Breast, "Group", n=660, Test=TCGA_Breast_Test, True=colnames(TCGA_Breast)[-ncol(TCGA_Breast)], h=100)





#TCGA Datasets
#Breast
XGBoost_Optimal(TCGA_Breast, "Group") #0.013
RF_Optimal(TCGA_Breast, "Group") #0.005

#Brain
XGBoost_Optimal(TCGA_Brain, "Group") #0.0001
RF_Optimal(TCGA_Brain, "Group") #0.0001


#Breast
XGBoost_Optimal(TCGA_Lung, "Group") #0.003
RF_Optimal(TCGA_Lung, "Group") #0.001

#IBD
XGBoost_Optimal(IBD, "Group") #0.005
RF_Optimal(IBD, "Group") #0.004

#PDAC
XGBoost_Optimal(PDAC, "Group") #0.006
RF_Optimal(PDAC, "Group") #0.004


#At original size
set.seed(2025)
XGBoost_Optimal(GBM_Sim[sample(1:nrow(GBM_Sim), 32),], "Group")
set.seed(2025)
RF_Optimal(GBM_Sim[sample(1:nrow(GBM_Sim), 32),], "Group")

NN_Curve_Evaluate(GBM_Sim, "Group", n=32*(4/5), True=colnames(GBM_Sim)[-ncol(GBM_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(HCC_Sim[sample(1:nrow(HCC_Sim), 63),], "Group")
set.seed(2025)
RF_Optimal(HCC_Sim[sample(1:nrow(HCC_Sim), 63),], "Group")

NN_Curve_Evaluate(HCC_Sim, "Group", n=63*(4/5), True=colnames(HCC_Sim)[-ncol(HCC_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(Kidney_Sim[sample(1:nrow(Kidney_Sim), 82),], "Group")
set.seed(2025)
RF_Optimal(Kidney_Sim[sample(1:nrow(Kidney_Sim), 82),], "Group")

NN_Curve_Evaluate(Kidney_Sim, "Group", n=82*(4/5), True=colnames(Kidney_Sim)[-ncol(Kidney_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(IPF_Sim[sample(1:nrow(IPF_Sim), 84),], "Group")
set.seed(2025)
RF_Optimal(IPF_Sim[sample(1:nrow(IPF_Sim), 84),], "Group")

NN_Curve_Evaluate(IPF_Sim, "Group", n=84*(4/5), True=colnames(IPF_Sim)[-ncol(IPF_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(COVID_Sim[sample(1:nrow(COVID_Sim), 84),], "Group")
set.seed(2025)
RF_Optimal(COVID_Sim[sample(1:nrow(COVID_Sim), 84),], "Group")

NN_Curve_Evaluate(COVID_Sim, "Group", n=84*(4/5), True=colnames(COVID_Sim)[-ncol(COVID_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(HBV_Sim[sample(1:nrow(HBV_Sim), 95),], "Group")
set.seed(2025)
RF_Optimal(HBV_Sim[sample(1:nrow(HBV_Sim), 95),], "Group")

NN_Curve_Evaluate(HBV_Sim, "Group", n=95*(4/5), True=colnames(HBV_Sim)[-ncol(HBV_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(NAFLD_Sim[sample(1:nrow(Hypertension_Sim), 97),], "Group")
set.seed(2025)
RF_Optimal(NAFLD_Sim[sample(1:nrow(Hypertension_Sim), 97),], "Group")

NN_Curve_Evaluate(Hypertension_Sim, "Group", n=97*(4/5), True=colnames(Hypertension_Sim)[-ncol(Hypertension_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(NAFLD_Sim[sample(1:nrow(NAFLD_Sim), 98),], "Group")
set.seed(2025)
RF_Optimal(NAFLD_Sim[sample(1:nrow(NAFLD_Sim), 98),], "Group")

NN_Curve_Evaluate(NAFLD_Sim, "Group", n=98*(4/5), True=colnames(NAFLD_Sim)[-ncol(NAFLD_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(PrePost_Sim[sample(1:nrow(PrePost_Sim), 99),], "Group")
set.seed(2025)
RF_Optimal(PrePost_Sim[sample(1:nrow(PrePost_Sim), 99),], "Group")

NN_Curve_Evaluate(PrePost_Sim, "Group", n=99*(4/5), True=colnames(PrePost_Sim)[-ncol(PrePost_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(RA_Sim[sample(1:nrow(RA_Sim), 100),], "Group")
set.seed(2025)
RF_Optimal(RA_Sim[sample(1:nrow(RA_Sim), 100),], "Group")

NN_Curve_Evaluate(RA_Sim, "Group", n=100*(4/5), True=colnames(RA_Sim)[-ncol(RA_Sim)], h=100)

XGBoost_Optimal(AVSC_Sim[sample(1:nrow(AVSC_Sim), 101),], "Group")
RF_Optimal(AVSC_Sim[sample(1:nrow(AVSC_Sim), 101),], "Group")

NN_Curve_Evaluate(AVSC_Sim, "Group", n=101*(4/5), True=colnames(AVSC_Sim)[-ncol(AVSC_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(Crohn_Sim[sample(1:nrow(Crohn_Sim), 125),], "Group")
set.seed(2025)
RF_Optimal(Crohn_Sim[sample(1:nrow(Crohn_Sim), 125),], "Group")

NN_Curve_Evaluate(Crohn_Sim, "Group", n=125*(4/5), True=colnames(Crohn_Sim)[-ncol(Crohn_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(MDD_Sim[sample(1:nrow(MDD_Sim), 139),], "Group")
set.seed(2025)
RF_Optimal(MDD_Sim[sample(1:nrow(MDD_Sim), 139),], "Group")

NN_Curve_Evaluate(MDD_Sim, "Group", n=139*(4/5), True=colnames(MDD_Sim)[-ncol(MDD_Sim)], h=100)

set.seed(2025)
XGBoost_Optimal(Ovarian_Sim[sample(1:nrow(Ovarian_Sim), 148),], "Group")
set.seed(2025)
RF_Optimal(Ovarian_Sim[sample(1:nrow(Ovarian_Sim), 148),], "Group")

NN_Curve_Evaluate(Ovarian_Sim, "Group", n=148*(4/5), True=colnames(Ovarian_Sim)[-ncol(Ovarian_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(ALS_Sim[sample(1:nrow(ALS_Sim), 162),], "Group")
set.seed(2025)
RF_Optimal(ALS_Sim[sample(1:nrow(ALS_Sim), 162),], "Group")

NN_Curve_Evaluate(ALS_Sim, "Group", n=162*(4/5), True=colnames(ALS_Sim)[-ncol(ALS_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(MISC_Sim[sample(1:nrow(MISC_Sim), 165),], "Group")
set.seed(2025)
RF_Optimal(MISC_Sim[sample(1:nrow(MISC_Sim), 165),], "Group")

NN_Curve_Evaluate(MISC_Sim, "Group", n=165*(4/5), True=colnames(MISC_Sim)[-ncol(MISC_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(CCA_Sim[sample(1:nrow(CCA_Sim), 170),], "Group")
set.seed(2025)
RF_Optimal(CCA_Sim[sample(1:nrow(CCA_Sim), 170),], "Group")

NN_Curve_Evaluate(CCA_Sim, "Group", n=170*(4/5), True=colnames(CCA_Sim)[-ncol(CCA_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(Glioma_Sim[sample(1:nrow(Glioma_Sim), 176),], "Group")
set.seed(2025)
RF_Optimal(Glioma_Sim[sample(1:nrow(Glioma_Sim), 176),], "Group")

NN_Curve_Evaluate(Glioma_Sim, "Group", n=176*(4/5), True=colnames(Glioma_Sim)[-ncol(Glioma_Sim)], h=100)



set.seed(2025)
XGBoost_Optimal(TB_Sim[sample(1:nrow(TB_Sim), 177),], "Group")
set.seed(2025)
RF_Optimal(TB_Sim[sample(1:nrow(TB_Sim), 177),], "Group")

NN_Curve_Evaluate(TB_Sim, "Group", n=177*(4/5), True=colnames(TB_Sim)[-ncol(TB_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(Tuberculosis_Sim[sample(1:nrow(Tuberculosis_Sim), 182),], "Group")
set.seed(2025)
RF_Optimal(Tuberculosis_Sim[sample(1:nrow(Tuberculosis_Sim), 182),], "Group")

NN_Curve_Evaluate(Tuberculosis_Sim, "Group", n=182*(4/5), True=colnames(Tuberculosis_Sim)[-ncol(Tuberculosis_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(EDS_Sim[sample(1:nrow(EDS_Sim), 200),], "Group")
set.seed(2025)
RF_Optimal(EDS_Sim[sample(1:nrow(EDS_Sim), 200),], "Group")


NN_Curve_Evaluate(EDS_Sim, "Group", n=200*(4/5), True=colnames(EDS_Sim)[-ncol(EDS_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(NSCLC_Sim[sample(1:nrow(NSCLC_Sim), 218),], "Group")
set.seed(2025)
RF_Optimal(NSCLC_Sim[sample(1:nrow(NSCLC_Sim), 218),], "Group")

NN_Curve_Evaluate(NSCLC_Sim, "Group", n=218*(4/5), True=colnames(NSCLC_Sim)[-ncol(NSCLC_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(Bipolar_Sim[sample(1:nrow(Bipolar_Sim), 239),], "Group")
set.seed(2025)
RF_Optimal(Bipolar_Sim[sample(1:nrow(Bipolar_Sim), 239),], "Group")

NN_Curve_Evaluate(Bipolar_Sim, "Group", n=239*(4/5), True=colnames(Bipolar_Sim)[-ncol(Bipolar_Sim)], h=100)


set.seed(2025)
XGBoost_Optimal(MS_Sim[sample(1:nrow(MS_Sim), 474),], "Group")
set.seed(2025)
RF_Optimal(MS_Sim[sample(1:nrow(MS_Sim), 474),], "Group")

NN_Curve_Evaluate(MS_Sim, "Group", n=474*(4/5), True=colnames(MS_Sim)[-ncol(MS_Sim)], h=100)



#XGBOOST CURVES
set.seed(2024)
GBM_XGB <- XGBoost_Pipeline_Curve(GBM_Sim, "Group", 50, 250, steps=10, True=colnames(GBM_Sim)[-ncol(GBM_Sim)])
set.seed(2024)
GBM_RF <- RF_Pipeline_Curve(GBM_Sim, "Group", 50, 250, steps=10, True=colnames(GBM_Sim)[-ncol(GBM_Sim)])

set.seed(2024)
HCC_XGB <- XGBoost_Pipeline_Curve(HCC_Sim, "Group", 50, 250, steps=10, True=colnames(HCC_Sim)[-ncol(HCC_Sim)])
set.seed(2024)
HCC_RF <- RF_Pipeline_Curve(HCC_Sim, "Group", 50, 250, steps=10, True=colnames(HCC_Sim)[-ncol(HCC_Sim)])

set.seed(2024)
Kidney_XGB <- XGBoost_Pipeline_Curve(Kidney_Sim, "Group", 50, 1000, steps=10, True=colnames(Kidney_Sim)[-ncol(Kidney_Sim)])
set.seed(2024)
Kidney_RF <- RF_Pipeline_Curve(Kidney_Sim, "Group", 50, 1000, steps=10, True=colnames(Kidney_Sim)[-ncol(Kidney_Sim)])

set.seed(2024)
IPF_XGB <- XGBoost_Pipeline_Curve(IPF_Sim, "Group", 25, 150, steps=10, True=colnames(IPF_Sim)[-ncol(IPF_Sim)])
set.seed(2024)
IPF_RF <- RF_Pipeline_Curve(IPF_Sim, "Group", 50, 150, steps=10, True=colnames(IPF_Sim)[-ncol(IPF_Sim)])

set.seed(2024)
COVID_XGB <- XGBoost_Pipeline_Curve(COVID_Sim, "Group", 50, 2000, steps=10, True=colnames(COVID_Sim)[-ncol(COVID_Sim)])
set.seed(2024)
COVID_RF <- RF_Pipeline_Curve(COVID_Sim, "Group", 50, 2000, steps=10, True=colnames(COVID_Sim)[-ncol(COVID_Sim)])

set.seed(2024)
HBV_XGB <- XGBoost_Pipeline_Curve(HBV_Sim, "Group", 50, 1000, steps=10, True=colnames(HBV_Sim)[-ncol(HBV_Sim)])
set.seed(2024)
HBV_RF <- RF_Pipeline_Curve(HBV_Sim, "Group", 50, 1000, steps=10, True=colnames(HBV_Sim)[-ncol(HBV_Sim)])

set.seed(2024)
Hypertension_XGB <- XGBoost_Pipeline_Curve(Hypertension_Sim, "Group", 50, 250, steps=10, True=colnames(Hypertension_Sim)[-ncol(Hypertension_Sim)])
set.seed(2024)
Hypertension_RF <- RF_Pipeline_Curve(Hypertension_Sim, "Group", 50, 250, steps=10, True=colnames(Hypertension_Sim)[-ncol(Hypertension_Sim)])

set.seed(2024)
NAFLD_XGB <- XGBoost_Pipeline_Curve(NAFLD_Sim, "Group", 50, 250, steps=10, True=colnames(NAFLD_Sim)[-ncol(NAFLD_Sim)])
set.seed(2024)
NAFLD_RF <- RF_Pipeline_Curve(NAFLD_Sim, "Group", 50, 250, steps=10, True=colnames(NAFLD_Sim)[-ncol(NAFLD_Sim)])

set.seed(2024)
PrePost_XGB <- XGBoost_Pipeline_Curve(PrePost_Sim, "Group", 50, 500, steps=10, True=colnames(PrePost_Sim)[-ncol(PrePost_Sim)])
set.seed(2024)
PrePost_RF <- RF_Pipeline_Curve(PrePost_Sim, "Group", 50, 500, steps=10, True=colnames(PrePost_Sim)[-ncol(PrePost_Sim)])

set.seed(2024)
RA_XGB <- XGBoost_Pipeline_Curve(RA_Sim, "Group", 50, 750, steps=10, True=colnames(RA_Sim)[-ncol(RA_Sim)])
set.seed(2024)
RA_RF <- RF_Pipeline_Curve(RA_Sim, "Group", 50, 750, steps=10, True=colnames(RA_Sim)[-ncol(RA_Sim)])

set.seed(2024)
AVSC_XGB <- XGBoost_Pipeline_Curve(AVSC_Sim, "Group", 50, 1500, steps=10, True=colnames(AVSC_Sim)[-ncol(AVSC_Sim)])
set.seed(2024)
AVSC_RF <- RF_Pipeline_Curve(AVSC_Sim, "Group", 50, 1500, steps=10, True=colnames(AVSC_Sim)[-ncol(AVSC_Sim)])

set.seed(2024)
Crohn_XGB <- XGBoost_Pipeline_Curve(Crohn_Sim, "Group", 50, 2000, steps=10, True=colnames(Crohn_Sim)[-ncol(Crohn_Sim)])
set.seed(2024)
Crohn_RF <- RF_Pipeline_Curve(Crohn_Sim, "Group", 50, 2000, steps=10, True=colnames(Crohn_Sim)[-ncol(Crohn_Sim)])

set.seed(2025)
MDD_XGB <- XGBoost_Pipeline_Curve(MDD_Sim, "Group", 50, 2000, steps=10, True=colnames(MDD_Sim)[-ncol(MDD_Sim)])
set.seed(2024)
MDD_RF <- RF_Pipeline_Curve(MDD_Sim, "Group", 50, 2000, steps=10, True=colnames(MDD_Sim)[-ncol(MDD_Sim)])

set.seed(2025)
Ovarian_XGB <- XGBoost_Pipeline_Curve(Ovarian_Sim, "Group", 50, 2000, steps=10, True=colnames(Ovarian_Sim)[-ncol(Ovarian_Sim)])
set.seed(2024)
Ovarian_RF <- RF_Pipeline_Curve(Ovarian_Sim, "Group", 50, 2000, steps=10, True=colnames(Ovarian_Sim)[-ncol(Ovarian_Sim)])

set.seed(2)
ALS_XGB <- XGBoost_Pipeline_Curve(ALS_Sim, "Group", 50, 2000, steps=10, True=colnames(ALS_Sim)[-ncol(ALS_Sim)])
set.seed(2)
ALS_RF <- RF_Pipeline_Curve(ALS_Sim, "Group", 50, 2000, steps=10, True=colnames(ALS_Sim)[-ncol(ALS_Sim)])

set.seed(2024)
MISC_XGB <- XGBoost_Pipeline_Curve(MISC_Sim, "Group", 50, 250, steps=10, True=colnames(MISC_Sim)[-ncol(MISC_Sim)])
set.seed(2024)
MISC_RF <- RF_Pipeline_Curve(MISC_Sim, "Group", 50, 250, steps=10, True=colnames(MISC_Sim)[-ncol(MISC_Sim)])

set.seed(2024)
CCA_XGB <- XGBoost_Pipeline_Curve(CCA_Sim, "Group", 50, 1000, steps=10, True=colnames(CCA_Sim)[-ncol(CCA_Sim)])
set.seed(2024)
CCA_RF <- RF_Pipeline_Curve(CCA_Sim, "Group", 50, 1000, steps=10, True=colnames(CCA_Sim)[-ncol(CCA_Sim)])

set.seed(2024)
Glioma_XGB <- XGBoost_Pipeline_Curve(Glioma_Sim, "Group", 50, 1000, steps=10, True=colnames(Glioma_Sim)[-ncol(Glioma_Sim)])
set.seed(2024)
Glioma_RF <- RF_Pipeline_Curve(Glioma_Sim, "Group", 50, 1000, steps=10, True=colnames(Glioma_Sim)[-ncol(Glioma_Sim)])

set.seed(2025)
TBProg_XGB <- XGBoost_Pipeline_Curve(TB_Sim, "Group", 50, 2000, steps=10, True=colnames(TB_Sim)[-ncol(TB_Sim)])
set.seed(2024)
TBProg_RF <- RF_Pipeline_Curve(TB_Sim, "Group", 50, 2000, steps=10, True=colnames(TB_Sim)[-ncol(TB_Sim)])

set.seed(2025)
Tuberculosis_XGB <- XGBoost_Pipeline_Curve(Tuberculosis_Sim, "Group", 75, 500, steps=10, True=colnames(Tuberculosis_Sim)[-ncol(Tuberculosis_Sim)])
set.seed(2025)
Tuberculosis_RF <- RF_Pipeline_Curve(Tuberculosis_Sim, "Group", 75, 500, steps=10, True=colnames(Tuberculosis_Sim)[-ncol(Tuberculosis_Sim)])

set.seed(2024)
EDS_XGB <- XGBoost_Pipeline_Curve(EDS_Sim, "Group", 50, 250, steps=10, True=colnames(EDS_Sim)[-ncol(EDS_Sim)])
set.seed(2024)
EDS_RF <- RF_Pipeline_Curve(EDS_Sim, "Group", 50, 250, steps=10, True=colnames(EDS_Sim)[-ncol(EDS_Sim)])

set.seed(2024)
NSCLC_XGB <- XGBoost_Pipeline_Curve(NSCLC_Sim, "Group", 50, 500, steps=10, True=colnames(NSCLC_Sim)[-ncol(NSCLC_Sim)])
set.seed(2024)
NSCLC_RF <- RF_Pipeline_Curve(NSCLC_Sim, "Group", 50, 500, steps=10, True=colnames(NSCLC_Sim)[-ncol(NSCLC_Sim)])

set.seed(2024)
Bipolar_XGB <- XGBoost_Pipeline_Curve(Bipolar_Sim, "Group", 50, 2000, steps=10, True=colnames(Bipolar_Sim)[-ncol(Bipolar_Sim)])
set.seed(2024)
Bipolar_RF <- RF_Pipeline_Curve(Bipolar_Sim, "Group", 50, 2000, steps=10, True=colnames(Bipolar_Sim)[-ncol(Bipolar_Sim)])

set.seed(2024)
MS_XGB <- XGBoost_Pipeline_Curve(MS_Sim, "Group", 50, 2000, steps=10, True=colnames(MS_Sim)[-ncol(MS_Sim)])
set.seed(2024)
MS_RF <- RF_Pipeline_Curve(MS_Sim, "Group", 50, 2000, steps=10, True=colnames(MS_Sim)[-ncol(MS_Sim)])

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(GBM_XGB, "GBM_XGB.csv") 
write.csv(HCC_XGB, "HCC_XGB.csv") 
write.csv(Kidney_XGB, "Kidney_XGB.csv")
write.csv(IPF_XGB, "IPF_XGB.csv") 
write.csv(COVID_XGB, "COVID_XGB.csv") 
write.csv(HBV_XGB, "HBV_XGB.csv") 
write.csv(Hypertension_XGB, "Hypertension_XGB.csv")
write.csv(NAFLD_XGB, "NAFLD_XGB.csv")
write.csv(PrePost_XGB, "PrePost_XGB.csv") 
write.csv(RA_XGB, "RA_XGB.csv") 
write.csv(AVSC_XGB, "AVSC_XGB.csv")  
write.csv(Crohn_XGB, "Crohn_XGB.csv") 
write.csv(MDD_XGB, "MDD_XGB.csv")
write.csv(Ovarian_XGB, "Ovarian_XGB.csv") 
write.csv(ALS_XGB, "ALS_XGB.csv") 
write.csv(MISC_XGB, "MISC_XGB.csv") 
write.csv(CCA_XGB, "CCA_XGB.csv")
write.csv(Glioma_XGB, "Glioma_XGB.csv")
write.csv(TBProg_XGB, "TBProg_XGB.csv")
write.csv(Tuberculosis_XGB, "Tuberculosis_XGB.csv") 
write.csv(EDS_XGB, "EDS_XGB.csv") 
write.csv(NSCLC_XGB, "NSCLC_XGB.csv") 
write.csv(Bipolar_XGB, "Bipolar_XGB.csv")
write.csv(MS_XGB, "MS_XGB.csv") 

write.csv(GBM_RF, "GBM_RF.csv") 
write.csv(HCC_RF, "HCC_RF.csv") 
write.csv(Kidney_RF, "Kidney_RF.csv")
write.csv(IPF_RF, "IPF_RF.csv") 
write.csv(COVID_RF, "COVID_RF.csv") 
write.csv(HBV_RF, "HBV_RF.csv") 
write.csv(Hypertension_RF, "Hypertension_RF.csv")
write.csv(NAFLD_RF, "NAFLD_RF.csv")
write.csv(PrePost_RF, "PrePost_RF.csv") 
write.csv(RA_RF, "RA_RF.csv") 
write.csv(AVSC_RF, "AVSC_RF.csv")  
write.csv(Crohn_RF, "Crohn_RF.csv") 
write.csv(MDD_RF, "MDD_RF.csv")
write.csv(Ovarian_RF, "Ovarian_RF.csv") 
write.csv(ALS_RF, "ALS_RF.csv") 
write.csv(MISC_RF, "MISC_RF.csv") 
write.csv(CCA_RF, "CCA_RF.csv")
write.csv(Glioma_RF, "Glioma_RF.csv")
write.csv(TBProg_RF, "TBProg_RF.csv")
write.csv(Tuberculosis_RF, "Tuberculosis_RF.csv") 
write.csv(EDS_RF, "EDS_RF.csv") 
write.csv(NSCLC_RF, "NSCLC_RF.csv") 
write.csv(Bipolar_RF, "Bipolar_RF.csv")
write.csv(MS_RF, "MS_RF.csv")


AUCs <- c(1, 0.995,0.848,1.000,0.936,0.883,0.999,0.999,0.998,0.981,0.852,0.903,0.862,0.845,0.949,0.999,0.924,0.964,0.893,0.995,1,0.995,0.933,0.999,0.883)

AUCs_RF <- c(1, 0.994,0.850,1.000,0.936,0.884,0.999,0.999,0.996,0.977,0.849,0.896,0.857,0.846,0.943,0.998,0.922,0.963,0.895,0.992,1,0.993,0.930,0.999,0.883)

TCGA_Brain_XGB <- read.csv("TCGA_Brain_XGB.csv")[,-1]
TCGA_Lung_XGB <- read.csv("TCGA_Lung_XGB.csv")[,-1]
TCGA_Breast_XGB <- read.csv("TCGA_Breast_XGB.csv")[,-1]

TCGA_Brain_RF <- read.csv("TCGA_Brain_RF.csv")[,-1]
TCGA_Lung_RF <- read.csv("TCGA_Lung_RF.csv")[,-1]
TCGA_Breast_RF <- read.csv("TCGA_Breast_RF.csv")[,-1]

N_fromcurve <- c(
  Get_N(GBM_XGB, AUCs[1], 2000),
  Get_N(HCC_XGB, AUCs[2], 2000),
  Get_N_LOG(Kidney_XGB, AUCs[3], 2000),
  Get_N_Fixed(IPF_XGB, AUCs[4], 2000),
  Get_N(COVID_XGB, AUCs[5], 2000),
  Get_N_LOG(HBV_XGB, AUCs[6], 2000),
  Get_N_Fixed(Hypertension_XGB, AUCs[7], 2000),
  Get_N(NAFLD_XGB, AUCs[8], 2000),
  Get_N(PrePost_XGB, AUCs[9], 2000),
  Get_N(RA_XGB, AUCs[10], 2000),
  Get_N(AVSC_XGB, AUCs[11], 2000),
  Get_N(Crohn_XGB, AUCs[12], 2000),
  Get_N_LOG(MDD_XGB, AUCs[13], 2000),
  Get_N(Ovarian_XGB, AUCs[14], 2000),
  Get_N(ALS_XGB, AUCs[15], 2000),
  Get_N_Fixed(MISC_XGB, AUCs[16], 2000),
  Get_N(CCA_XGB, AUCs[17], 2000),
  Get_N(Glioma_XGB, AUCs[18], 2000),
  Get_N(TBProg_XGB, AUCs[19], 2000),
  Get_N_Fixed(Tuberculosis_XGB, AUCs[20], 2000),
  Get_N(EDS_XGB, AUCs[21], 2000),
  Get_N(NSCLC_XGB, AUCs[22], 2000),
  Get_N(Bipolar_XGB, AUCs[23], 2000),
  Get_N_Fixed(PDAC_XGB, AUCs[24], 2000),
  Get_N(MS_XGB, AUCs[25], 2000),
  Get_N(TCGA_Brain_XGB, 1, 2000),
  Get_N_Fixed(TCGA_Breast_XGB, 0.718, 2000),
  Get_N(TCGA_Lung_XGB, 0.985, 2000)
)

N_fromcurve_0.01 <- c(
  Get_N_Th(GBM_XGB, AUCs[1], 2000, 0.01),
  Get_N_Th(HCC_XGB, AUCs[2], 2000, 0.01),
  Get_N_LOG_Th(Kidney_XGB, AUCs[3], 2000, 0.01),
  Get_N_Fixed_Th(IPF_XGB, AUCs[4], 2000, 0.01),
  Get_N_Th(COVID_XGB, AUCs[5], 2000, 0.01),
  Get_N_LOG_Th(HBV_XGB, AUCs[6], 2000, 0.01),
  Get_N_Fixed_Th(Hypertension_XGB, AUCs[7], 2000, 0.01),
  Get_N_Th(NAFLD_XGB, AUCs[8], 2000, 0.01),
  Get_N_Th(PrePost_XGB, AUCs[9], 2000, 0.01),
  Get_N_Th(RA_XGB, AUCs[10], 2000, 0.01),
  Get_N_Th(AVSC_XGB, AUCs[11], 2500, 0.01),
  Get_N_Th(Crohn_XGB, AUCs[12], 2000, 0.01),
  Get_N_LOG_Th(MDD_XGB, AUCs[13], 2000, 0.01),
  Get_N_Th(Ovarian_XGB, AUCs[14], 2000, 0.01),
  Get_N_Th(ALS_XGB, AUCs[15], 2000, 0.01),
  Get_N_Fixed_Th(MISC_XGB, AUCs[16], 2000, 0.01),
  Get_N_Th(CCA_XGB, AUCs[17], 2000, 0.01),
  Get_N_Th(Glioma_XGB, AUCs[18], 2000, 0.01),
  Get_N_Th(TBProg_XGB, AUCs[19], 2000, 0.01),
  Get_N_Fixed_Th(Tuberculosis_XGB, AUCs[20], 2000, 0.01),
  Get_N_Th(EDS_XGB, AUCs[21], 2000, 0.01),
  Get_N_Th(NSCLC_XGB, AUCs[22], 2000, 0.01),
  Get_N_Th(Bipolar_XGB, AUCs[23], 2000, 0.01),
  Get_N_Fixed_Th(PDAC_XGB, AUCs[24], 2000, 0.01),
  Get_N_Th(MS_XGB, AUCs[25], 2000, 0.01),
  Get_N_Th(TCGA_Brain_XGB, 1, 2000, 0.01),
  Get_N_Fixed_Th(TCGA_Breast_XGB, 0.718, 2500, 0.01),
  Get_N_Th(TCGA_Lung_XGB, 0.985, 2000, 0.01)
)

N_fromcurve_0.01
mean(c(87,303 ,1619 , 175, 1376 ,1138 , 279  ,134,  360 , 907 ,2500, 1311, 1536, 1377 ,1763,  120, 1034, 841 ,1591 , 483,  209 , 450,
         1707 , 1780 , 112  ,2500, 195))
sd(c(87,303 ,1619 , 175, 1376 ,1138 , 279  ,134,  360 , 907 ,2500, 1311, 1536, 1377 ,1763,  120, 1034, 841 ,1591 , 483,  209 , 450,
        1707 , 1780 , 112  ,2500, 195))


  median(c(51 ,2500, 706   , 25, 1792,  404,  777  , 84,  2500, 1542, 1361 ,1325 , 354, 1810 ,2265  , 92 ,2012, 2148  ,2500,  695 , 180 , 806,
 2500,2068   ,25 , 485  ,167))
range(c(51 ,2500, 706   , 25, 1792,  404,  777  , 84,  2500, 1542, 1361 ,1325 , 354, 1810 ,2265  , 92 ,2012, 2148  ,2500,  695 , 180 , 806,
        2500,2068   ,25 , 485  ,167))

N_fromcurve_0.05 <- c(
  Get_N_Th(GBM_XGB, AUCs[1], 2000, 0.05),
  Get_N_Th(HCC_XGB, AUCs[2], 2000, 0.05),
  Get_N_LOG_Th(Kidney_XGB, AUCs[3], 2000, 0.05),
  Get_N_Fixed_Th(IPF_XGB, AUCs[4], 2000, 0.05),
  Get_N_Th(COVID_XGB, AUCs[5], 2000, 0.05),
  Get_N_LOG_Th(HBV_XGB, AUCs[6], 2000, 0.05),
  Get_N_Fixed_Th(Hypertension_XGB, AUCs[7], 2000, 0.05),
  Get_N_Th(NAFLD_XGB, AUCs[8], 2000, 0.05),
  Get_N_Th(PrePost_XGB, AUCs[9], 2000, 0.05),
  Get_N_Th(RA_XGB, AUCs[10], 2000, 0.05),
  Get_N_Th(AVSC_XGB, AUCs[11], 2500, 0.05),
  Get_N_Th(Crohn_XGB, AUCs[12], 2000, 0.05),
  Get_N_LOG_Th(MDD_XGB, AUCs[13], 2000, 0.05),
  Get_N_Th(Ovarian_XGB, AUCs[14], 2000, 0.05),
  Get_N_Th(ALS_XGB, AUCs[15], 2000, 0.05),
  Get_N_Fixed_Th(MISC_XGB, AUCs[16], 2000, 0.05),
  Get_N_Th(CCA_XGB, AUCs[17], 2000, 0.05),
  Get_N_Th(Glioma_XGB, AUCs[18], 2000, 0.05),
  Get_N_Th(TBProg_XGB, AUCs[19], 2000, 0.05),
  Get_N_Fixed_Th(Tuberculosis_XGB, AUCs[20], 2000, 0.05),
  Get_N_Th(EDS_XGB, AUCs[21], 2000, 0.05),
  Get_N_Th(NSCLC_XGB, AUCs[22], 2000, 0.05),
  Get_N_Th(Bipolar_XGB, AUCs[23], 2000, 0.05),
  Get_N_Fixed_Th(PDAC_XGB, AUCs[24], 2000, 0.05),
  Get_N_Th(MS_XGB, AUCs[25], 2000, 0.05),
  Get_N_Th(TCGA_Brain_XGB, 1, 2000, 0.05),
  Get_N_Fixed_Th(TCGA_Breast_XGB, 0.718, 2500, 0.05),
  Get_N_Th(TCGA_Lung_XGB, 0.985, 2000, 0.05)
)

N_fromcurve_0.05

mean(c(
30 ,45, 103,  29 ,157,  53 , 44  ,28,  53, 122, 196, 148, 225, 164, 360 , 25, 202, 137, 286, 125,  58,  71, 247,  60, 477,  30,  69,  40)[-24])
sd(c(
  30 ,45, 103,  29 ,157,  53 , 44  ,28,  53, 122, 196, 148, 225, 164, 360 , 25, 202, 137, 286, 125,  58,  71, 247,  60, 477,  30,  69,  40)[-24])

median(c(29,  54 , 51 ,  25, 619,  37,  38 , 25,  48,  92, 383,  93  ,53 , 25, 621,  38 ,178, 122, 482,  58 , 44 , 96 ,369, 968,   25,  69,  25))
range(c(29,  54 , 51 ,  25, 619,  37,  38 , 25,  48,  92, 383,  93  ,53 , 25, 621,  38 ,178, 122, 482,  58 , 44 , 96 ,369, 968,   25,  69,  25))
  
N_RF_fromcurve <- c(
  Get_N(GBM_RF, AUCs_RF[1], 2000),
  Get_N(HCC_RF, AUCs_RF[2], 2000),
  Get_N_LOG(Kidney_RF, AUCs_RF[3], 2000),
  Get_N_LOG(IPF_RF, AUCs_RF[4], 2000),
  Get_N(COVID_RF, AUCs_RF[5], 2000),
  Get_N_LOG(HBV_RF, AUCs_RF[6], 2000),
  Get_N_Fixed(Hypertension_RF, AUCs_RF[7], 2000),
  Get_N_Fixed(NAFLD_RF, AUCs_RF[8], 2000),
  Get_N(PrePost_RF, AUCs_RF[9], 2000),
  Get_N(RA_RF, AUCs_RF[10], 2000),
  Get_N(AVSC_RF, AUCs_RF[11], 2000),
  Get_N(Crohn_RF, AUCs_RF[12], 2000),
  Get_N(MDD_RF, AUCs_RF[13], 2000),
  Get_N(Ovarian_RF, AUCs_RF[14], 2000),
  Get_N(ALS_RF, AUCs_RF[15], 2000),
  Get_N_Fixed(MISC_RF, AUCs_RF[16], 2000),
  Get_N(CCA_RF, AUCs_RF[17], 2000),
  Get_N(Glioma_RF, AUCs_RF[18], 2000),
  Get_N(TBProg_RF, AUCs_RF[19], 2000),
  Get_N(Tuberculosis_RF, AUCs_RF[20], 2000),
  Get_N(EDS_RF, AUCs_RF[21], 2000),
  Get_N(NSCLC_RF, AUCs_RF[22], 2000),
  Get_N(Bipolar_RF, AUCs_RF[23], 2000),
  Get_N_Fixed(PDAC_RF, AUCs_RF[24], 2000),
  Get_N(MS_RF, AUCs_RF[25], 2000),
  Get_N_Fixed(TCGA_Brain_RF, 1, 2000),
  Get_N_LOG(TCGA_Breast_RF, 0.75, 2000),
  Get_N_Fixed(TCGA_Lung_RF, 0.983, 2000)
)

N_RF_fromcurve_0.01 <- c(
  Get_N_Th(GBM_RF, AUCs_RF[1], 2000, 0.01),
  Get_N_Th(HCC_RF, AUCs_RF[2], 2000, 0.01),
  Get_N_LOG_Th(Kidney_RF, AUCs_RF[3], 2000, 0.01),
  Get_N_LOG_Th(IPF_RF, AUCs_RF[4], 2000, 0.01),
  Get_N_Th(COVID_RF, AUCs_RF[5], 2000, 0.01),
  Get_N_LOG_Th(HBV_RF, AUCs_RF[6], 2000, 0.01),
  Get_N_Fixed_Th(Hypertension_RF, AUCs_RF[7], 2000, 0.01),
  Get_N_Fixed_Th(NAFLD_RF, AUCs_RF[8], 2000, 0.01),
  Get_N_Th(PrePost_RF, AUCs_RF[9], 2000, 0.01),
  Get_N_Th(RA_RF, AUCs_RF[10], 2000, 0.01),
  Get_N_Th(AVSC_RF, AUCs_RF[11], 2000, 0.01),
  Get_N_Th(Crohn_RF, AUCs_RF[12], 2000, 0.01),
  Get_N_Th(MDD_RF, AUCs_RF[13], 2000, 0.01),
  Get_N_Th(Ovarian_RF, AUCs_RF[14], 2000, 0.01),
  Get_N_Th(ALS_RF, AUCs_RF[15], 2000, 0.01),
  Get_N_Fixed_Th(MISC_RF, AUCs_RF[16], 2000, 0.01),
  Get_N_Th(CCA_RF, AUCs_RF[17], 2000, 0.01),
  Get_N_Th(Glioma_RF, AUCs_RF[18], 2000, 0.01),
  Get_N_Th(TBProg_RF, AUCs_RF[19], 2000, 0.01),
  Get_N_Th(Tuberculosis_RF, AUCs_RF[20], 2000, 0.01),
  Get_N_Th(EDS_RF, AUCs_RF[21], 2000, 0.01),
  Get_N_Th(NSCLC_RF, AUCs_RF[22], 2000, 0.01),
  Get_N_Th(Bipolar_RF, AUCs_RF[23], 2000, 0.01),
  Get_N_Fixed_Th(PDAC_RF, AUCs_RF[24], 2000, 0.01),
  Get_N_Th(MS_RF, AUCs_RF[25], 2000, 0.01),
  Get_N_Fixed_Th(TCGA_Brain_RF, 1, 2000, 0.01),
  Get_N_LOG_Th(TCGA_Breast_RF, 0.75, 2000, 0.01),
  Get_N_Fixed_Th(TCGA_Lung_RF, 0.983, 2000, 0.01)
)

N_RF_fromcurve_0.01

median(c(42 , 531 , 918 ,   25, 881 , 463,  142 ,  52,  222 , 795 ,1275 , 404 , 654  ,639, 1752 ,  33 , 943,  772 ,1797,  831,  109 , 463
         ,1518  , 2500,   25, 686  , 83))

N_RF_fromcurve_0.05 <- c(
  Get_N_Th(GBM_RF, AUCs_RF[1], 2000, 0.05),
  Get_N_Th(HCC_RF, AUCs_RF[2], 2000, 0.05),
  Get_N_LOG_Th(Kidney_RF, AUCs_RF[3], 2000, 0.05),
  Get_N_LOG_Th(IPF_RF, AUCs_RF[4], 2000, 0.05),
  Get_N_Th(COVID_RF, AUCs_RF[5], 2000, 0.05),
  Get_N_LOG_Th(HBV_RF, AUCs_RF[6], 2000, 0.05),
  Get_N_Fixed_Th(Hypertension_RF, AUCs_RF[7], 2000, 0.05),
  Get_N_Fixed_Th(NAFLD_RF, AUCs_RF[8], 2000, 0.05),
  Get_N_Th(PrePost_RF, AUCs_RF[9], 2000, 0.05),
  Get_N_Th(RA_RF, AUCs_RF[10], 2000, 0.05),
  Get_N_Th(AVSC_RF, AUCs_RF[11], 2000, 0.05),
  Get_N_Th(Crohn_RF, AUCs_RF[12], 2000, 0.05),
  Get_N_Th(MDD_RF, AUCs_RF[13], 2000, 0.05),
  Get_N_Th(Ovarian_RF, AUCs_RF[14], 2000, 0.05),
  Get_N_Th(ALS_RF, AUCs_RF[15], 2000, 0.05),
  Get_N_Fixed_Th(MISC_RF, AUCs_RF[16], 2000, 0.05),
  Get_N_Th(CCA_RF, AUCs_RF[17], 2000, 0.05),
  Get_N_Th(Glioma_RF, AUCs_RF[18], 2000, 0.05),
  Get_N_Th(TBProg_RF, AUCs_RF[19], 2000, 0.05),
  Get_N_Th(Tuberculosis_RF, AUCs_RF[20], 2000, 0.05),
  Get_N_Th(EDS_RF, AUCs_RF[21], 2000, 0.05),
  Get_N_Th(NSCLC_RF, AUCs_RF[22], 2000, 0.05),
  Get_N_Th(Bipolar_RF, AUCs_RF[23], 2000, 0.05),
  Get_N_Fixed_Th(PDAC_RF, AUCs_RF[24], 2000, 0.05),
  Get_N_Th(MS_RF, AUCs_RF[25], 2000, 0.05),
  Get_N_Fixed_Th(TCGA_Brain_RF, 1, 2000, 0.05),
  Get_N_LOG_Th(TCGA_Breast_RF, 0.75, 2000, 0.05),
  Get_N_Fixed_Th(TCGA_Lung_RF, 0.983, 2000, 0.05)
)

N_RF_fromcurve_0.05
mean(c(25,  25,   25,25  ,94 ,  25,25 ,  25  ,35 , 64 ,136 , 71 , 40  ,25, 280,   25, 102,  53, 301,  96,  25,  32 ,125  ,365 ,  25, 96 ,  25))
sd(c(25,  25,   25,25  ,94 ,  25,25 ,  25  ,35 , 64 ,136 , 71 , 40  ,25, 280,   25, 102,  53, 301,  96,  25,  32 ,125  ,365 ,  25, 96 ,  25))

#CCA AVSC Bipolar

AUC25 <- c(
  AUC_N(GBM_XGB, 25),
  AUC_N(HCC_XGB, 25),
  AUC_N_LOG(Kidney_XGB, 25),
  AUC_N_Fixed(IPF_XGB, 1, 25),
  AUC_N(COVID_XGB, 25),
  AUC_N_LOG(HBV_XGB, 25),
  AUC_N_Fixed(Hypertension_XGB, 0.999, 25),
  AUC_N(NAFLD_XGB, 25),
  AUC_N(PrePost_XGB, 25),
  AUC_N(RA_XGB, 25),
  AUC_N(AVSC_XGB, 25),
  AUC_N(Crohn_XGB, 25),
  AUC_N_LOG(MDD_XGB, 25),
  AUC_N(Ovarian_XGB, 25),
  AUC_N(ALS_XGB, 25),
  AUC_N_Fixed(MISC_XGB, 0.999, 25),
  AUC_N(CCA_XGB, 25),
  AUC_N(Glioma_XGB, 25),
  AUC_N(TBProg_XGB, 25),
  AUC_N_Fixed(Tuberculosis_XGB, 0.995, 25),
  AUC_N(EDS_XGB, 25),
  AUC_N(NSCLC_XGB, 25),
  AUC_N(Bipolar_XGB, 25),
  AUC_N(MS_XGB, 25),
  AUC_N(TCGA_BRAIN_XGB, 25),
  AUC_N_LOG(TCGA_Breast_RF, 25),
  AUC_N(TCGA_LUNG_XGB, 25)
)


AUC25 <- round(AUC25, 3)

AUC_DF <- data.frame("AUC25"=AUC25*100, "AUC_Full"=AUC_Full, "High_Seperability"=as.factor(ifelse(AUC_Full>99, 1, 0)))

ggplot(data=AUC_DF, aes(x=AUC25, y=AUC_Full_NN, color=High_Seperability)) + geom_point(size=3) + labs(x="AUC at n=25 (Best-Performing Algorithm)", y="Maximum Achievable AUC (Best-Performing Algorithm)") + theme_bw()

#RANDOM FOREST CURVES ------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------

#Dispersions
Dispersion <- c(median(GetDispersion(GBM, F, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)])),
                          median(GetDispersion(HCC, F, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)])),
                          median(GetDispersion(Kidney, F, colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)])),
                          median(GetDispersion(IPF, F, colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)])),
                          median(GetDispersion(COVID, T, colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)])),
                          median(GetDispersion(HBV, F, colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)])),
                          median(GetDispersion(Hypertension, F, colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)])),
                          median(GetDispersion(NAFLD, F, colnames(NAFLD_Sim)[1:(length(colnames(NAFLD_Sim))-1)])),
                          median(GetDispersion(PrePost, F, colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)])),
                          median(GetDispersion(RA, F, colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)])),
                          median(GetDispersion(AVSC, T, colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)])),
                          median(GetDispersion(Crohn, F, colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)])),
                          median(GetDispersion(MDD, F, colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)])),
                          median(GetDispersion(Ovarian, F, colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)])),
                          median(GetDispersion(ALS, F, colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)])),
                          median(GetDispersion(MISC, F, colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)])),
                          median(GetDispersion(CCA, F, colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)])),
                          median(GetDispersion(Glioma, F, colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)])),
                          median(GetDispersion(TB, F, colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)])),
                          median(GetDispersion(Tub, F, colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)])),
                          median(GetDispersion(EDS, F, colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)])),
                          median(GetDispersion(NSCLC, F, colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)])),
                          median(GetDispersion(Bipolar, F, colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)])),
                          median(GetDispersion(MS, F, colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])))


Dispersion
# [1] 0.25418221 0.74206836 0.76029293 0.02470538 2.52212131 0.13857525 0.38745665 0.09856993 0.31207956 0.07373518
#[11] 0.85258905 0.17540939 0.16950090 0.26494593 0.12460711 0.32217798 1.00853593 0.43995221 2.36158105 0.06238073
#[21] 0.38176059 0.20305474 0.06723777 0.43745899


#Correlations
median(GetCorrelation(NAFLD, colnames(NAFLD_Sim)[1:52])) #0.54
median(GetCorrelation(HCC, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)])) #0.55
median(GetCorrelation(COVID, colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)])) #0.54
median(GetCorrelation(Kidney, colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)])) #0.67
median(GetCorrelation(IPF, colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)])) #0.13
median(GetCorrelation(Fibrosis, colnames(Fibrosis_Sim)[1:(length(colnames(Fibrosis_Sim))-1)])) #0.48
median(GetCorrelation(RA, colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)])) #0.44
median(GetCorrelation(AVSC, colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)])) #0.33
median(GetCorrelation(Crohn, colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)])) #0.06
median(GetCorrelation(Tub, colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)])) #0.26
median(GetCorrelation(MISC, colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)])) #0.63
median(GetCorrelation(MS, colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])) #0.33
median(GetCorrelation(Hypertension, colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)])) #0.57
median(GetCorrelation(CCA, colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)])) #0.24
median(GetCorrelation(Glioma, colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)])) #0.26
median(GetCorrelation(MDD, colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)])) #0.46
median(GetCorrelation(Bipolar, colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)])) #0.66
median(GetCorrelation(GC6, colnames(GC6_Sim)[1:(length(colnames(GC6_Sim))-1)])) #0.24
median(GetCorrelation(ALS, colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)])) #0.14
median(GetCorrelation(PrePost, colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)])) #0.32
median(GetCorrelation(PDAC, colnames(PDAC_Sim)[1:(length(colnames(PDAC_Sim))-1)])) #0.32
median(GetCorrelation(Gastric, colnames(Gastric_Sim)[1:(length(colnames(Gastric_Sim))-1)])) #0.16
median(GetCorrelation(HBV, colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)])) #0.08
median(GetCorrelation(NSCLC, colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)])) #0.43
median(GetCorrelation(GBM, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)])) #0.36
median(GetCorrelation(IBD, colnames(IBD_Sim)[1:(length(colnames(IBD_Sim))-1)])) #0.25
median(GetCorrelation(TB, colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)])) #0.17
median(GetCorrelation(Ovarian, colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)])) #0.22
median(GetCorrelation(EDS, colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)])) #0.32
median(GetCorrelation(BPD, colnames(BPD_Sim)[1:(length(colnames(BPD_Sim))-1)])) #0.24


#Characteristics
N <- c(54, 168, 813, 81,631,529,126, 72, 154,455, 1083, 613,950, 713, 1048,60,642,480,904,270,132, 235, 886, 1200, 68, 607, 86)
N_RF <- c(50, 128, 152,25,400,106,56,50,93,344,533,191,263,215,970,50,473,353,1005,357,50,189,644,1307, 25, 419, 25)
N_NN <- c(50,324,227,25,1374,144,212,50,199,761,888,681,156,565,1537,57,960,1047,1444,239,94,371,1471,1711, 25, 298, 69)

#N_DF <- data.frame("Sample Size"=c(N,N_RF,N_NN), "Algorithm"=c(rep("XGBoost",27), rep("Random Forest",27), rep("Neural Network",27)))
#ggplot(data=N_DF, aes(x=Algorithm, y=Sample.Size, color=Algorithm)) + geom_boxplot() + ylab("Sample Size") + xlab("")


Imbalance <- c(43.8, 28.6,50.0,41.7,43.5,30.5,26.8,48.0,49.5,49.0,39.6,17.6,49.6,39.9,10.5,17.6,50.0,25.0,17.5,7.7,50.0,23.4,36.4,18.1,41.0, 48.5, 39.6)

#Dispersion <- c(0.47, 0.02, 0.25, 0.38, 0.06, 0.06, 0.61, 0.14, 0.22, 0.80, 0.43, 0.06, 0.31, 0.07, 0.44)

Correlation <- c(0.36,0.55,0.13, 0.67, 0.54, 0.08, 0.57, 0.54, 0.32, 0.44, 0.33, 0.06,0.46, 0.22,0.14, 0.63, 0.24, 0.26, 0.17,0.26, 0.32, 0.43,0.66,0.33,0.60,0.25,0.56)
#Correlation <- ifelse(Correlation > 0.25, 1, 0)

#Correlation <- c(1.36, 1.02, 20.09, 1.23, 4.23, 1.33, 1.08, 1.06, 2.85, 1.18, 1.05, 2.27, 1.45, 2.30, 1.15)
#Correlation <- ifelse(Correlation > 3, 1, 0)

AUC_Full <- c(1, 0.995,0.848,1.000,0.936,0.883,0.999,0.999,0.998,0.981,0.852,0.903,0.862,0.845,0.949,0.999,0.924,0.964,0.893,0.995,1,0.995,0.933,0.883,1,0.718,0.985)*100
AUC_Full_RF <- c(1, 0.994,0.850,1.000,0.936,0.884,0.999,0.999,0.996,0.977,0.849,0.896,0.857,0.846,0.943,0.998,0.922,0.963,0.895,0.992,1,0.993,0.930,0.883,1,0.75,0.983)*100
AUC_Full_NN <- c(1.000,0.994,0.862,1.000,0.947,0.896,0.999,1.000,0.997,0.988,0.858,0.913,0.874,0.861,0.964,1.000,0.935,0.975,0.922,0.998,1.000,0.996,0.951,0.907,1,0.724,0.981)*100


#Change methods figure
#Add figure showing data set characteristics
#Figure 2 make black star a dot, add dots showing original sample size to TCGA ones. Also re-order 1st panel
#to order by cases where original sample is sufficient
#Maybe make fitted model formulas (for dissertaiton)

set.seed(2024)
LR_AUCs <- c(LR_Optimal(cbind(log2(GBM_Sim[,-ncol(GBM_Sim)]+0.01), "Group"=GBM_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(HCC_Sim[,-ncol(HCC_Sim)]+0.01), "Group"=HCC_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Kidney_Sim[,-ncol(Kidney_Sim)]+0.01), "Group"=Kidney_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(IPF_Sim[,-ncol(IPF_Sim)]+0.01), "Group"=IPF_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(COVID_Sim[,-ncol(COVID_Sim)]+0.01), "Group"=COVID_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(HBV_Sim[,-ncol(HBV_Sim)]+0.01), "Group"=HBV_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Hypertension_Sim[,-ncol(Hypertension_Sim)]+0.01), "Group"=Hypertension_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(NAFLD_Sim[,-ncol(NAFLD_Sim)]+0.01), "Group"=NAFLD_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(PrePost_Sim[,-ncol(PrePost_Sim)]+0.01), "Group"=PrePost_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(RA_Sim[,-ncol(RA_Sim)]+0.01), "Group"=RA_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(AVSC_Sim[,-ncol(AVSC_Sim)]+0.01), "Group"=AVSC_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Crohn_Sim[,-ncol(Crohn_Sim)]+0.01), "Group"=Crohn_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(MDD_Sim[,-ncol(MDD_Sim)]+0.01), "Group"=MDD_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Ovarian_Sim[,-ncol(Ovarian_Sim)]+0.01), "Group"=Ovarian_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(ALS_Sim[,-ncol(ALS_Sim)]+0.01), "Group"=ALS_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(MISC_Sim[,-ncol(MISC_Sim)]+0.01), "Group"=MISC_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(CCA_Sim[,-ncol(CCA_Sim)]+0.01), "Group"=CCA_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Glioma_Sim[,-ncol(Glioma_Sim)]+0.01), "Group"=Glioma_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(TB_Sim[,-ncol(TB_Sim)]+0.01), "Group"=TB_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Tuberculosis_Sim[,-ncol(Tuberculosis_Sim)]+0.01), "Group"=Tuberculosis_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(EDS_Sim[,-ncol(EDS_Sim)]+0.01), "Group"=EDS_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(NSCLC_Sim[,-ncol(NSCLC_Sim)]+0.01), "Group"=NSCLC_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(Bipolar_Sim[,-ncol(Bipolar_Sim)]+0.01), "Group"=Bipolar_Sim$Group), "Group")[1],
             LR_Optimal(cbind(log2(MS_Sim[,-ncol(MS_Sim)]+0.01), "Group"=MS_Sim$Group), "Group")[1],
             0.998,
             0.711,
             0.949
            
)*100

MaxAUC <- data.frame("AUC_Full"=AUC_Full, "AUC_Full_RF"=AUC_Full_RF, "AUC_Full_NN"=AUC_Full_NN)
MaxAUC$MAX <- apply(MaxAUC, 1, max)

Complexity <- MaxAUC$MAX - LR_AUCs
plot(Complexity, N)
#Complexity <- ifelse(Complexity > 2.5, 1, 0)

TEST <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "Complexity"=Complexity)

N_Plot <- ggplot(data=TEST, aes(x=Complexity, y=N)) + geom_point(color="blue", size=2) + geom_vline(xintercept=4.5 ,linetype = 'dotted', color="red", size=1.5) + theme_bw() + labs(title="XGB", y="Sample Size", x="Nonlinearity") 
N_Plot_RF <- ggplot(data=TEST, aes(x=Complexity, y=N_RF)) + geom_point(color="blue", size=2) + geom_vline(xintercept=4.5 ,linetype = 'dotted', color="red", size=1.5) + theme_bw() + labs(title="RF", y="Sample Size", x="Nonlinearity") 
N_Plot_NN <- ggplot(data=TEST, aes(x=Complexity, y=N_NN)) + geom_point(color="blue", size=2) + geom_vline(xintercept=4.5 ,linetype = 'dotted', color="red", size=1.5) + theme_bw() + labs(title="NN", y="Sample Size", x="Nonlinearity") 

library(ggpubr)
ggarrange(N_Plot, N_Plot_RF, N_Plot_NN, ncol=3, nrow=1)



#TO derive 4.5 cutoff for nonlinearity
plot(seq(2,8, by=0.5), sapply(seq(2,8, by=0.5), function(i) {
sqrt(cv.glm(glm.nb(N ~ I(Complexity>i)), data=TEST)$delta[1])}),xlab="Cutoff", ylab="LOO-CV MSE", main="XGB")

plot(seq(2,8, by=0.5), sapply(seq(2,8, by=0.5), function(i) {
  sqrt(cv.glm(glm.nb(N_RF ~ I(Complexity>i)), data=TEST)$delta[1]) }), xlab="Cutoff", ylab="LOO-CV MSE", main="RF")

plot(seq(2,8, by=0.5), sapply(seq(2,8, by=0.5), function(i) {
  sqrt(cv.glm(glm.nb(N_NN ~ I(Complexity>i)), data=TEST)$delta[1]) }),xlab="Cutoff", ylab="LOO-CV MSE", main="NN")



#----END

N_Features <- c(8,9,3,58,5,3,13,52,16,12,3,4,4,5,8,41,9,11,7,13,71,37,16,12,102,18,103)

Med_LFC <- c(1.12,2.32, 0.43, 2.96, 3.71, 0.86, 1.64, 0.59, 1.10, 0.34, 1.15, 1.09,0.45, 0.60,0.73,1.48, 1.52, 1.37, 1.19, 0.42, 0.86, 1.25,0.64,0.39, 3.39, 0.83, 4.42)

Min_LFC <- c(0.23,0.88,0.32, 0.74, 3.17, 0.56, 0.42, 0.27, 0.39, 0.21, 0.92, 0.19, 0.33,0.34, 0.15, 0.16, 0.46, 0.52, 0.57, 0.20, 0.20,0.67, 0.21,0.09, 0.43, 0.18, 0.21)

Max_LFC <- c(2.13,5.64, 1.44, 7.63, 4.92, 2.32, 4.70, 2.32, 2.61, 1.06, 2.23, 2.00, 0.66, 2.24,2.50, 5.33, 2.99, 2.84, 3.45, 1.43, 5.16,2.76, 1.95, 1.18, 8.62, 1.90, 11.19)

AvgMedRead <- c(round(c(median(AvgReadCount(GBM_Sim)),
                median(AvgReadCount(HCC_Sim)),
                median(AvgReadCount(Kidney_Sim)),
                median(AvgReadCount(IPF_Sim)),
                median(AvgReadCount(COVID_Sim)),
                median(AvgReadCount(HBV_Sim)),
                median(AvgReadCount(Hypertension_Sim)),
                median(AvgReadCount(NAFLD_Sim)),
                median(AvgReadCount(PrePost_Sim)),
                median(AvgReadCount(RA_Sim)),
                median(AvgReadCount(AVSC_Sim)),
                median(AvgReadCount(Crohn_Sim)),
                median(AvgReadCount(MDD_Sim)),
                median(AvgReadCount(Ovarian_Sim)),
                median(AvgReadCount(ALS_Sim)),
                median(AvgReadCount(MISC_Sim)),
                median(AvgReadCount(CCA_Sim)),
                median(AvgReadCount(Glioma_Sim)),
                median(AvgReadCount(TB_Sim)),
                median(AvgReadCount(Tuberculosis_Sim)),
                median(AvgReadCount(EDS_Sim)),
                median(AvgReadCount(NSCLC_Sim)),
                median(AvgReadCount(Bipolar_Sim)),
                median(AvgReadCount(MS_Sim))),2),
                35.49,
                615.84,
                1128.37)

Dispersion_Simulated <- c(round(c(0.32241168, 1.3329274, 0.87762952, 0.02871027, 9.64928220, 0.14182734, 1.31269917, 0.11748022, 0.41902412,
                                0.07555905, 2.04585115, 0.55446814, 0.16079212, 0.35024109, 2.86843758, 0.35780604, 1.50789033, 1.01346870,
                                9.73620391, 0.07331635, 0.38790421, 0.23388442, 0.07473904, 0.8003593), 2), 1.49, 0.84, 2.86)



#Univariable
#Imbalance
summary(glm.nb(N ~ Imbalance))
exp(coef(glm.nb(N ~ Imbalance)))
exp(confint(glm.nb(N ~ Imbalance)))

summary(glm.nb(N_RF ~ Imbalance))
exp(coef(glm.nb(N_RF ~ Imbalance)))
exp(confint(glm.nb(N_RF ~ Imbalance)))

summary(glm.nb(N_NN ~ Imbalance))
exp(coef(glm.nb(N_NN ~ Imbalance)))
exp(confint(glm.nb(N_NN ~ Imbalance)))

#AUC FULL
summary(glm.nb(N ~ AUC_Full))
exp(coef(glm.nb(N ~ AUC_Full)))
exp(confint(glm.nb(N ~ AUC_Full)))

summary(glm.nb(N_RF ~ AUC_Full_RF))
exp(coef(glm.nb(N_RF ~ AUC_Full_RF)))
exp(confint(glm.nb(N_RF ~ AUC_Full_RF)))

summary(glm.nb(N_NN ~ AUC_Full_NN))
exp(coef(glm.nb(N_NN ~ AUC_Full_NN)))
exp(confint(glm.nb(N_NN ~ AUC_Full_NN)))

#NUMBER of FEATURES 
summary(glm.nb(N ~ N_Features))
exp(coef(glm.nb(N ~ N_Features)))
exp(confint(glm.nb(N ~ N_Features)))

summary(glm.nb(N_RF ~ N_Features))
exp(coef(glm.nb(N_RF ~ N_Features)))
exp(confint(glm.nb(N_RF ~ N_Features)))

summary(glm.nb(N_NN ~ N_Features))
exp(coef(glm.nb(N_NN ~ N_Features)))
exp(confint(glm.nb(N_NN ~ N_Features)))

#Median LFC
summary(glm.nb(N ~ Med_LFC))
exp(coef(glm.nb(N ~ Med_LFC)))
exp(confint(glm.nb(N ~ Med_LFC)))

summary(glm.nb(N_RF ~ Med_LFC))
exp(coef(glm.nb(N_RF ~ Med_LFC)))
exp(confint(glm.nb(N_RF ~ Med_LFC)))

summary(glm.nb(N_NN ~ Med_LFC))
exp(coef(glm.nb(N_NN ~ Med_LFC)))
exp(confint(glm.nb(N_NN ~ Med_LFC)))

#Minimum LFC
summary(glm.nb(N ~ Min_LFC))
exp(coef(glm.nb(N ~ Min_LFC)))
exp(confint(glm.nb(N ~ Min_LFC)))

summary(glm.nb(N_RF ~ Min_LFC))
exp(coef(glm.nb(N_RF ~ Min_LFC)))
exp(confint(glm.nb(N_RF ~ Min_LFC)))

summary(glm.nb(N_NN ~ Min_LFC))
exp(coef(glm.nb(N_NN ~ Min_LFC)))
exp(confint(glm.nb(N_NN ~ Min_LFC)))

#Maximum LFC
summary(glm.nb(N ~ Max_LFC))
exp(coef(glm.nb(N ~ Max_LFC)))
exp(confint(glm.nb(N ~ Max_LFC)))

summary(glm.nb(N_RF ~ Max_LFC))
exp(coef(glm.nb(N_RF ~ Max_LFC)))
exp(confint(glm.nb(N_RF ~ Max_LFC)))

summary(glm.nb(N_NN ~ Max_LFC))
exp(coef(glm.nb(N_NN ~ Max_LFC)))
exp(confint(glm.nb(N_NN ~ Max_LFC)))

#Dispersion 
summary(glm.nb(N ~ I(Dispersion_Simulated>1)))
exp(coef(glm.nb(N ~ I(Dispersion_Simulated>1))))
exp(confint(glm.nb(N ~ I(Dispersion_Simulated>1))))

summary(glm.nb(N_RF ~ I(Dispersion_Simulated>1)))
exp(coef(glm.nb(N_RF ~ I(Dispersion_Simulated>1))))
exp(confint(glm.nb(N_RF ~ I(Dispersion_Simulated>1))))

summary(glm.nb(N_NN ~ I(Dispersion_Simulated>1)))
exp(coef(glm.nb(N_NN ~ I(Dispersion_Simulated>1))))
exp(confint(glm.nb(N_NN ~ I(Dispersion_Simulated>1))))

#Correlation
summary(glm.nb(N ~ Correlation))
exp(coef(glm.nb(N ~ Correlation)))
exp(confint(glm.nb(N ~ Correlation)))

summary(glm.nb(N_RF ~ Correlation))
exp(coef(glm.nb(N_RF ~ Correlation)))
exp(confint(glm.nb(N_RF ~ Correlation)))

summary(glm.nb(N_NN ~ Correlation))
exp(coef(glm.nb(N_NN ~ Correlation)))
exp(confint(glm.nb(N_NN ~ Correlation)))

#Dataset Nonlinearity
summary(glm.nb(N ~ Complexity))
exp(coef(glm.nb(N ~ Complexity)))
exp(confint(glm.nb(N ~ Complexity)))

summary(glm.nb(N_RF ~ Complexity))
exp(coef(glm.nb(N_RF ~ Complexity)))
exp(confint(glm.nb(N_RF ~ Complexity)))

summary(glm.nb(N_NN ~ Complexity))
exp(coef(glm.nb(N_NN ~ Complexity)))
exp(confint(glm.nb(N_NN ~ Complexity)))

#AvgMedRead
summary(glm.nb(N ~ AvgMedRead))
exp(coef(glm.nb(N ~ AvgMedRead)))
exp(confint(glm.nb(N ~ AvgMedRead)))

summary(glm.nb(N_RF ~ AvgMedRead))
exp(coef(glm.nb(N_RF ~ AvgMedRead)))
exp(confint(glm.nb(N_RF ~ AvgMedRead)))

summary(glm.nb(N_NN ~ AvgMedRead))
exp(coef(glm.nb(N_NN ~ AvgMedRead)))
exp(confint(glm.nb(N_NN ~ AvgMedRead)))

plot(log(N) ~ AvgMedRead)
plot(log(N_RF) ~ AvgMedRead)



plot(log(N) ~ Dispersion)
plot(log(N_RF) ~ Dispersion)

summary(glm(Complexity ~ Healthy, family=binomial(logit)))

Complexity <- ifelse(Complexity > 4.5, 1, 0)

Correlation <- ifelse(Correlation > 0.5, 1, 0)

AUC_Full <- ifelse(AUC_Full > 99, 1, 0)
AUC_Full_RF <- ifelse(AUC_Full_RF > 99, 1, 0)
AUC_Full_NN <- ifelse(AUC_Full_NN > 99, 1, 0)



Cancer <- c(1,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,1,0,1,1,1)

Healthy <- c(0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,1,0,1,1,0,0,0)

Disease <- c('Cancer', "Cancer", "Kidney", "Pulmonary", "Pulmonary", "Cancer", "Cardiovascular", "Liver",
             "Pulmonary", "Autoimmune", "Cardiovascular", "Gastro", "Psychiatric", "Cancer", "Neurological", "Autoimmune", "Cancer", "Cancer",
             "Pulmonary", "Pulmonary", "Autoimmune", "Cancer", "Psychiatric", "Cancer", "Autoimmune", "Cancer", "Cancer", "Cancer")

Grp1 <- ifelse(Disease %in% c("Kidney", "Liver", "Gastro"), 1, 0)
Grp2 <- ifelse(Disease %in% c("Autoimmune", "Neurological", "Psychiatric"), 1, 0)
Grp3 <- ifelse(Disease %in% c("Cardiovascular", "Pulmonary"), 1, 0)


Dataset_TEST <- data.frame("N" = N, "N_Features"=N_Features, "Max_LFC"=Max_LFC, "AUC_Full"=AUC_Full, "Min_LFC"=Min_LFC, "Med_LFC"=Med_LFC, "Imbalance"=Imbalance, 
                      "Correlation"=Correlation, "Complexity"=Complexity, "AvgMedRead"=AvgMedRead, "Dispersion"=Dispersion_Simulated)

Dataset_RF_TEST <- data.frame("N_RF" = N_RF, "N_Features"=N_Features,"Max_LFC"=Max_LFC,"AUC_Full_RF"=AUC_Full_RF, "Min_LFC"=Min_LFC, "Med_LFC"=Med_LFC, "Imbalance"=Imbalance,
                         "Correlation"=Correlation, "Complexity"=Complexity, "AvgMedRead"=AvgMedRead, "Dispersion"=Dispersion_Simulated)

Dataset_NN_TEST <- data.frame("N_NN" = N_NN, "N_Features"=N_Features,"Max_LFC"=Max_LFC,"AUC_Full_NN"=AUC_Full_NN, "Min_LFC"=Min_LFC, "Med_LFC"=Med_LFC, "Imbalance"=Imbalance,
                              "Correlation"=Correlation, "Complexity"=Complexity, "AvgMedRead"=AvgMedRead, "Dispersion"=Dispersion_Simulated)


#If needed, CV for other vars
TEST <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "AUC"=Dispersion_Simulated)

plot(seq(0.25, 1.5, by=0.05), sapply(seq(0.25, 1.5, by=0.05), function(i) {
  sqrt(cv.glm(glm.nb(N ~ I(AUC>i), data=TEST), data=TEST)$delta[1])}),xlab="Cutoff", ylab="LOO-CV MSE", main="XGB")

TEST <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "AUC"=Dispersion_Simulated)
plot(seq(0.25, 1.5, by=0.05), sapply(seq(0.25, 1.5, by=0.05), function(i) {
sqrt(cv.glm(glm.nb(N_RF ~ I(AUC>i), data=TEST), data=TEST)$delta[1])}),xlab="Cutoff", ylab="LOO-CV MSE", main="RF")

TEST <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "AUC"=Dispersion_Simulated)
plot(seq(0.25, 1.5, by=0.05), sapply(seq(0.25, 1.5, by=0.05), function(i) {
  sqrt(cv.glm(glm.nb(N_NN ~ I(AUC>i), data=TEST), data=TEST)$delta[1])}),xlab="Cutoff", ylab="LOO-CV MSE", main="NN")


#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
#write.csv(Dataset_TEST, "XGB_Vars.csv")
#write.csv(Dataset_RF_TEST, "RF_Vars.csv")
#write.csv(Dataset_NN_TEST, "NN_Vars.csv")


FindModel_XGB <- function(DS_Name) {
  Results <- NULL
  
  Full <- DS_Name
  
  Full$N_Features <- ifelse(Full$N_Features>10, 1, 0)
  
  Full$Dispersion <- ifelse(Full$Dispersion>1, 1, 0)
  
  for (i in 1:5000) {
  Variables <- sample(colnames(Full)[-1], 4)
  
  Reduced <- Full[,c(Variables, "N")]
  
  Model <- glm.nb(data=Reduced, N ~ .)
  
  AIC <- round(Model$aic, 2)
  
  R2 <- adjR2(Model)
  
  Results <- rbind(Results, c(AIC, R2, c(Variables))
  )
  
  }
  
  Results
}

FindModel_RF <- function(DS_Name) {
  Results <- NULL
  
  Full <- DS_Name
  
  Full$N_Features <- ifelse(Full$N_Features>10, 1, 0)
  
  Full$Dispersion <- ifelse(Full$Dispersion>1, 1, 0)
  
  for (i in 1:5000) {
    Variables <- sample(colnames(Full)[-1], 4)
    
    Reduced <- Full[,c(Variables, "N_RF")]
    
    Model <- glm.nb(data=Reduced, N_RF ~ .)
    
    AIC <- round(Model$aic, 2)
    
    R2 <- adjR2(Model)
    
    Results <- rbind(Results, c(AIC, R2, c(Variables)))
                     
  }
  
  Results
}

FindModel_NN <- function(DS_Name) {
  Results <- NULL
  
  Full <- DS_Name
  
  Full$N_Features <- ifelse(Full$N_Features>10, 1, 0)
  
  Full$Dispersion <- ifelse(Full$Dispersion>1, 1, 0)
  
  for (i in 1:5000) {
    Variables <- sample(colnames(Full)[-1], 4)
    
    Reduced <- Full[,c(Variables, "N_NN")]
    
    Model <- glm.nb(data=Reduced, N_NN ~ .)
    
    AIC <- round(Model$aic, 2)
    
    R2 <- adjR2(Model)
    
    Results <- rbind(Results, c(AIC, R2, c(Variables))
    )
    
  }
  
  Results
}

FindModel_XGB2 <- function(DS_Name) {
  Results <- NULL
  
  IBD_DF <- data.frame("N"=836, "AUC_Full"=0, "Complexity"=0, "AUC25"=0.737, "Imbalance"=18.1,
                       "Dispersion"=1.07, "AvgMedRead"=180.74, "Correlation"=0, "N_Features"=57, "Med_LFC"=0.90,
                       "Min_LFC"=0.17, "Max_LFC"=2.62)
  
  PDAC_DF <- data.frame("N"=176, "AUC_Full"=1, "Complexity"=0, "AUC25"=0.785, "Imbalance"=30.1,
                        "Dispersion"=0.89, "AvgMedRead"=5465.65, "Correlation"=0, "N_Features"=67, "Med_LFC"=0.87,
                        "Min_LFC"=0.01, "Max_LFC"=5.58)
  
  Full <- DS_Name
  
  Full$N_Features <- ifelse(Full$N_Features>10, 1, 0)
  
  Full$Dispersion <- ifelse(Full$Dispersion>1, 1, 0)
  
  for (i in 1:5000) {
    Variables <- sample(colnames(Full)[-1], 4)
    
    Reduced <- Full[,c(Variables, "N")]
    
    Model <- glm.nb(data=Reduced, N ~ .)
    
    IBD <- abs(836 - exp(predict(Model, newdata=IBD_DF[,Variables])))
    
    PDAC <- abs(176 - exp(predict(Model, newdata=PDAC_DF[,Variables])))
    
    Results <- rbind(Results, c(c(Variables), "IBD_error"=IBD, "PDAC_Error"=PDAC))
    
  }
  Results
}

FindModel_RF2 <- function(DS_Name) {
  Results <- NULL
  
  IBD_DF <- data.frame("N_RF"=367, "AUC_Full_RF"=0, "Complexity"=0, "AUC25_RF"=0.797, "Imbalance"=18.1,
                       "Dispersion"=1.07, "AvgMedRead"=180.74, "Correlation"=0, "N_Features"=57, "Med_LFC"=0.90,
                       "Min_LFC"=0.17, "Max_LFC"=2.62)
  
  PDAC_DF <- data.frame("N_RF"=80, "AUC_Full_RF"=1, "Complexity"=0, "AUC25_RF"=0.910, "Imbalance"=30.1,
                        "Dispersion"=0.89, "AvgMedRead"=5465.65, "Correlation"=0, "N_Features"=67, "Med_LFC"=0.87,
                        "Min_LFC"=0.01, "Max_LFC"=5.58)
  
  Full <- DS_Name
  
  #Full$AUC_LOG <- log(Full$AUC_Full) 
  #Full$N_Features_LOG <- log(Full$N_Features)
  #Full$Imbalance_LOG <- log(Full$Imbalance)
  #Full$Max_LFC_LOG <- log(Full$Max_LFC)
  #Full$Min_LFC_LOG <- log(Full$Min_LFC)
  #Full$Med_LFC_LOG <- log(Full$Med_LFC)
  #Full$AvgMedRead_LOG <- log(Full$AvgMedRead)
  #Full$Dispersion_LOG <- log(Full$Dispersion)
  
  Full$N_Features <- ifelse(Full$N_Features>10, 1, 0)
  
  Full$Dispersion <- ifelse(Full$Dispersion>1, 1, 0)
  
  for (i in 1:5000) {
    Variables <- sample(colnames(Full)[-1], 4)
    
    Reduced <- Full[,c(Variables, "N_RF")]
    
    Model <- glm.nb(data=Reduced, N_RF ~ .)
    
    IBD <- abs(367 - exp(predict(Model, newdata=IBD_DF[,Variables])))
    
    PDAC <- abs(80 - exp(predict(Model, newdata=PDAC_DF[,Variables])))
    
    Results <- rbind(Results, c(c(Variables), "IBD_error"=IBD, "PDAC_Error"=PDAC))
    
  }
  Results
}


FindModel_NN2 <- function(DS_Name) {
  Results <- NULL
  
  IBD_DF <- data.frame("N_NN"=905, "AUC_Full_NN"=0, "Complexity"=0, "Imbalance"=18.1,
                       "Dispersion"=1.07, "AvgMedRead"=180.74, "Correlation"=0, "N_Features"=57, "Med_LFC"=0.90,
                       "Min_LFC"=0.17, "Max_LFC"=2.62)
  
  PDAC_DF <- data.frame("N_NN"=113, "AUC_Full_NN"=1, "Complexity"=0, "Imbalance"=30.1,
                        "Dispersion"=0.89, "AvgMedRead"=5465.65, "Correlation"=0, "N_Features"=67, "Med_LFC"=0.87,
                        "Min_LFC"=0.01, "Max_LFC"=5.58)
  
  Full <- DS_Name
  
  Full$N_Features <- ifelse(Full$N_Features>10, 1, 0)
  
  Full$Dispersion <- ifelse(Full$Dispersion>1, 1, 0)
  
  for (i in 1:5000) {
    Variables <- sample(colnames(Full)[-1], 4)
    
    Reduced <- Full[,c(Variables, "N_NN")]
    
    Model <- glm.nb(data=Reduced, N_NN ~ .)
    
    IBD <- abs(905 - exp(predict(Model, newdata=IBD_DF[,Variables])))
    
    PDAC <- abs(113 - exp(predict(Model, newdata=PDAC_DF[,Variables])))
    
    Results <- rbind(Results, c(c(Variables), "IBD_error"=IBD, "PDAC_Error"=PDAC))
    
  }
  Results
}


XGB_Search <- FindModel_XGB2(Dataset_TEST[-c(24),])
XGB_Search <- as.data.frame(XGB_Search)
XGB_Search[order(as.numeric(XGB_Search[,5]), decreasing = F),]

RF_Search <- FindModel_RF2(Dataset_RF_TEST[-c(24),])
RF_Search <- as.data.frame(RF_Search)
View(RF_Search[order(as.numeric(RF_Search[,6]), decreasing = F),])

RF_Search$Flag <- ifelse(RF_Search[,1] == "AUC25_RF" | RF_Search[,2] == "AUC25_RF" | RF_Search[,3] == "AUC25_RF" | RF_Search[,4] == "AUC25_RF", 1, 0)

NN_Search <- FindModel_NN2(Dataset_NN_TEST[-c(24),])
NN_Search <- as.data.frame(NN_Search)
View(NN_Search[order(as.numeric(NN_Search[,6]), decreasing = F),])

#-----------------------------


XGB_Search <- FindModel_XGB(Dataset_TEST[-c(24, 28),])
XGB_Search <- as.data.frame(XGB_Search)
XGB_Search[order(as.numeric(XGB_Search[,1]), decreasing = F),]

XGB_Search$Flag <- ifelse(XGB_Search[,3] == "Complexity" | XGB_Search[,4] == "Complexity" | XGB_Search[,5] == "Complexity" | XGB_Search[,6] == "Complexity", 1, 0)
View(XGB_Search[order(as.numeric(XGB_Search[,2]), decreasing = T),])

RF_Search <- FindModel_RF(Dataset_RF_TEST[-24])
RF_Search <- as.data.frame(RF_Search)
RF_Search[order(as.numeric(RF_Search[,1]), decreasing = F),]

RF_Search$Flag <- ifelse(RF_Search[,3] == "Complexity" | RF_Search[,4] == "Complexity" | RF_Search[,5] == "Complexity" | RF_Search[,6] == "Complexity", 1, 0)
View(RF_Search[order(as.numeric(RF_Search[,2]), decreasing = T),])


NN_Search <- FindModel_NN(Dataset_NN_TEST[-24,])
NN_Search <- as.data.frame(NN_Search)
NN_Search[order(as.numeric(NN_Search[,1]), decreasing = F),]

NN_Search$Flag <- ifelse(NN_Search[,3] == "Complexity" | NN_Search[,4] == "Complexity" | NN_Search[,5] == "Complexity" | NN_Search[,6] == "Complexity", 1, 0)
View(NN_Search[order(as.numeric(NN_Search[,2]), decreasing = T),])

#MODEL XGB
Mod_XGB <- glm.nb(data=Dataset_TEST, N ~ I(Dispersion>1) + AUC_Full + Min_LFC + Med_LFC)
summary(Mod_XGB)
adjR2(Mod_XGB)

#IBD
exp(predict(Mod_XGB, newdata=data.frame("Dispersion"=1.07, "Min_LFC"=0.17, "Med_LFC"=0.90, "AUC_Full"=0))) #809

DF_new <- Dataset_TEST[, c("Dispersion", "Min_LFC","Med_LFC","AUC_Full","N")]

DF_new <- rbind(DF_new, data.frame("Dispersion"=1.07, "Min_LFC"=0.17, "Med_LFC"=0.90, "AUC_Full"=0, "N"=exp(predict(Mod_XGB, newdata=data.frame("Dispersion"=1.07, "Min_LFC"=0.17, "Med_LFC"=0.90, "AUC_Full"=0)))))

CI <- add_ci(DF_new, Mod_XGB)
CI[nrow(CI),c(6,7,8)] 

#PDAC
exp(predict(Mod_XGB, newdata=data.frame("Dispersion"=0.89, "Min_LFC"=0.01, "Med_LFC"=0.87, "AUC_Full"=1))) #136

DF_new <- Dataset_TEST[, c("Dispersion", "Min_LFC","Med_LFC","AUC_Full","N")]

DF_new <- rbind(DF_new, data.frame("Dispersion"=0.89, "Min_LFC"=0.01, "Med_LFC"=0.87, "AUC_Full"=1, "N"=exp(predict(Mod_XGB, newdata=data.frame("Dispersion"=0.89, "Min_LFC"=0.01, "Med_LFC"=0.87, "AUC_Full"=1)))))

CI <- add_ci(DF_new, Mod_XGB)
CI[nrow(CI),c(6,7,8)] 



Mod_RF <- glm.nb(data=Dataset_RF_TEST, N_RF ~ Max_LFC + I(Dispersion>1) + Imbalance + Complexity)
summary(Mod_RF)
adjR2(Mod_RF)

#IBD
exp(predict(Mod_RF, newdata=data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "Imbalance"=18.1, "Complexity"=0))) #302

DF_new <- Dataset_RF_TEST[, c("Dispersion", "Max_LFC","Imbalance","Complexity","N_RF")]

DF_new <- rbind(DF_new, data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "Imbalance"=18.1, "Complexity"=0, "N_RF"=exp(predict(Mod_RF, newdata=data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "Imbalance"=18.1, "Complexity"=0)))))

CI <- add_ci(DF_new, Mod_RF)
CI[nrow(CI),c(6,7,8)]

#PDAC
exp(predict(Mod_RF, newdata=data.frame("Dispersion"=0.89, "Max_LFC"=5.58, "Imbalance"=30.1, "Complexity"=0))) #63

DF_new <- Dataset_RF_TEST[, c("Dispersion", "Max_LFC","Imbalance","Complexity","N_RF")]

DF_new <- rbind(DF_new, data.frame("Dispersion"=0.89, "Max_LFC"=5.58, "Imbalance"=30.1, "Complexity"=0, "N_RF"=exp(predict(Mod_RF, newdata=data.frame("Dispersion"=0.89, "Max_LFC"=5.58, "Imbalance"=30.1, "Complexity"=0)))))

CI <- add_ci(DF_new, Mod_RF)
CI[nrow(CI),c(6,7,8)]


#NN NEXT!
Mod_NN <- glm.nb(data=Dataset_NN_TEST, N_NN ~ Complexity + AUC_Full_NN + Max_LFC + I(Dispersion>1))
summary(Mod_NN)
adjR2(Mod_NN)

#IBD
exp(predict(Mod_NN, newdata=data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "AUC_Full_NN"=0, "Complexity"=0))) #780

DF_new <- Dataset_NN_TEST[-c(24), c("Dispersion", "Max_LFC","AUC_Full_NN","Complexity","N_NN")]

DF_new <- rbind(DF_new, data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "AUC_Full_NN"=0, "Complexity"=0, "N_NN"=exp(predict(Mod_NN, newdata=data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "AUC_Full_NN"=0, "Complexity"=0)))))

CI <- add_ci(DF_new, Mod_NN)
CI[nrow(CI),c(6,7,8)]

#PDAC
exp(predict(Mod_NN, newdata=data.frame("Dispersion"=0.89, "Max_LFC"=5.58, "AUC_Full_NN"=1, "Complexity"=0))) #71

DF_new <- Dataset_NN_TEST[, c("Dispersion", "Max_LFC","AUC_Full_NN","Complexity","N_NN")]

DF_new <- rbind(DF_new, data.frame("Dispersion"=0.89, "Max_LFC"=5.58, "AUC_Full_NN"=1, "Complexity"=0, "N_NN"=exp(predict(Mod_NN, newdata=data.frame("Dispersion"=0.89, "Max_LFC"=5.58, "AUC_Full_NN"=1, "Complexity"=0)))))

CI <- add_ci(DF_new, Mod_NN)
CI[nrow(CI),c(6,7,8)]

#Cancer Subgroup
quantile(Dataset_TEST[which(Cancer==1), "N"])
quantile(Dataset_TEST[which(Cancer==0), "N"])
wilcox.test(Dataset_TEST$N ~ Cancer)

quantile(Dataset_RF_TEST[which(Cancer==1), "N_RF"])
quantile(Dataset_RF_TEST[which(Cancer==0), "N_RF"])
wilcox.test(Dataset_RF_TEST$N_RF ~ Cancer)

quantile(Dataset_NN_TEST[which(Cancer==1), "N_NN"])
quantile(Dataset_NN_TEST[which(Cancer==0), "N_NN"])
wilcox.test(Dataset_NN_TEST$N_NN ~ Cancer)

summary(glm.nb(Dataset_TEST$N[-24] ~ Cancer[-24]*(Dataset_TEST$Med_LFC[-24] + Dataset_TEST$AUC_Full[-24] + Dataset_TEST$Min_LFC[-24] + I(Dataset_TEST$Dispersion>1)[-24])))
summary(glm.nb(Dataset_RF_TEST$N_RF ~ Cancer*(Dataset_RF_TEST$Complexity + Dataset_RF_TEST$Max_LFC + Dataset_RF_TEST$Imbalance + I(Dataset_TEST$Dispersion>1))))
summary(glm.nb(Dataset_NN_TEST$N_NN[-24] ~ Cancer[-24]*(Dataset_NN_TEST$Complexity[-24] + Dataset_NN_TEST$AUC_Full_NN[-24] + I(Dataset_NN_TEST$Dispersion>1)[-24] + Dataset_NN_TEST$Max_LFC[-24])))

summary(glm.nb(Dataset_TEST$N ~ Grp1*(Dataset_TEST$Complexity + Dataset_TEST$AUC_Full + log(Dataset_TEST$Max_LFC))))
summary(glm.nb(Dataset_RF_TEST$N_RF ~ Grp1*(Dataset_RF_TEST$Complexity + Dataset_RF_TEST$Max_LFC + Dataset_RF_TEST$Imbalance)))
summary(glm.nb(Dataset_NN_TEST$N_NN ~ Grp1*(Dataset_NN_TEST$Complexity + Dataset_NN_TEST$AUC_Full_NN + Dataset_NN_TEST$Imbalance)))

summary(glm.nb(Dataset_TEST$N ~ Grp2*(Dataset_TEST$Complexity + Dataset_TEST$AUC_Full + log(Dataset_TEST$Max_LFC))))
summary(glm.nb(Dataset_RF_TEST$N_RF ~ Grp2*(Dataset_RF_TEST$Complexity + Dataset_RF_TEST$Max_LFC + Dataset_RF_TEST$Imbalance)))
summary(glm.nb(Dataset_NN_TEST$N_NN ~ Grp2*(Dataset_NN_TEST$Complexity + Dataset_NN_TEST$AUC_Full_NN + Dataset_NN_TEST$Imbalance)))

summary(glm.nb(Dataset_TEST$N ~ Grp3*(Dataset_TEST$Complexity + Dataset_TEST$AUC_Full + log(Dataset_TEST$Max_LFC))))
summary(glm.nb(Dataset_RF_TEST$N_RF ~ Grp3*(Dataset_RF_TEST$Complexity + Dataset_RF_TEST$Max_LFC + Dataset_RF_TEST$Imbalance)))
summary(glm.nb(Dataset_NN_TEST$N_NN ~ Grp3*(Dataset_NN_TEST$Complexity + Dataset_NN_TEST$AUC_Full_NN + Dataset_NN_TEST$Imbalance)))

Dataset_RF_TEST_CANCER <- Dataset_RF_TEST
Dataset_RF_TEST_CANCER$Cancer <- Cancer

Mod_RF_CANCER <- glm.nb(data=Dataset_RF_TEST_CANCER[-24,], N_RF ~ Cancer*(Complexity + Max_LFC + Imbalance + I(Dispersion>1)))
adjR2(Mod_RF_CANCER)
summary(Mod_RF_CANCER)

exp(predict(Mod_RF_CANCER, newdata=data.frame("Dispersion"=1.07, "Max_LFC"=2.62, "Imbalance"=18.1, "Complexity"=0, "Cancer"=0))) #301

summary(Mod_XGB)
summary(Mod_RF)
summary(Mod_NN)


#--------------------------------------------------------------------


adjR2(Mod_RF)
adjR2(Mod_NN)

sqrt(cv.glm(Mod_XGB, data=Dataset_TEST)$delta[1])

sqrt(cv.glm(Mod_RF, data=Dataset_RF_TEST)$delta[1])

sqrt(cv.glm(Mod_NN, data=Dataset_NN_TEST)$delta[1])





#PLOTS
DS_Names <- c("GBM", "HCC", "Kidney", "IPF", "COVID", "HBV",
              "Hypertension", "NAFLD", "PrePost", "RA", "AVSC", "Crohn",
              "MDD", "Ovarian", "ALS", "Kawasaki", "CCA", "Glioma", "TBProg",
              "Tuberculosis", "EDS", "NSCLC",
              "Bipolar","MS")


Original_N <- c(32,63,82,84,84,95,97,98,99,100,101,125,139,148,162,165,170,176,177,182,200,218,239,474)

Spread <- NULL
for (i in 1:24) {
  Spread <- append(Spread, max(N[i], N_RF[i], N_NN[i])-min(N[i], N_RF[i], N_NN[i]))
}

Spread_andN <- data.frame(DS_Names, Spread)
Labels <- c(Spread_andN[order(Spread_andN$Spread),"DS_Names"])
New_Labels <- c("4"="GBM",
            "10"="Kawasaki",
            "22"="NAFLD",
            "56"="IPF",
            "82"="EDS",
            "106"="PrePost",
            "118"="Tuberculosis",
            "156"="Hypertension",
            "182"="NSCLC",
            "196"="HCC",
            "417"="RA",
            "423"="HBV",
            "487"="CCA",
            "490"="Crohn",
            "498"="Ovarian",
            "511"="MS",
            "540"="TBProg",
            "550"="AVSC",
            "567"="ALS",
            "661"="Kidney",
            "694"="Glioma",
            "794"="MDD",
            "827"="Bipolar",
            "974"="COVID"
)
            

#PLOTS
All_Ns <- data.frame("N"=c(N[1:24], N_RF[1:24], N_NN[1:24]), "Dataset"=as.factor(rep(DS_Names, 3)),
                     "Algorithm"=as.factor(c(rep("XGB", 24), rep("RF", 24), rep("NN", 24))),
                     "Original_N"=rep(Original_N, 3), "Order"=factor(rep(Spread, 3)))


A <- ggplot(data=All_Ns, aes(y=N, colour=Algorithm, x="")) + geom_point(height=0, width=0.1, size=2) + facet_wrap(All_Ns$Order, nrow=2, ncol=12, labeller = as_labeller(New_Labels)) + theme_bw() + 
  labs(y="Required sample size") 


A

Spread2 <- NULL
for (i in 25:27) {
  Spread2 <- append(Spread2, max(N[i], N_RF[i], N_NN[i])-min(N[i], N_RF[i], N_NN[i]))
}

Spread_andN2 <- data.frame("DS_Names"=c("TCGA-Lung", "TCGA-Breast", "TCGA_Brain"), Spread2)
Labels <- c(Spread_andN2[order(Spread_andN2$Spread2),"DS_Names"])
New_Labels2 <- c("43"="TCGA-Lung",
                "61"="TCGA-Brain",
                "309"="TCGA-Breast"
)




NonSimulated <- data.frame("N"=c(N[25:27], N_RF[25:27], N_NN[25:27]), "Dataset"=as.factor(rep(c("TCGA-Lung", "TCGA-Breast", "TCGA_Brain"), 3)),
                     "Algorithm"=as.factor(c(rep("XGB", 3), rep("RF", 3), rep("NN", 3))), "Order"=factor(rep(Spread2, 3)),
                     "Original_N"=rep(c(1162,925,875), 3))

B <- ggplot(data=NonSimulated, aes(y=N, colour=Algorithm, x="")) + geom_point(height=0, width=0.1, size=2) + facet_wrap(NonSimulated$Order, nrow=2, ncol=12, labeller = as_labeller(New_Labels2)) + theme_bw() + 
  labs(x="",y="Required sample size") 

B

ggarrange(A, B, common.legend = T, ncol=1, nrow=2)


DS_Names <- c("GBM", "HCC", "Kidney", "IPF", "COVID", "HBV",
              "Hypertension", "NAFLD", "PrePost", "RA", "AVSC", "Crohn",
              "MDD", "Ovarian", "ALS", "Kawasaki", "CCA", "Glioma", "TBProg",
              "Tuberculosis", "EDS", "NSCLC",
              "Bipolar","MS", "TCGA-Lung", "TCGA-Breast", "TCGA-Brain")
 #AUC_DF <- data.frame("AUC"=Dataset_TEST$AUC_Full,
#                     "N"=Dataset_TEST$N)





#Figure comparing original DS AUCs vs Simulated AUCs-----------------------------------------------------------

AUC_Full <- c(1, 0.995,0.848,1.000,0.936,0.883,0.999,0.999,0.998,0.981,0.852,0.903,0.862,0.845,0.949,0.999,0.924,0.964,0.893,0.995,1,0.995,0.933,0.883,1,0.718,0.985)
AUC_Full_RF <- c(1, 0.994,0.850,1.000,0.936,0.884,0.999,0.999,0.996,0.977,0.849,0.896,0.857,0.846,0.943,0.998,0.922,0.963,0.895,0.992,1,0.993,0.930,0.883,1,0.75,0.983)
AUC_Full_NN <- c(1.000,0.994,0.862,1.000,0.947,0.896,0.999,1.000,0.997,0.988,0.858,0.913,0.874,0.861,0.964,1.000,0.935,0.975,0.922,0.998,1.000,0.996,0.951,0.907,1,0.724,0.981)


Original_DS_AUC_XGB <- c(0.942, 0.936, 0.843, 0.959, 0.806, 0.811, 0.983, 0.982, 0.984, 0.936, 0.718, 0.809, 
                         0.827, 0.763, 0.746, 0.994, 0.836, 0.926, 0.870, 0.971, 0.978, 0.971, 0.864, 0.839)

Original_DS_AUC_XGB_se <- c(0.058, 0.048, 0.043, 0.028, 0.061, 0.039, 0.010, 0.011, 0.014, 0.028, 0.073, 0.046,
                        0.014, 0.055, 0.055, 0.004, 0.004, 0.017, 0.020, 0.029, 0.008, 0.008, 0.030, 0.024)

Original_DS_AUC_RF <- c(0.983, 0.989, 0.830, 1, 0.827, 0.835, 0.995, 0.977, 0.986, 0.957, 0.766,
                        0.830, 0.839, 0.770, 0.766, 0.996, 0.857, 0.940, 0.849, 0.962, 0.993, 0.965,
                        0.884, 0.837)

Original_DS_AUC_RF_se <- c(0.017, 0.011, 0.029, 0, 0.054, 0.033, 0.005, 0.014, 0.010, 0.021, 0.030,
                           0.053, 0.027, 0.049, 0.065, 0.003, 0.019, 0.016, 0.030, 0.013, 0.006,
                           0.008, 0.019, 0.029)


Original_DS_AUC_NN <- c(0.981, 0.949, 0.800, 1, 0.834, 0.892, 0.962, 0.981, 0.958, 0.924, 0.810, 0.891,
                        0.844, 0.819, 0.865, 0.990, 0.878, 0.935, 0.810, 0.981, 0.994, 0.961, 0.885, 0.810)

Original_DS_AUC_NN_se <- c(0.017, 0.016, 0.059, 0.000, 0.032, 0.022, 0.016, 0.013, 0.014, 0.008, 0.021, 0.033,
                           0.017, 0.025, 0.017, 0.010, 0.028, 0.011, 0.046, 0.009, 0.009, 0.014, 0.005, 0.036)


DF_XGB_AUC_figure <- data.frame("AUC"=c(AUC_Full[1:24], Original_DS_AUC_XGB),
           "se"=c(rep(0.0015, 24), Original_DS_AUC_XGB_se),
           "Dataset"=rep(DS_Names, 2),
           "Group"=c(rep("Simulated", 24), rep("Original", 24)))

DF_XGB_AUC_figure$Lower <- DF_XGB_AUC_figure$AUC - 1.96*DF_XGB_AUC_figure$se
DF_XGB_AUC_figure$Upper <- ifelse(DF_XGB_AUC_figure$AUC + 1.96*DF_XGB_AUC_figure$se > 1, 1, 
                                  DF_XGB_AUC_figure$AUC + 1.96*DF_XGB_AUC_figure$se)

DF_RF_AUC_figure <- data.frame("AUC"=c(AUC_Full_RF[1:24], Original_DS_AUC_RF),
                                "se"=c(rep(0.0015, 24), Original_DS_AUC_RF_se),
                                "Dataset"=rep(DS_Names, 2),
                                "Group"=c(rep("Simulated", 24), rep("Original", 24)))

DF_RF_AUC_figure$Lower <- DF_RF_AUC_figure$AUC - 1.96*DF_RF_AUC_figure$se
DF_RF_AUC_figure$Upper <- ifelse(DF_RF_AUC_figure$AUC + 1.96*DF_RF_AUC_figure$se > 1, 1, 
                                  DF_RF_AUC_figure$AUC + 1.96*DF_RF_AUC_figure$se)

DF_NN_AUC_figure <- data.frame("AUC"=c(AUC_Full_NN[1:24], Original_DS_AUC_NN),
                               "se"=c(rep(0.0015, 24), Original_DS_AUC_NN_se),
                               "Dataset"=rep(DS_Names, 2),
                               "Group"=c(rep("Simulated", 24), rep("Original", 24)))

DF_NN_AUC_figure$Lower <- DF_NN_AUC_figure$AUC - 1.96*DF_NN_AUC_figure$se
DF_NN_AUC_figure$Upper <- ifelse(DF_NN_AUC_figure$AUC + 1.96*DF_NN_AUC_figure$se > 1, 1, 
                                 DF_NN_AUC_figure$AUC + 1.96*DF_NN_AUC_figure$se)



A <- ggplot(data=DF_XGB_AUC_figure, aes(y=AUC, x=factor(Group), color=Dataset)) + geom_point(size=2, position = position_dodge(width=0.1)) + geom_line(aes(group = Dataset), position = position_dodge(width = 0.1)) + labs(x="", title = "XGBoost") +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 1, alpha=0.8, linetype="dashed",
    position = position_dodge(width = 0.1)
  ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

B <- ggplot(data=DF_RF_AUC_figure, aes(y=AUC, x=factor(Group), color=Dataset)) + geom_point(size=2, position = position_dodge(width=0.1)) + geom_line(aes(group = Dataset), position = position_dodge(width = 0.1)) + labs(x="", y="", title = "Random Forest") +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 1, alpha=0.8, linetype="dashed",
    position = position_dodge(width = 0.1)
  ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

C <- ggplot(data=DF_NN_AUC_figure, aes(y=AUC, x=factor(Group), color=Dataset)) + geom_point(size=2, position = position_dodge(width=0.1)) + geom_line(aes(group = Dataset), position = position_dodge(width = 0.1)) + labs(x="", y="", title = "Neural Network") +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 1, alpha=0.8, linetype="dashed",
    position = position_dodge(width = 0.1)
  ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggarrange(A,B,C, ncol=3,nrow=1, common.legend = T)

#----------------------------------------------------------------------------------------------------------------------








Imbalance_DF <- data.frame("Imbalance"=c(Dataset_TEST$Imbalance, Dataset_RF_TEST$Imbalance, Dataset_NN_TEST$Imbalance),
                           "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                           "N"=c(N,N_RF,N_NN))

Nonlinearity_DF <- data.frame("Nonlinearity"=c(Complexity, Complexity, Complexity),
                              "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                              "N"=c(N,N_RF,N_NN))

N_Features_DF <- data.frame("N_Features"=c(Dataset_TEST$N_Features, Dataset_RF_TEST$N_Features, Dataset_NN_TEST$N_Features),
                            "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                            "N"=c(N,N_RF,N_NN))

Max_LFC_DF <- data.frame("MaxLFC"=c(Dataset_TEST$Max_LFC, Dataset_RF_TEST$Max_LFC, Dataset_NN_TEST$Max_LFC),
                         "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                         "N"=c(N,N_RF,N_NN))

Med_LFC_DF <- data.frame("MedLFC"=c(Dataset_TEST$Med_LFC, Dataset_RF_TEST$Med_LFC, Dataset_NN_TEST$Med_LFC),
                         "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                         "N"=c(N,N_RF,N_NN))

Min_LFC_DF <- data.frame("MinLFC"=c(Dataset_TEST$Min_LFC, Dataset_RF_TEST$Min_LFC, Dataset_NN_TEST$Min_LFC),
                         "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                         "N"=c(N,N_RF,N_NN))

Dispersion_DF <- data.frame("Dispersion"=c(ifelse(Dataset_TEST$Dispersion>1,1,0), ifelse(Dataset_RF_TEST$Dispersion>1,1,0), ifelse(Dataset_NN_TEST$Dispersion>1,1,0)),
                         "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                         "N"=c(N,N_RF,N_NN))

Correlation_DF <- data.frame("Correlation"=c(Dataset_TEST$Correlation, Dataset_RF_TEST$Correlation, Dataset_NN_TEST$Correlation),
                            "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                            "N"=c(N,N_RF,N_NN))

AvgMedRead_DF <- data.frame("AvgMedRead"=c(Dataset_TEST$AvgMedRead, Dataset_RF_TEST$AvgMedRead, Dataset_NN_TEST$AvgMedRead),
                             "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                             "N"=c(N,N_RF,N_NN))



#Do all of them...

AUC_DF <- data.frame("AUC"=c(ifelse(AUC_Full==1, "Full-Dataset AUC >= 0.99", "Full-Dataset AUC < 0.99"), ifelse(AUC_Full_RF==1, "Full-Dataset AUC >= 0.99", "Full-Dataset AUC < 0.99"), ifelse(AUC_Full_NN==1, "Full-Dataset AUC >= 0.99", "Full-Dataset AUC < 0.99")),
                     "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                     "N"=c(N,N_RF,N_NN))

A <- ggplot(data=AUC_DF, aes(y=log(N), colour=Algorithm)) + facet_wrap(facets=as.factor(AUC_DF$AUC)) + geom_boxplot(size=1) + theme_bw() + theme(strip.text = element_text(size=15), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x = "", y="Natural-Log Required Sample Size") 
A

B <- ggplot(data=Imbalance_DF, aes(x = Imbalance, y=log(N), colour=Algorithm)) + geom_point(size=2) + geom_smooth(method="lm", se=F) + theme_bw() + labs(x = "Minority class proportion", y="Natural-Log Required Sample Size")
B

Nonlinearity_DF <- data.frame("Nonlinearity"=c(ifelse(Dataset_TEST$Complexity>0,"Dataset Nonlinearity >= 4.5", "Dataset Nonlinearity < 4.5"), ifelse(Dataset_RF_TEST$Complexity>0,"Dataset Nonlinearity >= 4.5","Dataset Nonlinearity < 4.5"), ifelse(Dataset_NN_TEST$Complexity>0,"Dataset Nonlinearity >= 4.5", "Dataset Nonlinearity < 4.5")),
                            "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                            "N"=c(N,N_RF,N_NN))
C <- ggplot(data=Nonlinearity_DF, aes(y=log(N), colour=Algorithm)) + facet_wrap(facets=as.factor(Nonlinearity_DF$Nonlinearity)) + geom_boxplot(size=1) + theme_bw() + theme(strip.text = element_text(size=15), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x = "", y="") 
C


D <- ggplot(data=N_Features_DF, aes(x = N_Features, y=log(N), colour=Algorithm)) + geom_point(size=2) + geom_smooth(method="lm", se=F) + theme_bw() + labs(x = "Number of features", y="Natural-Log Required Sample Size")
D

E <- ggplot(data=Max_LFC_DF, aes(x = MaxLFC, y=log(N), colour=Algorithm)) + geom_point(size=2) + geom_smooth(method="lm", se=F) + theme_bw() + labs(x = "Max LFC", y="Natural-Log Required Sample Size")
E

G <- ggplot(data=Med_LFC_DF, aes(x = MedLFC, y=log(N), colour=Algorithm)) + geom_point(size=2) + geom_smooth(method="lm", se=F) + theme_bw() + labs(x = "Median LFC", y="")
G

H <- ggplot(data=Min_LFC_DF, aes(x = MinLFC, y=log(N), colour=Algorithm)) + geom_point(size=2) + geom_smooth(method="lm", se=F) + theme_bw() + labs(x = "Min LFC", y="")
H

Dispersion_DF <- data.frame("Dispersion"=c(ifelse(Dataset_TEST$Dispersion>1,"Median Dispersion >= 1", "Median Dispersion < 1"), ifelse(Dataset_RF_TEST$Dispersion>1,"Median Dispersion >= 1","Median Dispersion < 1"), ifelse(Dataset_NN_TEST$Dispersion>1,"Median Dispersion >= 1", "Median Dispersion < 1")),
                     "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                     "N"=c(N,N_RF,N_NN))

I <- ggplot(data=Dispersion_DF, aes(y=log(N), colour=Algorithm)) + facet_wrap(facets=as.factor(Dispersion_DF$Dispersion)) + geom_boxplot(size=1) + theme_bw() + theme(strip.text = element_text(size=15), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x = "", y="Natural-Log Required Sample Size") 
I

J <- ggplot(data=AvgMedRead_DF, aes(x = AvgMedRead, y=log(N), colour=Algorithm)) + geom_point(size=2) + geom_smooth(method="lm", se=F) + theme_bw() + labs(x = "Median Avg. Read Count", y="")
J

Correlation_DF <- data.frame("Correlation"=c(ifelse(Dataset_TEST$Correlation>0,"Avg. Between-Gene Correlation >= 0.5", "Avg. Between-Gene Correlation < 0.5"), ifelse(Dataset_RF_TEST$Correlation>0,"Avg. Between-Gene Correlation >= 0.5","Avg. Between-Gene Correlation < 0.5"), ifelse(Dataset_NN_TEST$Correlation>0,"Avg. Between-Gene Correlation >= 0.5", "Avg. Between-Gene Correlation < 0.5")),
                              "Algorithm"=as.factor(c(rep("XGB", 27), rep("RF", 27), rep("NN", 27))),
                              "N"=c(N,N_RF,N_NN))
K <- ggplot(data=Correlation_DF, aes(y=log(N), colour=Algorithm)) + facet_wrap(facets=as.factor(Correlation_DF$Correlation)) + geom_boxplot(size=1) + theme_bw() + theme(strip.text = element_text(size=15), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x = "", y="Natural-Log Required Sample Size") 
K


ggarrange(A, B, I, B, H, G, E, common.legend = T, nrow=2, ncol=4)

Theme <- theme(legend.text=
                 element_text(size=14), legend.title = element_text(size=14))
E <- E + Theme
G <- G + Theme
A <- A + Theme
C <- C + Theme

BigA <- ggarrange(E, G, common.legend = T, nrow=1, ncol=2)
BigB <- ggarrange(A, C, common.legend = T, nrow=1, ncol=2)

ggarrange(BigA, BigB, common.legend = T, nrow=2, ncol=1) 

#GGpredict on model estimates
library(ggeffects)
Mod_XGB <- (glm.nb(data=Dataset_TEST, N ~ factor(I(Dispersion>1)) + AUC_Full + Min_LFC + Med_LFC))
Mod_RF <- (glm.nb(data=Dataset_RF_TEST, N_RF ~ Max_LFC + factor(I(Dispersion>1)) + Imbalance + Complexity))
Mod_NN <- (glm.nb(data=Dataset_NN_TEST, N_NN ~ Complexity + AUC_Full_NN + Max_LFC + factor(I(Dispersion>1))))




#Check Agreement Dispersions-----------------------------------------------------------------------

set.seed(2024)
Dispersion_Simulated <- c(median(GetDispersion_Sim(GBM, F, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)])),
                median(GetDispersion_Sim(HCC, F, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)])),
                median(GetDispersion_Sim(Kidney, F, colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)])),
                median(GetDispersion_Sim(IPF, F, colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)])),
                median(GetDispersion_Sim(COVID, T, colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)])),
                median(GetDispersion_Sim(HBV, F, colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)])),
                median(GetDispersion_Sim(Hypertension, F, colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)])),
                median(GetDispersion_Sim(NAFLD, F, colnames(NAFLD_Sim)[1:(length(colnames(NAFLD_Sim))-1)])),
                median(GetDispersion_Sim(PrePost, F, colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)])),
                median(GetDispersion_Sim(RA, F, colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)])),
                median(GetDispersion_Sim(AVSC, T, colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)])),
                median(GetDispersion_Sim(Crohn, F, colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)])),
                median(GetDispersion_Sim(MDD, F, colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)])),
                median(GetDispersion_Sim(Ovarian, F, colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)])),
                median(GetDispersion_Sim(ALS, F, colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)])),
                median(GetDispersion_Sim(MISC, F, colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)])),
                median(GetDispersion_Sim(CCA, F, colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)])),
                median(GetDispersion_Sim(Glioma, F, colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)])),
                median(GetDispersion_Sim(TB, F, colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)])),
                median(GetDispersion_Sim(Tub, F, colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)])),
                median(GetDispersion_Sim(EDS, F, colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)])),
                median(GetDispersion_Sim(NSCLC, F, colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)])),
                median(GetDispersion_Sim(Bipolar, F, colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)])),
                median(GetDispersion_Sim(PDAC, F, colnames(PDAC_Sim)[1:(length(colnames(PDAC_Sim))-1)])),
                median(GetDispersion_Sim(MS, F, colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])))

Dispersion_Simulated
#[1] 0.32241168 1.33292748 0.87762952 0.02871027 9.64928220 0.14182734 1.31269917 0.11748022 0.41902412
#[10] 0.07555905 2.04585115 0.55446814 0.16079212 0.35024109 2.86843758 0.35780604 1.50789033 1.01346870
#[19] 9.73620391 0.07331635 0.38790421 0.23388442 0.07473904 0.84840942 0.80035934

set.seed(2024)
c(
median(GetDispersion_Sim(IBD, F, colnames(IBD_Sim)[1:(length(colnames(IBD_Sim))-1)])),
median(GetDispersion_Sim(Virus, F, colnames(Virus_Sim)[1:(length(colnames(Virus_Sim))-1)])),
median(GetDispersion_Sim(Gastric, F, colnames(Gastric_Sim)[1:(length(colnames(Gastric_Sim))-1)])))


plot(Dispersion, Dispersion_Simulated)
cor(Dispersion, Dispersion_Simulated)

summary(glm.nb(Dataset_TEST$N ~ Dispersion_Simulated))
summary(glm.nb(Dataset_RF_TEST$N_RF ~ Dispersion_Simulated))


#NEURAL NETWORKS
S0 <- seq(25, 250, by=225/10)
S1 <- seq(50, 250, by=200/10)
S2 <- seq(50, 500, by=450/10)
S3 <- seq(50, 2000, by=1950/10)

S4 <- seq(100, 2000, by=1900/10)

S5 <- seq(50, 1000, by=950/10)


NN_Optimal(GBM_Sim, "Group", 100)
NN_Optimal(HCC_Sim, "Group", 100)
NN_Optimal(Kidney_Sim, "Group", 100)
NN_Optimal(IPF_Sim, "Group", 100)
NN_Optimal(Hypertension_Sim, "Group", 100)
NN_Optimal(HBV_Sim, "Group", 100)
NN_Optimal(COVID_Sim, "Group", 100)
NN_Optimal(NAFLD_Sim, "Group", 100)
NN_Optimal(PrePost_Sim, "Group", 100)
NN_Optimal(RA_Sim, "Group", 100)
NN_Optimal(AVSC_Sim, "Group", 100)
NN_Optimal(Crohn_Sim, "Group", 100)
NN_Optimal(MDD_Sim, "Group", 100)
NN_Optimal(Ovarian_Sim, "Group", 100)
NN_Optimal(ALS_Sim, "Group", 100)
NN_Optimal(MISC_Sim, "Group", 100)
NN_Optimal(CCA_Sim, "Group", 100)
NN_Optimal(Glioma_Sim, "Group", 100)
NN_Optimal(TB_Sim, "Group", 100)
NN_Optimal(Tuberculosis_Sim, "Group", 100)
NN_Optimal(EDS_Sim, "Group", 100)
NN_Optimal(NSCLC_Sim, "Group", 100)
NN_Optimal(Bipolar_Sim, "Group", 100)
NN_Optimal(PDAC_Sim, "Group", 100)
NN_Optimal(MS_Sim, "Group", 100)

#GBM
GBM_NN <- NULL
for (i in S1) {
  GBM_NN <- rbind(GBM_NN, NN_Curve_Evaluate(GBM_Sim, "Group", n=i, True=colnames(GBM_Sim)[-ncol(GBM_Sim)], h=100))
}
GBM_NN

plot_curve_full(GBM_NN, "PowerLaw_Fixed", 1, 250)
Get_N_Fixed(GBM_NN, 1, 250) 

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(GBM_NN, "GBM_NN.csv")

#HCC
HCC_NN <- NULL
for (i in S2) {
  HCC_NN <- rbind(HCC_NN, NN_Curve_Evaluate(HCC_Sim, "Group", n=i, True=colnames(HCC_Sim)[-ncol(HCC_Sim)], h=100))
}
HCC_NN

plot_curve_full(HCC_NN, "PowerLaw", 0.994, 500)
Get_N(HCC_NN, 0.994, 500) #324

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(HCC_NN, "HCC_NN.csv")

#Kidney
Kidney_NN <- NULL
for (i in S3) {
  Kidney_NN <- rbind(Kidney_NN, NN_Curve_Evaluate(Kidney_Sim, "Group", n=i, True=colnames(Kidney_Sim)[-ncol(Kidney_Sim)], h=100))
}
Kidney_NN

plot_curve_full(Kidney_NN, "PowerLaw", 0.862, 2000)
Get_N_Fixed(Kidney_NN, 0.862, 2000) #227

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Kidney_NN, "Kidney_NN.csv")


#IPF
IPF_NN <- NULL
for (i in S0) {
  IPF_NN <- rbind(IPF_NN, NN_Curve_Evaluate(IPF_Sim, "Group", n=i, True=colnames(IPF_Sim)[-ncol(IPF_Sim)], h=100))
}
IPF_NN

plot_curve_full(IPF_NN, "PowerLaw_Fixed", 1, 250)
Get_N_Fixed(IPF_NN, 1, 250)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(IPF_NN, "IPF_NN.csv")


#COVID
COVID_NN <- NULL
for (i in S3) {
  COVID_NN <- rbind(COVID_NN, NN_Curve_Evaluate(COVID_Sim, "Group", n=i, True=colnames(COVID_Sim)[-ncol(COVID_Sim)], h=100))
}
COVID_NN

plot_curve_full(COVID_NN, "LOG", 0.947, 2000)
Get_N_LOG(COVID_NN, 0.947, 2000) #1374

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(COVID_NN, "COVID_NN.csv")


#HBV
HBV_NN <- NULL
for (i in S4) {
  HBV_NN <- rbind(HBV_NN, NN_Curve_Evaluate(HBV_Sim, "Group", n=i, True=colnames(HBV_Sim)[-ncol(HBV_Sim)], h=100))
}
HBV_NN

plot_curve_full(HBV_NN[-1,], "PowerLaw_Fixed", 0.896, 2000)
Get_N_Fixed(HBV_NN[-1,], 0.896, 2000) #144

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(HBV_NN, "HBV_NN.csv")


#Hypertension
Hypertension_NN <- NULL
for (i in S2) {
  Hypertension_NN <- rbind(Hypertension_NN, NN_Curve_Evaluate(Hypertension_Sim, "Group", n=i, True=colnames(Hypertension_Sim)[-ncol(Hypertension_Sim)], h=100))
}
Hypertension_NN

plot_curve_full(Hypertension_NN, "PowerLaw_Fixed", 0.999, 500)
Get_N_Fixed(Hypertension_NN, 0.999, 500) #212

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Hypertension_NN, "Hypertension_NN.csv")

#NAFLD
NAFLD_NN <- NULL
for (i in S2) {
  NAFLD_NN <- rbind(NAFLD_NN, NN_Curve_Evaluate(NAFLD_Sim, "Group", n=i, True=colnames(NAFLD_Sim)[-ncol(NAFLD_Sim)], h=100))
}
NAFLD_NN

plot_curve_full(NAFLD_NN, "PowerLaw", 1, 500)
Get_N(NAFLD_NN, 1, 500) #41

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(NAFLD_NN, "NAFLD_NN.csv")


#PrePost
PrePost_NN <- NULL
for (i in S2) {
  PrePost_NN <- rbind(PrePost_NN, NN_Curve_Evaluate(PrePost_Sim, "Group", n=i, True=colnames(PrePost_Sim)[-ncol(PrePost_Sim)], h=100))
}
PrePost_NN

plot_curve_full(PrePost_NN, "PowerLaw", 0.997, 500)
Get_N(PrePost_NN, 0.997, 500) #199

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(PrePost_NN, "PrePost_NN.csv")


#RA
RA_NN <- NULL
for (i in S3) {
  RA_NN <- rbind(RA_NN, NN_Curve_Evaluate(RA_Sim, "Group", n=i, True=colnames(RA_Sim)[-ncol(RA_Sim)], h=100))
}
RA_NN

plot_curve_full(RA_NN, "LOG", 0.988, 2000)
Get_N_LOG(RA_NN, 0.988, 2000) #761

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(RA_NN, "RA_NN.csv")


#AVSC
AVSC_NN <- NULL
for (i in S3) {
  AVSC_NN <- rbind(AVSC_NN, NN_Curve_Evaluate(AVSC_Sim, "Group", n=i, True=colnames(AVSC_Sim)[-ncol(AVSC_Sim)], h=100))
}
AVSC_NN

plot_curve_full(AVSC_NN[-1,], "PowerLaw", 0.858, 2000)
Get_N(AVSC_NN[-1,], 0.858, 2000) #888

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(AVSC_NN, "AVSC_NN.csv")

#Crohn
Crohn_NN <- NULL
for (i in S3) {
  Crohn_NN <- rbind(Crohn_NN, NN_Curve_Evaluate(Crohn_Sim, "Group", n=i, True=colnames(Crohn_Sim)[-ncol(Crohn_Sim)], h=100))
}
Crohn_NN

plot_curve_full(Crohn_NN, "LOG", 0.913, 2000)
Get_N_LOG(Crohn_NN, 0.913, 2000) #681

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Crohn_NN, "Crohn_NN.csv")

#MDD
MDD_NN <- NULL
for (i in S3) {
  MDD_NN <- rbind(MDD_NN, NN_Curve_Evaluate(MDD_Sim, "Group", n=i, True=colnames(MDD_Sim)[-ncol(MDD_Sim)], h=100))
}
MDD_NN

plot_curve_full(MDD_NN, "PowerLaw_Fixed", 0.874, 2000)
Get_N_Fixed(MDD_NN, 0.874, 2000) #156

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(MDD_NN, "MDD_NN.csv")

#Ovarian
Ovarian_NN <- NULL
for (i in S3) {
  Ovarian_NN <- rbind(Ovarian_NN, NN_Curve_Evaluate(Ovarian_Sim, "Group", n=i, True=colnames(Ovarian_Sim)[-ncol(Ovarian_Sim)], h=100))
}
Ovarian_NN

plot_curve_full(Ovarian_NN, "LOG", 0.861, 2000)
Get_N_LOG(Ovarian_NN, 0.861, 2000) #565

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Ovarian_NN, "Ovarian_NN.csv")

#ALS
ALS_NN <- NULL
for (i in S3) {
  ALS_NN <- rbind(ALS_NN, NN_Curve_Evaluate(ALS_Sim, "Group", n=i, True=colnames(ALS_Sim)[-ncol(ALS_Sim)], h=100))
}
ALS_NN

plot_curve_full(ALS_NN, "PowerLaw", 0.964, 2000)
Get_N(ALS_NN, 0.964, 2000) #1537

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(ALS_NN, "ALS_NN.csv")

#MISC
MISC_NN <- NULL
for (i in S1) {
  MISC_NN <- rbind(MISC_NN, NN_Curve_Evaluate(MISC_Sim, "Group", n=i, True=colnames(MISC_Sim)[-ncol(MISC_Sim)], h=100))
}
MISC_NN

plot_curve_full(MISC_NN, "PowerLaw", 1, 250)
Get_N(MISC_NN, 1, 250) #57

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(MISC_NN, "MISC_NN.csv")


#CCA
CCA_NN <- NULL
for (i in S3) {
  CCA_NN <- rbind(CCA_NN, NN_Curve_Evaluate(CCA_Sim, "Group", n=i, True=colnames(CCA_Sim)[-ncol(CCA_Sim)], h=100))
}
CCA_NN

plot_curve_full(CCA_NN, "PowerLaw", 0.935, 2000)
Get_N(CCA_NN, 0.935, 2000) #960

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(CCA_NN, "CCA_NN.csv")


#Glioma
Glioma_NN <- NULL
for (i in S3) {
  Glioma_NN <- rbind(Glioma_NN, NN_Curve_Evaluate(Glioma_Sim, "Group", n=i, True=colnames(Glioma_Sim)[-ncol(Glioma_Sim)], h=100))
}
Glioma_NN

plot_curve_full(Glioma_NN, "LOG", 0.975, 2000)
Get_N_LOG(Glioma_NN, 0.975, 2000) #1047

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Glioma_NN, "Glioma_NN.csv")

#TB_Prog
TB_NN <- NULL
for (i in S3) {
  TB_NN <- rbind(TB_NN, NN_Curve_Evaluate(TB_Sim, "Group", n=i, True=colnames(TB_Sim)[-ncol(TB_Sim)], h=100))
}
TB_NN

plot_curve_full(TB_NN, "PowerLaw", 0.922, 2000)
Get_N(TB_NN, 0.922, 2000) #1444

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(TB_NN, "TB_NN.csv")

#Tuberculosus
Tub_NN <- NULL
for (i in S2) {
  Tub_NN <- rbind(Tub_NN, NN_Curve_Evaluate(Tuberculosis_Sim, "Group", n=i, True=colnames(Tuberculosis_Sim)[-ncol(Tuberculosis_Sim)], h=100))
}
Tub_NN

plot_curve_full(Tub_NN, "PowerLaw_Fixed", 0.998, 500)
Get_N_Fixed(Tub_NN, 0.998, 500) #239

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Tub_NN, "Tub_NN.csv")

#EDS
EDS_NN <- NULL
for (i in S1) {
  EDS_NN <- rbind(EDS_NN, NN_Curve_Evaluate(EDS_Sim, "Group", n=i, True=colnames(EDS_Sim)[-ncol(EDS_Sim)], h=100))
}
EDS_NN

plot_curve_full(EDS_NN, "PowerLaw", 1, 250)
Get_N(EDS_NN, 1, 250) #94

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(EDS_NN, "EDS_NN.csv")

#NSCLC
NSCLC_NN <- NULL
for (i in S5) {
  NSCLC_NN <- rbind(NSCLC_NN, NN_Curve_Evaluate(NSCLC_Sim, "Group", n=i, True=colnames(NSCLC_Sim)[-ncol(NSCLC_Sim)], h=100))
}
NSCLC_NN

plot_curve_full(NSCLC_NN, "PowerLaw", 0.996, 1000)
Get_N(NSCLC_NN, 0.996, 1000) #371

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(NSCLC_NN, "NSCLC_NN.csv")


#Bipolar
Bipolar_NN <- NULL
for (i in S3) {
  Bipolar_NN <- rbind(Bipolar_NN, NN_Curve_Evaluate(Bipolar_Sim, "Group", n=i, True=colnames(Bipolar_Sim)[-ncol(Bipolar_Sim)], h=100))
}
Bipolar_NN

plot_curve_full(Bipolar_NN, "PowerLaw", 0.951, 2000)
Get_N(Bipolar_NN, 0.951, 2000) #1471

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Bipolar_NN, "Bipolar_NN.csv")


#PDAC
PDAC_NN <- NULL
for (i in S1) {
  PDAC_NN <- rbind(PDAC_NN, NN_Curve_Evaluate(PDAC_Sim, "Group", n=i, True=colnames(PDAC_Sim)[-ncol(PDAC_Sim)], h=100))
}
PDAC_NN

plot_curve_full(PDAC_NN, "PowerLaw", 0.999, 250)
Get_N(PDAC_NN, 0.999, 250) #81

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(PDAC_NN, "PDAC_NN.csv")

#MS
MS_NN <- NULL
for (i in S3) {
  MS_NN <- rbind(MS_NN, NN_Curve_Evaluate(MS_Sim, "Group", n=i, True=colnames(MS_Sim)[-ncol(MS_Sim)], h=100))
}
MS_NN

plot_curve_full(MS_NN, "LOG", 0.907, 2000)
Get_N_LOG(MS_NN, 0.907, 2000) #1711

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(MS_NN, "MS_NN.csv")

#VALDIATION
#IBD
IBD_NN <- NULL
for (i in S3) {
  IBD_NN <- rbind(IBD_NN, NN_Curve_Evaluate(IBD_Sim, "Group", n=i, True=colnames(IBD_Sim)[-ncol(IBD_Sim)], h=100))
}
IBD_NN

plot_curve_full(IBD_NN, "PowerLaw", 0.989, 2000)
Get_N(IBD_NN, 0.989, 2000) #1307

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(IBD_NN, "IBD_NN.csv")

#Virus
Virus_NN <- NULL
for (i in S0) {
  Virus_NN <- rbind(Virus_NN, NN_Curve_Evaluate(Virus_Sim, "Group", n=i, True=colnames(Virus_Sim)[-ncol(Virus_Sim)], h=100))
}
Virus_NN

plot_curve_full(Virus_NN, "LOG", 1, 250)
Get_N(Virus_NN, 1, 250) #25

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Virus_NN, "Virus_NN.csv")

#Gastric
Gastric_NN <- NULL
for (i in S3) {
  Gastric_NN <- rbind(Gastric_NN, NN_Curve_Evaluate(Gastric_Sim, "Group", n=i, True=colnames(Gastric_Sim)[-ncol(Gastric_Sim)], h=100))
}
Gastric_NN

plot_curve_full(Gastric_NN, "LOG", 0.827, 2000)
Get_N_LOG(Gastric_NN, 0.827, 2000) #688

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
write.csv(Gastric_NN, "Gastric_NN.csv")

NN_Optimal(Gastric_Sim, "Group", 100) #0.827
LR_Optimal(cbind(log2(Gastric_Sim[,-ncol(Gastric_Sim)]+0.01), "Group"=Gastric_Sim$Group), "Group")[1]
#LOW

#max 0.24


#LOAD NN AND GET AUC25
-----------------------------------

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
  
GBM_NN <- read.csv("GBM_NN.csv")[,-1]  
HCC_NN <- read.csv("HCC_NN.csv")[,-1]  
Kidney_NN <- read.csv("Kidney_NN.csv")[,-1]  
IPF_NN <- read.csv("IPF_NN.csv")[,-1]  
COVID_NN <- read.csv("COVID_NN.csv")[,-1]  
HBV_NN <- read.csv("HBV_NN.csv")[,-1]  
Hypertension_NN <- read.csv("Hypertension_NN.csv")[,-1]  
NAFLD_NN <- read.csv("NAFLD_NN.csv")[,-1]  
PrePost_NN <- read.csv("PrePost_NN.csv")[,-1]  
RA_NN <- read.csv("RA_NN.csv")[,-1]  
AVSC_NN <- read.csv("AVSC_NN.csv")[,-1]  
Crohn_NN <- read.csv("Crohn_NN.csv")[,-1]  
MDD_NN <- read.csv("MDD_NN.csv")[,-1]  
Ovarian_NN <- read.csv("Ovarian_NN.csv")[,-1]  
ALS_NN <- read.csv("ALS_NN.csv")[,-1]  
MISC_NN <- read.csv("MISC_NN.csv")[,-1]  
CCA_NN <- read.csv("CCA_NN.csv")[,-1]  
Glioma_NN <- read.csv("Glioma_NN.csv")[,-1]  
TB_NN <- read.csv("TB_NN.csv")[,-1]  
Tub_NN <- read.csv("Tub_NN.csv")[,-1]  
EDS_NN <- read.csv("EDS_NN.csv")[,-1]  
NSCLC_NN <- read.csv("NSCLC_NN.csv")[,-1]  
Bipolar_NN <- read.csv("Bipolar_NN.csv")[,-1]  
MS_NN <- read.csv("MS_NN.csv")[,-1]  
TCGA_BRCA_NN <- read.csv("TCGA_BRCA_NN.csv")[,-1]  
TCGA_BRAIN_NN <- read.csv("TCGA_BRAIN_NN.csv")[,-1]  
TCGA_LUNG_NN <- read.csv("TCGA_LUNG_NN.csv")[,-1]  

AUCs_NN <- c(1.000,0.994,0.862,1.000,0.947,0.896,0.999,1.000,0.997,0.988,0.858,0.913,0.874,0.861,0.964,1.000,0.935,0.975,0.922,0.998,1.000,0.996,0.951,0.907,1,0.724,0.981)


N_fromcurve_NN <- c(
  50,
  Get_N(HCC_NN, AUCs_NN[2], 2000),
  Get_N_Fixed(Kidney_NN, AUCs_NN[3], 2000),
  25,
  Get_N_LOG(COVID_NN, AUCs_NN[5], 2000),
  Get_N_Fixed(HBV_NN[-1,], AUCs_NN[6], 2000),
  Get_N_Fixed(Hypertension_NN, AUCs_NN[7], 2000),
  50,
  Get_N(PrePost_NN, AUCs_NN[9], 2000),
  Get_N_LOG(RA_NN, AUCs_NN[10], 2000),
  Get_N(AVSC_NN[-1,], AUCs_NN[11], 2000),
  Get_N_LOG(Crohn_NN, AUCs_NN[12], 2000),
  Get_N_Fixed(MDD_NN, AUCs_NN[13], 2000),
  Get_N_LOG(Ovarian_NN, AUCs_NN[14], 2000),
  Get_N(ALS_NN, AUCs_NN[15], 2000),
  Get_N(MISC_NN, AUCs_NN[16], 2000),
  Get_N(CCA_NN, AUCs_NN[17], 2000),
  Get_N_LOG(Glioma_NN, AUCs_NN[18], 2000),
  Get_N(TB_NN, AUCs_NN[19], 2000),
  Get_N_Fixed(Tub_NN, AUCs_NN[20], 2000),
  Get_N(EDS_NN, AUCs_NN[21], 2000),
  Get_N(NSCLC_NN, AUCs_NN[22], 2000),
  Get_N(Bipolar_NN, AUCs_NN[23], 2000),
  Get_N_LOG(MS_NN, AUCs_NN[24], 2000),
  25,
  Get_N_LOG(TCGA_BRCA_NN, AUCs_NN[26], 2000),
  Get_N(TCGA_LUNG_NN, AUCs_NN[27], 2000)
)

N_fromcurve_NN_0.01 <- c(
  Get_N_Fixed_Th(GBM_NN, AUCs_NN[1], 2500, 0.01),
  Get_N_Th(HCC_NN[-1], AUCs_NN[2], 2500, 0.01),
  Get_N_Fixed_Th(Kidney_NN, AUCs_NN[3], 2500, 0.01),
  Get_N_Fixed_Th(IPF_NN, AUCs_NN[4], 2500, 0.01),
  Get_N_LOG_Th(COVID_NN, AUCs_NN[5], 2500, 0.01),
  Get_N_Fixed_Th(HBV_NN[-1,], AUCs_NN[6], 2500, 0.01),
  Get_N_Fixed_Th(Hypertension_NN, AUCs_NN[7], 2500, 0.01),
  Get_N_Th(NAFLD_NN, AUCs_NN[8], 2500, 0.01),
  Get_N_Th(PrePost_NN, AUCs_NN[9], 2500, 0.01),
  Get_N_LOG_Th(RA_NN, AUCs_NN[10], 2500, 0.01),
  Get_N_Th(AVSC_NN[-1,], AUCs_NN[11], 2500, 0.01),
  Get_N_LOG_Th(Crohn_NN, AUCs_NN[12], 2500, 0.01),
  Get_N_Fixed_Th(MDD_NN, AUCs_NN[13], 2500, 0.01),
  Get_N_LOG_Th(Ovarian_NN, AUCs_NN[14], 2500, 0.01),
  Get_N_Th(ALS_NN, AUCs_NN[15], 2500, 0.01),
  Get_N_Th(MISC_NN, AUCs_NN[16], 2500, 0.01),
  Get_N_Th(CCA_NN, AUCs_NN[17], 2500, 0.01),
  Get_N_LOG_Th(Glioma_NN, AUCs_NN[18], 2500, 0.01),
  Get_N_Th(TB_NN, AUCs_NN[19], 2500, 0.01),
  Get_N_Fixed_Th(Tub_NN, AUCs_NN[20], 2500, 0.01),
  Get_N_Th(EDS_NN, AUCs_NN[21], 2500, 0.01),
  Get_N_Th(NSCLC_NN, AUCs_NN[22], 2500, 0.01),
  Get_N_Th(Bipolar_NN, AUCs_NN[23], 2500, 0.01),
  Get_N_LOG_Th(MS_NN, AUCs_NN[24], 2500, 0.01),
  Get_N_Th(TCGA_BRAIN_NN, AUCs_NN[25], 2500, 0.01),
  Get_N_LOG_Th(TCGA_BRCA_NN, AUCs_NN[26], 2500, 0.01),
  Get_N_Th(TCGA_LUNG_NN, AUCs_NN[27], 2500, 0.01)
)

N_fromcurve_NN_0.01

N_fromcurve_NN_0.05 <- c(
  Get_N_Fixed_Th(GBM_NN, AUCs_NN[1], 2000, 0.05),
  Get_N_Th(HCC_NN, AUCs_NN[2], 2000, 0.05),
  Get_N_Fixed_Th(Kidney_NN, AUCs_NN[3], 2000, 0.05),
  Get_N_Fixed_Th(IPF_NN, AUCs_NN[4], 2000, 0.05),
  Get_N_LOG_Th(COVID_NN, AUCs_NN[5], 2000, 0.05),
  Get_N_Fixed_Th(HBV_NN[-1,], AUCs_NN[6], 2000, 0.05),
  Get_N_Fixed_Th(Hypertension_NN, AUCs_NN[7], 2000, 0.05),
  Get_N_Th(NAFLD_NN, AUCs_NN[8], 2000, 0.05),
  Get_N_Th(PrePost_NN, AUCs_NN[9], 2000, 0.05),
  Get_N_LOG_Th(RA_NN, AUCs_NN[10], 2000, 0.05),
  Get_N_Th(AVSC_NN[-1,], AUCs_NN[11], 2000, 0.05),
  Get_N_LOG_Th(Crohn_NN, AUCs_NN[12], 2000, 0.05),
  Get_N_Fixed_Th(MDD_NN, AUCs_NN[13], 2000, 0.05),
  Get_N_LOG_Th(Ovarian_NN, AUCs_NN[14], 2000, 0.05),
  Get_N_Th(ALS_NN, AUCs_NN[15], 2000, 0.05),
  Get_N_Th(MISC_NN, AUCs_NN[16], 2000, 0.05),
  Get_N_Th(CCA_NN, AUCs_NN[17], 2000, 0.05),
  Get_N_LOG_Th(Glioma_NN, AUCs_NN[18], 2000, 0.05),
  Get_N_Th(TB_NN, AUCs_NN[19], 2000, 0.05),
  Get_N_Fixed_Th(Tub_NN, AUCs_NN[20], 2000, 0.05),
  Get_N_Th(EDS_NN, AUCs_NN[21], 2000, 0.05),
  Get_N_Th(NSCLC_NN, AUCs_NN[22], 2000, 0.05),
  Get_N_Th(Bipolar_NN, AUCs_NN[23], 2000, 0.05),
  Get_N_LOG_Th(MS_NN, AUCs_NN[24], 2000, 0.05),
  Get_N_Th(TCGA_BRAIN_NN, AUCs_NN[25], 2000, 0.05),
  Get_N_LOG_Th(TCGA_BRCA_NN, AUCs_NN[26], 2000, 0.05),
  Get_N_Th(TCGA_LUNG_NN, AUCs_NN[27], 2000, 0.05)
)

N_fromcurve_NN_0.05


#-------------------------------------------------------------------------------------------------------------------------------------

#SAMPLE SIZES NEEDED FOR DIFF EX FROM PACKAGES
library(ssizeRNA)
library(RNASeqPower)
library(RnaSeqSampleSize)

GetCV <- function(DF, Keep) {
  
  Sim <- Simulate_OMICs_Data_BNG_Only(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                                      p=ncol(DF)-1, n=5000)
  
  
  #Back-transform to neg-bin distribution
  Sim <- cbind((2^Sim[,-ncol(Sim)])-0.01, "Group"=Sim$Group)
  
  Sim[,1:(ncol(Sim)-1)] <- round(Sim[,1:(ncol(Sim)-1)], 0)
  
  median(apply(Sim[,-ncol(Sim)], 2, sd) / apply(Sim[,-ncol(Sim)], 2, mean))
}

GetMeanRdCt <- function(DF, Keep) {
  
  Sim <- Simulate_OMICs_Data_BNG_Only(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                                      p=ncol(DF)-1, n=5000)
  
  
  #Back-transform to neg-bin distribution
  Sim <- cbind((2^Sim[,-ncol(Sim)])-0.01, "Group"=Sim$Group)
  
  Sim[,1:(ncol(Sim)-1)] <- round(Sim[,1:(ncol(Sim)-1)], 0)
  
  median(apply(Sim[,-ncol(Sim)], 2, mean))
}



CVs <- c(
GetCV(GBM, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)]),
GetCV(HCC, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)]),
GetCV(IPF,colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)]),
GetCV(Kidney,colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)]),
GetCV(COVID,colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)]),
GetCV(HBV,colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)]),
GetCV(Hypertension,colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)]),
GetCV(NAFLD,colnames(NAFLD_Sim)[1:(length(colnames(NAFLD_Sim))-1)]),
GetCV(PrePost,colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)]),
GetCV(RA,colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)]),
GetCV(AVSC,colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)]),
GetCV(Crohn,colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)]),
GetCV(MDD,colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)]),
GetCV(Ovarian,colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)]),
GetCV(ALS,colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)]),
GetCV(MISC,colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)]),
GetCV(CCA,colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)]),
GetCV(Glioma,colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)]),
GetCV(TB,colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)]),
GetCV(Tub,colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)]),
GetCV(EDS,colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)]),
GetCV(NSCLC,colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)]),
GetCV(Bipolar,colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)]),
GetCV(PDAC,colnames(PDAC_Sim)[1:(length(colnames(PDAC_Sim))-1)]),
GetCV(MS,colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])
)

RawReadCts <-  c(
    GetMeanRdCt(GBM, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)]),
    GetMeanRdCt(HCC, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)]),
    GetMeanRdCt(IPF,colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)]),
    GetMeanRdCt(Kidney,colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)]),
    GetMeanRdCt(COVID,colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)]),
    GetMeanRdCt(HBV,colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)]),
    GetMeanRdCt(Hypertension,colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)]),
    GetMeanRdCt(NAFLD,colnames(NAFLD_Sim)[1:(length(colnames(NAFLD_Sim))-1)]),
    GetMeanRdCt(PrePost,colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)]),
    GetMeanRdCt(RA,colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)]),
    GetMeanRdCt(AVSC,colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)]),
    GetMeanRdCt(Crohn,colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)]),
    GetMeanRdCt(MDD,colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)]),
    GetMeanRdCt(Ovarian,colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)]),
    GetMeanRdCt(ALS,colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)]),
    GetMeanRdCt(MISC,colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)]),
    GetMeanRdCt(CCA,colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)]),
    GetMeanRdCt(Glioma,colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)]),
    GetMeanRdCt(TB,colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)]),
    GetMeanRdCt(Tub,colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)]),
    GetMeanRdCt(EDS,colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)]),
    GetMeanRdCt(NSCLC,colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)]),
    GetMeanRdCt(Bipolar,colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)]),
    GetMeanRdCt(PDAC,colnames(PDAC_Sim)[1:(length(colnames(PDAC_Sim))-1)]),
    GetMeanRdCt(MS,colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])
)

#CVs <- 0.7009680  1.6937870  0.3492100  1.4366484 22.6857659  0.6517944  2.4434466  0.3799816  1.0862418  0.4242705  3.1731267  2.2724914
#0.6556912  0.8484268  4.1434765  1.2048906  3.8957280  1.8343512 16.0055238  0.3069531  0.9337846  1.1075895  0.6365479  2.1515552
#1.8316428

#Raw Read Means <- 
#1274.2212   271.2234    83.5392   300.7956 13147.5010    53.4326    89.3072   401.6335    68.0118   350.9078    28.8212    36.1783
#157.2697   429.5918    32.0327   404.7818   712.0524    90.7632   313.7218  1380.6476   340.7968   361.7378   121.1902  6252.7318
#97.3696
#

TotalGenes <- c(33445,37826,34202,29434,26769,39054,38985,33377,
                                         38972,37288,34428,33258,39333,37901,38787,39308,
                                         39110,39042,37056,34825,23239,39230,39342,35922, 39290
)

FC <- 2^Med_LFC

CVs <- ifelse(CVs>5, 5, CVs)

2*round(sapply(1:25, function(i)
{rnapower(depth=RawReadCts[i],
         cv=CVs[i],
         effect=FC[i],
         alpha=0.05,
         power=0.8)}
),0)



set.seed(2024)
ssizeRNA_sizes <- 2*sapply(1:25, function(i) {
ssizeRNA_single(nGenes = TotalGenes[i],
                pi0= 1 - (N_Features[i] / TotalGenes[i]), 
                m=200,
                mu=RawReadCts[i],
                disp=Dispersion_Simulated[i],
                fc=FC[i],
                up=1,
                fdr=0.05,
                power=0.8,
                maxN = 2000,
                side="two-sided")$ssize[2]
})


ssizeRNA_sizes

ssizeRNA_single(nGenes = TotalGenes[19],
                pi0= 1 - (N_Features[19] / TotalGenes[19]), 
                m=200,
                mu=RawReadCts[19],
                disp=Dispersion_Simulated[19],
                fc=FC[19],
                up=1,
                fdr=0.05,
                power=0.8,
                maxN = 10000,
                side="two-sided")$ssize[2] #18082

ssizeRNA_single(nGenes = TotalGenes[15],
                pi0= 1 - (N_Features[15] / TotalGenes[15]), 
                m=200,
                mu=RawReadCts[15],
                disp=Dispersion_Simulated[15],
                fc=FC[15],
                up=1,
                fdr=0.05,
                power=0.8,
                maxN = 10000,
                side="two-sided")$ssize[2] #2866




#EXPLORE IBD DAtatset
setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2")

IBD <- read.csv("GSE193677_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE193677")
IBD <- IBD[,-1]
IBD <- data.frame(t(IBD))

Meta <- data.frame("ID"=Meta$GSE193677_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE193677_series_matrix.txt.gz$`ibd_disease:ch1`=="Control", 0, 1))

IBD$ID <- rownames(IBD)

IBD <- left_join(IBD, Meta, "ID")
IBD <- IBD[,-39377]

Vars <- apply(IBD[,-ncol(IBD)], 2, var)
IBD <- IBD[,-which(colnames(IBD) %in% names(Vars[Vars==0]))]
IBD <- data.frame(apply(IBD, 2, as.numeric))

#Predictors

ExpDesign <- data.frame(row.names=rownames(IBD),
                        Group = IBD[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(IBD[,-which(colnames(IBD) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

IBD[,1:(ncol(IBD)-1)] <- DDS

#Find Predictors in IBDing Set + Simulate BNG Set
set.seed(2024)
B <- Boruta(x=IBD[,-which(colnames(IBD) %in% "Group")], y=as.factor(IBD[,which(colnames(IBD) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Keep

# "X447"   "X940"   "X1284"  "X2783"  "X2849"  "X2851"  "X3217"  "X3268"  "X3912"  "X4019"  "X4481"  "X5289"  "X5294"  "X5549"  "X5550" 
# "X5794"  "X5795"  "X5796"  "X5797"  "X5798"  "X6035"  "X7581"  "X7622"  "X7819"  "X8866"  "X9301"  "X9448"  "X9463"  "X9490"  "X9732" 
# "X9743"  "X10148" "X10149" "X10426" "X10749" "X12467" "X12978" "X13008" "X13451" "X14080" "X14741" "X15021" "X15203" "X15732" "X16349"
# "X17957" "X18100" "X18470" "X18595" "X18600" "X19078" "X19857" "X23138" "X23171" "X23825" "X25581" "X26663" "X26769" "X28159" "X28160"
# "X28161" "X28349" "X30000" "X31298" "X31362" "X32342" "X32463" "X32810" "X32834" "X34032" "X34876" "X35859" "X36142" "X36488" "X36712"
# "X37025" "X37026" "X37144" "X37192" "X38336" "X38821" "X38890"


#Analysis
IBD_Keep <- IBD[,c(Keep, "Group")]


IBD_Keep <- cbind(log2(IBD_Keep[,-ncol(IBD_Keep)]+0.01), "Group"=IBD_Keep$Group)

LFC_True <- sapply(1:(ncol(IBD_Keep)-1), function(i) {abs(mean(IBD_Keep[IBD_Keep$Group==1,i]) - mean(IBD_Keep[IBD_Keep$Group==0,i]))})

set.seed(2024)
IBD_Sim <- Simulate_OMICs_Data_BNG_Only(IBD_Keep, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(IBD_Keep[which(IBD_Keep$Group==1),])), 
                                     p=ncol(IBD_Keep)-1, n=5000)

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(IBD_Keep, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(IBD_Keep[which(IBD_Keep$Group==1),])), 
                                             p=ncol(IBD_Keep)-1, n=1000)


LFC_Sim <- sapply(1:(ncol(IBD_Sim)-1), function(i) {abs(mean(IBD_Sim[IBD_Sim$Group==1,i]) - mean(IBD_Sim[IBD_Sim$Group==0,i]))})
plot(LFC_Sim, LFC_True)
LFC_True
LFC_Sim

#Learning Curves for Each Dataset
#IBD_Keep is true
#IBD_Sim_2000 is simulated
#both evaluated on 471 observation test set

XGB_True_IBD <- XGBoost_Pipeline_Curve_2(IBD_Keep, 'Group', start=50, size=1500, steps=10, True=Keep, Test=Test[c(Keep, "Group")])

XGB_Sim_IBD <- XGBoost_Pipeline_Curve_2(IBD_Sim, 'Group', start=50, size=1500, steps=10, True=Keep, Test=Test[c(Keep, "Group")])

XGBoost_Optimal(IBD_Keep, "Group", Test)
XGBoost_Optimal(IBD_Sim_2000, "Group", Test)

plot_curve_full(XGB_True_IBD, "PowerLaw", 0.869, 1500)
Get_N(XGB_True_IBD, 0.869, 1500)

plot_curve_full(XGB_Sim_IBD, "PowerLaw", 0.981, 1500)
Get_N(XGB_Sim_IBD, 0.981, 1500)







ExpDesign <- data.frame(row.names=rownames(AVSC),
                        Group = AVSC[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(AVSC[,-which(colnames(AVSC) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

AVSC[,1:(ncol(AVSC)-1)] <- DDS

set.seed(2024)
B <- Boruta(x=AVSC[,-which(colnames(AVSC) %in% "Group")], y=as.factor(AVSC[,which(colnames(AVSC) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])


AVSC_Keep <- AVSC[,c(Keep, "Group")]


Vars_True <- apply(AVSC_Keep[,-ncol(AVSC_Keep)], 2, var)

AVSC_Keep <- cbind(log2(AVSC_Keep[,-ncol(AVSC_Keep)]+0.01), "Group"=AVSC_Keep$Group)

LFC_True <- sapply(1:(ncol(AVSC_Keep)-1), function(i) {abs(mean(AVSC_Keep[AVSC_Keep$Group==1,i]) - mean(AVSC_Keep[AVSC_Keep$Group==0,i]))})

set.seed(2024)
Test <- Simulate_OMICs_Data_BNG_Only(AVSC, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(AVSC[which(AVSC$Group==1),])), 
                                     p=ncol(AVSC)-1, n=5000)




setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Simulated Datasets")
write.csv(Test, "AVSC_Sim.csv")

LFC_Sim <- sapply(1:(ncol(Test)-1), function(i) {abs(mean(Test[Test$Group==1,i]) - mean(Test[Test$Group==0,i]))})
plot(LFC_Sim, LFC_True)
LFC_Sim

Test <- cbind((2^Test[,-ncol(Test)])-0.01, "Group"=Test$Group)

Vars_Sim <- apply(Test[,-ncol(Test)], 2, var)

plot(Vars_Sim, Vars_True)

Test

set.seed(2024)
XGBoost_Optimal(Test, "Group")

set.seed(2024)
RF_Optimal(Test, "Group")

NN_Optimal(Test, "Group", 3*3)

set.seed(2024)
AVSC_XGB <- XGBoost_Pipeline_Curve(Test, "Group", 50, 1500, steps=10, True=Keep)

set.seed(2024)
AVSC_RF <- RF_Pipeline_Curve(Test, "Group", 50, 1500, steps=10, True=Keep)

plot_curve_full(AVSC_XGB, "PowerLaw", 0.852, 1500)

plot_curve_full(AVSC_RF, "PowerLaw", 0.849, 1500)

Get_N(AVSC_RF, 0.852, 1500)
Get_N(AVSC_XGB, 0.852, 1500)

#Analysis of Covariates
#-------------------------------------------------------------------------------------------------------------------------

cor.test(AUC_Full, Max_LFC) #p=0.009
cor.test(AUC_Full, Min_LFC)
cor.test(AUC_Full, Med_LFC)

cor.test(AUC_Full, Complexity) #0.048
cor.test(AUC_Full, Dispersion_Simulated)
cor.test(AUC_Full, AvgMedRead)
cor.test(AUC_Full, Correlation) #p=0.013
cor.test(AUC_Full, N_Features) #p=0.001
cor.test(AUC_Full, Imbalance)

#COMPLEXITY
t.test(Max_LFC ~ Complexity) 
t.test(Min_LFC ~ Complexity) 
t.test(Med_LFC ~ Complexity) 

t.test(Dispersion_Simulated ~ Complexity)
t.test(AvgMedRead ~ Complexity) 
chisq.test(Correlation,Complexity)  
t.test(Dataset_TEST[Dataset_TEST$AUC_Full==0, "AvgMedRead"] ~ Dataset_TEST[Dataset_TEST$AUC_Full==0, "Complexity"])   #p=0.017
chisq.test(AUC_Full, Complexity) #0.0206

ComplexityMod <- glm(Complexity ~ AUC_Full + Dispersion_Simulated, family=binomial(logit))
summary(ComplexityMod)

table(unname(ifelse(predict(ComplexityMod, type="response") > 0.5, 1, 0)), Complexity)

predict(ComplexityMod, newdata=data.frame("N_Features"=c(67),
                                          "AUC_Full"=c(0),
                                          "Dispersion_Simulated"=c(1.07)), type="response")

Dispersion_25 <- function(DF, UseWhole, Keep) {
  Sim <- Simulate_OMICs_Data_BNG_Only(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                                      p=ncol(DF)-1, n=25)
  
  #Simulate DF from raw counts
  if(UseWhole==T) {
    
    Sim <- Simulate_OMICs_Data_SPLIT(DF, "Group", transform=T, fast=F, Keep=Keep, Case_Rows=as.numeric(rownames(DF[which(DF$Group==1),])), 
                                     p=ncol(DF)-1, n=25)
    
  }
  
  
  #Back-transform to neg-bin distribution
  Sim <- cbind((2^Sim[,-ncol(Sim)])-0.01, "Group"=Sim$Group)
  
  Sim[,1:(ncol(Sim)-1)] <- round(Sim[,1:(ncol(Sim)-1)], 0)
  
  ExpDesign <- data.frame(row.names=rownames(Sim))
  
  DiffEx <- estimateSizeFactors(DESeqDataSetFromMatrix(countData = as.matrix(t(Sim[,-which(colnames(Sim) %in% "Group")])), 
                                                       colData=ExpDesign, design=~1))
  
  
  DIS <- dispersions(estimateDispersions(DiffEx))
  
  DIS <- DIS[which(colnames(Sim) %in% Keep)]
  DIS
  
}

set.seed(2024)
median(Dispersion_25(Gastric, F, colnames(Gastric_Sim)[1:(length(colnames(Gastric_Sim))-1)])) #<1, 0.02
set.seed(2024)
median(Dispersion_25(IBD, F, colnames(IBD_Sim)[1:(length(colnames(IBD_Sim))-1)])) #>1, 1.52
set.seed(2024)
median(Dispersion_25(Virus, F, colnames(Virus_Sim)[1:(length(colnames(Virus_Sim))-1)])) #<1, 0.51

table(unname(ifelse(predict(ComplexityMod, newdata=data.frame("N_Features"=c(20,20,0),
                                                              "AUC_Full_RF"=c(0,0,1),
                                                              "Dispersion_Simulated"=c(0,1.07,0)), type="response") > 0.5, 1, 0)), c(0,1,0))

#-------------------------------------------------------------------------------------------------------------------------

#IBD Real
setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2")

IBD <- read.csv("GSE193677_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE193677")
IBD <- IBD[,-1]
IBD <- data.frame(t(IBD))

Meta <- data.frame("ID"=Meta$GSE193677_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE193677_series_matrix.txt.gz$`ibd_disease:ch1`=="Control", 0, 1))

IBD$ID <- rownames(IBD)

IBD <- left_join(IBD, Meta, "ID")
IBD <- IBD[,-39377]

Vars <- apply(IBD[,-ncol(IBD)], 2, var)
IBD <- IBD[,-which(colnames(IBD) %in% names(Vars[Vars==0]))]
IBD <- data.frame(apply(IBD, 2, as.numeric))

set.seed(2025)
Test <- IBD[sample(1:nrow(IBD), 471), ]
IBD <- IBD[-which(rownames(IBD) %in% rownames(Test)),]

#Predictors
ExpDesign <- data.frame(row.names=rownames(IBD),
                        Group = IBD[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(IBD[,-which(colnames(IBD) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

IBD[,1:(ncol(IBD)-1)] <- DDS

#Test Transformation
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test[,-which(colnames(Test) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))


DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Test[,1:(ncol(Test)-1)] <- DDS


#Find Predictors in IBD Set + Simulate BNG Set
set.seed(2024)
B <- Boruta(x=IBD[,-which(colnames(IBD) %in% "Group")], y=as.factor(IBD[,which(colnames(IBD) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

Keep
# "X447"   "X940"   "X1284"  "X2783"  "X2849"  "X2851"  "X3217"  "X3268"  "X3912"  "X4019"  "X4481"  "X5289"  "X5294"  "X5549"  "X5550" 
# "X5794"  "X5795"  "X5796"  "X5797"  "X5798"  "X6035"  "X7581"  "X7622"  "X7819"  "X8866"  "X9301"  "X9448"  "X9463"  "X9490"  "X9732" 
# "X9743"  "X10148" "X10149" "X10426" "X10749" "X12467" "X12978" "X13008" "X13451" "X14080" "X14741" "X15021" "X15203" "X15732" "X16349"
# "X17957" "X18100" "X18470" "X18595" "X18600" "X19078" "X19857" "X23138" "X23171" "X23825" "X25581" "X26663" "X26769" "X28159" "X28160"
# "X28161" "X28349" "X30000" "X31298" "X31362" "X32342" "X32463" "X32810" "X32834" "X34032" "X34876" "X35859" "X36142" "X36488" "X36712"
# "X37025" "X37026" "X37144" "X37192" "X38336" "X38821" "X38890"

#AT 2000
#"X447"   "X940"   "X1284"  "X2783"  "X2849"  "X3217"  "X3912"  "X4019"  "X4481"  "X5289"  "X5294"  "X5549"  "X5794"  "X5795"  "X5797" 
#"X5798"  "X7581"  "X7622"  "X7819"  "X8161"  "X8936"  "X9448"  "X9463"  "X10148" "X10149" "X10749" "X13451" "X15046" "X15203" "X15674"
#"X15732" "X15742" "X16349" "X17957" "X18100" "X19857" "X23138" "X25581" "X26663" "X26769" "X28159" "X28160" "X28161" "X28349" "X30000"
#"X31298" "X32463" "X32810" "X32834" "X34032" "X34876" "X34878" "X36488" "X36712" "X37026" "X38336" "X38821"


IBD_Keep <- IBD[,c("X447"  , "X940"   ,"X1284"  ,"X2783"  ,"X2849" , "X3217" , "X3912" , "X4019",  "X4481" , "X5289"  ,"X5294" , "X5549" , "X5794"  ,"X5795" , "X5797" ,
"X5798" , "X7581"  ,"X7622"  ,"X7819"  ,"X8161" , "X8936" , "X9448",  "X9463",  "X10148", "X10149" ,"X10749" ,"X13451" ,"X15046", "X15203", "X15674",
"X15732" ,"X15742" ,"X16349" ,"X17957", "X18100", "X19857", "X23138" ,"X25581" ,"X26663", "X26769" ,"X28159", "X28160", "X28161", "X28349", "X30000",
"X31298" ,"X32463","X32810" ,"X32834", "X34032", "X34876", "X34878" ,"X36488", "X36712" ,"X37026", "X38336" ,"X38821", "Group")]

#median(AvgReadCount(IBD_Keep))

#Analysis
IBD_Keep <- IBD[,c(Keep, "Group")]

IBD_Keep <- cbind(log2(IBD_Keep[,-ncol(IBD_Keep)]+0.01), "Group"=IBD_Keep$Group)

LFC_True <- sapply(1:(ncol(IBD_Keep)-1), function(i) {abs(mean(IBD_Keep[IBD_Keep$Group==1,i]) - mean(IBD_Keep[IBD_Keep$Group==0,i]))})
max(LFC_True) #2.63
median(LFC_True)
min(LFC_True)
1-mean(IBD_Keep$Group)

colnames(IBD_Keep)[48]

median(AvgReadCount(IBD_Keep))

#write.csv(Test_Keep, "IBD_Test.csv")
#write.csv(IBD_Keep, "IBD_Keep.csv")

LFC_True

Test_Keep <- Test[,c(Keep, "Group")]
Test_Keep <- Test[,c(colnames(IBD_Keep)[-58], "Group")]

Test_Keep <- cbind(log2(Test_Keep[,-ncol(Test_Keep)]+0.01), "Group"=Test_Keep$Group)

XGBoost_Optimal(IBD_Keep, "Group", Test_Keep) #0.917

RF_Optimal_2(IBD_Keep, "Group", Test_Keep) #0.911

set.seed(2025)
XGB_True_IBD <- XGBoost_Pipeline_Curve_2(IBD_Keep, 'Group', start=50, size=1500, steps=10, True=Keep, Test=Test_Keep)

set.seed(2025)
RF_True_IBD <- RF_Pipeline_Curve_2(IBD_Keep, 'Group', start=50, size=1500, steps=10, True=Keep, Test=Test_Keep)

LR_Optimal_2(IBD_Keep, "Group", Test_Keep)

NN_Optimal_2(IBD_Keep, "Group", Test_Keep, h=100)
mean(c(0.935, 0.920, 0.931, 0.931, 0.931, 0.930, 0.931, 0.919, 0.928, 0.933)) #0.929

IBD_NN <- NULL
for (i in round(seq(50, 1500, by=1450/10),0)) {
  IBD_NN <- rbind(IBD_NN, NN_Curve_Evaluate_2(IBD_Keep, "Group", n=i, True=colnames(IBD_Keep)[-ncol(IBD_Keep)], Test=Test_Keep, h=100))
}

IBD_NN<-read.csv("IBD_NN_True.csv")[,-1]
AUC_N_Fixed(IBD_NN, 0.929, 779)

plot_curve_full(IBD_NN_True, "PowerLaw_Fixed", 0.929, 1500)
Get_N_Fixed(IBD_NN_True, 0.929, 1500) #905

Simulate_OMICs_Data_BNG_Only(IBD_Keep, "Group", Keep=colnames(IBD_Keep)[-58], transform=F, fast=F, Case_Rows = as.numeric(rownames(IBD_Keep[which(IBD_Keep$Group==1),])),
                             p=57, n=2000)


#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
#write.csv(IBD_NN, "IBD_NN_True.csv")

#Learning Curves for Each Dataset
#IBD_Keep is true
#IBD_Sim_2000 is simulated
#both evaluated on 471 observation test set

plot_curve_full(XGB_True_IBD, "PowerLaw", 0.917, 1500)
Get_N(XGB_True_IBD, 0.917, 1500) #836
AUC_N_Fixed(read.csv("IBD_NN_True.csv")[,-1], 0.929, 736)
AUC_N_Fixed(XGB_True_IBD, 0.917, 25) #0.737


plot_curve_full(RF_True_IBD, "PowerLaw", 0.911, 1500)
Get_N(RF_True_IBD, 0.911, 1500) #367
AUC_N(RF_True_IBD, 291)
AUC_N(RF_True_IBD, 25) #0.797

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
#write.csv(XGB_True_IBD, "IBD_XGB.csv")
#write.csv(RF_True_IBD, "IBD_RF.csv")
library(ciTools)

XGB_True_IBD <- read.csv("IBD_XGB.csv")[,-1]
RF_True_IBD <- read.csv("IBD_RF.csv")[,-1]
IBD_NN_True <- read.csv("IBD_NN_True.csv")[,-1]

Mod_XGB <- glm.nb(data=Dataset_TEST[-24,], N ~ Med_LFC + Min_LFC + AUC_Full + log(Imbalance))
summary(Mod_XGB)
adjR2(Mod_XGB)
exp(predict(Mod_XGB, newdata=data.frame("Min_LFC"=0.17, "Med_LFC"=0.90, "AUC_Full"=0, "Imbalance"=18.1)))

#IBD XGB
DF_new <- Dataset_TEST[-24, c("Min_LFC", "Med_LFC", "AUC_Full","Imbalance","N")]

DF_new <- rbind(DF_new, data.frame("Min_LFC"=0.17, "Med_LFC"=0.90, "AUC_Full"=0, "Imbalance"=18.1, "N"=exp(predict(Mod_XGB, newdata=data.frame("Min_LFC"=0.17, "Med_LFC"=0.90, "AUC_Full"=0, "Imbalance"=18.1)))))

CI <- add_ci(DF_new, Mod_XGB)
CI[nrow(CI),c(6,7,8)] 

#PDAC XGB
DF_new <- Dataset_TEST[-24, c("Min_LFC", "Med_LFC", "AUC_Full","Imbalance","N")]

DF_new <- rbind(DF_new, data.frame("Min_LFC"=0.01, "Med_LFC"=0.87, "AUC_Full"=1, "Imbalance"=30.1, "N"=exp(predict(Mod_XGB, newdata=data.frame("Min_LFC"=0.01, "Med_LFC"=0.87, "AUC_Full"=1, "Imbalance"=30.1)))))

CI <- add_ci(DF_new, Mod_XGB)
CI[nrow(CI),c(6,7,8)] 



Mod_RF <- glm.nb(data=Dataset_RF_TEST[-24,], N_RF ~ log(Max_LFC) + Imbalance + Complexity + I(Dispersion>1))
summary(Mod_RF)
adjR2(Mod_RF)

Mod_RF <- glm.nb(data=Dataset_RF_TEST, N_RF ~ log(Max_LFC) + Imbalance + Complexity + I(Dispersion>1))
exp(predict(Mod_RF, newdata=data.frame("Max_LFC"=2.63, "Imbalance"=18.1, "Complexity"=0, "Dispersion"=1.07)))
exp(predict(glm.nb(N_RF ~ 1), newdata=data.frame("Complexity"=0, "Imbalance"=18.1, "AUC_Full_NN"=0))) #191 using median


#IBD RF
DF_new <- Dataset_RF_TEST[-24, c("Max_LFC", "Complexity", "Dispersion","Imbalance","N_RF")]

DF_new <- rbind(DF_new, data.frame("Max_LFC"=2.63, "Imbalance"=18.1, "Complexity"=0, "Dispersion"=1.07, "N_RF"=exp(predict(Mod_RF, newdata=data.frame("Max_LFC"=2.63, "Imbalance"=18.1, "Complexity"=0, "Dispersion"=1.07)))))

CI <- add_ci(DF_new, Mod_RF)
CI[nrow(CI),c(6,7,8)] 

#PDAC RF
DF_new <- Dataset_RF_TEST[-24, c("Max_LFC", "Complexity", "Dispersion","Imbalance","N_RF")]

DF_new <- rbind(DF_new, data.frame("Max_LFC"=5.58, "Imbalance"=30.1, "Complexity"=0, "Dispersion"=0.89, "N_RF"=exp(predict(Mod_RF, newdata=data.frame("Max_LFC"=5.58, "Imbalance"=30.1, "Complexity"=0, "Dispersion"=0.89)))))

CI <- add_ci(DF_new, Mod_RF)
CI[nrow(CI),c(6,7,8)] 


#IBD NN
DF_new <- Dataset_NN_TEST[-24, c("AUC_Full_NN", "Complexity","Imbalance","N_NN")]

DF_new <- rbind(DF_new, data.frame("Complexity"=0, "Imbalance"=18.1, "AUC_Full_NN"=0, "N_NN"=exp(predict(Mod_NN, newdata=data.frame("Complexity"=0, "Imbalance"=18.1, "AUC_Full_NN"=0)))))

CI <- add_ci(DF_new, Mod_NN)
CI[nrow(CI),c(5,6,7)] 

#PDAC NN
DF_new <- Dataset_NN_TEST[-24, c("AUC_Full_NN", "Complexity","Imbalance","N_NN")]

DF_new <- rbind(DF_new, data.frame("AUC_Full_NN"=1, "Complexity"=0, "Imbalance"=30.1, "N_NN"=exp(predict(Mod_NN, newdata=data.frame("AUC_Full_NN"=1, "Complexity"=0, "Imbalance"=30.1)))))

CI <- add_ci(DF_new, Mod_NN)
CI[nrow(CI),c(5,6,7)] 





Mod_NN <- glm.nb(data=Dataset_NN_TEST[-24,], N_NN ~ Complexity + AUC_Full_NN + Imbalance)
summary(Mod_NN)
adjR2(Mod_NN)
exp(predict(Mod_NN, newdata=data.frame("Complexity"=0, "Imbalance"=18.1, "AUC_Full_NN"=0))) #731
exp(predict(glm.nb(N_NN ~ 1), newdata=data.frame("Complexity"=0, "Imbalance"=18.1, "AUC_Full_NN"=0))) #324 using median




#MODEL ESTIMATES PLOTS
Dataset_TEST <- data.frame("N" = N, "N_Features"=N_Features, "Max_LFC"=Max_LFC, "AUC_Full"=AUC_Full, "Min_LFC"=Min_LFC, "Med_LFC"=Med_LFC, "Imbalance"=Imbalance, 
                           "Correlation"=Correlation, "Complexity"=Complexity, "AvgMedRead"=AvgMedRead, "Dispersion"=Dispersion_Simulated)

Dataset_RF_TEST <- data.frame("N_RF" = N_RF, "N_Features"=N_Features,"Max_LFC"=Max_LFC,"AUC_Full_RF"=AUC_Full_RF, "Min_LFC"=Min_LFC, "Med_LFC"=Med_LFC, "Imbalance"=Imbalance,
                              "Correlation"=Correlation, "Complexity"=Complexity, "AvgMedRead"=AvgMedRead, "Dispersion"=Dispersion_Simulated)

Dataset_NN_TEST <- data.frame("N_NN" = N_NN, "N_Features"=N_Features,"Max_LFC"=Max_LFC,"AUC_Full_NN"=as.factor(ifelse(AUC_Full_NN==1, ">=0.99", "<0.99")), "Min_LFC"=Min_LFC, "Med_LFC"=Med_LFC, "Imbalance"=Imbalance,
                              "Correlation"=Correlation, "Complexity"=Complexity, "AvgMedRead"=AvgMedRead, "Dispersion"=Dispersion_Simulated)



Dataset_TEST$Dispersion <- as.factor(ifelse(Dataset_TEST$Dispersion>1, "Dispersion >= 1", "Dispersion < 1"))
Dataset_RF_TEST$Dispersion <- as.factor(ifelse(Dataset_RF_TEST$Dispersion>1, "Dispersion >= 1", "Dispersion < 1"))
Dataset_NN_TEST$Dispersion <- as.factor(ifelse(Dataset_NN_TEST$Dispersion>1, "Dispersion >= 1", "Dispersion < 1"))
Mod_XGB <- glm.nb(data=Dataset_TEST, N ~ Med_LFC + Min_LFC + Dispersion + AUC_Full)
Mod_RF <- glm.nb(data=Dataset_RF_TEST, N_RF ~ Complexity + Dispersion + Imbalance + Max_LFC)
Mod_NN <- glm.nb(data=Dataset_NN_TEST, N_NN ~ Complexity + Dispersion + AUC_Full_NN + Max_LFC)

library(ggeffects)
library(ggpubr)
A <- ggpredict(Mod_XGB, terms=c("Med_LFC [all]", "AUC_Full [all]", "Dispersion [all]")) %>% plot() + labs(title="XGBoost Model", x="Median LFC", y="Predicted Sample Size") + scale_color_manual(labels=c("Full Dataset AUC < 0.99", "Full Dataset AUC >= 0.99"), values=c("red", "blue")) + theme(legend.title = element_blank())
B <- ggpredict(Mod_XGB, terms=c("Min_LFC [all]", "AUC_Full [all]", "Dispersion [all]")) %>% plot() + labs(title="", x="Minimum LFC", y="Predicted Sample Size") + scale_color_manual(labels=c("Full Dataset AUC < 0.99", "Full Dataset AUC >= 0.99"), values=c("red", "blue")) + theme(legend.title = element_blank())

C <- ggpredict(Mod_RF, terms=c("Max_LFC [all]", "Complexity [all]", "Dispersion [all]")) %>% plot() + labs(title="Random Forest Model", x="Maximum LFC", y="Predicted Sample Size") + scale_color_manual(labels=c("Dataset Nonlinearity < 4.5", "Dataset Nonlinearity >= 4.5"), values=c("red", "blue")) + theme(legend.title = element_blank())
D <- ggpredict(Mod_RF, terms=c("Imbalance [all]", "Complexity [all]", "Dispersion [all]")) %>% plot() + labs(title="", x="Minority Class Proportion", y="Predicted Sample Size") + scale_color_manual(labels=c("Dataset Nonlinearity < 4.5", "Dataset Nonlinearity >= 4.5"), values=c("red", "blue")) + theme(legend.title = element_blank())

E <- ggpredict(Mod_NN, terms=c("Max_LFC [all]", "Complexity [all]", "Dispersion [all]")) %>% plot() + labs(title="Neural Network Model", x="Maximum LFC", y="Predicted Sample Size") + scale_color_manual(labels=c("Dataset Nonlinearity < 4.5", "Dataset Nonlinearity >= 4.5"), values=c("red", "blue")) + theme(legend.title = element_blank())
EF <- ggpredict(Mod_NN, terms=c("AUC_Full_NN [all]", "Complexity [all]", "Dispersion [all]")) %>% plot() + labs(title="", x="Separability", y="Predicted Sample Size") + scale_color_manual(labels=c("Dataset Nonlinearity < 4.5", "Dataset Nonlinearity >= 4.5"), values=c("red", "blue")) + theme(legend.title = element_blank()) 


BigA <- ggarrange(A,B, common.legend = T) + theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))
BigB <-  ggarrange(C,D, common.legend = T) + theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))
BigC <- ggarrange(E,EF, common.legend = T) + theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))

ggarrange(BigA, BigB, BigC, nrow=3, ncol=1, common.legend = F)



#GET DISPERSIONS
setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2")

IBD <- read.csv("GSE193677_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE193677")
IBD <- IBD[,-1]
IBD <- data.frame(t(IBD))

Meta <- data.frame("ID"=Meta$GSE193677_series_matrix.txt.gz$geo_accession,
                   "Group"=ifelse(Meta$GSE193677_series_matrix.txt.gz$`ibd_disease:ch1`=="Control", 0, 1))

IBD$ID <- rownames(IBD)

IBD <- left_join(IBD, Meta, "ID")
IBD <- IBD[,-39377]

Vars <- apply(IBD[,-ncol(IBD)], 2, var)
IBD <- IBD[,-which(colnames(IBD) %in% names(Vars[Vars==0]))]
IBD <- data.frame(apply(IBD, 2, as.numeric))

median(GetDispersion(IBD, F, Keep=Keep)) #1.07
min(GetDispersion(IBD, F, Keep=Keep))

Keep <- c("X447","X940"  , "X1284" , "X2783"  ,"X2849",  "X3217" , "X3912", "X4019" , "X4481" , "X5289" , "X5294",  "X5549" , "X5794" , "X5795" , "X5797" ,
"X5798" , "X7581",  "X7622",  "X7819" , "X8161" , "X8936" , "X9448" , "X9463",  "X10148" ,"X10149", "X10749" ,"X13451" ,"X15046" ,"X15203", "X15674",
"X15732" ,"X15742", "X16349", "X17957", "X18100", "X19857" ,"X23138" ,"X25581", "X26663", "X26769", "X28159" ,"X28160" ,"X28161", "X28349" ,"X30000",
"X31298", "X32463", "X32810" ,"X32834" ,"X34032", "X34876" ,"X34878", "X36488", "X36712" ,"X37026" ,"X38336", "X38821")

  
#Try this out
setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2")

PDAC <- read.csv("GSE133684_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")
Meta <- getGEO("GSE133684")

PDAC  <- PDAC[,-1]
PDAC  <- data.frame(t(PDAC))

Meta$GSE133684_series_matrix.txt.gz$`disease state:ch1`

Meta <- data.frame("ID"=c(Meta$GSE133684_series_matrix.txt.gz$geo_accession),
                   "Group"=ifelse(Meta$GSE133684_series_matrix.txt.gz$`disease state:ch1`=="healthy",0,1))

PDAC$ID <- rownames(PDAC)

PDAC <- left_join(PDAC, Meta, "ID")
PDAC <- PDAC[,-39377]

PDAC <- na.omit(PDAC)

Vars <- apply(PDAC[,-ncol(PDAC)], 2, var)
PDAC <- PDAC[,-which(colnames(PDAC) %in% names(Vars[Vars==0]))]
PDAC <- data.frame(apply(PDAC, 2, as.numeric))




set.seed(2024)
Test <- PDAC[sample(1:nrow(PDAC), 39), ]
PDAC <- PDAC[-which(rownames(PDAC) %in% rownames(Test)),]


ExpDesign <- data.frame(row.names=rownames(PDAC),
                        Group = PDAC[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(PDAC[,-which(colnames(PDAC) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

PDAC[,1:(ncol(PDAC)-1)] <- DDS

#Test Counts
ExpDesign <- data.frame(row.names=rownames(Test),
                        Group = Test[,"Group"])

DiffEx <- DESeq(DESeqDataSetFromMatrix(countData = as.matrix(t(Test[,-which(colnames(Test) %in% "Group")])), 
                                       colData=ExpDesign, design=~Group))



DDS <- estimateSizeFactors(DiffEx)
DDS <- data.frame(t(counts(DDS, normalized=TRUE)))

Test[,1:(ncol(Test)-1)] <- DDS




set.seed(2024)
B <- Boruta(x=PDAC[,-which(colnames(PDAC) %in% "Group")], y=as.factor(PDAC[,which(colnames(PDAC) %in% "Group")]))

B <- data.frame(B$finalDecision)

B$names <- rownames(B)

B[B$B.finalDecision=="Confirmed", "names"]

DiffEx <- results(DiffEx, alpha=0.05)
DiffEx <- na.omit(DiffEx)

LFC <- abs(DiffEx[DiffEx$padj<0.05,"log2FoldChange"])
rownames(DiffEx[DiffEx$padj<0.05,])

Keep <- intersect(rownames(DiffEx[DiffEx$padj<0.05,]), B[B$B.finalDecision=="Confirmed", "names"])

PDAC_Keep <- PDAC[,c(Keep, "Group")]

PDAC_Keep <- cbind(log2(PDAC_Keep[,-ncol(PDAC_Keep)]+0.01), "Group"=PDAC_Keep$Group)

Test_Keep <- Test[,c(Keep, "Group")]

Test_Keep <- cbind(log2(Test_Keep[,-ncol(Test_Keep)]+0.01), "Group"=Test_Keep$Group)

XGBoost_Optimal(PDAC_Keep, "Group", Test_Keep) #>0.99

RF_Optimal_2(PDAC_Keep, "Group", Test_Keep) #>0.99

LR_Optimal_2(PDAC_Keep, "Group", Test_Keep) #0.987

1-mean(PDAC_Keep$Group) #30.3%
LFC_True <- sapply(1:(ncol(PDAC_Keep)-1), function(i) {abs(mean(PDAC_Keep[PDAC_Keep$Group==1,i]) - mean(PDAC_Keep[PDAC_Keep$Group==0,i]))})
max(LFC_True) #5.58
median(LFC_True) #0.87
min(LFC_True) #0.01

median(LFC_True[-which(LFC_True<0.5)])
min(LFC_True[-which(LFC_True<0.5)])
max(LFC_True[-which(LFC_True<0.5)])


quantile(GetDispersion(PDAC, F, colnames(PDAC_Keep)[1:(length(colnames(PDAC_Keep))-1)]))

median(GetDispersion(PDAC, F, colnames(PDAC_Keep)[1:(length(colnames(PDAC_Keep))-1)]), LFC_True)

ssizeRNA_single(nGenes = 39354,
                pi0= 1 - (57/39354), 
                m=200,
                mu=39.48,
                disp=0.958,
                fc=2^0.90,
                up=1,
                fdr=0.05,
                power=0.8,
                maxN = 2000,
                side="two-sided")$ssize[2] #192 for DIFFEX

set.seed(1)
XGB_True_PDAC <- XGBoost_Pipeline_Curve_2(PDAC_Keep, 'Group', start=50, size=290, steps=10, True=Keep, Test=Test_Keep)

XGB_True_PDAC

#PLOTS - MODEL ESTIMATES
#PLOTS of AUC 25 vs AUC FULL

#setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
#write.csv(XGB_True_PDAC, "PDAC_XGB.csv")
#write.csv(RF_True_PDAC, "PDAC_RF.csv")

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")
XGB_True_PDAC <- read.csv("PDAC_XGB.csv")[,-1]
RF_True_PDAC <- read.csv("PDAC_RF.csv")[,-1]
PDAC_NN_True <- read.csv("PDAC_NN_True.csv")[,-1]

Get_N_Fixed(XGB_True_PDAC, 1, 290) #176
plot_curve_full(XGB_True_PDAC, "PowerLaw", FullAUC = 1, 300)
AUC_N(XGB_True_PDAC, 136)
AUC_N(XGB_True_PDAC, 25) #0.785

Mod_XGB <- glm.nb(data=Dataset_TEST[-24,], N ~ Med_LFC + Min_LFC + AUC_Full + log(Imbalance))
summary(Mod_XGB)
exp(predict(Mod_XGB, newdata=data.frame("Min_LFC"=0.01, "Med_LFC"=0.87, "AUC_Full"=1, "Imbalance"=30.1)))
exp(predict(glm.nb(N ~ 1), newdata=data.frame("Complexity"=0, "Imbalance"=18.1, "AUC_Full_NN"=0))) #480 using median

set.seed(2025)
RF_True_PDAC <- RF_Pipeline_Curve_2(PDAC_Keep, 'Group', start=50, size=290, steps=10, True=Keep, Test=Test_Keep)

RF_True_PDAC

Get_N(RF_True_PDAC, 1, 290) #80
plot_curve_full(RF_True_PDAC, "PowerLaw", FullAUC = 1, 290)
AUC_N(RF_True_PDAC, 63)
AUC_N(RF_True_PDAC, 25)

Mod_RF <- glm.nb(data=Dataset_RF_TEST[-24,], N_RF ~ log(Max_LFC) + Imbalance + Complexity + I(Dispersion>1))
exp(predict(Mod_RF, newdata=data.frame("Max_LFC"=5.58, "Imbalance"=30.1, "Complexity"=0, "Dispersion"=0.89))) #55

median(GetDispersion(PDAC, F, colnames(PDAC_Keep)[1:(length(colnames(PDAC_Keep))-1)])) #0.89

exp(predict(Mod_NN, newdata=data.frame("AUC_Full_NN"=1, "Complexity"=0, "Imbalance"=30.1))) #164


NN_Optimal_2(PDAC_Keep, "Group", Test_Keep, h=100)
#1

PDAC_NN <- NULL
for (i in round(seq(50, 290, by=290/10),0)) {
  PDAC_NN <- rbind(PDAC_NN, NN_Curve_Evaluate_2(PDAC_Keep, "Group", n=i, True=colnames(PDAC_Keep)[-ncol(PDAC_Keep)], Test=Test_Keep, h=100))
}
PDAC_NN

plot_curve_full(PDAC_NN, "PowerLaw_Fixed", 1, 300)
Get_N_Fixed(PDAC_NN, 1, 300) #113

AUC_N_Fixed(PDAC_NN_True, 1, 70)
AUC_N_Fixed(PDAC_XGB, 1, 70)
AUC_N(RF_True_PDAC, 70)

AUC_N(XGB_True_IBD, 384)
AUC_N(RF_True_IBD, 384)
AUC_N_Fixed(IBD_NN, 0.929, 384)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Aim 2/Learning Curve Backup Data")

median(LFC_True)

#write.csv(PDAC_NN, "PDAC_NN_True.csv")
PDAC_NN <- read.csv("PDAC_NN_True.csv")

#At the median value
ssizeRNA_single(nGenes = 35923,
                pi0= 1 - (67/35923), 
                m=200,
                mu=13069.69,
                disp=0.209,
                fc=2^0.87,
                up=1,
                fdr=0.05,
                power=0.8,
                maxN = 2000,
                side="two-sided")$ssize[2] #67 for DIFFEX

#5465.65

LFC_True[29]
GetDispersion(PDAC, F, colnames(PDAC_Keep)[1:(length(colnames(PDAC_Keep))-1)])[29]

colnames(PDAC_Keep)[29]

GetDispersion(IBD, F, colnames(IBD_Keep)[1:(length(colnames(IBD_Keep))-1)])[48]

AvgReadCount(cbind(2^(IBD_Keep[,-ncol(IBD_Keep)]-0.01), "Group"=IBD_Keep$Group))[48]

AvgReadCount(cbind(2^(PDAC_Keep[,-ncol(PDAC_Keep)]-0.01), "Group"=PDAC_Keep$Group))[29]

LFCs_True <- rbind(
GetLFC(GBM, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)]),
GetLFC(HCC, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)]),
GetLFC(IPF, colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)]),
GetLFC(Kidney, colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)]),
GetLFC(COVID, colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)]),
GetLFC(HBV, colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)]),
GetLFC(Hypertension, colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)]),
GetLFC(NAFLD, colnames(NAFLD_Sim)[1:(length(colnames(NAFLD_Sim))-1)]),
GetLFC(PrePost, colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)]),
GetLFC(RA, colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)]) ,
GetLFC(AVSC, colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)]) ,
GetLFC(Crohn, colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)]),
GetLFC(MDD, colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)]),
GetLFC(Ovarian, colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)]),
GetLFC(ALS, colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)]),
GetLFC(MISC, colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)]),
GetLFC(CCA, colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)]),
GetLFC(Glioma, colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)]),
GetLFC(TB, colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)]),
GetLFC(Tub, colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)]),
GetLFC(EDS, colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)]),
GetLFC(NSCLC, colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)]),
GetLFC(Bipolar, colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)]),
GetLFC(MS, colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])
)

plot(LFCs_True[,1], Med_LFC[-24])

Med_LFC_DF <- data.frame("Dataset"=DS_Names, "Median_LFC_True"=LFCs_True[,1], "Median_LFC_Sim"=Med_LFC[-24])
Min_LFC_DF <- data.frame("Dataset"=DS_Names, "Min_LFC_True"=LFCs_True[,2], "Min_LFC_Sim"=Min_LFC[-24])
Max_LFC_DF <- data.frame("Dataset"=DS_Names, "Max_LFC_True"=LFCs_True[,3], "Max_LFC_Sim"=Max_LFC[-24])
Dispersion_DF <- data.frame("Dataset"=DS_Names, "Dispersion_True"=Dispersion, "Dispersion_Sim"=Dispersion_Simulated[1:24])
Correlation_DF <- data.frame("Dataset"=DS_Names, "Correlation_True"=Cors_True, "Correlation_Sim"=Correlation[1:24])


ggplot(data=Med_LFC_DF, aes(x=Median_LFC_True, y=Median_LFC_Sim)) + geom_point(aes(color=Dataset), size=3, shape="triangle") + geom_smooth(method="lm", se=F, linetype="dashed", color="grey") + theme_bw() + theme(panel.grid.major = element_blank()) + labs(x="True Median Log-Fold Change", y="Simulated Median Log-Fold Change")

ggplot(data=Max_LFC_DF, aes(x=Max_LFC_True, y=Max_LFC_Sim)) + geom_point(aes(color=Dataset), size=3, shape="triangle") + geom_smooth(method="lm", se=F, linetype="dashed", color="grey") + theme_bw() + theme(panel.grid.major = element_blank()) + labs(x="True Max. Log-Fold Change", y="Simulated Max. Log-Fold Change")

ggplot(data=Min_LFC_DF, aes(x=Min_LFC_True, y=Min_LFC_Sim)) + geom_point(aes(color=Dataset), size=3, shape="triangle") + geom_smooth(method="lm", se=F, linetype="dashed", color="grey") + theme_bw() + theme(panel.grid.major = element_blank()) + labs(x="True Min. Log-Fold Change", y="Simulated Min. Log-Fold Change")

ggplot(data=Correlation_DF, aes(x=Correlation_True, y=Correlation_Sim)) + geom_point(aes(color=Dataset), size=3, shape="triangle") + geom_smooth(method="lm", se=F, linetype="dashed", color="grey") + theme_bw() + theme(panel.grid.major = element_blank()) + labs(x="True Median Between-Gene Correlation", y="Simulated Median Between-Gene Correlation")

ggplot(data=Dispersion_DF, aes(x=Dispersion_True, y=Dispersion_Sim)) + geom_point(aes(color=Dataset), size=3, shape="triangle") + geom_smooth(method="lm", se=F, linetype="dashed", color="grey") + theme_bw() + theme(panel.grid.major = element_blank()) + labs(x="True Median Dispersion", y="Simulated Median Dispersion")
cor(Dispersion_DF$Dispersion_True, Dispersion_DF$Dispersion_Sim) #r=0.939
cor(Correlation_DF$Correlation_True, Correlation_DF$Correlation_Sim) #r=1

#Things I need to do
#Before submitting change fig 2 back to remove those stupid stars
#Also remove sentences about stars BS from paper
#Training AUC in original DS as comparitor? Try this out for actual dissertation
#Add example compariosn of correlations between BNG data and original to supplementary

Cors_True <- 
c(
median(GetCorrelation(GBM_Sim, colnames(GBM_Sim)[1:(length(colnames(GBM_Sim))-1)])) ,
median(GetCorrelation(HCC_Sim, colnames(HCC_Sim)[1:(length(colnames(HCC_Sim))-1)])) ,
median(GetCorrelation(Kidney_Sim, colnames(Kidney_Sim)[1:(length(colnames(Kidney_Sim))-1)])) ,
median(GetCorrelation(IPF_Sim, colnames(IPF_Sim)[1:(length(colnames(IPF_Sim))-1)])) ,
median(GetCorrelation(COVID_Sim, colnames(COVID_Sim)[1:(length(colnames(COVID_Sim))-1)])),
median(GetCorrelation(HBV_Sim, colnames(HBV_Sim)[1:(length(colnames(HBV_Sim))-1)])) ,
median(GetCorrelation(Hypertension_Sim, colnames(Hypertension_Sim)[1:(length(colnames(Hypertension_Sim))-1)])),
median(GetCorrelation(NAFLD_Sim, colnames(NAFLD_Sim)[1:52])) ,
median(GetCorrelation(PrePost_Sim, colnames(PrePost_Sim)[1:(length(colnames(PrePost_Sim))-1)])) ,
median(GetCorrelation(RA_Sim, colnames(RA_Sim)[1:(length(colnames(RA_Sim))-1)])) ,
median(GetCorrelation(AVSC_Sim, colnames(AVSC_Sim)[1:(length(colnames(AVSC_Sim))-1)])) ,
median(GetCorrelation(Crohn_Sim, colnames(Crohn_Sim)[1:(length(colnames(Crohn_Sim))-1)])) ,
median(GetCorrelation(MDD_Sim, colnames(MDD_Sim)[1:(length(colnames(MDD_Sim))-1)])) ,
median(GetCorrelation(Ovarian_Sim, colnames(Ovarian_Sim)[1:(length(colnames(Ovarian_Sim))-1)])) ,
median(GetCorrelation(ALS_Sim, colnames(ALS_Sim)[1:(length(colnames(ALS_Sim))-1)])) ,
median(GetCorrelation(MISC_Sim, colnames(MISC_Sim)[1:(length(colnames(MISC_Sim))-1)])) ,
median(GetCorrelation(CCA_Sim, colnames(CCA_Sim)[1:(length(colnames(CCA_Sim))-1)])),
median(GetCorrelation(Glioma_Sim, colnames(Glioma_Sim)[1:(length(colnames(Glioma_Sim))-1)])) ,
median(GetCorrelation(TB_Sim, colnames(TB_Sim)[1:(length(colnames(TB_Sim))-1)])) ,
median(GetCorrelation(Tuberculosis_Sim, colnames(Tuberculosis_Sim)[1:(length(colnames(Tuberculosis_Sim))-1)])), 
median(GetCorrelation(EDS_Sim, colnames(EDS_Sim)[1:(length(colnames(EDS_Sim))-1)])),
median(GetCorrelation(NSCLC_Sim, colnames(NSCLC_Sim)[1:(length(colnames(NSCLC_Sim))-1)])) ,
median(GetCorrelation(Bipolar_Sim, colnames(Bipolar_Sim)[1:(length(colnames(Bipolar_Sim))-1)])) ,
median(GetCorrelation(MS_Sim, colnames(MS_Sim)[1:(length(colnames(MS_Sim))-1)])) )

cor.test(Cors_True, Correlation[-24])

#0.01
Get_N_Th(XGB_True_PDAC, 1, 2000, 0.01) #217
#App predicts 258 [198,335]

Get_N_Th(RF_True_PDAC, 1, 2000, 0.01) #123
#App predicts 132 [78,224]

Get_N_Fixed_Th(PDAC_NN, 1, 2000, 0.01) #174
#App predicts 122 [69,217]

Get_N_Th(XGB_True_IBD, 0.917, 2000, 0.01) #1549
#App predicts 1454 [1036,2041]

Get_N_Th(RF_True_IBD, 0.911, 2000, 0.01) #709
#App predicts 877 [407,1887]

Get_N_Fixed_Th(IBD_NN_True, 0.929, 2500, 0.01) #2108
#App predicts 1944 [860,4396]

#0.05
Get_N_Th(XGB_True_PDAC, 1, 2000, 0.05) #100
#App predicts 52 [38,72]

Get_N_Th(RF_True_PDAC, 1, 2000, 0.05) #41
#App predicts 35 [25,48]

Get_N_Fixed_Th(PDAC_NN, 1, 2000, 0.05) #64
#App predicts 35 [25,50]

Get_N_Th(XGB_True_IBD, 0.917, 2000, 0.05) #238
#App predicts 55 [159,361]

Get_N_Th(RF_True_IBD, 0.911, 2000, 0.05) #105
#App predicts 1454 [35,87]

Get_N_Fixed_Th(IBD_NN_True, 0.929, 2500, 0.05) #296
#App predicts 72 [44,118]













