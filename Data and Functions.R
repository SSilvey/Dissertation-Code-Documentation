#Set WD

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Data")

library(pathviewr)
library(randomForest)
library(rpart)
library(dplyr)
library(ggplot2)
library(splitTools)
library(xgboost)
library(pROC)
library(stats)
library(ranger)
library(glmnet)
library(caret)

#Data Cleaning and Preparation

#Diabetes
Diabetes <- read.csv("diabetes_binary_health_indicators_BRFSS2015.csv")

str(Diabetes)

mean(Diabetes$Diabetes_binary) #Outcome : 13.9% 

sum(Diabetes$Diabetes_binary)

#CDC Heart Disease 2020

CDC_2020 <- read.csv("heart_2020_cleaned.csv")

CDC_2020$HeartDisease <- ifelse(CDC_2020$HeartDisease=="Yes", 1, 0)

CDC_2020$Smoking <- ifelse(CDC_2020$Smoking=="Yes", 1, 0)
CDC_2020$AlcoholDrinking <- ifelse(CDC_2020$AlcoholDrinking=="Yes", 1, 0)
CDC_2020$Stroke <- ifelse(CDC_2020$Stroke=="Yes", 1, 0)
CDC_2020$DiffWalking <- ifelse(CDC_2020$DiffWalking=="Yes", 1, 0)
CDC_2020$Sex <- ifelse(CDC_2020$Sex=="Male", 1, 0)
CDC_2020$AgeCategory <- ifelse(CDC_2020$AgeCategory %in% c("50-54", "55-59", "60-64", "65-69", "70-74","75-79","80 or older"), 1, 0)
CDC_2020$Race <- ifelse(CDC_2020$Race=="White", 1, 0)
CDC_2020$Diabetic <- ifelse(CDC_2020$Diabetic=="Yes", 1, 0)
CDC_2020$PhysicalActivity <- ifelse(CDC_2020$PhysicalActivity=="Yes", 1, 0)
CDC_2020$GenHealth <- ifelse(CDC_2020$GenHealth %in% c("Good", "Very Good", "Excellent"), 1, 0)
CDC_2020$Asthma <- ifelse(CDC_2020$Asthma=="Yes", 1, 0)
CDC_2020$KidneyDisease <- ifelse(CDC_2020$KidneyDisease=="Yes", 1, 0)
CDC_2020$SkinCancer <- ifelse(CDC_2020$SkinCancer=="Yes", 1, 0)

mean(CDC_2020$HeartDisease) #8.6%

#Artifical 1

HeartSynthetic <- read.csv("csv_result-BNG_heart-statlog.csv")

HeartSynthetic_Y <- ifelse(HeartSynthetic$class=='present',1,0)

HeartSynthetic$class <- ifelse(HeartSynthetic$class=='present',1,0)

HeartSynthetic <- HeartSynthetic[,-which(colnames(HeartSynthetic) %in% c("id"))]

mean(HeartSynthetic_Y) #44.4%

#Artificial 2

Hepatitis <- read.csv("csv_result-BNG_hepatitis.csv")

Hepatitis_Y <- ifelse(Hepatitis$Class=="LIVE",0,1)

Hepatitis$Class <- ifelse(Hepatitis$Class=="LIVE",0,1)

Hepatitis <- Hepatitis[,-which(colnames(Hepatitis) %in% c("id"))]

Hepatitis$SEX <- ifelse(Hepatitis$SEX=="male",1,0)
Hepatitis[Hepatitis=="yes"] <- 1
Hepatitis[Hepatitis=="no"] <- 0

Hepatitis <- data.frame(apply(Hepatitis, 2, as.numeric))

mean(Hepatitis$Class) #20.8%

#Artificial 3

Lymph <- read.csv("csv_result-BNG_lymph.csv")

Lymph <- Lymph[,-which(colnames(Lymph) %in% c("id"))]

Lymph$class <- ifelse(Lymph$class == "metastases", 0, 1)

Lymph$lymphatics <- ifelse(Lymph$lymphatics == "normal", 0, 1)
Lymph$changes_in_lym <- ifelse(Lymph$changes_in_lym == "round", 0, 1)
Lymph$defect_in_node <- ifelse(Lymph$defect_in_node=="no",0,1)
Lymph$changes_in_node <- ifelse(Lymph$changes_in_node=="no",0,1)
Lymph$changes_in_stru <- ifelse(Lymph$changes_in_stru=="no",0,1)
Lymph$special_forms <- ifelse(Lymph$special_forms=="no",0,1)

Lymph[Lymph=="yes"] <- 1
Lymph[Lymph=="no"] <- 0

Lymph <- data.frame(apply(Lymph, 2, as.numeric))

mean(Lymph$class) #45.7


#Artificial 4

Pharynx <- read.csv('csv_result-BNG_pharynx.csv')

Pharynx <- Pharynx[,-which(colnames(Pharynx) %in% c("id"))]

Pharynx$class <- ifelse(Pharynx$class > 750, 1, 0)

mean(Pharynx$class) #25.6

Pharynx$sex <- Pharynx$sex - 1
Pharynx$Treatment <- Pharynx$Treatment - 1
Pharynx$Grade <- ifelse(Pharynx$Grade==1, 0, 1)
Pharynx$Condition <- ifelse(Pharynx$Condition == 1 | Pharynx$Condition == 0, 0, 1)
Pharynx$Site <- ifelse(Pharynx$Site == 4, 1, 0)
Pharynx$T <- ifelse(Pharynx$T == 1 | Pharynx$T == 2, 1, 0)
Pharynx$N <- ifelse(Pharynx$N == 0 | Pharynx$N == 1, 0, 1)

Pharynx <- data.frame(apply(Pharynx, 2, as.numeric))

#Validation Set 1
Chol <- read.csv("csv_result-BNG_Cholesterol.csv")

Chol <- Chol[,-which(colnames(Chol) %in% "id")]

Chol$chol <- ifelse(Chol$chol > 200, 0, 1)

mean(Chol$chol)


#Validation Set 2

CDC_2022 <- read.csv("heart_2022_with_nans.csv")

CDC_2022$State <- ifelse(CDC_2022$State %in% c("Connecticut", "Delaware", "District of Columbia",
                                               "Maine", "Maryland", "Massachusetts", "New Hampshire",
                                               "New Jersey", "New York", "Pennsylvania", "Rhode Island",
                                               "Virginia", "Vermont"), 1, 0)

CDC_2022$Sex <- ifelse(CDC_2022$Sex=="Male", 1, 0)
CDC_2022$GeneralHealth <- ifelse(CDC_2022$GeneralHealth == "Poor", 0, 1)
CDC_2022$LastCheckupTime <- ifelse(CDC_2022$LastCheckupTime == "Within past year (anytime less than 12 months ago)", 0, 1)
CDC_2022$RemovedTeeth <- ifelse(CDC_2022$RemovedTeeth=="None of them", 0, 1)
CDC_2022$SmokerStatus <- ifelse(CDC_2022$SmokerStatus=="Never smoked", 0, 1)
CDC_2022$ECigaretteUsage <- ifelse(CDC_2022$ECigaretteUsage == "Never used e-cigarettes in my entire life", 0, 1)
CDC_2022$RaceEthnicityCategory <- ifelse(CDC_2022$RaceEthnicityCategory=="White only, Non-Hispanic", 0,1)
CDC_2022$AgeCategory <- ifelse(CDC_2022$AgeCategory %in% c("Age 65 to 69", "Age 70 to 74", "Age 80 or older"),1,0)
CDC_2022$TetanusLast10Tdap <- ifelse(CDC_2022$TetanusLast10Tdap=="Yes, received Tdap",1,0)

CDC_2022[CDC_2022=="Yes"] <- 1
CDC_2022[CDC_2022=="No"] <- 0

CDC_2022 <- data.frame(apply(CDC_2022, 2, as.numeric))
CDC_2022 <- CDC_2022[-which(is.na(CDC_2022$HighRiskLastYear)==T),]

CDC_2022 <- na.roughfix(CDC_2022)

#Cardio
Cardio <- read.csv("Cardio.csv")

mean(Cardio$CVD) #50.0%

#Diabetes130
Diabetes130 <- read.csv("Diabetes130.csv")
Diabetes130 <- Diabetes130[,-which(colnames(Diabetes130) %in% c("id", "diag_1", "diag_2", "diag_3", "max_glu_serum", "A1Cresult",
                                                                "examide", "citaglipton"))]

Diabetes130[Diabetes130=="?"] <- NA

Diabetes130$race <- ifelse(Diabetes130$race=="Caucasian", 0, 1)
Diabetes130$gender <- ifelse(Diabetes130$gender=="Male", 1, 0)
Diabetes130$age <- ifelse(Diabetes130$age %in% c("[10-20)", "[0-10)"), 0, 1)
Diabetes130$insulin <- ifelse(Diabetes130$insulin=="Steady",0,1)
Diabetes130$change <- ifelse(Diabetes130$change=="Ch", 1, 0)
Diabetes130$diabetesMed <- ifelse(Diabetes130$diabetesMed=="Yes",1,0)

Diabetes130[Diabetes130=="Steady"] <- 1
Diabetes130[Diabetes130=="No"] <- 0

Diabetes130$readmitted <- ifelse(Diabetes130$readmitted=="NO",0,1)

Diabetes130 <- data.frame(apply(Diabetes130, 2, as.numeric))

length(Diabetes130[is.na(Diabetes130)==T]) / (nrow(Diabetes130)*ncol(Diabetes130))

length(Diabetes130[is.na(Diabetes130)==T]) / (nrow(Diabetes130)*ncol(Diabetes130))

length(CDC_2022[is.na(CDC_2022)==T]) / (nrow(CDC_2022)*ncol(CDC_2022))

Diabetes130 <- na.roughfix(Diabetes130)

mean(Diabetes130$readmitted) #46.1%

#NoShow

NoShow <- read.csv("NoShow.csv")

NoShow$Gender <- ifelse(NoShow$Gender=="M", 1, 0)
NoShow$NoShow <- ifelse(NoShow$NoShow=="No", 0, 1)

mean(NoShow$NoShow) #20.2%


#PBC
PBC <- read.csv('csv_result-BNG_pbc.csv')

PBC <- PBC[,-which(colnames(PBC) %in% "id")]

PBC$class <- ifelse(PBC$class > 3000, 1, 0)

mean(PBC$class) #17.8%


#Dermatology
Dermatology <- read.csv('csv_result-BNG_dermatology.csv')

Dermatology <- Dermatology[,-which(colnames(Dermatology) %in% c("id", "class"))]

mean(Dermatology$family_history) #13.2%


#BreastTumor
BreastTumor <- read.csv('csv_result-BNG_breastTumor.csv')

BreastTumor <- BreastTumor[,-which(colnames(BreastTumor) %in% "id")]

BreastTumor$menopause <- ifelse(BreastTumor$menopause=="premenopausal", 1, 0)
BreastTumor$irradiation <- ifelse(BreastTumor$irradiation=="yes", 1, 0)
BreastTumor$breast.quad <- ifelse(BreastTumor$breast.quad=="central", 1, 0)
BreastTumor$node.caps <- ifelse(BreastTumor$node.caps=="no",0,1)
BreastTumor$breast <- ifelse(BreastTumor$breast=="right", 1,0)
BreastTumor$recurrence <- ifelse(BreastTumor$recurrence=="r", 1, 0)

mean(BreastTumor$recurrence) #34.6%

#COVID Mexico

COVID <- read.csv("COVID19_Mexico.csv")

COVID[COVID==1] <- 0
COVID[COVID==2] <- 1

COVID$RESULTADO <- ifelse(COVID$RESULTADO==1, 0, 1)

mean(COVID$RESULTADO) #39.0%

COVID[COVID=="99"] <- NA
COVID[COVID=="97"] <- NA

length(COVID[is.na(COVID)==T]) / (nrow(COVID)*ncol(COVID))

COVID <- na.roughfix(COVID)

#LOS

LOS <- read.csv("COVID_LOS.csv")
LOS[LOS=="?"] <- NA

length(LOS[is.na(LOS)==T])/nrow(LOS)

LOS$HospitalType <- ifelse(LOS$HospitalType==1, 0, 1)
LOS$Department <- ifelse(LOS$Department=="gynecology", 0, 1)
LOS$WardType <- ifelse(LOS$WardType=="R", 0, 1)
LOS$WardFacility <- ifelse(LOS$WardFacility=="F", 0, 1)
LOS$BedGrade <- ifelse(LOS$BedGrade==2, 0, 1)
LOS$TypeOfAdmission <- ifelse(LOS$TypeOfAdmission=="Trauma", 0, 1)
LOS$Severity <- ifelse(LOS$Severity=="Extreme",1,0)
LOS$Age <- ifelse(LOS$Age %in% c("81-90", "91-100"), 1, 0)

LOS <- na.roughfix(LOS)

mean(LOS$LOS) #2.1%

#---------------------------------------------------------------------------------------------



#Learning curves

Subsample_DF_random <- function(DF, DF_Response, n) {
  
  Reduced_Rows <- NULL
  
  for (i in 1:10) {
    Reduced_Rows <- append(Reduced_Rows, list(c(sample(rownames(DF), size=n, replace = F))
    ))
  }
  Reduced_Rows
}


#Figure out :
#Increase n in the learning curves all the way up to the optimal one
#Maybe do random subsample with no class balance criteria 

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

XGBoost_Curve_Random <- function(DF, DF_Outcome) {
  
  a <- 500
  b <- 50000
  
  stepsize <- round((b - a) / 10)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > 10) {samps <- samps[1:10]}
  
  Curve <- NULL
  
  for (n in samps) {
    
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
    
    
  }
  Curve
}

XGBoost_Curve_Random_2 <- function(DF, DF_Outcome, size) {
  
  a <- 500
  b <- size
  
  stepsize <- round((b - a) / 10)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > 10) {samps <- samps[1:10]}
  
  Curve <- NULL
  
  for (n in samps) {
    
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
    
    
  }
  Curve
}

XGBoost_Optimal <- function(DF, DF_Outcome) {
  
  AUC <- NULL
  
  params <- list(booster="gbtree",
                 eta=0.3,
                 gamma=0,
                 max_depth=6,
                 lambda=0,
                 subsample=1)
  
  Folds_inv <- create_folds(as.factor(DF[,colnames(DF) %in% DF_Outcome]), k=5, type="stratified", invert = T)
  
  for (j in 1:length(Folds_inv)) {
    HoldOut <- Folds_inv[[j]]
    
    XGB_Example <- xgboost(data=as.matrix(DF[-HoldOut, -which(colnames(DF) %in% DF_Outcome)]),
                           label=DF[-HoldOut,colnames(DF) %in% DF_Outcome],
                           missing=NA, nrounds=25,
                           eval_metric="auc",verbose=F,
                           objective = "binary:logistic")
    
    Output <- predict(XGB_Example, newdata = as.matrix(DF[HoldOut, -which(colnames(DF) %in% DF_Outcome)]),
                      type="response")
    
    AUC <- append(AUC, as.numeric(roc(as.factor(DF[HoldOut,colnames(DF) %in% DF_Outcome]), Output, quiet=T)$auc))
  }
  c(mean(AUC), sd(AUC) / sqrt(5))
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

Get_Fitted_Values_PowerLaw <- function(DF, Size) {
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5, c=0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  c <- summary(test)$coef[3,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:Size)
  
  Values
  
}

Get_Fitted_Values_PowerLaw_Fixed <- function(DF, Size, FullAUC) {
  c <- FullAUC
  
  test <- nls(formula = Curve ~ a*(1/n^-b) + c, data = DF, start = c(a=-0.5,b=-0.5), control = nls.control(maxiter = 1000))
  
  a <- summary(test)$coef[1,1]
  b <- summary(test)$coef[2,1]
  Model <- function(a, b, c, n) {a*(1/n^-b) + c}
  
  Values <- Model(a=a, b=b, c=c, n=1:Size)
  
  Values
  
}

Get_Fitted_Values_LOG <- function(DF, Size) {
  test <- lm(data=DF, Curve ~ log(n))
  
  B0 <- test$coefficients[1]
  
  B1 <- test$coefficients[2]
  
  Model <- function(B0, B1, n) {B0 + B1*log(n)}
  
  Values <- data.frame("Size"=1:Size, "Fitted"=Model(B0=B0, B1=B1, n=1:Size))
  
  Values$Fitted
  
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


#Fix this

plot_curve <- function(Curve, Type, FullAUC, Size) {
  
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
  
  Plot <- ggplot(data=Data_forPlot_Raw, aes(x=N, y=Raw)) + geom_point(color="grey") +
    labs(x="Sample Size", y="AUC") + ylim(c(min(Data_forPlot_Raw$Raw)-0.01, FullAUC+0.01))
  
  
  Plot + geom_line(data=Values, aes(y=Fitted, x=Size), col="black") + annotate("Text", x=Marker$N,
                                                                               y=Marker$AUC, label="X", size=6)
  
  
}


LASSO_Optimal_CoreFeatures <- function(DF, DF_Outcome) {
  
  CoreFeatures <- NULL
  
  LASSO_Example <- glmnet(x=DF[,-which(colnames(DF) %in% DF_Outcome)], 
                          y=DF[,colnames(DF) %in% DF_Outcome], 
                          family = "binomial", alpha = 1, lambda = cv.glmnet(x=as.matrix(DF[,-which(colnames(DF) %in% DF_Outcome)]), 
                                                                             y=DF[,colnames(DF) %in% DF_Outcome], nfolds = 5, type.measure="auc",
                                                                             family = "binomial", alpha = 1)$lambda.1se)
  
  CoreFeatures <- sum(LASSO_Example$beta != 0)
  
  CoreFeatures
}

LASSO_Optimal_CoreFeatures_2 <- function(DF, DF_Outcome) {
  
  CoreFeatures <- NULL
  
  LASSO_Example <- glmnet(x=DF[,-which(colnames(DF) %in% DF_Outcome)], 
                          y=DF[,colnames(DF) %in% DF_Outcome], 
                          family = "binomial", alpha = 0.999, lambda = cv.glmnet(x=as.matrix(DF[,-which(colnames(DF) %in% DF_Outcome)]), 
                                                                             y=DF[,colnames(DF) %in% DF_Outcome], nfolds = 5, type.measure="auc",
                                                                             family = "binomial", alpha = 1)$lambda.1se)
  
  CoreFeatures <- sum(LASSO_Example$beta != 0)
  
  CoreFeatures
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

RF_Curve_Evaluate <- function(DF, DF_Outcome, n) {
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
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


#----------------------------------

RF_Curve_Random <- function(DF, DF_Outcome, size) {
  
  a <- 500
  b <- size
  
  stepsize <- round((b - a) / 10)
  
  S <- seq(a, b, by=stepsize)
  
  samps <- S
  
  if (length(samps) > 10) {samps <- samps[1:10]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
    
    ComboAUC <- NULL
    
    for (i in 1:10) {
      
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
    
    
  }
  Curve
}

LASSO_Optimal <- function(DF, DF_Outcome) {
  
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

LASSO_Optimal_2 <- function() {
  
  AUC <- NULL
  
  Folds_inv <- create_folds(as.factor(LOS[,colnames(LOS) %in% "LOS"]), k=5, type="stratified", invert = T)
  
  for (j in 1:length(Folds_inv)) {
    HoldOut <- Folds_inv[[j]]
    
    LASSO_Example <- glm(data=LOS[-HoldOut,], LOS ~ ., family=binomial(logit))
    
    Output <- predict(LASSO_Example, newdata = LOS[HoldOut, -which(colnames(LOS) %in% "LOS")],type="response")
    
    AUC <- append(AUC, as.numeric(roc(as.factor(LOS[HoldOut,colnames(LOS) %in% "LOS"]), Output, quiet=T)$auc))
  }
  c(mean(AUC), sd(AUC) / sqrt(5))

}

LASSO_Curve_Evaluate <- function(DF, DF_Outcome, n) {
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
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
      
      LASSO_Example <- glmnet(x=Reduced[-HoldOut,-which(colnames(Reduced) %in% DF_Outcome)], 
                              y=Reduced[-HoldOut,colnames(Reduced) %in% DF_Outcome], 
                              family = "binomial", alpha = 1, lambda = 0)
      
      Output <- unname(predict(LASSO_Example, newx = as.matrix(Reduced[HoldOut, -which(colnames(Reduced) %in% DF_Outcome)]),
                               type="response")[,1])
      
      AUC <- append(AUC, as.numeric(roc(as.factor(Reduced[HoldOut,colnames(Reduced) %in% DF_Outcome]), Output, quiet=T)$auc))
    }
    
    ComboAUC <- append(ComboAUC, mean(AUC))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}

LASSO_Curve_Random <- function(DF, DF_Outcome) {
  
  if(mean(DF[,colnames(DF) %in% DF_Outcome]) - 0.15 < 0) {
    a <- 500
    b <- 10000
    
    stepsize <- round((b - a) / 25)
    
    S <- seq(a, b, by=stepsize)
    
  }
  
  else {
    
    a <- 100
    b <- 5000
    
    stepsize <- round((b - a) / 25)
    
    S <- seq(a, b, by=stepsize)
    
  }
  
  samps <- S
  
  if (length(samps) > 25) {samps <- samps[1:25]}
  
  Curve <- NULL
  
  for (n in samps) {
    
    Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
    
    ComboAUC <- NULL
    
    for (i in 1:10) {
      
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
        
        LASSO_Example <- glmnet(x=Reduced[-HoldOut,-which(colnames(Reduced) %in% DF_Outcome)], 
                                y=Reduced[-HoldOut,colnames(Reduced) %in% DF_Outcome], 
                                family = "binomial", alpha = 1, lambda = 0)
        
        Output <- unname(predict(LASSO_Example, newx = as.matrix(Reduced[HoldOut, -which(colnames(Reduced) %in% DF_Outcome)]),
                                 type="response")[,1])
        
        AUC <- append(AUC, as.numeric(roc(as.factor(Reduced[HoldOut,colnames(Reduced) %in% DF_Outcome]), Output, quiet=T)$auc))
      }
      
      ComboAUC <- append(ComboAUC, mean(AUC))
      
    }
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
    )))
    
    
  }
  Curve
}

NN_Optimal <- function(DF, DF_Outcome) {
  
  DF <- data.frame(cbind(apply(DF[,-which(sapply(1:ncol(DF), function(i) {length(table(DF[,i]))==2}))], 
                               2, scale), DF[,which(sapply(1:ncol(DF), function(i) {length(table(DF[,i]))==2}))]))
  
  DF[,which(colnames(DF) %in% DF_Outcome)] <- as.factor(DF[,which(colnames(DF) %in% DF_Outcome)])
  
  DF <- as.h2o(DF)
  
  NN_Example <- h2o.deeplearning(y=DF_Outcome,
                                 training_frame = DF, nfolds=5, standardize = F, fold_assignment="Stratified",
                                 keep_cross_validation_predictions=T,  reproducible = TRUE, hidden = c(20), epochs = 10,
                                 seed = 1)
  
  h2o.auc(NN_Example, xval = T)
}

NN_Optimal_2 <- function(DF, DF_Outcome) {
  
  DF[,which(colnames(DF) %in% DF_Outcome)] <- as.factor(DF[,which(colnames(DF) %in% DF_Outcome)])
  
  DF <- as.h2o(DF)
  
  NN_Example <- h2o.deeplearning(y=DF_Outcome,
                                 training_frame = DF, nfolds=5, standardize = F, fold_assignment="Stratified",
                                 keep_cross_validation_predictions=T,  reproducible = TRUE, hidden = c(20), epochs = 10,
                                 seed = 1)
  
  h2o.auc(NN_Example, xval = T)
}

NN_Curve_Evaluate <- function(DF, DF_Outcome, n) {
  
  DF <- data.frame(cbind(apply(DF[,-which(sapply(1:ncol(DF), function(i) {length(table(DF[,i]))==2}))], 
                               2, scale), DF[,which(sapply(1:ncol(DF), function(i) {length(table(DF[,i]))==2}))]))
  
  DF[,which(colnames(DF) %in% DF_Outcome)] <- as.factor(DF[,which(colnames(DF) %in% DF_Outcome)])
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
    Reduced <- DF[Reduced_Rows[[i]],]
    
    #if (n < 100) {
    
    #Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=2, type="stratified", invert = T)
    
    #}
    
    #else {
    
    Reduced <- as.h2o(Reduced)
    
    NN_Example <- h2o.deeplearning(y=DF_Outcome,
                                   training_frame = Reduced, nfolds=5, standardize = F, fold_assignment="Stratified",
                                   keep_cross_validation_predictions=T,  reproducible = FALSE, hidden = c(20), epochs = 10)
    
    ComboAUC <- append(ComboAUC, h2o.auc(NN_Example, xval = T))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}


NN_Curve_Evaluate_2 <- function(DF, DF_Outcome, n) {
  
  DF[,which(colnames(DF) %in% DF_Outcome)] <- as.factor(DF[,which(colnames(DF) %in% DF_Outcome)])
  
  Curve <- NULL
  
  Reduced_Rows <- Subsample_DF_random(DF=DF, DF_Response=DF_Outcome, n=n)
  
  ComboAUC <- NULL
  
  for (i in 1:10) {
    
    Reduced <- DF[Reduced_Rows[[i]],]
    
    #if (n < 100) {
    
    #Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=2, type="stratified", invert = T)
    
    #}
    
    #else {
    
    Reduced <- as.h2o(Reduced)
    
    NN_Example <- h2o.deeplearning(y=DF_Outcome,
                                   training_frame = Reduced, nfolds=5, standardize = F, fold_assignment="Stratified",
                                   keep_cross_validation_predictions=T,  reproducible = FALSE, hidden = c(20), epochs = 10)
    
    ComboAUC <- append(ComboAUC, h2o.auc(NN_Example, xval = T))
    
  }
  
  Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
  )))
  
  Curve
}


LASSO_Curve_Random_2 <- function() {
  
  S <- seq(500, 10000, by=round(10000-500)/25)
  
  samps <- S[1:25]
  
  Curve <- NULL
  
  for (n in samps) {
    
    Reduced_Rows <- Subsample_DF_random(DF=LOS, DF_Response="LOS", n=n)
    
    ComboAUC <- NULL
    
    for (i in 1:10) {
      
      Reduced <- LOS[Reduced_Rows[[i]],]
      
      AUC <- NULL
      
      #if (n < 100) {
      
      #Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% DF_Outcome]), k=2, type="stratified", invert = T)
      
      #}
      
      #else {
      
      Folds_inv <- create_folds(as.factor(Reduced[,colnames(Reduced) %in% "LOS"]), k=5, type="stratified", invert = T)
      
      
      #}
      
      for (j in 1:length(Folds_inv)) {
        HoldOut <- Folds_inv[[j]]
        
        LASSO_Example <- glm(data=Reduced[-HoldOut,], LOS ~ ., family=binomial(logit))
        
        Output <- predict(LASSO_Example, newdata = Reduced[HoldOut, -which(colnames(Reduced) %in% "LOS")],
                                 type="response")
        
        AUC <- append(AUC, as.numeric(roc(as.factor(Reduced[HoldOut,colnames(Reduced) %in% "LOS"]), Output, quiet=T)$auc))
      }
      
      ComboAUC <- append(ComboAUC, mean(AUC))
      
    }
    
    Curve <- data.frame(rbind(Curve, c("n" = n, "Curve" = mean(ComboAUC)
    )))
    
    
  }
  Curve
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
  
  Plot <- ggplot(data=Data_forPlot_Raw, aes(x=N, y=Raw)) + geom_point(color="blue", size=2.5) + ylim(c(0.75, 0.829)) + xlim(c(0,50000)) +
  geom_line(data=Values, aes(y=Fitted, x=Size), col="Blue", size=1) + 
  labs(x="Sample Size", y="AUC") + theme_bw() + theme(axis.title=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #+ annotate("Text", x=9306, y=0.804, label="*", size=12, col="red") 
  #geom_segment(x=5500, xend=50000, y=0.8, yend=0.8, linetype="dotted", col="red", size=1.2) + 
    theme(
                                                                         axis.text.x=element_blank(),
                                                                         axis.ticks.x=element_blank(),
    
                                                                        axis.text.y=element_blank(),
                                                                        axis.ticks.y=element_blank())
  
  Plot
  
  #Plot <- ggplot(data=Data_forPlot_Raw, aes(x=N, y=Raw)) + geom_point(color="grey") +
  #  labs(x="Sample Size", y="AUC") + ylim(c(min(Data_forPlot_Raw$Raw)-0.01, FullAUC+0.01))
  
}

ggarrange(plot_curve_full(Diabetes_XGB, Type="PowerLaw", FullAUC=0.829, Size=50000),
          plot_curve_full(Diabetes_XGB, Type="PowerLaw_Fixed", FullAUC=0.829, Size=50000),
          plot_curve_full(Diabetes_XGB, Type="LOG", FullAUC=0.829, Size=50000))

#Actually use this one below
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
  
  #Plot <- ggplot(data=Data_forPlot_Raw, aes(x=N, y=Raw)) + geom_point(color="grey") +
  #  labs(x="Sample Size", y="AUC") + ylim(c(min(Data_forPlot_Raw$Raw)-0.01, FullAUC+0.01))
  
  
  #Plot + 
  #ggplot(data=Data_forPlot_Raw, aes(x=N, y=Raw)) + theme_bw() +  theme(
  #                                                                     axis.text.x=element_blank(),
  #                                                                     axis.ticks.x=element_blank(),
  
  #                                                                     axis.text.y=element_blank(),
  #                                                                     axis.ticks.y=element_blank()) +
  #labs(x="Sample Size", y="AUC", title="The Limiting Performance of a Classifier") + theme(axis.title=element_text(size=14)) + ylim(c(min(Data_forPlot_Raw$Raw)-0.01, FullAUC+0.01)) + geom_line(data=Values, aes(y=Fitted, x=Size), col="Blue", size=1.5) 
  
  
}

plot_curve_full(Diabetes_XGB, Type="PowerLaw", 0.829, 50000)





