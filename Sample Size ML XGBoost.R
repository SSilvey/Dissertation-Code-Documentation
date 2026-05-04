
#XGBoost
set.seed(2023)
XGBoost_Optimal(Cardio, "CVD")
XGBoost_Optimal(Cardio, "CVD")[2]*2.25

set.seed(2023)
XGBoost_Optimal(BreastTumor, "recurrence")
XGBoost_Optimal(BreastTumor, "recurrence")[2]*2.25

set.seed(2023)
XGBoost_Optimal(NoShow, "NoShow")
XGBoost_Optimal(NoShow, "NoShow")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Diabetes130, "readmitted")
XGBoost_Optimal(Diabetes130, "readmitted")[2]*2.25

set.seed(2023)
XGBoost_Optimal(COVID, "RESULTADO")
XGBoost_Optimal(COVID, "RESULTADO")[2]*2.25

set.seed(2023)
XGBoost_Optimal(LOS, "LOS")
XGBoost_Optimal(LOS, "LOS")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Diabetes, "Diabetes_binary")
XGBoost_Optimal(Diabetes, "Diabetes_binary")[2]*2.25

set.seed(2023)
XGBoost_Optimal(CDC_2020, "HeartDisease")
XGBoost_Optimal(CDC_2020, "HeartDisease")[2]*2.25

set.seed(2023)
XGBoost_Optimal(CDC_2022, "HighRiskLastYear")
XGBoost_Optimal(CDC_2022, "HighRiskLastYear")[2]*2.25

set.seed(2023)
XGBoost_Optimal(HeartSynthetic, "class")
XGBoost_Optimal(HeartSynthetic, "class")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Hepatitis, "Class")
XGBoost_Optimal(Hepatitis, "Class")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Lymph, "class")
XGBoost_Optimal(Lymph, "class")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Pharynx, "class")
XGBoost_Optimal(Pharynx, "class")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Chol, "chol")
XGBoost_Optimal(Chol, "chol")[2]*2.25

set.seed(2023)
XGBoost_Optimal(Dermatology, "family_history")
XGBoost_Optimal(Dermatology, "family_history")[2]*2.25

set.seed(2023)
XGBoost_Optimal(PBC, "class")
XGBoost_Optimal(PBC, "class")[2]*2.25

#Actual Curve Fit ---------------------------------------------------------

plot_curve(Curve=Hepatitis_XGB, T, Title="Diabetes130 Dataset, All Features", Size=50000)


#-------------------------------------------------------------------------

#------------------------------------------------------------------------
#More Refined Curves 
  
set.seed(2023)
Cardio_XGB <- XGBoost_Curve_Random(Cardio, "CVD")

set.seed(2023)
BreastTumor_XGB <- XGBoost_Curve_Random(BreastTumor, "recurrence")

set.seed(2023)
NoShow_XGB <- XGBoost_Curve_Random(NoShow, "NoShow")

set.seed(2023)
Diabetes130_XGB <- XGBoost_Curve_Random(Diabetes130, "readmitted")

set.seed(2023)
Diabetes_XGB <- XGBoost_Curve_Random(Diabetes, "Diabetes_binary")

set.seed(2023)
COVID_XGB <- XGBoost_Curve_Random(COVID, "RESULTADO")

set.seed(2024)
LOS_XGB <- XGBoost_Curve_Random(LOS, "LOS")

set.seed(2023)
CDC_2020_XGB <- XGBoost_Curve_Random(CDC_2020, "HeartDisease")

set.seed(2023)
CDC_2022_XGB <- XGBoost_Curve_Random(CDC_2022, "HighRiskLastYear")

set.seed(2023)
Heart_XGB <- XGBoost_Curve_Random(HeartSynthetic, "class")

set.seed(2023)
Hepatitis_XGB <- XGBoost_Curve_Random(Hepatitis, "Class")

set.seed(2023)
Lymph_XGB <- XGBoost_Curve_Random(Lymph, "class")

set.seed(2023)
Pharynx_XGB <- XGBoost_Curve_Random(Pharynx, "class")

set.seed(2023)
Chol_XGB <- XGBoost_Curve_Random_2(Chol, "chol", 100000)

set.seed(2023)
Dermatology_XGB <- XGBoost_Curve_Random(Dermatology, "family_history")

set.seed(2023)
PBC_XGB <- XGBoost_Curve_Random(PBC, "class")

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Learning Curve Backup Data")
#write.csv(Diabetes_XGB, "Diabetes_XGB.csv")
#write.csv(Cardio_XGB, "Cardio_XGB.csv")
#write.csv(Diabetes130_XGB, "Diabetes130_XGB.csv")
#write.csv(NoShow_XGB, "NoShow_XGB.csv")
#write.csv(BreastTumor_XGB, "BreastTumor_XGB.csv")
#write.csv(Heart_XGB, "Heart_XGB.csv")
#write.csv(Hepatitis_XGB, "Hepatitis_XGB.csv")
#write.csv(Lymph_XGB, "Lymph_XGB.csv")
#write.csv(Pharynx_XGB, "Pharynx_XGB.csv")
#write.csv(Chol_XGB, "Chol_XGB.csv")
#write.csv(Dermatology_XGB, "Dermatology_XGB.csv")
#write.csv(PBC_XGB, "PBC_XGB.csv")
#write.csv(CDC_2020_XGB, "CDC_2020_XGB.csv")
#write.csv(CDC_2022_XGB, "CDC_2022_XGB.csv")
#write.csv(COVID_XGB, "COVID_XGB.csv")
#write.csv(LOS_XGB, "LOS_XGB.csv")

Diabetes_XGB <- read.csv("Diabetes_XGB.csv")
Cardio_XGB <- read.csv("Cardio_XGB.csv")
Diabetes130_XGB <- read.csv("Diabetes130_XGB.csv")
NoShow_XGB <- read.csv("NoShow_XGB.csv")
BreastTumor_XGB <- read.csv("BreastTumor_XGB.csv")
Heart_XGB <- read.csv("Heart_XGB.csv")
Hepatitis_XGB <- read.csv("Hepatitis_XGB.csv")
Lymph_XGB <- read.csv("Lymph_XGB.csv")
Pharynx_XGB <- read.csv("Pharynx_XGB.csv")
Chol_XGB <- read.csv("Chol_XGB.csv")
PBC_XGB <- read.csv("PBC_XGB.csv")
CDC_2020_XGB <- read.csv("CDC_2020_XGB.csv")
CDC_2022_XGB <- read.csv("CDC_2022_XGB.csv")
Dermatology_XGB <- read.csv("Dermatology_XGB.csv")
LOS_XGB <- read.csv("LOS_XGB.csv")
COVID_XGB <- read.csv("COVID_XGB.csv")


#Plots and Results ----------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------
#PART 2 - INFERENCE


#XGB
N <- c(
Get_N(Cardio_XGB, 0.802, 50000),
Get_N(Diabetes130_XGB, 0.662, 50000),
Get_N_LOG(NoShow_XGB, 0.608, 50000),
Get_N_LOG(BreastTumor_XGB, 0.777, 50000),
Get_N(Diabetes_XGB, 0.829, 50000),
Get_N(COVID_XGB, 0.664, 50000),
Get_N_LOG(LOS_XGB, 0.917, 50000),
Get_N(CDC_2020_XGB, 0.815, 50000),
Get_N_LOG(CDC_2022_XGB, 0.815, 50000),
Get_N(Heart_XGB, 0.965, 50000),
Get_N(Hepatitis_XGB, 0.979, 50000),
Get_N(Lymph_XGB, 0.957, 50000),
Get_N_LOG(Pharynx_XGB, 0.858, 50000),
Get_N_LOG(Chol_XGB, 0.736, 100000),
Get_N(Dermatology_XGB, 0.859, 50000),
Get_N_LOG(PBC_XGB, 0.850, 50000))

#RF
N_RF <- c(
  Get_N(Cardio_RF, 0.796, 50000),
  Get_N(Diabetes130_RF, 0.661, 50000),
  Get_N(NoShow_RF, 0.609, 50000),
  Get_N_LOG(BreastTumor_RF, 0.780, 50000),
  Get_N(Diabetes_RF, 0.822, 50000),
  Get_N(COVID_RF, 0.661, 50000),
  Get_N_LOG(LOS_RF, 0.915, 50000),
  Get_N_LOG(CDC_2020_RF, 0.810, 50000),
  Get_N_LOG(CDC_2022_RF, 0.801, 50000),
  max(Get_N_LOG(Heart_RF, 0.963, 50000), 250),
  Get_N(Hepatitis_RF, 0.976, 50000),
  Get_N_LOG(Lymph_RF, 0.957, 50000),
  Get_N(Pharynx_RF, 0.858, 50000),
  Get_N(Chol_RF, 0.728, 200000),
  Get_N(Dermatology_RF, 0.857, 50000),
  Get_N_LOG(PBC_RF, 0.850, 100000))

#LR
N_LR <- c(
  Get_N_Fixed(Cardio_LR, 0.784, 5000),
  Get_N_LOG(Diabetes130_LR, 0.654, 5000),
  Get_N_LOG(NoShow_LR, 0.596, 5000),
  Get_N_Fixed(BreastTumor_LR, 0.682, 5000),
  Get_N_Fixed(Diabetes_LR, 0.822, 10000),
  Get_N_Fixed(COVID_LR, 0.642, 5000),
  Get_N_Fixed(LOS_LR, 0.898, 5000),
  Get_N_Fixed(CDC_2020_LR, 0.810, 10000),
  Get_N_Fixed(CDC_2022_LR, 0.809, 10000),
  Get_N_Fixed(Heart_LR, 0.949, 5000),
  Get_N_Fixed(Hepatitis_LR, 0.939, 5000),
  Get_N_Fixed(Lymph_LR, 0.934, 5000),
  Get_N_Fixed(Pharynx_LR, 0.843, 5000),
  Get_N_LOG(Chol_LR, 0.669, 5000),
  Get_N_Fixed(Dermatology_LR, 0.823, 10000),
  Get_N_Fixed(PBC_LR, 0.788, 5000))

#NN
N_NN <- c(
  Get_N_Fixed(Cardio_NN, 0.795, 50000),
  Get_N_LOG(Diabetes130_NN, 0.661, 50000),
  Get_N_LOG(NoShow_NN, 0.603, 50000),
  Get_N_Fixed(BreastTumor_NN, 0.730, 50000),
  Get_N_Fixed(Diabetes_NN, 0.826, 50000),
  Get_N(COVID_NN, 0.661, 50000),
  Get_N_Fixed(LOS_NN, 0.901, 50000),
  Get_N_Fixed(CDC_2020_NN, 0.809, 50000),
  Get_N_Fixed(CDC_2022_NN, 0.801, 50000),
  Get_N(Heart_NN, 0.962, 50000),
  Get_N_Fixed(Hepatitis_NN, 0.974, 50000),
  Get_N(Lymph_NN, 0.956, 50000),
  Get_N(Pharynx_NN, 0.856, 50000),
  Get_N_LOG(Chol_NN, 0.714, 250000),
  Get_N_LOG(Dermatology_NN, 0.859, 50000),
  Get_N_LOG(PBC_NN, 0.823, 100000))


#XGB
AUC_Full <- c(0.802, 0.662, 0.608, 0.777, 0.829, 0.664, 0.917, 0.815, 0.815, 0.965, 0.979, 0.957, 0.858, 0.736, 0.859, 0.850)*100

#RF
AUC_Full <- c(0.796, 0.661, 0.609, 0.780, 0.822, 0.661, 0.915, 0.810, 0.801, 0.963, 0.976, 0.957, 0.858, 0.728, 0.857, 0.850)*100

#LR
AUC_Full <- c(0.784, 0.654, 0.596, 0.682, 0.822, 0.642, 0.898, 0.810, 0.809, 0.949, 0.939, 0.934, 0.843, 0.669, 0.823, 0.788)*100

#NN
AUC_Full <- c(0.795, 0.661, 0.603, 0.730, 0.826, 0.661, 0.901, 0.809, 0.801, 0.962, 0.974, 0.956, 0.856, 0.714, 0.852, 0.823)*100

Imbalance <- c(50.0, 46.1, 20.2, 34.6, 13.9, 39.0, 2.1, 8.6, 4.4, 44.4, 20.8, 45.7, 25.6, 16.5, 17.8, 13.2)

N_Features <- c(11,35,8,9,21,16,11,17,39,13,19,18,11,13,33,18)

Complexity <- c(1.8, 0.8, 1.2, 9.5, 0.7, 2.2, 1.9, 0.5, 0.6, 1.6, 4.0, 2.3, 1.5, 6.7, 3.6, 6.2)

#Justification for Nonlinearity Threshold--------------------------------------------------------------------------------------------------------------------

AUC_Full <- c(0.802, 0.662, 0.608, 0.777, 0.829, 0.664, 0.917, 0.815, 0.815, 0.965, 0.979, 0.957, 0.858, 0.736, 0.859, 0.850)*100
TEST <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "Complexity"=Complexity)

TEST_XGB <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "AUC_Full"=AUC_Full)

AUC_Full_RF <- c(0.796, 0.661, 0.609, 0.780, 0.822, 0.661, 0.915, 0.810, 0.801, 0.963, 0.976, 0.957, 0.858, 0.728, 0.857, 0.850)*100
TEST_RF <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "AUC_Full"=AUC_Full_RF)

AUC_Full_NN <- c(0.795, 0.661, 0.603, 0.730, 0.826, 0.661, 0.901, 0.809, 0.801, 0.962, 0.974, 0.956, 0.856, 0.714, 0.852, 0.823)*100
TEST_NN <- data.frame("N"=N, "N_RF"=N_RF, "N_NN"=N_NN, "AUC_Full"=AUC_Full_NN)

N_Plot <- ggplot(data=TEST, aes(x=Complexity, y=N)) + geom_point(color="blue", size=2) + geom_vline(xintercept=5 ,linetype = 'dotted', color="red", size=1.5) + theme_bw() + labs(title="XGB", y="Sample Size", x="Nonlinearity") 
N_Plot_RF <- ggplot(data=TEST, aes(x=Complexity, y=N_RF)) + geom_point(color="blue", size=2) + geom_vline(xintercept=5 ,linetype = 'dotted', color="red", size=1.5) + theme_bw() + labs(title="RF", y="Sample Size", x="Nonlinearity") 
N_Plot_NN <- ggplot(data=TEST, aes(x=Complexity, y=N_NN)) + geom_point(color="blue", size=2) + geom_vline(xintercept=5 ,linetype = 'dotted', color="red", size=1.5) + theme_bw() + labs(title="NN", y="Sample Size", x="Nonlinearity") 

library(ggpubr)
ggarrange(N_Plot, N_Plot_RF, N_Plot_NN, ncol=3, nrow=1)




library(boot)
#TO derive 5 cutoff for nonlinearity
plot(seq(2,8, by=0.5), sapply(seq(2,8, by=0.5), function(i) {
  sqrt(cv.glm(glm.nb(N ~ I(Complexity>i)), data=TEST)$delta[1])}),xlab="Cutoff", ylab="LOO-CV MSE", main="XGB")

plot(seq(2,8, by=0.5), sapply(seq(2,8, by=0.5), function(i) {
  sqrt(cv.glm(glm.nb(N_RF ~ I(Complexity>i)), data=TEST)$delta[1]) }), xlab="Cutoff", ylab="LOO-CV MSE", main="RF")

plot(seq(2,8, by=0.5), sapply(seq(2,8, by=0.5), function(i) {
  sqrt(cv.glm(glm.nb(N_NN ~ I(Complexity>i)), data=TEST)$delta[1]) }),xlab="Cutoff", ylab="LOO-CV MSE", main="NN")






#End--------------------------------------------------------------------------------------------------------------------


Complexity2 <- ifelse(Complexity >= 5, 1, 0)

PercentageBinary <- c(5, 7, 1, 3, 4, 1, 3, 4, 6, 6, 6, 1, 2, 4, 1, 10) / N_Features

CoreFeatures <- c(11/11, 16/35, 2/8, 4/9, 18/21, 6/16, 4/11, 13/17, 23/39, 12/13, 18/19, 15/18, 8/11, 8/13, 28/33, 15/18)*100
CoreFeatures <- c(95.4, 41.4, 26.3, 44.4, 67.1, 37.5, 46.4, 65.9, 55.9, 92.3, 93.7, 83.3, 68.2, 53.1, 85.2, 83.3)

Event_Per_Var <- N / N_Features

plot(Imbalance, log(N))
plot(AUC_Full, log(N))
plot(N_Features, log(N))
plot(Complexity, log(N))
plot(CoreFeatures, log(N))
plot(PercentageBinary, log(N))

summary(glm.nb(N ~ Imbalance))
exp(coef(glm.nb(N ~ Imbalance)))
exp(confint(glm.nb(N ~ Imbalance)))

summary(glm.nb(N ~ I(AUC_Full-50)))
exp(coef(glm.nb(N ~ I(AUC_Full-50))))
exp(confint(glm.nb(N ~ I(AUC_Full-50))))

summary(glm.nb(N ~ N_Features))
exp(coef(glm.nb(N ~ N_Features)))
exp(confint(glm.nb(N ~ N_Features)))

summary(glm.nb(N ~ CoreFeatures))
exp(coef(glm.nb(N ~ CoreFeatures)))
exp(confint(glm.nb(N ~ CoreFeatures)))

summary(glm.nb(N ~ Complexity2))
exp(coef(glm.nb(N ~ Complexity2)))
exp(confint(glm.nb(N ~ Complexity2)))

summary(glm.nb(N ~ I(PercentageBinary*100)))
exp(coef(glm.nb(N ~ I(PercentageBinary*100))))
exp(confint(glm.nb(N ~ I(PercentageBinary*100))))


Dataset <- data.frame("N"= N, "Imbalance"=Imbalance, "Complexity"=Complexity, "Complexity2"=Complexity2, "AUCFull"=c(0.802, 0.662, 0.608, 0.777, 0.829, 0.664, 0.917, 0.815, 0.815, 0.965, 0.979, 0.957, 0.858, 0.736, 0.859, 0.850)*100 - 50, "PercentageBinary"=PercentageBinary*100, "N_Features"=N_Features, "CoreFeatures"=CoreFeatures)

Dataset_2 <- data.frame("N"= N_RF, "Imbalance"=Imbalance, "Complexity"=Complexity, "Complexity2"=Complexity2, "AUC_Full"=c(0.796, 0.661, 0.609, 0.780, 0.822, 0.661, 0.915, 0.810, 0.801, 0.963, 0.976, 0.957, 0.858, 0.728, 0.857, 0.850)*100 - 50,
                        "PercentageBinary"=PercentageBinary, "N_Features"=N_Features, "CoreFeatures"=CoreFeatures)

Dataset_4 <- data.frame("N"= N_LR, "Imbalance"=Imbalance, "N_Features"=N_Features, "CoreFeatures"=CoreFeatures, "PercentageBinary"=PercentageBinary, "N_Features"=N_Features, "CoreFeatures"=CoreFeatures,"Complexity2"=Complexity2,
                        "AUC_Full"=c(0.784, 0.654, 0.596, 0.682, 0.822, 0.642, 0.898, 0.810, 0.809, 0.949, 0.939, 0.934, 0.843, 0.669, 0.823, 0.788)*100)

Dataset_3 <- data.frame("N"= N_NN, "Imbalance"=Imbalance, "N_Features"=N_Features, "Complexity"=Complexity, "Complexity2"=Complexity2, "PercentageBinary"=PercentageBinary*100, "N_Features"=N_Features, "CoreFeatures"=CoreFeatures,
                        "AUC_Full"=c(0.795, 0.661, 0.603, 0.730, 0.826, 0.661, 0.901, 0.809, 0.801, 0.962, 0.974, 0.956, 0.856, 0.714, 0.852, 0.823)*100)

library(qpcR)
library(MASS)
library(modEvA)
library(glmtoolbox)
library(boot)

adjR2(Mod_XGB)
adjR2(Mod_RF)
adjR2(Mod_LR)
adjR2(Mod_NN) #Need to discuss this in manuscript
                                                                                                                                             
Mod_XGB <- (glm.nb(data=Dataset, N ~ Imbalance + AUCFull + Complexity2))
summary(Mod_XGB)
sqrt(cv.glm(Mod_XGB, data=Dataset)$delta[1])
exp(coef(Mod_XGB))
exp(confint(Mod_XGB))

Mod_RF <- (glm.nb(data=Dataset_2, N ~ Imbalance + Complexity2 + AUC_Full))
summary(Mod_RF)
sqrt(cv.glm(Mod_RF, data=Dataset_2)$delta[1])
exp(coef(Mod_RF))
exp(confint(Mod_RF))

Mod_LR <- (glm.nb(data=Dataset_4, N ~ Imbalance + N_Features + CoreFeatures))
summary(Mod_LR)
sqrt(cv.glm(Mod_LR, data=Dataset_4)$delta[1])
exp(coef(Mod_LR))
exp(confint(Mod_LR))

Mod_NN <- (glm.nb(data=Dataset_3, N ~ Imbalance + PercentageBinary + Complexity2))
summary(Mod_NN)
sqrt(cv.glm(Mod_NN, data=Dataset_3)$delta[1])
exp(coef(Mod_NN))
exp(confint(Mod_NN))

exp(predict(Mod_XGB, newdata=data.frame("Imbalance"=15, "Complexity2"=0, "AUCFull"=25), se.fit=T)$fit)

library(ciTools)

DF_new <- Dataset[,c("Imbalance", "Complexity2", "AUCFull", "N")]

DF_new <- rbind(DF_new, data.frame("Imbalance"=15, "Complexity2"=0, "AUCFull"=25, "N"=exp(predict(Mod_XGB, newdata=data.frame("Imbalance"=15, "Complexity2"=0, "AUCFull"=25)))))

PI <- add_pi(DF_new, Mod_XGB)
PI[nrow(PI),c(5,6,7)] 

#Model Estimates Plots -----------------------------------
A <- ggpredict(Mod_XGB, terms=c("AUCFull [all]", "Complexity2 [all]", "Imbalance [5,10,20,30,40,50]")) %>% plot() + scale_color_manual(labels=c("Low", "High"), values=c("red", "blue")) + scale_fill_manual(labels=c("Low", "High"), values=c("red", "blue")) + labs(x="Separability", colour = "Nonlinearity", y="Predicted sample size", title="XGB") +
  scale_x_continuous(breaks = c(10,20,30,40),
                     labels = paste0(c("60", "70", "80", "90"))) + theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))

B <- ggpredict(Mod_RF, terms=c("AUC_Full [all]", "Complexity2 [all]", "Imbalance [5,10,20,30,40,50]")) %>% plot() + scale_color_manual(labels=c("Low", "High"), values=c("red", "blue")) + scale_fill_manual(labels=c("Low", "High"), values=c("red", "blue")) + labs(x="Separability", colour = "Nonlinearity", y="Predicted sample size", title="RF") +
  scale_x_continuous(breaks = c(10,20,30,40),
                     labels = paste0(c("60", "70", "80", "90")))+ theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))

C <- ggpredict(Mod_NN, terms=c("PercentageBinary [0,20,40,60,80,100]", "Complexity2 [all]","Imbalance [5,10,20,30,40,50]")) %>% plot() + scale_color_manual(labels=c("Low", "High"), values=c("red", "blue")) + scale_fill_manual(labels=c("Low", "High"), values=c("red", "blue")) + labs(x="Percentage of continuous numeric features", colour = "Nonlinearity", y="Predicted sample size", title="NN") + theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))

D <- ggpredict(Mod_LR, terms=c("CoreFeatures [all]", "N_Features [10,20,30]", "Imbalance [5,10,20,30,40,50]")) %>% plot() + scale_color_manual(labels=c("10", "20", "30"), values=c("green", "red", "blue")) + scale_fill_manual(labels=c("10", "20", "30"), values=c("green", "red", "blue")) + labs(x="Percentage of core linear features", colour = "Number of features", y="Predicted sample size", title="LR")+ theme(plot.background = element_rect(color = "black")) + theme(plot.title = element_text(hjust = 0.5))


#------------------------------------------------------------

exp(coef(Mod_XGB))
exp(confint(Mod_XGB))

exp(coef(Mod_RF))

exp(coef(Mod_LR))

exp(coef(Mod_NN))

summary(Mod_NN)
predict(Mod_XGB, newdata = data.frame("Imbalance"=50, "Complexity"=1, "AUC_Full"=75))
#8.148

exp(11.901644)*exp(-0.04827*50)*exp(0.128074*1)*exp(-0.058694*25)

exp(11.90)*(0.95^50)

exp(2*3)

exp(2)^3

library(boot)





#--------------Inference

#The difference between the fullDSAUC using XGBoost and using LASSO.
#Higher difference indicates a higher degree of data complexity.



#Strength of the top linear 3 predictors (average)


#Binary features line
sum(sapply(1:ncol(LOS), function(i) {length(table(LOS[,i]))==2}))-1

sum(sapply(1:ncol(PBC), function(i) {length(table(PBC[,i]))>10}))



#Plots and Predictions


set.seed(2023)
LASSO_Optimal_CoreFeatures(Diabetes, "Diabetes_binary") #18/21
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Diabetes, "Diabetes_binary")))/21 #67.1

set.seed(2023)
LASSO_Optimal_CoreFeatures(CDC_2020, "HeartDisease") #13/17
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(CDC_2020, "HeartDisease")))/17 #65.9

set.seed(2023)
LASSO_Optimal_CoreFeatures(CDC_2022, "HighRiskLastYear") #23/39
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(CDC_2022, "HighRiskLastYear")))/39 #55.9

set.seed(2023)
LASSO_Optimal_CoreFeatures(HeartSynthetic, "class") #12/13
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(HeartSynthetic, "class")))/13 #92.3

set.seed(2023)
LASSO_Optimal_CoreFeatures(Hepatitis, "Class") #18/19
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Hepatitis, "Class")))/19 #93.7

set.seed(2023)
LASSO_Optimal_CoreFeatures(Lymph, "class") #15/18
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Lymph, "class")))/18 #83.3

set.seed(2023)
LASSO_Optimal_CoreFeatures(Pharynx, "class") #8/11
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Pharynx, "class")))/11 #68.2

set.seed(2023)
LASSO_Optimal_CoreFeatures(Chol, "chol") #8/13
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Chol, "chol")))/13 #53.1%

set.seed(2023)
LASSO_Optimal_CoreFeatures(Cardio, "CVD") #11/11
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Cardio, "CVD")))/11 #95.4%

set.seed(2023)
LASSO_Optimal_CoreFeatures(BreastTumor, "recurrence") #4/9
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(BreastTumor, "recurrence")))/9 #44.4%

set.seed(2023)
LASSO_Optimal_CoreFeatures(NoShow, "NoShow") #2/8
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(NoShow, "NoShow")))/8 #26.3%

set.seed(2023)
LASSO_Optimal_CoreFeatures(Diabetes130, "readmitted") #16/35
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Diabetes130, "readmitted")))/35 #41.4%

set.seed(2023)
LASSO_Optimal_CoreFeatures(COVID, "RESULTADO") #6/16
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(COVID, "RESULTADO")))/16 #37.5%

set.seed(2023)
LASSO_Optimal_CoreFeatures_2(LOS, "LOS") #4/11
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures_2(LOS, "LOS")))/11 #46.4

set.seed(2023)
LASSO_Optimal_CoreFeatures(Dermatology, "family_history") #28/33
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(Dermatology, "family_history")))/33 #85.2%

set.seed(2023)
LASSO_Optimal_CoreFeatures(PBC, "class") #15/18
set.seed(2023)
mean(replicate(10, LASSO_Optimal_CoreFeatures(PBC, "class")))/18 #83.3%











  


                  