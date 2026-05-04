
#Logistic Regression

set.seed(2023)
LASSO_Optimal(Cardio, "CVD")
LASSO_Optimal(Cardio, "CVD")[2]*2.25

set.seed(2023)
LASSO_Optimal(BreastTumor, "recurrence")
LASSO_Optimal(BreastTumor, "recurrence")[2]*2.25

set.seed(2023)
LASSO_Optimal(NoShow, "NoShow")
LASSO_Optimal(NoShow, "NoShow")[2]*2.25

set.seed(2023)
LASSO_Optimal(Diabetes130, "readmitted")
LASSO_Optimal(Diabetes130, "readmitted")[2]*2.25

set.seed(2023)
LASSO_Optimal(Diabetes, "Diabetes_binary")
LASSO_Optimal(Diabetes, "Diabetes_binary")[2]*2.25

set.seed(2023)
LASSO_Optimal(COVID, "RESULTADO")
LASSO_Optimal(COVID, "RESULTADO")[2]*2.25

set.seed(2023)
LASSO_Optimal_2() 
LASSO_Optimal_2()[2]*2.25

set.seed(2023)
LASSO_Optimal(CDC_2020, "HeartDisease")
LASSO_Optimal(CDC_2020, "HeartDisease")[2]*2.25

set.seed(2023)
LASSO_Optimal(CDC_2022, "HighRiskLastYear")
LASSO_Optimal(CDC_2022, "HighRiskLastYear")[2]*2.25

set.seed(2023)
LASSO_Optimal(HeartSynthetic, "class")
LASSO_Optimal(HeartSynthetic, "class")[2]*2.25

set.seed(2023)
LASSO_Optimal(Hepatitis, "Class")
LASSO_Optimal(Hepatitis, "Class")[2]*2.25

set.seed(2023)
LASSO_Optimal(Lymph, "class")
LASSO_Optimal(Lymph, "class")[2]*2.25

set.seed(2023)
LASSO_Optimal(Pharynx, "class")
LASSO_Optimal(Pharynx, "class")[2]*2.25

set.seed(2023)
LASSO_Optimal(Chol, "chol")
LASSO_Optimal(Chol, "chol")[2]*2.25

set.seed(2023)
LASSO_Optimal(Dermatology, "family_history")
LASSO_Optimal(Dermatology, "family_history")[2]*2.25

set.seed(2023)
LASSO_Optimal(PBC, "class")
LASSO_Optimal(PBC, "class")[2]*2.25
#Functions


set.seed(2023)
Cardio_LR <- LASSO_Curve_Random(Cardio, "CVD")

set.seed(2023)
BreastTumor_LR <- LASSO_Curve_Random(BreastTumor, "recurrence")

set.seed(2023)
NoShow_LR <- LASSO_Curve_Random(NoShow, "NoShow")
NoShow_LR <- NoShow_LR[-1,]

set.seed(2023)
Diabetes130_LR <- LASSO_Curve_Random(Diabetes130, "readmitted")

set.seed(2023)
Diabetes_LR <- LASSO_Curve_Random(Diabetes, "Diabetes_binary")

set.seed(2023)
COVID_LR <- LASSO_Curve_Random(COVID, "RESULTADO")

set.seed(2024)
LOS_LR <- LASSO_Curve_Random_2()

set.seed(2023)
CDC_2020_LR <- LASSO_Curve_Random(CDC_2020, "HeartDisease")

set.seed(2023)
CDC_2022_LR <- LASSO_Curve_Random(CDC_2022, "HighRiskLastYear")

set.seed(2023)
Heart_LR <- LASSO_Curve_Random(HeartSynthetic, "class")

set.seed(2023)
Hepatitis_LR <- LASSO_Curve_Random(Hepatitis, "Class")

set.seed(2023)
Lymph_LR <- LASSO_Curve_Random(Lymph, "class")

set.seed(2023)
Pharynx_LR <- LASSO_Curve_Random(Pharynx, "class")

set.seed(2023)
Chol_LR <- LASSO_Curve_Random(Chol, "chol")

set.seed(2023)
Dermatology_LR <- LASSO_Curve_Random(Dermatology, "family_history")

set.seed(2023)
PBC_LR <- LASSO_Curve_Random(PBC, "class")

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Learning Curve Backup Data")
write.csv(Diabetes_LR, "Diabetes_LR.csv")
write.csv(Cardio_LR, "Cardio_LR.csv")
write.csv(Diabetes130_LR, "Diabetes130_LR.csv")
write.csv(NoShow_LR, "NoShow_LR.csv")
write.csv(BreastTumor_LR, "BreastTumor_LR.csv")
write.csv(Heart_LR, "Heart_LR.csv")
write.csv(Hepatitis_LR, "Hepatitis_LR.csv")
write.csv(Lymph_LR, "Lymph_LR.csv")
write.csv(Pharynx_LR, "Pharynx_LR.csv")
write.csv(Chol_LR, "Chol_LR.csv")
write.csv(Dermatology_LR, "Dermatology_LR.csv")
write.csv(PBC_LR, "PBC_LR.csv")
write.csv(CDC_2020_LR, "CDC_2020_LR.csv")
write.csv(CDC_2022_LR, "CDC_2022_LR.csv")
write.csv(COVID_LR, "COVID_LR.csv")
write.csv(LOS_LR, "LOS_LR.csv")

Diabetes_LR <- read.csv("Diabetes_LR.csv")
Cardio_LR <- read.csv("Cardio_LR.csv")
Diabetes130_LR <- read.csv("Diabetes130_LR.csv")
NoShow_LR <- read.csv("NoShow_LR.csv")
BreastTumor_LR <- read.csv("BreastTumor_LR.csv")
Heart_LR <- read.csv("Heart_LR.csv")
Hepatitis_LR <- read.csv("Hepatitis_LR.csv")
Lymph_LR <- read.csv("Lymph_LR.csv")
Pharynx_LR <- read.csv("Pharynx_LR.csv")
Chol_LR <- read.csv("Chol_LR.csv")
PBC_LR <- read.csv("PBC_LR.csv")
CDC_2020_LR <- read.csv("CDC_2020_LR.csv")
CDC_2022_LR <- read.csv("CDC_2022_LR.csv")
Dermatology_LR <- read.csv("Dermatology_LR.csv")
LOS_LR <- read.csv("LOS_LR.csv")
COVID_LR <- read.csv("COVID_LR.csv")




#Plots and Results ----------------------------------------------------------------------

#------------------------------

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

AUC_Full <- c(0.784, 0.654, 0.596, 0.682, 0.822, 0.642, 0.898, 0.810, 0.809, 0.949, 0.939, 0.934, 0.843, 0.669, 0.823, 0.788)*100

plot(Imbalance, log(N_LR))
plot(AUC_Full, log(N_LR))
plot(N_Features, log(N_LR))
plot(Complexity, log(N_LR))
plot(CoreFeatures, log(N_LR))
plot(PercentageBinary, log(N_LR))

summary(glm.nb(N_LR ~ Imbalance))
exp(coef(glm.nb(N_LR ~ Imbalance)))
exp(confint(glm.nb(N_LR ~ Imbalance)))

summary(glm.nb(N_LR ~ I(AUC_Full-50)))
exp(coef(glm.nb(N_LR ~ I(AUC_Full-50))))
exp(confint(glm.nb(N_LR ~ I(AUC_Full-50))))

summary(glm.nb(N_LR ~ N_Features))
exp(coef(glm.nb(N_LR ~ N_Features)))
exp(confint(glm.nb(N_LR ~ N_Features)))

summary(glm.nb(N_LR ~ CoreFeatures))
exp(coef(glm.nb(N_LR ~ CoreFeatures)))
exp(confint(glm.nb(N_LR ~ CoreFeatures)))

summary(glm.nb(N_LR ~ Complexity2))
exp(coef(glm.nb(N_LR ~ Complexity2)))
exp(confint(glm.nb(N_LR ~ Complexity2)))

summary(glm.nb(N_LR ~ I(PercentageBinary*100)))
exp(coef(glm.nb(N_LR ~ I(PercentageBinary*100))))
exp(confint(glm.nb(N_LR ~ I(PercentageBinary*100))))

summary(glm.nb(N_LR ~ Imbalance + CoreFeatures + N_Features))
 
E <- round(N_LR / N_Features, 0)

summary(glm.nb(E ~ Imbalance + CoreFeatures))

summary(glm.nb(E ~ Imbalance + AUC_Full + CoreFeatures + PercentageBinary))

#--------------------------------











