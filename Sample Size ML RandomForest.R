#Random Forest

set.seed(2023)
RF_Optimal(Cardio, "CVD")
RF_Optimal(Cardio, "CVD")[2]*2.25

set.seed(2023)
RF_Optimal(BreastTumor, "recurrence")
RF_Optimal(BreastTumor, "recurrence")[2]*2.25

set.seed(2023)
RF_Optimal(NoShow, "NoShow")
RF_Optimal(NoShow, "NoShow")[2]*2.25

set.seed(2023)
RF_Optimal(Diabetes130, "readmitted")
RF_Optimal(Diabetes130, "readmitted")[2]*2.25

set.seed(2023)
RF_Optimal(Diabetes, "Diabetes_binary")
RF_Optimal(Diabetes, "Diabetes_binary")[2]*2.25

set.seed(2023)
RF_Optimal(COVID, "RESULTADO")
RF_Optimal(COVID, "RESULTADO")[2]*2.25

set.seed(2023)
RF_Optimal(LOS, "LOS")
RF_Optimal(LOS, "LOS")[2]*2.25

set.seed(2023)
RF_Optimal(CDC_2020, "HeartDisease")
RF_Optimal(CDC_2020, "HeartDisease")[2]*2.25

set.seed(2023)
RF_Optimal(CDC_2022, "HighRiskLastYear")
RF_Optimal(CDC_2022, "HighRiskLastYear")[2]*2.25

set.seed(2023)
RF_Optimal(HeartSynthetic, "class")
RF_Optimal(HeartSynthetic, "class")[2]*2.25

set.seed(2023)
RF_Optimal(Hepatitis, "Class")
RF_Optimal(Hepatitis, "Class")[2]*2.25

set.seed(2023)
RF_Optimal(Lymph, "class")
RF_Optimal(Lymph, "class")[2]*2.25

set.seed(2023)
RF_Optimal(Pharynx, "class")
RF_Optimal(Pharynx, "class")[2]*2.25

set.seed(2023)
RF_Optimal(Chol, "chol")
RF_Optimal(Chol, "chol")[2]*2.25

set.seed(2023)
RF_Optimal(Dermatology, "family_history")
RF_Optimal(Dermatology, "family_history")[2]*2.25

set.seed(2023)
RF_Optimal(PBC, "class")
RF_Optimal(PBC, "class")[2]*2.25

#-----------------------------------------------------------------------------

set.seed(2023)
Cardio_RF <- RF_Curve_Random(Cardio, "CVD", 50000)

set.seed(2023)
BreastTumor_RF <- RF_Curve_Random(BreastTumor, "recurrence", 50000)

set.seed(2023)
NoShow_RF <- RF_Curve_Random(NoShow, "NoShow", 50000)

set.seed(2023)
Diabetes130_RF <- RF_Curve_Random(Diabetes130, "readmitted", 50000)

set.seed(2023)
Diabetes_RF <- RF_Curve_Random(Diabetes, "Diabetes_binary", 50000)

set.seed(2023)
COVID_RF <- RF_Curve_Random(COVID, "RESULTADO", 50000)

set.seed(2024)
LOS_RF <- RF_Curve_Random(LOS, "LOS", 50000)

set.seed(2023)
CDC_2020_RF <- RF_Curve_Random(CDC_2020, "HeartDisease", 50000)

set.seed(2023)
CDC_2022_RF <- RF_Curve_Random(CDC_2022, "HighRiskLastYear", 50000)

set.seed(2023)
Heart_RF <- RF_Curve_Random(HeartSynthetic, "class", 50000)

set.seed(2023)
Hepatitis_RF <- RF_Curve_Random(Hepatitis, "Class", 50000)

set.seed(2023)
Lymph_RF <- RF_Curve_Random(Lymph, "class", 50000)

set.seed(2023)
Pharynx_RF <- RF_Curve_Random(Pharynx, "class", 50000)

set.seed(2023)
Chol_RF <- RF_Curve_Random(Chol, "chol", 200000)

set.seed(2023)
Dermatology_RF <- RF_Curve_Random(Dermatology, "family_history", 50000)

set.seed(2023)
PBC_RF <- RF_Curve_Random(PBC, "class", 100000)

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Learning Curve Backup Data")
write.csv(Diabetes_RF, "Diabetes_RF.csv")
write.csv(Cardio_RF, "Cardio_RF.csv")
write.csv(Diabetes130_RF, "Diabetes130_RF.csv")
write.csv(NoShow_RF, "NoShow_RF.csv")
write.csv(BreastTumor_RF, "BreastTumor_RF.csv")
write.csv(Heart_RF, "Heart_RF.csv")
write.csv(Hepatitis_RF, "Hepatitis_RF.csv")
write.csv(Lymph_RF, "Lymph_RF.csv")
write.csv(Pharynx_RF, "Pharynx_RF.csv")
write.csv(Chol_RF, "Chol_RF.csv")
write.csv(Dermatology_RF, "Dermatology_RF.csv")
write.csv(PBC_RF, "PBC_RF.csv")
write.csv(CDC_2020_RF, "CDC_2020_RF.csv")
write.csv(CDC_2022_RF, "CDC_2022_RF.csv")
write.csv(COVID_RF, "COVID_RF.csv")
write.csv(LOS_RF, "LOS_RF.csv")


Diabetes_RF <- read.csv("Diabetes_RF.csv")
Cardio_RF <- read.csv("Cardio_RF.csv")
Diabetes130_RF <- read.csv("Diabetes130_RF.csv")
NoShow_RF <- read.csv("NoShow_RF.csv")
BreastTumor_RF <- read.csv("BreastTumor_RF.csv")
Heart_RF <- read.csv("Heart_RF.csv")
Hepatitis_RF <- read.csv("Hepatitis_RF.csv")
Lymph_RF <- read.csv("Lymph_RF.csv")
Pharynx_RF <- read.csv("Pharynx_RF.csv")
Chol_RF <- read.csv("Chol_RF.csv")
PBC_RF <- read.csv("PBC_RF.csv")
CDC_2020_RF <- read.csv("CDC_2020_RF.csv")
CDC_2022_RF <- read.csv("CDC_2022_RF.csv")
Dermatology_RF <- read.csv("Dermatology_RF.csv")
LOS_RF <- read.csv("LOS_RF.csv")
COVID_RF <- read.csv("COVID_RF.csv")



#------------------------------------------------------------


#----------------------------------------

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

AUC_Full <- c(0.796, 0.661, 0.609, 0.780, 0.822, 0.661, 0.915, 0.810, 0.801, 0.963, 0.976, 0.957, 0.858, 0.728, 0.857, 0.850)*100



plot(Imbalance, log(N_RF))
plot(AUC_Full, log(N_RF))
plot(N_Features, log(N_RF))
boxplot(Complexity, log(N_RF))
plot(CoreFeatures, log(N_RF))
plot(PercentageBinary, log(N_RF))

summary(glm.nb(N_RF ~ Imbalance))
exp(coef(glm.nb(N_RF ~ Imbalance)))
exp(confint(glm.nb(N_RF ~ Imbalance)))

summary(glm.nb(N_RF ~ I(AUC_Full-50)))
exp(coef(glm.nb(N_RF ~ I(AUC_Full-50))))
exp(confint(glm.nb(N_RF ~ I(AUC_Full-50))))

summary(glm.nb(N_RF ~ N_Features))
exp(coef(glm.nb(N_RF ~ N_Features)))
exp(confint(glm.nb(N_RF ~ N_Features)))

summary(glm.nb(N_RF ~ CoreFeatures))
exp(coef(glm.nb(N_RF ~ CoreFeatures)))
exp(confint(glm.nb(N_RF ~ CoreFeatures)))

summary(glm.nb(N_RF ~ Complexity2))
exp(coef(glm.nb(N_RF ~ Complexity2)))
exp(confint(glm.nb(N_RF ~ Complexity2)))

summary(glm.nb(N_RF ~ I(PercentageBinary*100)))
exp(coef(glm.nb(N_RF ~ I(PercentageBinary*100))))
exp(confint(glm.nb(N_RF ~ I(PercentageBinary*100))))

summary(glm.nb(N_RF ~ Imbalance + Complexity + I(AUC_Full-50)))

RsqGLM(glm.nb(N_RF ~ Imbalance + Complexity + I(AUC_Full-50)))


#Plots and Predictions




