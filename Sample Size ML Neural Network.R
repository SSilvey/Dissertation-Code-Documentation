library(h2o)
h2o.init()
h2o.shutdown()
h2o.no_progress()
h2o.clearLog()

set.seed(2023)
NN_Optimal(Cardio, "CVD")

set.seed(2023)
NN_Optimal(BreastTumor, "recurrence")

set.seed(2023)
NN_Optimal(NoShow, "NoShow")

set.seed(2023)
NN_Optimal(Diabetes130, "readmitted")

set.seed(2023)
NN_Optimal(Diabetes, "Diabetes_binary")

set.seed(2023)
NN_Optimal(COVID, "RESULTADO")

set.seed(2023)
NN_Optimal(LOS, "LOS")

set.seed(2023)
NN_Optimal(CDC_2020, "HeartDisease")

set.seed(2023)
NN_Optimal(CDC_2022, "HighRiskLastYear")

set.seed(2023)
NN_Optimal(HeartSynthetic, "class")

set.seed(2023)
NN_Optimal(Hepatitis, "Class")

set.seed(2023)
NN_Optimal(Lymph, "class")

set.seed(2023)
NN_Optimal(Pharynx, "class")

set.seed(2023)
NN_Optimal(Chol, "chol")

set.seed(2023)
NN_Optimal(Dermatology, "V34")

set.seed(2023)
NN_Optimal(PBC, "class")
#Functions


#-----------------------------------------

S <- seq(500, 50000, by=45000/10)
S2 <- seq(500, 200000, by=195000/10)
S3 <- seq(500, 250000, by=245000/10)

Diabetes_NN <- NULL
for (i in S) {
  Diabetes_NN <- rbind(Diabetes_NN, NN_Curve_Evaluate(Diabetes, "Diabetes_binary", i))
}
Diabetes_NN

Cardio_NN <- NULL
for (i in S) {
  Cardio_NN <- rbind(Cardio_NN, NN_Curve_Evaluate(Cardio, "CVD", i))
}
Cardio_NN

Diabetes130_NN <- NULL
for (i in S) {
  Diabetes130_NN <- rbind(Diabetes130_NN, NN_Curve_Evaluate(Diabetes130, "readmitted", i))
}
Diabetes130_NN

NoShow_NN <- NULL
for (i in S) {
  NoShow_NN <- rbind(NoShow_NN, NN_Curve_Evaluate(NoShow, "NoShow", i))
}
NoShow_NN

BreastTumor_NN <- NULL
for (i in S) {
  BreastTumor_NN <- rbind(BreastTumor_NN, NN_Curve_Evaluate(BreastTumor, "recurrence", i))
}
BreastTumor_NN

COVID_NN <- NULL
for (i in S) {
  COVID_NN <- rbind(COVID_NN, NN_Curve_Evaluate(COVID, "RESULTADO", i))
}
COVID_NN

CDC_2020_NN <- NULL
for (i in S) {
  CDC_2020_NN <- rbind(CDC_2020_NN, NN_Curve_Evaluate(CDC_2020, "HeartDisease", i))
}
CDC_2020_NN

CDC_2022_NN <- NULL
for (i in S) {
  CDC_2022_NN <- rbind(CDC_2022_NN, NN_Curve_Evaluate(CDC_2022, "HighRiskLastYear", i))
}
CDC_2022_NN

Heart_NN <- NULL
for (i in S) {
  Heart_NN <- rbind(Heart_NN, NN_Curve_Evaluate(HeartSynthetic, "class", i))
}
Heart_NN

Hepatitis_NN <- NULL
for (i in S) {
  Hepatitis_NN <- rbind(Hepatitis_NN, NN_Curve_Evaluate(Hepatitis, "Class", i))
}
Hepatitis_NN

Lymph_NN <- NULL
for (i in S) {
  Lymph_NN <- rbind(Lymph_NN, NN_Curve_Evaluate(Lymph, "class", i))
}
Lymph_NN

Pharynx_NN <- NULL
for (i in S) {
  Pharynx_NN <- rbind(Pharynx_NN, NN_Curve_Evaluate(Pharynx, "class", i))
}
Pharynx_NN

Dermatology_NN <- NULL
for (i in S) {
  Dermatology_NN <- rbind(Dermatology_NN, NN_Curve_Evaluate(Dermatology, "V34", i))
}
Dermatology_NN

PBC_NN <- NULL
for (i in S2) {
  PBC_NN <- rbind(PBC_NN, NN_Curve_Evaluate(PBC, "class", i))
}
PBC_NN

#Run these
LOS_NN <- NULL
for (i in S) {
  LOS_NN <- rbind(LOS_NN, NN_Curve_Evaluate(LOS, "LOS", i))
}
LOS_NN

Chol_NN <- NULL
for (i in S3) {
  Chol_NN <- rbind(Chol_NN, NN_Curve_Evaluate(Chol, "chol", i))
}
Chol_NN


#---------------------

setwd("C:/Users/scott/OneDrive/Desktop/Dissertation/Learning Curve Backup Data")
write.csv(Diabetes_NN, "Diabetes_NN.csv")
write.csv(Cardio_NN, "Cardio_NN.csv")
write.csv(Diabetes130_NN, "Diabetes130_NN.csv")
write.csv(NoShow_NN, "NoShow_NN.csv")
write.csv(BreastTumor_NN, "BreastTumor_NN.csv")
write.csv(Heart_NN, "Heart_NN.csv")
write.csv(Hepatitis_NN, "Hepatitis_NN.csv")
write.csv(Lymph_NN, "Lymph_NN.csv")
write.csv(Pharynx_NN, "Pharynx_NN.csv")
write.csv(Chol_NN, "Chol_NN.csv")
write.csv(Dermatology_NN, "Dermatology_NN.csv")
write.csv(PBC_NN, "PBC_NN.csv")
write.csv(CDC_2020_NN, "CDC_2020_NN.csv")
write.csv(CDC_2022_NN, "CDC_2022_NN.csv")
write.csv(COVID_NN, "COVID_NN.csv")
write.csv(LOS_NN, "LOS_NN.csv")

Diabetes_NN <- read.csv("Diabetes_NN.csv")
Cardio_NN <- read.csv("Cardio_NN.csv")
Diabetes130_NN <- read.csv("Diabetes130_NN.csv")
NoShow_NN <- read.csv("NoShow_NN.csv")
BreastTumor_NN <- read.csv("BreastTumor_NN.csv")
Heart_NN <- read.csv("Heart_NN.csv")
Hepatitis_NN <- read.csv("Hepatitis_NN.csv")
Lymph_NN <- read.csv("Lymph_NN.csv")
Pharynx_NN <- read.csv("Pharynx_NN.csv")
Chol_NN <- read.csv("Chol_NN.csv")
PBC_NN <- read.csv("PBC_NN.csv")
CDC_2020_NN <- read.csv("CDC_2020_NN.csv")
CDC_2022_NN <- read.csv("CDC_2022_NN.csv")
Dermatology_NN <- read.csv("Dermatology_NN.csv")
LOS_NN <- read.csv("LOS_NN.csv")
COVID_NN <- read.csv("COVID_NN.csv")

#-------------------------------------------------------

AUC_Full <- c(0.795, 0.661, 0.603, 0.730, 0.826, 0.661, 0.901, 0.809, 0.801, 0.962, 0.974, 0.956, 0.856, 0.714, 0.852, 0.823)*100

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

plot(Imbalance, log(N_NN))
plot(AUC_Full, log(N_NN))
plot(N_Features, log(N_NN))
plot(Complexity, log(N_NN))
plot(CoreFeatures, log(N_NN))
plot(PercentageBinary, log(N_NN))

summary(glm.nb(N_NN ~ Imbalance))
exp(coef(glm.nb(N_NN ~ Imbalance)))
exp(confint(glm.nb(N_NN ~ Imbalance)))

summary(glm.nb(N_NN ~ I(AUC_Full-50)))
exp(coef(glm.nb(N_NN ~ I(AUC_Full-50))))
exp(confint(glm.nb(N_NN ~ I(AUC_Full-50))))

summary(glm.nb(N_NN ~ N_Features))
exp(coef(glm.nb(N_NN ~ N_Features)))
exp(confint(glm.nb(N_NN ~ N_Features)))

summary(glm.nb(N_NN ~ CoreFeatures))
exp(coef(glm.nb(N_NN ~ CoreFeatures)))
exp(confint(glm.nb(N_NN ~ CoreFeatures)))

summary(glm.nb(N_NN ~ Complexity2))
exp(coef(glm.nb(N_NN ~ Complexity2)))
exp(confint(glm.nb(N_NN ~ Complexity2)))

summary(glm.nb(N_NN ~ I(PercentageBinary*100)))
exp(coef(glm.nb(N_NN ~ I(PercentageBinary*100))))
exp(confint(glm.nb(N_NN ~ I(PercentageBinary*100))))

summary(glm.nb(N_NN ~ Imbalance + I(AUC_Full-50) + Complexity))

summary(glm.nb(N_NN ~ Imbalance + N_Features + Complexity))

summary(glm.nb(N_NN ~ Imbalance + N_Features + Complexity))


