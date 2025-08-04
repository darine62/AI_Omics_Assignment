getwd()
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")
dir.create("Tasks")

data <- read.csv(file.choose())
View(data)
str(data)

data$gender_fac <- as.factor(data$gender)
data$diagnosis_fac <- as.factor(data$diagnosis)
data$smoker_fac <- as.factor(data$smoker)
str(data)
View(data)

data$smoker_num <- ifelse(data$smoker_fac == "Yes", 1, 0)
View(data)

write.csv(data, "clean_data/patient_info_clean.csv")
save.image(file = "class_Ib.RData")
