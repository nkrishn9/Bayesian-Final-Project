library(ridge)

adipose_train = as.data.frame(read.csv("processed_data/adipose_train.csv", row.names=1, header=TRUE))
muscle_train = as.data.frame(read.csv("processed_data/muscle_train.csv", row.names=1, header=TRUE))
thyroid_train = as.data.frame(read.csv("processed_data/thyroid_train.csv", row.names=1, header=TRUE))
whole_blood_train = as.data.frame(read.csv("processed_data/whole_blood_train.csv", row.names=1, header=TRUE))

adipose_test = as.data.frame(read.csv("processed_data/adipose_test.csv", row.names=1, header=TRUE))
muscle_test = as.data.frame(read.csv("processed_data/muscle_test.csv", row.names=1, header=TRUE))
thyroid_test = as.data.frame(read.csv("processed_data/thyroid_test.csv", row.names=1, header=TRUE))
whole_blood_test = as.data.frame(read.csv("processed_data/whole_blood_test.csv", row.names=1, header=TRUE))

adipose_model = linearRidge(AGE ~ ., data=data.frame(adipose_train), lambda=1)
muscle_model = linearRidge(AGE ~ ., data=data.frame(muscle_train), lambda=1)
thyroid_model = linearRidge(AGE ~ ., data=data.frame(thyroid_train), lambda=1)
whole_blood_model = linearRidge(AGE ~ ., data=data.frame(whole_blood_train), lambda=1)

results = data.frame(row.names=c("adipose", "muscle", "thyroid", "whole_blood"))

predictions = predict(adipose_model, adipose_test[, -dim(adipose_test)[2]])
results["adipose", "rmse"] = sqrt(sum((predictions - adipose_test[, dim(adipose_test)[2]])^2) / length(predictions))
box_adipose = sqrt(((predictions - adipose_test[, dim(adipose_test)[2]])^2))

predictions = predict(muscle_model, muscle_test[, -dim(muscle_test)[2]])
results["muscle", "rmse"] = sqrt(sum((predictions - muscle_test[, dim(muscle_test)[2]])^2) / length(predictions))
box_muscle = sqrt(((predictions - muscle_test[, dim(muscle_test)[2]])^2))

predictions = predict(thyroid_model, thyroid_test[, -dim(thyroid_test)[2]])
results["thyroid", "rmse"] = sqrt(sum((predictions - thyroid_test[, dim(thyroid_test)[2]])^2) / length(predictions))
box_thyroid = sqrt(((predictions - thyroid_test[, dim(thyroid_test)[2]])^2))

predictions = predict(whole_blood_model, whole_blood_test[, -dim(whole_blood_test)[2]])
results["whole_blood", "rmse"] = sqrt(sum((predictions - whole_blood_test[, dim(whole_blood_test)[2]])^2) / length(predictions))
box_whole_blood = sqrt(((predictions - whole_blood_test[, dim(whole_blood_test)[2]])^2))

png("figures/freq_results.png", width=6, height=4, units="in", res=300)
boxplot(box_adipose, box_muscle, box_thyroid, box_whole_blood, 
        names=c("adipose","muscle", "thyroid", "whole_blood"), ylab="Prediction Error (Years)", 
        main="Prediction Error per Tissue Using Ridge Regression")
dev.off()

