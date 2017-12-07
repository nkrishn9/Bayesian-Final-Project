library(ridge)
lambdas = c(.001, .01, 0.1, 1, 10, 100)

adipose_train = as.data.frame(read.csv("processed_data/adipose_train.csv", row.names=1, header=TRUE))
muscle_train = as.data.frame(read.csv("processed_data/muscle_train.csv", row.names=1, header=TRUE))
thyroid_train = as.data.frame(read.csv("processed_data/thyroid_train.csv", row.names=1, header=TRUE))
whole_blood_train = as.data.frame(read.csv("processed_data/whole_blood_train.csv", row.names=1, header=TRUE))

cv_results = data.frame(row.names=c("adipose", "muscle", "thyroid", "whole_blood"))

for (i in lambdas) {
  k=5
  folds = createFolds(adipose_train$AGE, k=5)
  error = 0
  for (j in 1:k) {
    cv_train = adipose_train[unlist(folds[-j]), ]
    cv_test = adipose_train[unlist(folds[j]), ]
    adipose_model = linearRidge(AGE ~ ., data=data.frame(cv_train), lambda=i)
    predictions = predict(adipose_model, cv_test[, -dim(cv_test)[2]])
    error = error + sqrt(sum((predictions - cv_test[, dim(cv_test)[2]])^2) / length(predictions))
  }
  error = error / k 
  cv_results["adipose", toString(i)] = error
  
  folds = createFolds(muscle_train$AGE, k=5)
  error = 0
  for (j in 1:k) {
    cv_train = muscle_train[unlist(folds[-j]), ]
    cv_test = muscle_train[unlist(folds[j]), ]
    muscle_model = linearRidge(AGE ~ ., data=data.frame(cv_train), lambda=i)
    predictions = predict(muscle_model, cv_test[, -dim(cv_test)[2]])
    error = error + sqrt(sum((predictions - cv_test[, dim(cv_test)[2]])^2) / length(predictions))
  }
  error = error / k 
  cv_results["muscle", toString(i)] = error
  
  folds = createFolds(thyroid_train$AGE, k=5)
  error = 0
  for (j in 1:k) {
    cv_train = thyroid_train[unlist(folds[-j]), ]
    cv_test = thyroid_train[unlist(folds[j]), ]
    thyroid_model = linearRidge(AGE ~ ., data=data.frame(cv_train), lambda=i)
    predictions = predict(thyroid_model, cv_test[, -dim(cv_test)[2]])
    error = error + sqrt(sum((predictions - cv_test[, dim(cv_test)[2]])^2) / length(predictions))
  }
  error = error / k 
  cv_results["thyroid", toString(i)] = error
  
  folds = createFolds(whole_blood_train$AGE, k=5)
  error = 0
  for (j in 1:k) {
    cv_train = whole_blood_train[unlist(folds[-j]), ]
    cv_test = whole_blood_train[unlist(folds[j]), ]
    whole_blood_model = linearRidge(AGE ~ ., data=data.frame(cv_train), lambda=i)
    predictions = predict(whole_blood_model, cv_test[, -dim(cv_test)[2]])
    error = error + sqrt(sum((predictions - cv_test[, dim(cv_test)[2]])^2) / length(predictions))
  }
  error = error / k 
  cv_results["whole_blood", toString(i)] = error
}

png("figures/cv.png", width=6, height=4, units="in", res=300)
plot(log(lambdas), cv_results[1, ], 
     type='l', main="Cross-Validation Results for Lambda in Ridge Regression", 
     xlab="Log(lambdas)", ylab="RMSE (Years)", col="blue", ylim=c(10, 25))
lines(log(lambdas), cv_results[2, ], col="darkorange")
lines(log(lambdas), cv_results[3, ], col="darkgreen")
lines(log(lambdas), cv_results[4, ], col="red")
legend("topright", c("adipose","muscle", "thyroid", "whole_blood"), 
       fill=c("blue","darkorange", "darkgreen", "red"))
dev.off()
