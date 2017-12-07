library(MASS)
tissue_files = c("processed_data/adipose.csv", "processed_data/muscle.csv", "processed_data/thyroid.csv", "processed_data/whole_blood.csv")
gs = seq(0, 10, 0.5)
bayesian_errors_across_tissues = matrix(0, length(tissue_files), length(gs))

for (t in 1:length(tissue_files)) {
  print(tissue_files[t])
  data_file = read.csv(tissue_files[t], row.names=1, header=TRUE)
  labels_column = dim(data_file)[2]
  train_include = -c(labels_column)
  label_include = c(labels_column)
  
  y = as.matrix(data_file[,label_include])
  X = as.matrix(cbind(rep(1, dim(y)[1], 1), data_file[,train_include]))
  X_no_intercept = as.matrix(data_file[, train_include])
  x_names = c("intercept", names(data_file)[train_include])
  
  n = dim(X)[1]
  p = dim(X)[2]
  #gs = seq(0, 10, 0.5)
  #g = n
  nu0 = 2
  s20 = 1
  
  S = 10000
  bayesian_errors = c()
  num_iterations = 10
  
  for (k in 1:length(gs)) {
    g = gs[k]
    bayesian_error = 0
    for (i in 1:num_iterations) {
      train_index = sample(n, size = n/2)
      train_data = X_no_intercept[train_index, ]
      train_labels = y[train_index, ]
      test_data = X_no_intercept[-train_index, ]
      test_labels = y[-train_index, ]
      test_data = as.matrix(cbind(rep(1, length(test_labels), 1), test_data))
      
      # Compute error for bayesian model
      train_data_with_intercept = as.matrix(cbind(rep(1, length(test_labels), 1), train_data))
      d = dim(train_data_with_intercept)[1]
      p = dim(train_data_with_intercept)[2]
      
      Hg = (g/(g + 1)) * train_data_with_intercept %*% ginv(t(train_data_with_intercept) %*%train_data_with_intercept)%*%t(train_data_with_intercept)
      SSRg = t(train_labels)%*%(diag(1,nrow=d) - Hg)%*%train_labels
      
      s2 = 1/rgamma(S, (nu0 + d)/2, (nu0*s20+SSRg)/2)
      Vb = g*ginv(t(train_data_with_intercept)%*%train_data_with_intercept)/(g+1)
      Eb = Vb%*%t(train_data_with_intercept)%*%train_labels
      
      E = matrix(rnorm(S*p, 0, sqrt(s2)), S, p)
      beta = t(t(E%*%chol(Vb, pivot = TRUE)) +c(Eb))
      
      posterior_means = c()
      for (j in 1:p) {
        post_mean = mean(beta[, j])
        posterior_means = c(posterior_means, post_mean)
      }
      
      betas = posterior_means
      predictions = test_data%*%betas
      error = 1/length(predictions) * sum((test_labels - predictions)^2)
      bayesian_error = bayesian_error + error
    }
    bayesian_error = bayesian_error/num_iterations
    cat("Bayesian error for G value ", g, ":  ",  sqrt(bayesian_error), "\n")
    bayesian_errors = c(bayesian_errors, bayesian_error)
  }
  bayesian_errors_across_tissues[t,] = bayesian_errors
}

png("figures/cv_bayes_gprior.png", width=8, height=4, units="in", res=300)
plot(gs, sqrt(bayesian_errors_across_tissues[1,]), 
     type='l', main="Cross-Validation Results for G Prior in Bayesian Regression", 
     xlab="G Value Prior", ylab="RMSE (Years)", col="blue")
lines(gs, sqrt(bayesian_errors_across_tissues[2,]), col="darkorange")
lines(gs, sqrt(bayesian_errors_across_tissues[3,]), col="darkgreen")
lines(gs, sqrt(bayesian_errors_across_tissues[4,]), col="red")
legend("topright", c("adipose","muscle", "thyroid", "whole_blood"),
       fill=c("blue","darkorange", "darkgreen", "red"))
dev.off()