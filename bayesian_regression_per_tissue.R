library(MASS)

data_file = read.csv("processed_data/for_emily.csv", row.names=1, header=TRUE)
labels_column = dim(data_file)[2]
train_include = -c(labels_column)
label_include = c(labels_column)

# data_file = read.table("crime.dat", header = TRUE)
# train_include = -c(1)
# label_include = c(1)

y = as.matrix(data_file[,label_include])
X = as.matrix(cbind(rep(1, dim(y)[1], 1), data_file[,train_include]))
X_no_intercept = as.matrix(data_file[, train_include])
x_names = c("intercept", names(data_file)[train_include])

n = dim(X)[1]
p = dim(X)[2]
gs = seq(1, n, 10)
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
  print(k)
  cat("Bayesian error: ", bayesian_error, "\n")
  bayesian_error = bayesian_error/num_iterations
  bayesian_errors = c(bayesian_errors, bayesian_error)
}