### Luca Bigonzi

# Import Libraries
library(boot)
library(bootstrap)

# Number of Eruptions
volcano <- c(0, 0, 4, 3, 0, 5, 2, 0, 1, 4, 3, 9, 0, 3, 0)

# A) Computing a 95% Confidence Interval for the eruptions of the Volcanoes

# Taking into account the nature of the data (positive or null) we opt for a
# coherent confidence interval (can't be negative). The available
# methods are percentile and BCa. The BCa is chosen due to a 
# better generalization

theta_mean <- function(data,i){
  mean(data[i])
}
set.seed(1)
boot_a <- boot(volcano, theta_mean, R = 10000)
conf_int <- boot.ci(boot_a , R = 10000 , conf = 0.95 , type = "bca")
conf_int # 95% (1.2 , 3.8)

# B) Reconstruct the distribution of the variance through a NON Parametric
# Bootstrap Analysis

theta_var <- function(data,i){
  var(data[i])
}
set.seed(1)
boot_b <- boot(volcano, theta_var, R = 10000)
stima_boot <- boot_b$t
#estimate and bias
est_b <- mean(stima_boot)
est_b
TF.var <- var(volcano)
bias <- est_b - TF.var
bias
## Estimate = 6.208524
## Bias = -0.4295714

# C) PARAMETRIC BOOTSTRAP -- POISSON

c_var <- function(data){
  var(data) 
}
set.seed(1)
rp <- function(data, lambda){
  out <- rpois(length(data), lambda = lambda)
  out
}

boot_pois <- boot(volcano, c_var , R = 10000, sim = "parametric" , ran.gen = rp, mle = mean(volcano)) # the mle of a Poisson is the mean

boot_c <- boot_pois$t
# Estimate and Bias
est_c <- mean(boot_c)
est_c
bias_c <- est_c - TF.var
bias_c
## Estimate = 2.264464
## Bias = -4.373631

# D) Comparison between NonParametric and Parametric Bootstrap

# The NonParametric is more reliable because its Bias is smaller 
# than the Parametric (-0.4295714 vs -4.373631).
# Hence, it seems better to rely on the NonParametric Bootstrap Distribution

# Moreover, the distribution of bootstrap estimated variances is better
# centered (as coherent with bias) with our sample value
par(mfrow = c(1,2))
plot(density(stima_boot), col = 'green', ylim = c(0,0.5), main = 'Bootstrap densities') 
lines(density(boot_c), col = 'red')
abline(v = var(volcano), col = 'blue', lty = 'dashed')

# In addition, the ecdf and theorical poisson distributions are different:
plot(ecdf(volcano), main = 'ecdf vs  theorical poisson')
curve(ppois(volcano, lambda = mean(volcano)), add = TRUE, col = 'red')

# E) Bootstrap Hypothesis Testing:
# Is the Poisson Distribution suitable in order to describe the data?

# H0: h(Y) = Var(Y) - E(Y) = 0

# h(Y) = Var(Y) - E(Y)
h_y <- var(volcano) - mean(volcano)

statistic_e <- function(data, i){
  sel <- data[i]
  var(sel) - mean(sel)
}

mean(volcano) # 2.266667
var(volcano) # 6.638095
volcano2 <- volcano + (var(volcano) - mean(volcano))
mean(volcano2) # 6.638095
var(volcano2) # 6.638095

set.seed(1)
bootstrap_e <- boot(volcano2, statistic_e, R = 10000)
b_stat_e <- bootstrap_e$t

# We are interested in a two tail test in order to detect
# any difference between mean and variance (in any direction)
n <- length(b_stat_e)
(pvalue <- sum(abs(b_stat_e) >= abs(h_y))/n)

# Due to the fact that the p-value is 0.0427 < 0.05, we can say that the 
# Poisson probabilistic model is NOT suitable for describing the data,
# hence we reject the Null Hypothesis.
# This also means that the mean is not equal to the variance