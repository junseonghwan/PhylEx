library(matrixStats)
library(TailRank)

options(digits = 10)

seq_err <- 0.001
phi <- 0.5
var_read_count <- 10
total_read_count <- 20
var_cp_prob <- 0.5

# Bulk log likelihood with total copy number.
log_lik_total_cn <- function(b, d, phi, seq_err, var_cp_prob, total_cn) {
  if (d == 0) {
    return(0)
  }
  if (total_cn == 0) {
    return(dbinom(b, d, prob = seq_err, log=TRUE))
  }
  log_vals <- rep(0, total_cn)
  log_genotype_prior <- rep(0, total_cn)
  log_genotype_norm <- log(1 - dbinom(0, total_cn, var_cp_prob))
  for (var_cn in 1:total_cn) {
    xi <- (1 - phi) * seq_err
    if (var_cn == total_cn) {
      xi <- xi + phi * (1 - seq_err)
    } else {
      xi <- xi + phi * var_cn/total_cn
    }
    log_vals[var_cn] <- dbinom(b, d, prob = xi, log=TRUE)
    log_genotype_prior[var_cn] = dbinom(var_cn, total_cn, prob = var_cp_prob, log=TRUE) - log_genotype_norm
  }

  logSumExp(log_vals + log_genotype_prior)
}

log_lik_genotype <- function(b, d, phi, seq_err, major_cn, minor_cn) {
  total_cn <- major_cn + minor_cn
  log_liks <- rep(0, total_cn)
  log_lik_prior <- rep(log(1/total_cn), total_cn)
  for (var_cn in 1:minor_cn) {
    xi <- (1 - phi) * seq_err
    if (var_cn == total_cn) {
      xi <- xi + phi * (1 - seq_err)
    } else {
      xi <- xi + phi * var_cn/total_cn
    }
    log_liks[var_cn] <- dbinom(b, d, prob = xi, log=TRUE)
  }
  for (var_cn in 1:major_cn) {
    xi <- (1 - phi) * seq_err
    if (var_cn == total_cn) {
      xi <- xi + phi * (1 - seq_err)
    } else {
      xi <- xi + phi * var_cn/total_cn
    }
    log_liks[var_cn+minor_cn] <- dbinom(b, d, prob = xi, log=TRUE)
  }
  logSumExp(log_liks + log_lik_prior)
}

log_lik_sc <- function(b, d, delta0, alpha, beta, alpha0, beta0, has_snv, seq_err = 0.001) {
  if (d == 0) {
    return(0);
  }
  if (has_snv) {
    log_val1 <- dbb(b, d, alpha, beta, log = TRUE) + log(1 - delta0)
    log_val2 <- dbb(b, d, alpha0, beta0, log = TRUE) + log(delta0)
    log_lik <- logSumExp(c(log_val1, log_val2))
    return(log_lik)
  } else {
    return(dbb(b, d, seq_err, 1 - seq_err, log = TRUE))
  }
}

total_cn <- 0
log_lik_total_cn(var_read_count, total_read_count, phi, seq_err, var_cp_prob, total_cn)
total_cn <- 2
log_lik_total_cn(var_read_count, total_read_count, phi, seq_err, var_cp_prob, total_cn)
total_cn <- 8
log_lik_total_cn(var_read_count, total_read_count, phi, seq_err, var_cp_prob, total_cn)

# Possible genotypes:
# (v, rrr), 3x(r, vrr), 3x(r, vvr), 1x(r, vvv).
major_cn <- 3
minor_cn <- 1
log_lik_genotype(var_read_count, total_read_count, phi, seq_err, major_cn, minor_cn)

major_cn <- 2
minor_cn <- 0
log_lik_genotype(var_read_count, total_read_count, phi, seq_err, major_cn, minor_cn)

b <- 10
d <- 20
alpha <- 3
beta <- 12
delta0 <- 0.8
alpha0 <- 0.01
beta0 <- 0.01
seq_err <- 0.001
log_lik_sc(b, d, delta0, alpha, beta, alpha0, beta0, TRUE, seq_err)
log_lik_sc(b, d, delta0, alpha, beta, alpha0, beta0, FALSE, seq_err)

log_lik_genotype(b = 277, d = 1032, phi = 0.56555, seq_err = 0.001, major_cn = 1, minor_cn = 1)
log_lik_genotype(b = 277, d = 1032, phi = 0.6553, seq_err = 0.001, major_cn = 1, minor_cn = 1)
log_lik_genotype(b = 277, d = 1032, phi = 0.718629, seq_err = 0.001, major_cn = 1, minor_cn = 1)
