#CHAPTER 6: Statistical Inference 

# Estimators

variance <- function(x, y, estimator) {
  mean <- sum(x*y)/sum(y)
  num <- sum((x-mean)**2)
  N <- sum(y)
  if (estimator) {
    N <- N-1
  }
  num/N
}



# When using class marks
breaks <- seq(0,100, by=10)
mids <- (breaks[-length(breaks)] + breaks[-1])/2 # xi but cool



# if sum(p) != 1 
normalise <- function(p) {
  if (abs(sum(p) - 1) >= 1e-6) {
    p <- p/sum(p)
  }
}


#CHAPTER 6: Statistical Inference 

# ------------------------------------------------------------------------------
# SECTION 1: 
# POPULATION MEAN (Mu)
# ------------------------------------------------------------------------------

#' Confidence Interval for Mean (Known Variance)
#' 
#' @param x_bar Sample mean
#' @param sigma Population standard deviation (Known)
#' @param n Sample size
#' @param conf Confidence level (default 0.95)
ci_mean_known_var <- function(x_bar, sigma, n, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2) # Critical value Z
  
  margin_error <- z * (sigma / sqrt(n))
  
  lower <- x_bar - margin_error
  upper <- x_bar + margin_error
  
  cat(sprintf("Confidence Level: %.2f%%\n", conf*100))
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

#' Confidence Interval for Mean (Unknown Variance - Large Sample n > 30)
#' 
#' @param x_bar Sample mean
#' @param s Sample standard deviation (Quasi-variance root -> sqrt(var()))
#' @param n Sample size (Must be > 30)
#' @param conf Confidence level
ci_mean_unknown_var_large <- function(x_bar, s, n, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2) # Uses Z distribution because n is large
  
  margin_error <- z * (s / sqrt(n))
  
  lower <- x_bar - margin_error
  upper <- x_bar + margin_error
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

#' Confidence Interval for Mean (Unknown Variance - Small Sample n <= 30)
#' 
#' @param x_bar Sample mean
#' @param s Sample standard deviation
#' @param n Sample size (Must be <= 30)
#' @param conf Confidence level
ci_mean_unknown_var_small <- function(x_bar, s, n, conf = 0.95) {
  alpha <- 1 - conf
  # Uses T-Student distribution with n-1 degrees of freedom
  t_val <- qt(1 - alpha/2, df = n - 1) 
  
  margin_error <- t_val * (s / sqrt(n))
  
  lower <- x_bar - margin_error
  upper <- x_bar + margin_error
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

# ------------------------------------------------------------------------------
# POPULATION VARIANCE (Sigma^2)
# ------------------------------------------------------------------------------

#' Confidence Interval for Variance
#' 
#' @param s2 Sample Quasivariance (s^2)
#' @param n Sample size
#' @param conf Confidence level
ci_variance <- function(s2, n, conf = 0.95) {
  alpha <- 1 - conf
  
  # Chi-square critical values
  chi_upper <- qchisq(1 - alpha/2, df = n - 1)
  chi_lower <- qchisq(alpha/2, df = n - 1)
  
  # Formula: [(n-1)s^2 / Chi_upper, (n-1)s^2 / Chi_lower]
  lower <- ((n - 1) * s2) / chi_upper
  upper <- ((n - 1) * s2) / chi_lower
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

# ------------------------------------------------------------------------------
# POPULATION PROPORTION (p)
# ------------------------------------------------------------------------------

#' Confidence Interval for Proportion
#' 
#' @param p_hat Sample proportion (successes / n)
#' @param n Sample size (Assumes n > 30)
#' @param conf Confidence level
ci_proportion <- function(p_hat, n, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  
  margin_error <- z * sqrt((p_hat * (1 - p_hat)) / n)
  
  lower <- p_hat - margin_error
  upper <- p_hat + margin_error
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

# ------------------------------------------------------------------------------
# DIFFERENCE OF MEANS (Mu1 - Mu2)
# ------------------------------------------------------------------------------

#' CI Difference of Means (Known Variances)
#' 
#' @param mean1 Sample mean 1
#' @param mean2 Sample mean 2
#' @param var1 Population variance 1 (Sigma1^2)
#' @param var2 Population variance 2 (Sigma2^2)
#' @param n1 Sample size 1
#' @param n2 Sample size 2
ci_diff_means_known_vars <- function(mean1, mean2, var1, var2, n1, n2, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  
  diff <- mean1 - mean2
  se <- sqrt((var1 / n1) + (var2 / n2))
  
  lower <- diff - z * se
  upper <- diff + z * se
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

#' CI Diff Means (Unknown Variances, Big Samples)
#' Requirement: n1 + n2 > 30
ci_diff_means_unknown_large <- function(mean1, mean2, s1_sq, s2_sq, n1, n2, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  
  diff <- mean1 - mean2
  se <- sqrt((s1_sq / n1) + (s2_sq / n2))
  
  lower <- diff - z * se
  upper <- diff + z * se
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

#' CI Diff Means (Unknown Equal Variances, Small Samples)
#' Requirement: n1 + n2 <= 30 AND variances assumed equal
ci_diff_means_unknown_equal <- function(mean1, mean2, s1_sq, s2_sq, n1, n2, conf = 0.95) {
  alpha <- 1 - conf
  
  # Calculate Pooled Variance (Sp^2)
  sp_sq <- ((n1 - 1) * s1_sq + (n2 - 1) * s2_sq) / (n1 + n2 - 2)
  sp <- sqrt(sp_sq)
  
  # T-statistic with n1 + n2 - 2 degrees of freedom
  t_val <- qt(1 - alpha/2, df = n1 + n2 - 2)
  
  diff <- mean1 - mean2
  margin <- t_val * sp * sqrt((1/n1) + (1/n2))
  
  lower <- diff - margin
  upper <- diff + margin
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

#' CI Diff Means (Unknown DIFFERENT Variances, Small Samples)
ci_diff_means_unknown_unequal <- function(mean1, mean2, s1_sq, s2_sq, n1, n2, conf = 0.95) {
  alpha <- 1 - conf
  
  # Calculate degrees of freedom 'f' (Welch formula from PDF page 21)
  num <- ( (s1_sq/n1) + (s2_sq/n2) )^2
  denom <- ( ((s1_sq/n1)^2)/(n1+1) ) + ( ((s2_sq/n2)^2)/(n2+1) )
  f <- (num / denom) - 2
  
  # T-statistic with f degrees of freedom
  t_val <- qt(1 - alpha/2, df = f)
  
  diff <- mean1 - mean2
  se <- sqrt((s1_sq / n1) + (s2_sq / n2))
  
  lower <- diff - t_val * se
  upper <- diff + t_val * se
  
  cat(sprintf("Degrees of Freedom (f): %.2f\n", f))
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

# ------------------------------------------------------------------------------
# RATIO OF VARIANCES (sigma1/sigma2) & DIFF PROPORTIONS (p1-p2)
# ------------------------------------------------------------------------------

#' CI Ratio of Variances (Sigma1^2 / Sigma2^2) 
#' Notes: 2 normal distributions
ci_ratio_variances <- function(s1_sq, s2_sq, n1, n2, conf = 0.95) {
  alpha <- 1 - conf
  
  ratio <- s1_sq / s2_sq
  
  # F Critical Values
  f_lower_crit <- qf(alpha/2, df1 = n1 - 1, df2 = n2 - 1)
  f_upper_crit <- qf(1 - alpha/2, df1 = n1 - 1, df2 = n2 - 1)
  
  # PDF Formula: [Ratio / F_upper, Ratio / F_lower]
  lower <- ratio / f_upper_crit
  upper <- ratio / f_lower_crit
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

#' CI Difference of Proportions (p1 - p2) 
#' Notes: 2 binomial distributions
#' Requirements: big samples n1 + n2 > 30 && n1 ≈ n2.
ci_diff_proportions <- function(p1, p2, n1, n2, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  
  diff <- p1 - p2
  
  # Standard Error
  se <- sqrt( (p1*(1-p1)/n1) + (p2*(1-p2)/n2) )
  
  lower <- diff - z * se
  upper <- diff + z * se
  
  cat(sprintf("Interval: [%.4f, %.4f]\n", lower, upper))
  return(c(lower, upper))
}

# ==============================================================================
# PART 2: HYPOTHESIS TESTING (PARAMETRIC)
# Based on PDF Pages 49 - 57
# ==============================================================================

# HELPER FUNCTION: DECISION MAKER
# calculates decision based on p-value and prints result
print_decision <- function(p_value, alpha) {
  cat(sprintf("   P-Value: %.5f\n", p_value))
  cat(sprintf("   Alpha:   %.2f\n", alpha))
  if(p_value <= alpha) {
    cat("   RESULT: REJECT Null Hypothesis (H0)\n")
    cat("   (There is significant evidence for the alternative)\n")
  } else {
    cat("   RESULT: ACCEPT (Fail to Reject) Null Hypothesis (H0)\n")
  }
}

#' HELPER FUNCTION: Universal P-Value Calculator
#' Calculates the p-value for any distribution based on the test type.
#' 
#' @param stat The calculated test statistic
#' @param dist The distribution name: "norm", "t", "chisq", or "f"
#' @param type The hypothesis type: "less", "greater", "two.sided"
#' @param df1 Degrees of freedom 1 (for t, chisq, f)
#' @param df2 Degrees of freedom 2 (only for f)
get_p_value <- function(stat, dist, type, df1=NULL, df2=NULL) {
  
  p_val <- 0
  
  # --- NORMAL DISTRIBUTION (Z) ---
  if (dist == "norm") {
    if (type == "less")      p_val <- pnorm(stat)
    if (type == "greater")   p_val <- 1 - pnorm(stat)
    if (type == "two.sided") p_val <- 2 * (1 - pnorm(abs(stat)))
  }
  
  # --- T-STUDENT DISTRIBUTION (t) ---
  else if (dist == "t") {
    if (type == "less")      p_val <- pt(stat, df1)
    if (type == "greater")   p_val <- 1 - pt(stat, df1)
    if (type == "two.sided") p_val <- 2 * (1 - pt(abs(stat), df1))
  }
  
  # --- CHI-SQUARE DISTRIBUTION (Variance Test) ---
  # Note: Chi-sq is asymmetric. Two-sided is 2 * min(tail).
  else if (dist == "chisq") {
    if (type == "less")      p_val <- pchisq(stat, df1)
    if (type == "greater")   p_val <- 1 - pchisq(stat, df1)
    if (type == "two.sided") {
      left <- pchisq(stat, df1)
      right <- 1 - pchisq(stat, df1)
      p_val <- 2 * min(left, right)
    }
  }
  
  # --- F-DISTRIBUTION (Ratio Variances) ---
  else if (dist == "f") {
    if (type == "less")      p_val <- pf(stat, df1, df2)
    if (type == "greater")   p_val <- 1 - pf(stat, df1, df2)
    if (type == "two.sided") {
      left <- pf(stat, df1, df2)
      right <- 1 - pf(stat, df1, df2)
      p_val <- 2 * min(left, right)
    }
  }
  
  return(p_val)
}

# ------------------------------------------------------------------------------
# 1. TEST FOR MEAN (Mu0) OF A NORMAL POPULATION
# ------------------------------------------------------------------------------

#' Hypothesis Test for Mean (One Population) - CRITICAL REGION METHOD
#' 
#' @param x_bar Sample mean
#' @param mu0 Hypothesized population mean (value in H0)
#' @param n Sample size
#' @param sigma Known population STD (use NULL if unknown)
#' @param s Sample STD (use NULL if sigma is known)
#' @param type "two.sided" (!=), "less" (<), or "greater" (>) (comparator of the H1)
#' @param alpha Significance level (default 0.05)
test_mean <- function(x_bar, mu0, n, sigma=NULL, s=NULL, type="two.sided", alpha=0.05) {
  
  dist_type <- ""
  df <- 0
  
  if (!is.null(sigma)) {
    # Case 1: Known Variance (Z-Test)
    se <- sigma / sqrt(n)
    stat <- (x_bar - mu0) / se
    dist_type <- "norm"
    
  } else if (!is.null(s) && n > 30) {
    # Case 2: Unknown Variance, Large Sample n > 30 (Z-Test)
    se <- s / sqrt(n)
    stat <- (x_bar - mu0) / se
    dist_type <- "norm"
    
  } else if (!is.null(s) && n <= 30) {
    # Case 3: Unknown Variance, Small Sample n <= 30 (T-Test)
    se <- s / sqrt(n)
    stat <- (x_bar - mu0) / se
    dist_type <- "t"
    df <- n - 1
  } else {
    stop("Error: Provide either sigma (known) or s (unknown).")
  }
  
  # 2. Get P-Value using our new function
  p_val <- get_p_value(stat, dist_type, type, df1=df)
  
  # --- STEP 2: DETERMINE REJECTION REGION & DECISION ---
  
  reject_null <- FALSE
  crit_val <- 0
  crit_lower <- 0
  crit_upper <- 0
  
  if (dist_type == "norm") {
    # --- Z-DISTRIBUTION ---
    if (type == "less") {
      # H1: mu < mu0. Reject if stat < -z_alpha
      crit_val <- qnorm(alpha) # -z_alpha
      if (stat < crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region: (-Inf, %.4f)\n", crit_val))
      
    } else if (type == "greater") {
      # H1: mu > mu0. Reject if stat > z_alpha
      crit_val <- qnorm(1 - alpha) # z_alpha
      if (stat > crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
      
    } else { 
      # H1: mu != mu0. Reject if |stat| > z_alpha/2
      crit_lower <- qnorm(alpha / 2)     # -z_alpha/2
      crit_upper <- qnorm(1 - alpha / 2) # +z_alpha/2
      if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
      cat(sprintf("Critical Region: (-Inf, %.4f) U (%.4f, Inf)\n", crit_lower, crit_upper))
    }
    
  } else {
    # --- T-DISTRIBUTION ---
    if (type == "less") {
      # H1: mu < mu0. Reject if stat < -t_alpha
      crit_val <- qt(alpha, df) 
      if (stat < crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region (df=%d): (-Inf, %.4f)\n", df, crit_val))
      
    } else if (type == "greater") {
      # H1: mu > mu0. Reject if stat > t_alpha
      crit_val <- qt(1 - alpha, df)
      if (stat > crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region (df=%d): (%.4f, Inf)\n", df, crit_val))
      
    } else { 
      # H1: mu != mu0. Reject if |stat| > t_alpha/2
      crit_lower <- qt(alpha / 2, df)
      crit_upper <- qt(1 - alpha / 2, df)
      if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
      cat(sprintf("Critical Region (df=%d): (-Inf, %.4f) U (%.4f, Inf)\n", df, crit_lower, crit_upper))
    }
  }
  
  if(reject_null) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ------------------------------------------------------------------------------
# 2. TEST FOR THE VARIANCE (Sigma^2_0 -> "sigma0_sq") of a Normal population N(µ, sigma(std))
# ------------------------------------------------------------------------------

#' Hypothesis Test for Variance
#' 
#' @param s2 Sample Quasivariance (s^2)
#' @param sigma0_sq Hypothesized variance (value in H0)
#' @param n Sample size
#' @param type "two.sided" (!=), "less" (<), or "greater" (>)
#' @param alpha Significance level
test_variance <- function(s2, sigma0_sq, n, type="two.sided", alpha=0.05) {
  
  # --- STEP 1: CALCULATE TEST STATISTIC ---
  # [cite_start]Statistic Formula: (n-1)s^2 / sigma0^2 [cite: 692]
  stat <- ((n - 1) * s2) / sigma0_sq
  df <- n - 1
  
  ## 2. Get P-Value
  p_val <- get_p_value(stat, "chisq", type, df1=df)
  cat(sprintf("P-Value:        %.5f\n", p_val))
  # --- STEP 2: CHECK REJECTION REGION ---
  
  reject_null <- FALSE
  
  if (type == "less") {
    # H1: sigma^2 < sigma0^2. Reject if stat < Chi^2_(1-alpha)
    # Note: PDF table uses 1-alpha for the left tail critical value notation
    crit_val <- qchisq(alpha, df) 
    if (stat < crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: [0, %.4f)\n", crit_val))
    
  } else if (type == "greater") {
    # H1: sigma^2 > sigma0^2. Reject if stat > Chi^2_alpha
    crit_val <- qchisq(1 - alpha, df)
    if (stat > crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
    
  } else { 
    # H1: !=. Reject if stat not in [Chi^2_(1-a/2), Chi^2_(a/2)]
    crit_lower <- qchisq(alpha / 2, df)
    crit_upper <- qchisq(1 - alpha / 2, df)
    if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
    cat(sprintf("Critical Region: [0, %.4f) U (%.4f, Inf)\n", crit_lower, crit_upper))
  }
  
  if(reject_null) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ------------------------------------------------------------------------------
# 3. TEST FOR PROPORTION (p) OF A BINOMIAL DISTTRIBUTION B(n,p)
# ------------------------------------------------------------------------------

#' Hypothesis Test for Proportion (Binomial Approx)
#' 
#' @param p_hat Sample proportion
#' @param p0 Hypothesized proportion (value in H0)
#' @param n Sample size
#' @param type "two.sided" (!=), "less" (<), or "greater" (>)
#' @param alpha Significance level
test_proportion <- function(p_hat, p0, n, type="two.sided", alpha=0.05) {
  
  se <- sqrt( (p0 * (1 - p0)) / n )
  stat <- (p_hat - p0) / se
  
  # 2. Get P-Value
  p_val <- get_p_value(stat, "norm", type)
  cat(sprintf("P-Value:        %.5f\n", p_val))
  # --- STEP 2: CHECK REJECTION REGION ---
  
  reject_null <- FALSE
  
  if (type == "less") {
    # H1: p < p0. Reject if stat < -z_alpha
    crit_val <- qnorm(alpha)
    if (stat < crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: (-Inf, %.4f)\n", crit_val))
    
  } else if (type == "greater") {
    # H1: p > p0. Reject if stat > z_alpha
    crit_val <- qnorm(1 - alpha)
    if (stat > crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
    
  } else { 
    # H1: p != p0. Reject if |stat| > z_alpha/2
    crit_lower <- qnorm(alpha / 2)
    crit_upper <- qnorm(1 - alpha / 2)
    if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
    cat(sprintf("Critical Region: (-Inf, %.4f) U (%.4f, Inf)\n", crit_lower, crit_upper))
  }
  
  if(reject_null) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ------------------------------------------------------------------------------
# 4. TEST FOR EQUALITY OF TWO MEANS (Mu1 = Mu2) OF 2 NORMAL POPULATIONS
# ------------------------------------------------------------------------------

#' Hypothesis Test for Difference of Means
#' 
#' @param m1 Sample Mean 1
#' @param m2 Sample Mean 2
#' @param n1 Sample Size 1
#' @param n2 Sample Size 2
#' @param v1 Known Var1 (NULL if unknown)
#' @param v2 Known Var2 (NULL if unknown)
#' @param s1 Sample Var1 (s^2) (NULL if known)
#' @param s2 Sample Var2 (s^2) (NULL if known)
#' @param var.equal TRUE if we assume variances are equal (for small samples)
test_diff_means <- function(m1, m2, n1, n2, v1=NULL, v2=NULL, s1=NULL, s2=NULL, 
                            var.equal=FALSE, type="two.sided", alpha=0.05) {
  
  stat <- 0
  dist_type <- ""
  df <- 0
  
  
  if (!is.null(v1) && !is.null(v2)) {
    # [cite_start]Case 1: Known Variances (Z-Test) [cite: 722]
    se <- sqrt( (v1/n1) + (v2/n2) )
    stat <- (m1 - m2) / se
    dist_type <- "norm"
    
  } else if ((n1 + n2 > 30)) {
    # [cite_start]Case 2: Unknown Variances, Large Samples (Z-Test) [cite: 736]
    se <- sqrt( (s1/n1) + (s2/n2) )
    stat <- (m1 - m2) / se
    dist_type <- "norm"
    
  } else if (var.equal == TRUE) {
    # [cite_start]Case 3: Small Samples, EQUAL Variances (Pooled T-Test) [cite: 755]
    sp_sq <- ((n1 - 1)*s1 + (n2 - 1)*s2) / (n1 + n2 - 2)
    sp <- sqrt(sp_sq)
    se <- sp * sqrt( (1/n1) + (1/n2) )
    stat <- (m1 - m2) / se
    dist_type <- "t"
    df <- n1 + n2 - 2
    
  } else {
    # [cite_start]Case 4: Small Samples, UNEQUAL Variances (Welch T-Test) [cite: 775]
    se <- sqrt( (s1/n1) + (s2/n2) )
    stat <- (m1 - m2) / se
    dist_type <- "t"
    # [cite_start]Welch Degrees of Freedom [cite: 776]
    num <- ( (s1/n1) + (s2/n2) )^2
    den <- ( ((s1/n1)^2)/(n1+1) ) + ( ((s2/n2)^2)/(n2+1) )
    df <- (num/den) - 2
    cat(sprintf("Welch Degrees of Freedom: %.2f\n", df))
  }
  
  # 2. Get P-Value
  p_val <- get_p_value(stat, dist_type, type, df1=df)
  cat(sprintf("P-Value:        %.5f\n", p_val))
  # --- STEP 2: CHECK REJECTION REGION ---
  
  reject_null <- FALSE
  crit_val <- 0
  
  if (dist_type == "norm") {
    # Z-Test Logic
    if (type == "less") {
      crit_val <- qnorm(alpha)
      if (stat < crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region: (-Inf, %.4f)\n", crit_val))
    } else if (type == "greater") {
      crit_val <- qnorm(1 - alpha)
      if (stat > crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
    } else {
      crit_lower <- qnorm(alpha / 2)
      crit_upper <- qnorm(1 - alpha / 2)
      if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
      cat(sprintf("Critical Region: (-Inf, %.4f) U (%.4f, Inf)\n", crit_lower, crit_upper))
    }
    
  } else {
    # T-Test Logic
    if (type == "less") {
      crit_val <- qt(alpha, df)
      if (stat < crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region (df=%.2f): (-Inf, %.4f)\n", df, crit_val))
    } else if (type == "greater") {
      crit_val <- qt(1 - alpha, df)
      if (stat > crit_val) reject_null <- TRUE
      cat(sprintf("Critical Region (df=%.2f): (%.4f, Inf)\n", df, crit_val))
    } else {
      crit_lower <- qt(alpha / 2, df)
      crit_upper <- qt(1 - alpha / 2, df)
      if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
      cat(sprintf("Critical Region (df=%.2f): (-Inf, %.4f) U (%.4f, Inf)\n", df, crit_lower, crit_upper))
    }
  }
  
  if(reject_null) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ------------------------------------------------------------------------------
# 5. TEST FOR EQUALITY OF VARIANCES (Sigma1^2 / Sigma2^2) OF NORMAL DISTRIBUTIONS
# ------------------------------------------------------------------------------

#' Hypothesis Test for Equality of Variances
#' 
#' @param s1_sq Sample Variance 1
#' @param s2_sq Sample Variance 2
test_ratio_variances <- function(s1_sq, s2_sq, n1, n2, type="two.sided", alpha=0.05) {
  
  
  stat <- s1_sq / s2_sq
  df1 <- n1 - 1
  df2 <- n2 - 1
  
  cat(sprintf("Test Statistic (F): %.4f\n", stat))
  
  # 2. Get P-Value
  p_val <- get_p_value(stat, "f", type, df1=df1, df2=df2)
  cat(sprintf("P-Value:        %.5f\n", p_val))
  # --- STEP 2: CHECK REJECTION REGION ---
  
  reject_null <- FALSE
  
  if (type == "less") {
    # H1: sigma1^2 < sigma2^2. Reject if F < F_(1-alpha)
    crit_val <- qf(alpha, df1, df2)
    if (stat < crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: [0, %.4f)\n", crit_val))
    
  } else if (type == "greater") {
    # H1: sigma1^2 > sigma2^2. Reject if F > F_alpha
    crit_val <- qf(1 - alpha, df1, df2)
    if (stat > crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
    
  } else { 
    # H1: !=. Reject if F not in [F_(1-a/2), F_(a/2)]
    crit_lower <- qf(alpha / 2, df1, df2)
    crit_upper <- qf(1 - alpha / 2, df1, df2)
    if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
    cat(sprintf("Critical Region: [0, %.4f) U (%.4f, Inf)\n", crit_lower, crit_upper))
  }
  
  if(reject_null) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ------------------------------------------------------------------------------
# 6. TEST FOR EQUALITY OF PROPORTIONS (p1 = p2)
# ------------------------------------------------------------------------------

#' Hypothesis Test for Difference of Proportions
#' 
#' @param p1 Sample Proportion 1
#' @param p2 Sample Proportion 2
test_diff_proportions <- function(p1, p2, n1, n2, type="two.sided", alpha=0.05) {
  
  # --- STEP 1: CALCULATE STATISTIC ---
  num <- p1 - p2
  den <- sqrt( (p1*(1-p1)/n1) + (p2*(1-p2)/n2) )
  stat <- num / den
  
  # 2. Get P-Value
  p_val <- get_p_value(stat, "norm", type)
  cat(sprintf("P-Value:        %.5f\n", p_val))
  # --- STEP 2: CHECK REJECTION REGION ---
  
  reject_null <- FALSE
  crit_val <- 0
  
  if (type == "less") {
    crit_val <- qnorm(alpha)
    if (stat < crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: (-Inf, %.4f)\n", crit_val))
    
  } else if (type == "greater") {
    crit_val <- qnorm(1 - alpha)
    if (stat > crit_val) reject_null <- TRUE
    cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
    
  } else { 
    crit_lower <- qnorm(alpha / 2)
    crit_upper <- qnorm(1 - alpha / 2)
    if (stat < crit_lower || stat > crit_upper) reject_null <- TRUE
    cat(sprintf("Critical Region: (-Inf, %.4f) U (%.4f, Inf)\n", crit_lower, crit_upper))
  }
  
  if(reject_null) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ==============================================================================
# PART 3: CHI-SQUARE TESTS (NON-PARAMETRIC)
# Based on PDF Pages 67 - 73
# ==============================================================================

# TO FIT ANY DISTRIBUTION (PLEASE READ)
#' Helper: Fit Distribution & Get Probabilities
#' Calculates expected probabilities for Poisson, Binomial, or Normal distributions.
#' 
#' @param dist Type of distribution: "poisson", "binomial", or "normal"
#' @param observed Vector of observed frequencies
#' @param x_vals (For Poisson/Binomial) Vector of discrete values (e.g., 0, 1, 2...)
#' @param breaks (For Normal) Vector of interval boundaries (e.g., c(0, 10, 20...))
#' @param n_trials (For Binomial) Total number of trials N
#' @return A vector of probabilities summing to 1
fit_probs <- function(dist, observed, x_vals=NULL, breaks=NULL, n_trials=NULL) {
  
  n <- sum(observed)
  probs <- c()
  m <- 0 # Number of estimated parameters
  
  # --- CASE 1: POISSON DISTRIBUTION ---
  # Estimates Lambda (mean)
  if (dist == "poisson") {
    if (is.null(x_vals)) stop("Error: 'poisson' requires 'x_vals'.")
    
    # Estimate Lambda [cite: 900]
    lambda <- sum(x_vals * observed) / n
    m <- 1
    cat(sprintf("Fitted Poisson: Lambda = %.4f (m=1)\n", lambda))
    
    # Calculate probs for all specific values except the last "tail"
    exact_vals <- x_vals[1:(length(x_vals) - 1)]
    probs <- dpois(exact_vals, lambda)
    # The last bin gets the remainder (1 - sum)
    probs <- c(probs, 1 - sum(probs))
  }
  
  # --- CASE 2: BINOMIAL DISTRIBUTION ---
  # Estimates p (proportion)
  else if (dist == "binomial") {
    if (is.null(x_vals) || is.null(n_trials)) stop("Error: 'binomial' requires 'x_vals' and 'n_trials'.")
    
    # Estimate p: Mean = n*p -> p = Mean / n [cite: 901]
    mean_val <- sum(x_vals * observed) / n
    p_hat <- mean_val / n_trials
    m <- 1
    cat(sprintf("Fitted Binomial: p = %.4f (m=1)\n", p_hat))
    
    # Calculate probs
    probs <- dbinom(x_vals, size = n_trials, prob = p_hat)
    
    # Normalize if necessary (if x_vals covers full range 0..n, sum is 1. If grouped, handle tail)
    if(abs(sum(probs) - 1) > 0.001) {
      probs[length(probs)] <- 1 - sum(probs[-length(probs)])
    }
  }
  
  # --- CASE 3: NORMAL DISTRIBUTION ---
  # Estimates Mu (mean) and S (standard deviation) from grouped data
  else if (dist == "normal") {
    if (is.null(breaks)) stop("Error: 'normal' requires 'breaks' (interval boundaries).")
    
    # Calculate midpoints (class marks)
    mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
    
    # Estimate Mean (Mu) [cite: 902]
    mu <- sum(mids * observed) / n
    
    # Estimate Quasi-Standard Deviation (S) [cite: 902]
    var_s <- sum((mids - mu)^2 * observed) / (n - 1)
    s <- sqrt(var_s)
    m <- 2
    cat(sprintf("Fitted Normal: Mu = %.4f, S = %.4f (m=2)\n", mu, s))
    
    # Calculate Z-scores for boundaries
    # To ensure sum=1, we treat the first lower bound as -Infinity and last upper as +Infinity
    cum_probs <- pnorm(breaks, mean = mu, sd = s)
    cum_probs[1] <- 0 
    cum_probs[length(cum_probs)] <- 1
    
    # Prob of interval = P(upper) - P(lower)
    probs <- diff(cum_probs)
  }
  
  cat("NOTE: Use the 'm' value printed above in the test_goodness_of_fit function.\n")
  return(probs)
}

# ------------------------------------------------------------------------------
# 1. GOODNESS OF FIT
# ------------------------------------------------------------------------------

#' Chi-Square Goodness of Fit Test
#' 
#' @param observed Vector of observed frequencies (O_i)
#' @param probs Vector of expected probabilities (p_i)
#' @param m Number of estimated parameters (e.g. if you estimated mean, m=1) 
test_goodness_of_fit <- function(observed, probs, m=0, alpha=0.05) {
  
  n <- sum(observed)
  expected <- n * probs # (e_i)
  
  # Check condition: Expected frequencies >= 5 [cite: 950]
  if (any(expected < 5)) {
    warning("WARNING: Some expected frequencies are < 5. Results may be inaccurate.")
  }
  
  # Statistic Formula: Sum( (O - E)^2 / E )
  chi_stat <- sum( (observed - expected)^2 / expected )
  
  # Degrees of freedom: k - m - 1 
  k <- length(observed)
  df <- k - m - 1
  
  cat(sprintf("Test Statistic (Chi-Sq): %.4f\n", chi_stat))
  cat(sprintf("Degrees of Freedom: %d\n", df))
  
  # --- STEP 2: CHECK REJECTION REGION ---
  crit_val <- qchisq(1 - alpha, df)
  
  cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
  
  if (chi_stat >= crit_val) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}

# ------------------------------------------------------------------------------
# 2. TEST OF INDEPENDENCE (Contingency Tables)
# ------------------------------------------------------------------------------

#' Chi-Square Test of Independence
#' IF H0 ACCEPTED = INDEPENDENT
#' IF H0 REJECTED = DEPENDENT
#' @param matrix_obs A matrix of observed frequencies (O_ij)
test_independence <- function(matrix_obs, alpha=0.05) {
  
  row_sums <- rowSums(matrix_obs)
  col_sums <- colSums(matrix_obs)
  n <- sum(matrix_obs)
  
  # Calculate Expected Matrix: E_ij = (RowSum * ColSum) / n 
  expected_matrix <- outer(row_sums, col_sums) / n
  
  cat("Expected Frequencies Matrix:\n")
  print(round(expected_matrix, 2))
  
  # Statistic Formula: Sum( (O - E)^2 / E ) 
  chi_stat <- sum( (matrix_obs - expected_matrix)^2 / expected_matrix )
  
  # Degrees of freedom: (rows - 1) * (cols - 1) 
  df <- (nrow(matrix_obs) - 1) * (ncol(matrix_obs) - 1)
  
  cat(sprintf("Test Statistic (Chi-Sq): %.4f\n", chi_stat))
  cat(sprintf("Degrees of Freedom: %d\n", df))
  
  # --- STEP 2: CHECK REJECTION REGION ---
  # Reject H0 if chi^2 >= chi^2_alpha (upper tail) 
  crit_val <- qchisq(1 - alpha, df)
  
  cat(sprintf("Critical Region: (%.4f, Inf)\n", crit_val))
  
  if (chi_stat >= crit_val) {
    cat("RESULT: REJECT H0 (Statistic is inside the Rejection Region)\n")
  } else {
    cat("RESULT: ACCEPT H0 (Statistic is NOT in the Rejection Region)\n")
  }
}


# ==============================================================================
# PART 4: SAMPLE SIZE CALCULATION (Finding n)
# ==============================================================================

#' 1. Sample Size for Single Mean
#' Formula: n = (Z * sigma / E)^2
#' @param sigma Population standard deviation (or estimate s)
#' @param E Desired Margin of Error
#' @param conf Confidence Level
n_for_mean <- function(sigma, E, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  n <- (z * sigma / E)^2
  cat(sprintf("Required Sample Size: %d\n", ceiling(n)))
  return(ceiling(n))
}

#' 2. Sample Size for Single Proportion
#' Formula: n = p(1-p) * (Z/E)^2
#' @param E Desired Margin of Error (e.g., 0.03)
#' @param p Estimated proportion (default 0.5 for max sample size)
n_for_proportion <- function(E, p = 0.5, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  n <- p * (1 - p) * (z / E)^2
  cat(sprintf("Required Sample Size: %d\n", ceiling(n)))
  return(ceiling(n))
}

#' 3. Sample Size for Difference of Means (n1 = n2)
#' Formula: n = (Z/E)^2 * (sigma1^2 + sigma2^2)
#' @param s1 Standard Deviation Group 1
#' @param s2 Standard Deviation Group 2
#' @param E Desired Margin of Error for the difference
n_for_diff_means <- function(s1, s2, E, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  n <- (z / E)^2 * (s1^2 + s2^2)
  cat(sprintf("Required Sample Size (per group): %d\n", ceiling(n)))
  return(ceiling(n))
}

#' 4. Sample Size for Difference of Proportions (n1 = n2)
#' Formula: n = (Z/E)^2 * (p1(1-p1) + p2(1-p2))
#' @param p1 Estimated Proportion Group 1
#' @param p2 Estimated Proportion Group 2
#' @param E Desired Margin of Error for the difference
n_for_diff_proportions <- function(p1, p2, E, conf = 0.95) {
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  n <- (z / E)^2 * (p1 * (1 - p1) + p2 * (1 - p2))
  cat(sprintf("Required Sample Size (per group): %d\n", ceiling(n)))
  return(ceiling(n))
}

