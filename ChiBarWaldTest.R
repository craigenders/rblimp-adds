chibar_test_slopes <- function(model, varnames = NULL) {
  
  # get all column names from the model's iterations slot.
  col_names <- names(model@iterations)
  
  # count number of level-2 random effect parameters
  raneff_pattern <- "(?i)level-2 intercept covariance with|level-2 covariance between|level-2 slope variance|level-2 intercept variance"
  raneff_matches <- grep(raneff_pattern, col_names, perl = TRUE, value = TRUE)
  total_raneff <- length(raneff_matches)
  
  # count total number of randonm slopes
  slp_pattern <- "(?i)level-2 slope variance"
  slp_matches <- grep(slp_pattern, col_names, perl = TRUE, value = TRUE)
  total_slp <- length(slp_matches)
  
  # identify column names to test
  if (!is.null(varnames)) {
    target_terms <- c("level-2 slope variance", "level-2 covariance between", "level-2 intercept covariance with")
    target_terms_pattern <- paste(target_terms, collapse = "|")
    varnames_pattern <- paste(varnames, collapse = "|")
    combined_pattern <- paste0("(?i)(", varnames_pattern, ").*(", target_terms_pattern, ")|(",target_terms_pattern, ").*(", varnames_pattern, ")")
    combined_matches <- grep(combined_pattern, col_names, perl = TRUE, value = TRUE)
  } else {
    combined_matches <- NULL
  }
  
  # subset the iterations data
  if (!is.null(combined_matches) && length(combined_matches) > 0) {
    est2test <- model@iterations[, combined_matches, drop = FALSE]
  } else {
    est2test <- model@iterations[, raneff_matches, drop = FALSE]
  }
  
  # compute wald
  wald <- colMeans(est2test) %*% solve(cov(est2test)) %*% colMeans(est2test)
  df <- length(combined_matches)
  pvalue_central <- pchisq(wald, df = df, lower.tail = FALSE) 
  
  # pvalue for a model with a single random slope
  if(total_raneff == 3){
    if (wald == 0) {
      pvalue_chibar <- 1
    } else {
      pvalue_chibar <- 0.50 * pchisq(wald, df = 1, lower.tail = FALSE) +
        0.25 * pchisq(wald, df = 2, lower.tail = FALSE)
    }
  }
  
  # pvalues for models with two random slopes
  if(total_raneff == 6){
    # test one of the slopes
    if(total_slp - length(varnames) == 1){
      if (wald == 0) {
        pvalue_chibar <- 1
      } else {
        pvalue_chibar <- (3/8) * pchisq(wald, df = 1, lower.tail = FALSE) +
          (3/8) * pchisq(wald, df = 2, lower.tail = FALSE) +
          (1/8) * pchisq(wald, df = 3, lower.tail = FALSE)
      }
    } 
    # test both slopes
    if(total_slp - length(varnames) == 0){

      if (wald == 0) {
        pvalue_chibar <- 1
      } else {
        pvalue_chibar <- (5/32) * pchisq(wald, df = 1, lower.tail = FALSE) +
          (10/32) * pchisq(wald, df = 2, lower.tail = FALSE) +
          (10/32) * pchisq(wald, df = 3, lower.tail = FALSE) +
          (5/32) * pchisq(wald, df = 4, lower.tail = FALSE) +
          (1/32) * pchisq(wald, df = 5, lower.tail = FALSE)
      }
    } 
  }
  
  # round printed values
  wald_r <- round(wald, 3)
  pvalue_chibar_r <- round(pvalue_chibar, 3)
  pvalue_central_r <- round(pvalue_central, 3)
  
  # create formatted lines: left-align the description in a 40-character field,
  # then right-align the numeric value in a 10-character field.
  line1 <- sprintf("%-40s %10.3f", "Wald Statistic", wald_r)
  line2 <- sprintf("%-40s %10d", "Number of Parameters Tested (df)", df)
  line3 <- sprintf("%-40s %10.3f", "Probability (Chi-Bar)", pvalue_chibar_r)
  line4 <- sprintf("%-40s %10.3f", "Probability (Central)", pvalue_central_r)
  
  # combine the lines into a single multi-line string
  result <- paste(line1, line2, line3, line4, sep = "\n")
  
  # print the formatted table to the console
  cat(result, "\n")
  
  # return table
  if(total_raneff > 6) {
    return("This function currently supports models with only two random slopes.")
  } else {
    return(invisible(result))
  }
  
}
