plot_posteriors <- function(model, var) {
  
  # Error handling for the variable input:
  if(missing(var) || !is.character(var) || length(var) != 1) {
    stop("The requested variable does not exist in the rblimp model object.")
  }
  
  # Load ggplot2, stopping with an error if it's not installed.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but is not installed. Please install it with install.packages('ggplot2').")
  }
  library(ggplot2)
  
  # Extract the iterations data and convert it to a data frame.
  iterations_df <- as.data.frame(model@iterations)
  
  # Subset columns:
  # If var equals "All" (case-insensitive), use all columns;
  # otherwise, select only columns whose names start with var (case-insensitive).
  if(tolower(var) == "all") {
    cols_to_keep <- names(iterations_df)
  } else {
    cols_to_keep <- grep(paste0("^", var), names(iterations_df), ignore.case = TRUE, value = TRUE)
    if (length(cols_to_keep) == 0) {
      stop("The requested variable does not exist in the rblimp model object.")
    }
  }
  iterations_df_subset <- iterations_df[, cols_to_keep, drop = FALSE]
  
  # Clean the column names for printing: replace periods with spaces.
  cleaned_names <- gsub("\\.", " ", colnames(iterations_df_subset))
  colnames(iterations_df_subset) <- cleaned_names
  
  # Reshape the data using base R's stack() function.
  params_long <- stack(iterations_df_subset)
  names(params_long) <- c("value", "variable")
  
  # Set the factor levels for 'variable' so that facets appear in the original column order.
  params_long$variable <- factor(params_long$variable, levels = cleaned_names)
  
  # Assign fill colors based on the parameter name.
  params_long$fill_color <- ifelse(grepl("standardized", params_long$variable, ignore.case = TRUE),
                                   "#439A9D",
                                   ifelse(grepl("r2", params_long$variable, ignore.case = TRUE),
                                          "#583BBF",
                                          "#D95C14"))
  
  # Compute quantiles (2.5%, 50%, and 97.5%) and density values for each parameter.
  quantiles_list <- lapply(split(params_long$value, params_long$variable), function(x) {
    qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
    d <- density(x)
    lower_y <- as.numeric(approx(x = d$x, y = d$y, xout = qs[1], rule = 2)$y)
    upper_y <- as.numeric(approx(x = d$x, y = d$y, xout = qs[3], rule = 2)$y)
    data.frame(lower = qs[1], median = qs[2], upper = qs[3],
               lower_y = lower_y, upper_y = upper_y)
  })
  quantiles_df <- do.call(rbind, quantiles_list)
  
  # Ensure the factor levels for quantiles_df$variable match the original order.
  quantiles_df$variable <- factor(gsub("\\.", " ", rownames(quantiles_df)), levels = cleaned_names)
  
  # Create a label for each parameter using fixed-width formatting.
  quantiles_df$label <- sprintf("Est. = %-7.3f, CI = (%-7.3f, %-7.3f)",
                                quantiles_df$median, quantiles_df$lower, quantiles_df$upper)
  
  # Create new facet labels that combine the cleaned parameter name with the annotation.
  new_labels <- setNames(
    paste0(cleaned_names, "\n", quantiles_df$label[match(cleaned_names, as.character(quantiles_df$variable))]),
    cleaned_names
  )
  
  # Create the base ggplot: density plot with facets for each parameter.
  p <- ggplot(params_long, aes(x = value)) +
    geom_density(aes(fill = fill_color), alpha = 0.5) +
    scale_fill_identity() +
    facet_wrap(~ variable, scales = "free", labeller = as_labeller(new_labels)) +
    labs(title = paste("Posterior Distributions for", var),
         x = "Parameter Estimate", y = "Density") +
    theme_minimal() +
    theme(axis.title.y = element_blank(),  # Remove y-axis labels and ticks
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # Add a solid vertical line at the median.
  p <- p + geom_vline(data = quantiles_df, aes(xintercept = median),
                      color = "black", linetype = "solid")
  
  # Add solid vertical lines at the 2.5% and 97.5% quantiles,
  # with each line extending from y = 0 to the corresponding density value.
  p <- p + geom_segment(data = quantiles_df,
                        aes(x = lower, xend = lower, y = 0, yend = lower_y),
                        color = "black", linetype = "solid") +
    geom_segment(data = quantiles_df,
                 aes(x = upper, xend = upper, y = 0, yend = upper_y),
                 color = "black", linetype = "solid")
  
  return(p)
}

chibar_test <- function(model, raneff = NULL) {
  
  # error handling: if raneff is "Intercept", exit with an appropriate message
  if (!is.null(raneff) && any(tolower(raneff) == "intercept")) {
    return("This function currently supports tests of random slopes only.")
  }
  
  # get all column names from the model's iterations slot.
  col_names <- names(model@iterations)
  
  # count number of level-2 random effect parameters
  raneff_pattern <- "(?i)level-2 intercept covariance with|level-2 covariance between|level-2 slope variance|level-2 intercept variance"
  raneff_matches <- grep(raneff_pattern, col_names, perl = TRUE, value = TRUE)
  total_raneff <- length(raneff_matches)
  
  # count total number of random slopes
  slp_pattern <- "(?i)level-2 slope variance"
  slp_matches <- grep(slp_pattern, col_names, perl = TRUE, value = TRUE)
  total_slp <- length(slp_matches)
  
  # check whether input list has a random slope
  missing_vars <- raneff[!sapply(raneff, function(v) {
    any(grepl(v, col_names, ignore.case = TRUE) & 
          grepl("level-2 slope variance", col_names, ignore.case = TRUE))
  })]
  
  # identify column names to test
  if (!is.null(raneff)) {
    target_terms <- c("level-2 slope variance", "level-2 covariance between", "level-2 intercept covariance with")
    target_terms_pattern <- paste(target_terms, collapse = "|")
    raneff_pattern <- paste(raneff, collapse = "|")
    combined_pattern <- paste0("(?i)(", raneff_pattern, ").*(", target_terms_pattern, ")|(", 
                               target_terms_pattern, ").*(", raneff_pattern, ")")
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
    if(total_slp - length(raneff) == 1){
      if (wald == 0) {
        pvalue_chibar <- 1
      } else {
        pvalue_chibar <- (3/8) * pchisq(wald, df = 1, lower.tail = FALSE) +
          (3/8) * pchisq(wald, df = 2, lower.tail = FALSE) +
          (1/8) * pchisq(wald, df = 3, lower.tail = FALSE)
      }
    } 
    # test both slopes
    if(total_slp - length(raneff) == 0){
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
  line3 <- sprintf("%-40s %10.3f", "Probability (Chi-Bar Approximate)", pvalue_chibar_r)
  line4 <- sprintf("%-40s %10.3f", "Probability (Central)", pvalue_central_r)
  
  # combine the lines into a single multi-line string
  result <- paste(line1, line2, line3, line4, sep = "\n")
  
  # print the formatted table to the console
  if(length(missing_vars) == 0){
    cat(result, "\n")
  }
  
  # return table or error message
  if(total_raneff > 6) {
    return("This function currently supports models with only two random slopes.")
  } else if(length(missing_vars) >= 1) {
    return("Error: One or more variables on the test list do not have a random slope in the fitted model.")
  } else {
    return(invisible(result))
  }
  
}