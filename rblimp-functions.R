################################################################################
# posterior distribution plotting function (one parameter)
################################################################################

plot_one_posterior <- function(model, parname) {
  
  # model = med_model1
  # parname = 'indirect_ACE'
  
  # Check that parname is provided as a single character string
  if(missing(parname) || !is.character(parname) || length(parname) != 1) {
    stop("The requested parameter name must be a single character string.")
  }
  
  # Check that the model object is provided
  if(missing(model) || is.null(model)) {
    stop("Error: Model input object not found. Please supply a valid model object.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but is not installed. Please install it with install.packages('ggplot2').")
  }
  library(ggplot2)
  
  # Convert the iterations slot to a data frame
  iterations_df <- as.data.frame(model@iterations)
  
  # Search for columns with names that contain parname (ignoring case)
  matching_cols <- grep(parname, names(iterations_df), ignore.case = TRUE, value = TRUE)
  
  if (length(matching_cols) == 0) {
    stop("The requested parameter was not found in the model object.")
  }
  if (length(matching_cols) > 1) {
    warning("More than one column matches the parameter name. The first match will be used.")
  }
  
  # Use the first matching column
  selected_col <- matching_cols[1]
  
  # Extract the data vector
  data_vec <- iterations_df[[selected_col]]
  
  # Define fill color depending on parname content (same logic as your original)
  fill_color <- "#D95C14"
  
  # Compute quantiles and the corresponding y-values from the density estimate
  qs <- quantile(data_vec, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  d_obj <- density(data_vec, na.rm = TRUE)
  lower_y <- as.numeric(approx(x = d_obj$x, y = d_obj$y, xout = qs[1], rule = 2)$y)
  upper_y <- as.numeric(approx(x = d_obj$x, y = d_obj$y, xout = qs[3], rule = 2)$y)
  
  # Create a data frame for the density curve
  density_df <- data.frame(value = d_obj$x, density = d_obj$y)
  
  # Create a label for quantiles
  label_text <- sprintf("Est. = %-7.3f, CI = (%-7.3f, %-7.3f)", qs[2], qs[1], qs[3])
  
  # Set plot title
  plot_title <- paste("Posterior Distribution for", parname, "\n", label_text)
  
  # Build the density plot with ggplot2
  p <- ggplot(density_df, aes(x = value, y = density)) +
    # Add a ribbon to shade under the density curve
    geom_ribbon(aes(ymin = 0, ymax = density), fill = fill_color, alpha = 0.15) +
    # Then plot the density line on top
    geom_line(color = fill_color, size = 1.0) +
    labs(title = plot_title, x = "Parameter Estimate", y = "Density") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          plot.title   = element_text(size = 18, face = "bold"),
          plot.subtitle= element_text(size = 14))
  
  # Add vertical line at the median
  p <- p + geom_vline(xintercept = qs[2], color = "black", linetype = "solid")
  
  # Add segments for the lower and upper quantile density values
  ( p + 
      geom_segment(x = qs[1], xend = qs[1], y = 0, yend = lower_y,
                   color = "black", linetype = "solid") +
      geom_segment(x = qs[3], xend = qs[3], y = 0, yend = upper_y,
                   color = "black", linetype = "solid")
  )
  
}

################################################################################
# posterior distribution plotting function (multiple parameters)
################################################################################

plot_posteriors <- function(model, outcome) {
  
  if(missing(outcome) || !is.character(outcome) || length(outcome) != 1) {
    stop("The requested outcome does not exist in the rblimp model object.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but is not installed. Please install it with install.packages('ggplot2').")
  }
  library(ggplot2)
  
  iterations_df <- as.data.frame(model@iterations)
  
  if(tolower(outcome) == "all") {
    cols_to_keep <- names(iterations_df)
  } else {
    cols_to_keep <- grep(paste0("^", outcome), names(iterations_df), ignore.case = TRUE, value = TRUE)
    if (length(cols_to_keep) == 0) {
      stop("The requested outcome does not exist in the rblimp model object.")
    }
  }
  iterations_df_subset <- iterations_df[, cols_to_keep, drop = FALSE]
  
  # Replace periods with spaces and then collapse multiple spaces into one.
  cleaned_names <- gsub("\\.", " ", colnames(iterations_df_subset))
  cleaned_names <- gsub("\\s+", " ", cleaned_names)
  colnames(iterations_df_subset) <- cleaned_names
  
  params_long <- stack(iterations_df_subset)
  names(params_long) <- c("value", "variable")
  
  params_long$variable <- factor(params_long$variable, levels = cleaned_names)
  
  # Create a group variable based on variable names.
  params_long$group <- ifelse(grepl("standardized", params_long$variable, ignore.case = TRUE),
                              "Standardized",
                              ifelse(grepl("r2", params_long$variable, ignore.case = TRUE),
                                     "R2",
                                     ifelse(grepl("variance|covariance|correlation", params_long$variable, ignore.case = TRUE),
                                            "VarCovCor",
                                            "Other")))
  params_long$group <- factor(params_long$group, levels = c("Standardized", "R2", "VarCovCor", "Other"))
  
  quantiles_list <- lapply(split(params_long$value, params_long$variable), function(x) {
    qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
    d <- density(x)
    lower_y <- as.numeric(approx(x = d$x, y = d$y, xout = qs[1], rule = 2)$y)
    upper_y <- as.numeric(approx(x = d$x, y = d$y, xout = qs[3], rule = 2)$y)
    data.frame(lower = qs[1], median = qs[2], upper = qs[3],
               lower_y = lower_y, upper_y = upper_y)
  })
  quantiles_df <- do.call(rbind, quantiles_list)
  
  quantiles_df$variable <- factor(gsub("\\.", " ", rownames(quantiles_df)), levels = cleaned_names)
  levels(quantiles_df$variable) <- gsub("\\s+", " ", levels(quantiles_df$variable))
  
  # Revised label formatting using %.3f to avoid extra blank spaces.
  quantiles_df$label <- sprintf("Est. = %.3f, CI = (%.3f, %.3f)",
                                quantiles_df$median, quantiles_df$lower, quantiles_df$upper)
  
  new_labels <- setNames(
    paste0(cleaned_names, "\n", quantiles_df$label[match(cleaned_names, as.character(quantiles_df$variable))]),
    cleaned_names
  )
  
  if(tolower(outcome) == "all") {
    plot_title <- "Posterior Distributions for all Estimated Parameters"
  } else {
    plot_title <- paste("Posterior Distributions for the", outcome, "Model Parameters")
  }
  
  # Set the global mapping to x only.
  p <- ggplot(params_long, aes(x = value)) +
    geom_density(aes(fill = group, color = group), alpha = 0.25, size = 1.0) +
    facet_wrap(~ variable, scales = "free", labeller = as_labeller(new_labels)) +
    labs(title = plot_title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.position = "none") +
    scale_fill_hue() +
    scale_color_hue()
  
  # The quantiles layers now use their own mapping.
  p <- p + geom_vline(data = quantiles_df, mapping = aes(xintercept = median),
                      color = "black", linetype = "solid")
  
  p <- p + geom_segment(data = quantiles_df, mapping = aes(x = lower, xend = lower, y = 0, yend = lower_y),
                        color = "black", linetype = "solid") +
    geom_segment(data = quantiles_df, mapping = aes(x = upper, xend = upper, y = 0, yend = upper_y),
                 color = "black", linetype = "solid")
  
  return(p)
}

################################################################################
# function to compute wald-based chibar test of random slopes
################################################################################

chibar_test <- function(model, raneff = NULL, print = TRUE) {
  
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
    pvalue_binom <- if (wald == 0) {0.25 * 1 + 0.50 * 1 + 0.25 * 1 # If LRT == 0, the 0-df component gives probability 1.
    } else {0.25 * 0 + 0.50 * (1 - pchisq(wald, df = 1)) + 0.25 * (1 - pchisq(wald, df = 2))}
    pvalue_mixture <- 0.5 * (1 - pchisq(wald, df = 1)) + 0.5 * (1 - pchisq(wald, df = 2))
  }
  
  # pvalues for models with two random slopes
  if(total_raneff == 6){
    # test one slope
    if(total_slp - length(raneff) == 1){
      if (wald == 0) {pvalue_binom <- 1
      } else {
        pvalue_binom <- (3/8) * pchisq(wald, df = 1, lower.tail = FALSE) + (3/8) * pchisq(wald, df = 2, lower.tail = FALSE) + (1/8) * pchisq(wald, df = 3, lower.tail = FALSE)
      }
      pvalue_mixture <- 0.5 * (1 - pchisq(wald, df = 2)) + 0.5 * (1 - pchisq(wald, df = 3))
    } 
    # test both slopes
    if(total_slp - length(raneff) == 0){
      if (wald == 0) {
        pvalue_binom <- 1
      } else {
        pvalue_binom <- (5/32) * pchisq(wald, df = 1, lower.tail = FALSE) +
          (10/32) * pchisq(wald, df = 2, lower.tail = FALSE) +
          (10/32) * pchisq(wald, df = 3, lower.tail = FALSE) +
          (5/32) * pchisq(wald, df = 4, lower.tail = FALSE) +
          (1/32) * pchisq(wald, df = 5, lower.tail = FALSE)
      }
      pvalue_mixture <- NA
    } 
  }
  
  # round printed values
  wald_r <- round(wald, 3)
  pvalue_mixture_r <- round(pvalue_mixture, 3)
  pvalue_binom_r <- round(pvalue_binom, 3)
  pvalue_central_r <- round(pvalue_central, 3)
  
  # create formatted lines: left-align the description in a 40-character field,
  # then right-align the numeric value in a 10-character field.
  line1 <- sprintf("%-40s %10.3f", "Wald Statistic", wald_r)
  line2 <- sprintf("%-40s %10d", "Number of Parameters Tested (df)", df)
  line3 <- sprintf("%-40s %10.3f", "Probability (Chi-Bar Mixture Method)", pvalue_mixture_r)
  line4 <- sprintf("%-40s %10.3f", "Probability (Chi-Bar Binomial Method)", pvalue_binom_r)
  line5 <- sprintf("%-40s %10.3f", "Probability (Central Chi-Square)", pvalue_central_r)
  
  # combine the lines into a single multi-line string
  results2print <- paste(line1, line2, line3, line4, line5, sep = "\n")
  
  # print the formatted table to the console
  if(length(missing_vars) == 0 & print == TRUE){
    cat(results2print, "\n")
  }
  
  # return table or error message
  if(total_raneff > 6) {
    return("This function currently supports models with only two random slopes.")
  } else if(length(missing_vars) >= 1) {
    return("Error: One or more variables on the test list do not have a random slope in the fitted model.")
  } else {
    return(invisible(list(wald = wald,pmixture = pvalue_mixture,pbinom = pvalue_binom,pcentral = pvalue_central)))
  }
  
}

################################################################################
# simple slope plotting function
################################################################################

plot_interaction <- function(model, outcome, focal, moderator, bands = T) {
  
  # Load ggplot2, stopping with an error if it's not installed.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but is not installed. Please install it with install.packages('ggplot2').")
  }
  library(ggplot2)

  # Identify estimates to select
  iter_names <- names(model@iterations)
  selected_info <- data.frame(col_name = character(0),
                              type = character(0),
                              stringsAsFactors = FALSE)
  
  # Check whether the interaction term (focal*moderator or moderator*focal) is present in any column names
  interaction_found <- any(
    grepl(
      paste0("(", focal, "(\\.\\d+)?\\*", moderator, "(\\.\\d+)?|",
             moderator, "(\\.\\d+)?\\*", focal, "(\\.\\d+)?", ")"
      ),
      iter_names,
      ignore.case = TRUE
    )
  )
  if (!interaction_found) {
    stop("Error: The model does not include an interaction term between the focal and moderator variables.")
  }
  
  warning(paste("This function currently requires a model with only one interaction effect.",
                "Focal predictors must be numeric (manifest or latent). Moderators can be numeric (manifest or latent) or nominal.",
                sep = "\n"))
  
  # Loop over each column name.
  for (col in iter_names) {
    # Only consider columns that start with the outcome string.
    if (!startsWith(col, outcome)) next
    
    # Check for intercept columns.
    if (grepl("intercept", col, ignore.case = TRUE)) {
      selected_info <- rbind(selected_info, 
                             data.frame(col_name = col, type = "intercept", 
                                        stringsAsFactors = FALSE))
      next
    }
    
    # Check for "regressed on." columns.
    if (grepl("regressed on\\.", col, perl = TRUE)) {
      # Split on "regressed on." (using perl = TRUE for regex).
      parts <- strsplit(col, "regressed on\\.", perl = TRUE)[[1]]
      if (length(parts) != 2) next
      
      # Clean the left part: trim and remove any trailing period.
      left_part <- trimws(parts[1])
      left_part <- sub("\\.$", "", left_part)
      
      # Only proceed if the left part equals the outcome.
      if (left_part != outcome) next
      
      # Clean the right part.
      right_part <- trimws(parts[2])
      
      # Determine if right part contains the focal and/or moderator.
      contains_focal <- grepl(focal, right_part, fixed = TRUE)
      contains_moderator <- grepl(moderator, right_part, fixed = TRUE)
      
      if (contains_focal && !contains_moderator) {
        selected_info <- rbind(selected_info, 
                               data.frame(col_name = col, type = "focal",
                                          stringsAsFactors = FALSE))
      } else if (contains_moderator && !contains_focal) {
        selected_info <- rbind(selected_info, 
                               data.frame(col_name = col, type = "moderator",
                                          stringsAsFactors = FALSE))
      } else if (contains_focal && contains_moderator) {
        selected_info <- rbind(selected_info, 
                               data.frame(col_name = col, type = "both",
                                          stringsAsFactors = FALSE))
      }
    }
  }
  
  # Order the selected columns: intercept, focal, moderator, then both.
  if(nrow(selected_info) > 0) {
    type_order <- factor(selected_info$type, levels = c("intercept", "focal", "moderator", "both"))
    selected_info <- selected_info[order(type_order), ]
    selected_data <- model@iterations[, selected_info$col_name, drop = FALSE]
  } else {
    selected_data <- data.frame()  # Empty data frame if no columns match.
  }
  
  # Extract parameter draws for regression coefficients
  mat_p <- as.matrix(selected_data)
  
  # Define syntax as an object
  syntax_text <- as.character(model@syntax)
  syntax_lines <- unlist(strsplit(syntax_text, "\n"))
  
  # Check whether variables are nominal
  nominal_line <- grep("(?i)^NOMINAL:", syntax_lines, perl = TRUE, value = TRUE)
  if (length(nominal_line) > 0) {
    nominal_text <- sub("(?i)^NOMINAL:\\s*([^;]+);.*", "\\1", nominal_line[1], perl = TRUE)
  } else {
    nominal_text <- ""
  }
  # Create T/F flags for focal and moderator appearing in NOMINAL
  focal_categorical_flag <- grepl(focal, nominal_text, ignore.case = TRUE)
  moderator_categorical_flag <- grepl(moderator, nominal_text, ignore.case = TRUE)
  
  if (focal_categorical_flag) {
    stop("This function currently allows categorical moderators but requires numeric focal predictors.")
  }

  # Check whether either variable is centered
  center_line <- grep("(?i)^CENTER:", syntax_lines, perl = TRUE, value = TRUE)
  if (length(center_line) > 0) {
    centered_text <- sub("(?i)^CENTER:\\s*([^;]+);.*", "\\1", center_line[1], perl = TRUE)
  } else {
    centered_text <- ""
  }
  moderator_center_flag <- grepl(moderator, centered_text, fixed = TRUE)
  focal_center_flag     <- grepl(focal, centered_text, fixed = TRUE)
  
  if (moderator_categorical_flag & moderator_center_flag) {
    stop("This function currently requires categorical moderators that are not centered.")
  }
  
  # Check whether either variable is latent
  latent_line <- grep("(?i)^LATENT:", syntax_lines, perl = TRUE, value = TRUE)
  if (length(latent_line) > 0) {
    latent_text <- sub("(?i)^CENTER:\\s*([^;]+);.*", "\\1", latent_line[1], perl = TRUE)
  } else {
    latent_text <- ""
  }
  moderator_latent_flag <- grepl(moderator, latent_text, fixed = TRUE)
  focal_latent_flag     <- grepl(focal, latent_text, fixed = TRUE)
  
  # Alter names of any latent variables
  if (moderator_latent_flag) {
    moderator <- paste0(moderator, ".latent")
  }
  
  if (focal_latent_flag) {
    focal <- paste0(focal, ".latent")
  }
  
  if (moderator_categorical_flag) {
    model@average_imp[[moderator]] <- round(model@average_imp[[moderator]])
    num_moderator_cats <- length(unique(model@average_imp[[moderator]]))
    unique_moderator_vals <- sort(unique(model@average_imp[[moderator]]))
  }
  
  # Compute means and standard deviations
  mean_focal <- mean(model@average_imp[[focal]], na.rm = TRUE)
  sd_focal <- sd(model@average_imp[[focal]], na.rm = TRUE)
  mean_moderator <- mean(model@average_imp[[moderator]], na.rm = TRUE)
  sd_moderator <- sd(model@average_imp[[moderator]], na.rm = TRUE)
  
  # Define variable ranges
  if(focal_center_flag){
    x_low <- -sd_focal*1.5
    x_mean <- 0
    x_high <- +sd_focal*1.5
  } else{
    x_low <- mean_focal - sd_focal*1.5
    x_mean <- mean_focal
    x_high <- mean_focal + sd_focal*1.5
  }
  if (moderator_categorical_flag == F) {
    if(moderator_center_flag){
      m_low <- -sd_moderator
      m_mean <- 0
      m_high <- +sd_moderator
    } else{
      m_low <- mean_moderator - sd_moderator
      m_mean <- mean_moderator
      m_high <- mean_moderator + sd_moderator
    }
  }
  if (moderator_categorical_flag == T) {
    m_cat <- matrix(0, nrow = num_moderator_cats, ncol = num_moderator_cats - 1)
    for(g in 2:num_moderator_cats){
      m_cat[g,g-1] <- 1
    }
  }
  
  # Function to Compute conditional Effects
  cond_eff <- function(x, m) {
    cbind(b0 = x[,1] + x[,3] * m, b1 = x[,2] + x[,4] * m, m = m )
  }
  cond_eff_cat <- function(x, i) {
    m <- m_cat[i, ]
    m_score <- unique_moderator_vals[i]
    cbind(b0 = x[,1] + rowSums(x[,3:(3 + num_moderator_cats - 2)] * matrix(m, nrow = nrow(x), ncol = length(m), byrow = T)), 
          b1 = x[,2] + rowSums(x[,(3 + num_moderator_cats - 1):ncol(x)] * matrix(m, nrow = nrow(x), ncol = length(m), byrow = T)), m = m_score )
  }
  
  # Compute all conditional effects into data.frame
  if (moderator_categorical_flag == F) {
    simple_data <- as.data.frame(rbind(
      cond_eff(mat_p, m = m_low),
      cond_eff(mat_p, m = m_mean),
      cond_eff(mat_p, m =  m_high)
    ))
  }
  if (moderator_categorical_flag == T) {
    simple_data <- as.data.frame(
      do.call(rbind, lapply(seq_len(nrow(m_cat)), function(i) {
        cond_eff_cat(mat_p, i = i)
      }))
    )
  }
  # freq(simple_data$m)
  # Create factor based on levels of m
  simple_data$mf <- as.factor(simple_data$m)
  
  # Create range of predictor scores
  xvals <- seq(x_low, x_high, by = 0.1)
  pred_score <- \(d) d[1] + d[2] * xvals
  
  # Compute predicted scores
  pred <- lapply(split(simple_data[,1:2], simple_data$mf), \(x) apply(x, 1, pred_score))
  
  # Compute quantiles (2.5%, 50%, 97.5%)
  quan <- lapply(pred, \(x) apply(x, 1, quantile, p = c(0.025, 0.5, 0.975)))
  
  # Combine into data.frame
  rib_data  <- do.call('rbind', lapply(names(quan), \(x) {
    data.frame( l = quan[[x]][1,], fit = quan[[x]][2,], h = quan[[x]][3,], x = xvals, m = x)
  }))
  
  # Create factor with labels
  if (moderator_categorical_flag == F) {
    rib_data$mf <- factor(
      rib_data$m,
      levels = c(m_high, m_mean, m_low),
      labels = c('+1 SD','Mean','-1 SD')
    )
  }
  if (moderator_categorical_flag == T) {
    rib_data$mf <- factor(
      rib_data$m,
      levels = unique_moderator_vals,
      labels = paste0('Group ',unique_moderator_vals)
    )
  }
  
  ## Make Conditional Effects Plot
  if(bands){
    cond_plot <- (
      ggplot(rib_data, aes(x, color = mf, fill = mf)) +
        # geom_ribbon(aes(ymin = l, ymax = h), alpha = 0, color = NA) +
        geom_ribbon(aes(ymin = l, ymax = h), alpha = .15) +
        geom_line(aes(y = fit), linewidth = 1.5)
    )
  }
  if(!bands){
    cond_plot <- (
      ggplot(rib_data, aes(x, color = mf, fill = mf)) +
        geom_ribbon(aes(ymin = l, ymax = h), alpha = 0, color = NA) +
        # geom_ribbon(aes(ymin = l, ymax = h), alpha = .15) +
        geom_line(aes(y = fit), linewidth = 1.5)
    )
  }

  
  ## Print Plot with labels
(
    cond_plot
    + scale_x_continuous(
      paste(focal, "Scores"),
      limits = c(x_low, x_high)
    )
    + scale_y_continuous(
      paste(outcome, " Scores")
    )
    # + scale_color_manual(values = c("#D95C14","#439A9D","#583BBF","#AE151D","#2E62EA"))
    # + scale_fill_manual(values = c("#D95C14","#439A9D","#583BBF","#AE151D","#2E62EA"))
    + labs(
      color = "Moderator",
      fill = "Moderator",
      title = 'Plot of Conditional Regressions',
      subtitle = paste(moderator, "as Moderator")
    )
    + theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14)
    )
  )
  
}
