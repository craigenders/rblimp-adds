plot_posteriors <- function(model, var) {
  
  # Error handling:
  # Check if the model has an 'iterations' slot.
  if(!("iterations" %in% slotNames(model))) {
    stop("The rblimp model object does not include the MCMC iterations slot")
  }
  # Check if 'var' is provided, is a character string, and is of length 1.
  if(missing(var) || !is.character(var) || length(var) != 1) {
    stop("The requested variable does not exist in the rblimp model object.")
  }
  
  # Load ggplot2, stopping with an error if it's not installed.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but is not installed. Please install it using install.packages('ggplot2').")
  }
  library(ggplot2)
  
  # Extract the iterations data and convert it to a data frame.
  iterations_df <- as.data.frame(model@iterations)
  
  # Subset columns whose names start with the provided variable (case-insensitive).
  cols_to_keep <- grep(paste0("^", var), names(iterations_df), ignore.case = TRUE, value = TRUE)
  if (length(cols_to_keep) == 0) {
    stop("The requested variable does not exist in the rblimp model object.")
  }
  iterations_df_subset <- iterations_df[, cols_to_keep, drop = FALSE]
  
  # Clean the column names for printing: replace periods with spaces.
  cleaned_names <- gsub("\\.", " ", colnames(iterations_df_subset))
  colnames(iterations_df_subset) <- cleaned_names
  
  # Reshape the data using base R's stack() function.
  params_long <- stack(iterations_df_subset)
  names(params_long) <- c("value", "variable")
  
  # Set the factor levels for 'variable' so that the facets appear in the original column order.
  params_long$variable <- factor(params_long$variable, levels = cleaned_names)
  
  # Assign fill colors based on the parameter name:
  # If the name contains "standardized", use "#439A9D";
  # if it contains "r2", use "#583BBF";
  # otherwise use "#D95C14".
  params_long$fill_color <- ifelse(grepl("standardized", params_long$variable, ignore.case = TRUE),
                                   "#439A9D",
                                   ifelse(grepl("r2", params_long$variable, ignore.case = TRUE),
                                          "#583BBF",
                                          "#D95C14"))
  
  # Compute quantiles (2.5%, 50% and 97.5%) and density values for each parameter.
  quantiles_list <- lapply(split(params_long$value, params_long$variable), function(x) {
    qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
    d <- density(x)
    lower_y <- as.numeric(approx(x = d$x, y = d$y, xout = qs[1], rule = 2)$y)
    upper_y <- as.numeric(approx(x = d$x, y = d$y, xout = qs[3], rule = 2)$y)
    data.frame(lower = qs[1], median = qs[2], upper = qs[3],
               lower_y = lower_y, upper_y = upper_y)
  })
  
  quantiles_df <- do.call(rbind, quantiles_list)
  # Ensure the factor levels match the original order.
  quantiles_df$variable <- factor(gsub("\\.", " ", rownames(quantiles_df)), levels = cleaned_names)
  
  # Create a label for each parameter using fixed-width formatting.
  quantiles_df$label <- sprintf("Est. = %-7.3f\nCI = (%-7.3f, %-7.3f)",
                                quantiles_df$median, quantiles_df$lower, quantiles_df$upper)
  
  # Create the base ggplot: density plot with facets for each parameter.
  p <- ggplot(params_long, aes(x = value)) +
    geom_density(aes(fill = fill_color), alpha = 0.5) +
    scale_fill_identity() +
    facet_wrap(~ variable, scales = "free") +
    labs(title = paste("Posterior Distributions for", var),
         x = "Parameter Estimate", y = "Density") +
    theme_minimal()
  
  # Add a solid vertical line at the median.
  p <- p + geom_vline(data = quantiles_df, aes(xintercept = median),
                      color = "black", linetype = "solid")
  
  # Add solid vertical lines at the 2.5% and 97.5% quantiles,
  # with each line extending from y=0 to the corresponding density value.
  p <- p + geom_segment(data = quantiles_df,
                        aes(x = lower, xend = lower, y = 0, yend = lower_y),
                        color = "black", linetype = "solid") +
    geom_segment(data = quantiles_df,
                 aes(x = upper, xend = upper, y = 0, yend = upper_y),
                 color = "black", linetype = "solid")
  
  # Add text annotation in the upper left corner of each facet.
  p <- p + geom_text(data = quantiles_df, 
                     aes(x = -Inf, y = Inf, label = label),
                     hjust = -0.1, vjust = 1.1, size = 3)
  
  return(p)
}