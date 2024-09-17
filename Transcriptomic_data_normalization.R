
## Transcriptomic data normalization

# Load transcriptomic data 
transcriptomic_data = readRDS('expr_sclc_ucologne_2015.rds')

# Load necessary library
library(dplyr)

# 1. Perform an appropriate log transformation on the data.

#constant to handle 0 values in matrix
pseudocount = 0.01
# Apply the log transformation to all columns (assuming all are numeric)
transformed_data = log2(transcriptomic_data + pseudocount)

#2. Implement a median polish algorithm from scratch.

median_polish = function(data, max_iter = 1000, tol = 1e-6) {
  
  data_matrix = as.matrix(data[, -1])  # Convert to matrix
  
  # Initialize row and column effects
  reff = rep(0, nrow(data_matrix))
  ceff = rep(0, ncol(data_matrix))
  
  for (iter in 1:max_iter) {
    # Store previous row effects for convergence check
    previous_reff = reff
    
    # Update row effects
    for (i in 1:nrow(data_matrix)) {
      reff[i] = median(data_matrix[i, ] - ceff)
    }
    
    # Update column effects
    for (j in 1:ncol(data_matrix)) {
      ceff[j] = median(data_matrix[, j] - reff)
    }
    
    # Check for convergence: if changes in row effects are small
    if (max(abs(reff - previous_reff)) < tol) {
      print("minimal changes detected")
      break
    }
  }
  
  # Calculate final residuals
  final_residuals = data_matrix - outer(reff, rep(1, ncol(data_matrix))) - outer(rep(1, nrow(data_matrix)), ceff)
  
  return(final_residuals)
}


# Run the custom median polish algorithm
custom_residuals = median_polish(transformed_data)

# Display the results
print(custom_residuals)


# 3. Compare the residuals of your algorithm and stats::medpolish

# Run median polish algo from stats::medpolish
library(stats)
medianpolish_result = medpolish(as.matrix(transformed_data[,-1]), trace = FALSE)
medianpolish_result_residuals = medianpolish_result$residuals

# compare
comparison = data.frame(custom_residuals = as.vector(custom_residuals), stats_medpolish_residuals = as.vector(medianpolish_result_residuals))
print(comparison)

# 4. Plot heatmaps of the results before and after median polish.

#install.packages("reshape2")
#install.packages("tidyverse")
library(tidyverse)

before_median_polish = as.matrix(transformed_data[,-1])
after_custom_median_polish = custom_residuals

before_median_polish_melted = reshape2::melt(before_median_polish)
after_custom_median_polish_melted = reshape2::melt(after_custom_median_polish)

ggplot(before_median_polish_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(before_median_polish_melted$value), limit = c(min(before_median_polish_melted$value), max(before_median_polish_melted$value)), name = "Expression") +
  labs(title = "Original Gene Expression Data", x = "Conditions", y = "Genes") +
  theme_minimal()

ggplot(after_custom_median_polish_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = median(after_custom_median_polish_melted$value), 
                       limit = c(min(after_custom_median_polish_melted$value), max(after_custom_median_polish_melted$value)), 
                       name = "Residuals") +
  labs(title = "Residuals After Median Polish", x = "Conditions", y = "Genes") +
  theme_minimal()


# 5.Output the median polished residual matrix as the normalized transcriptomic data.

normalized_data = as.data.frame(custom_residuals)
write.csv(normalized_data,"normalized_transcriptomic_data.csv", row.names = TRUE)
