## Mutation data analysis

# Install required packages, if not already installed
install.packages("dplyr")
install.packages("tidyr")

library(dplyr)
library(tidyr)

# Load the .rds file
data = readRDS("mutations_sclc_ucologne_2015.rds")

# Count mutations per gene
mutation_counts = data %>%
                  group_by(gene) %>%
                  summarise(count = n()) %>%
                  arrange(desc(count))


# Q1. Identify the top 10 most frequently mutated genes. Identify samples whose mutation count is in the 80 to 90 percentile.
# Get the top 10 most frequently mutated genes
top_genes = head(mutation_counts, 10)
print(top_genes)

# Count mutations per sample
sample_mutation_counts = data %>%
                         group_by(sample_id) %>%  
                         summarise(mutation_count = n())

# Calculate 80th and 90th percentiles
percentiles = quantile(sample_mutation_counts$mutation_count, probs = c(0.8, 0.9))

# Filter samples in the 80th to 90th percentile
samples_in_percentile = sample_mutation_counts %>%
  filter(mutation_count >= percentiles[1] & mutation_count <= percentiles[2])

print(samples_in_percentile)
# Output
# sample_id                 mutation_count

# 1 sclc_ucologne_2015_S00841            497
# 2 sclc_ucologne_2015_S01020            624
# 3 sclc_ucologne_2015_S01022            548
# 4 sclc_ucologne_2015_S01023            494
# 5 sclc_ucologne_2015_S01861            483
# 6 sclc_ucologne_2015_S02242            479
# 7 sclc_ucologne_2015_S02248            492
# 8 sclc_ucologne_2015_S02285            561
# 9 sclc_ucologne_2015_S02328            486
# 10 sclc_ucologne_2015_S02344           574
# 11 sclc_ucologne_2015_S02376           504
# 12 sclc_ucologne_2015_S02384           479

#2. Categorize variants based on their expected effects using the data/mutation_effects.tsv table. Generate a count matrix containing the numbers of loss-of-function and neutral mutations for each gene.

# Load the mutation effects data
mutation_effects = read.table("mutation_effects.tsv", header = TRUE, sep = "\t")

# Create a vector for categorization
mutation_effects = mutation_effects %>%
  mutate(effect_type = case_when(
    effect == "loss_of_function" ~ "loss_of_function",
    effect == "neutral" ~ "neutral",
    TRUE ~ "other"  # Optional, for any other categories
  ))

# Join the mutations with effects
mutations_with_effects = data %>%
  left_join(mutation_effects, by = "variant_class")

# Count the occurrences of each effect type per gene
count_matrix = mutations_with_effects %>%
  group_by(gene, effect_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = effect_type, values_from = count, values_fill = 0)

# View the count matrix
print(count_matrix)

#3.Implement a statistical test that determines whether a gene has a significantly higher proportion of loss-of-function mutations (excluding mutations with uncertain effects), compared to other genes.

# Filter for relevant effects
filtered_mutations = mutations_with_effects %>%
  filter(effect != "uncertain")  # Exclude uncertain effects

# group filtered data by genes and count the loss_of_function effect and non loss_of_function effects
gene_by_effect_table = filtered_mutations %>%
  group_by(gene) %>%
  summarise(
    loss_of_function = sum(effect == "loss_of_function"),
    non_loss_of_function = sum(effect != "loss_of_function"),
    .groups = 'drop'
  )

print(gene_by_effect_table)


# Total counts for other genes
total_loss_of_function = sum(gene_by_effect_table$loss_of_function)
total_non_loss_of_function = sum(gene_by_effect_table$non_loss_of_function)

# Initialize a function to perform Chi-squared test and calculate effect size
calculate_stats = function(row) {
  gene = row["gene"]
  
  # Extract counts for the current gene
  counts = c(
    as.numeric(row["loss_of_function"]),  # Loss-of-function for the gene
    total_loss_of_function - as.numeric(row["loss_of_function"]),  # Loss-of-function for others
    as.numeric(row["non_loss_of_function"]),  # Non-loss-of-function for the gene
    total_non_loss_of_function - as.numeric(row["non_loss_of_function"])  # Non-loss-of-function for others
  )
  
  # Create contingency matrix for Chi squared test
  contingency_matrix = matrix(counts, nrow = 2)
  
  # Perform Chi-squared tests for each gene
  chi_square_result = chisq.test(contingency_matrix)
  
  # Calculate effect size (odds ratio)
  odds_ratio = (counts[1] / counts[3]) / (counts[2] / counts[4])
  
  # Return results as a list
  return(c(gene_symbol = gene, 
           effect_size = odds_ratio, 
           p_value = chi_square_result$p.value))
}

# Apply the function to each row of the gene_by_effect_table
results_matrix = apply(gene_by_effect_table, 1, calculate_stats)

#4. Identify candidate tumour suppressor genes using this statistical test, adjusting for multiple hypothesis testing
# Convert the results into a data frame
results = as.data.frame(t(results_matrix), stringsAsFactors = FALSE)

# Convert numeric columns to appropriate types
results$effect_size = as.numeric(results$effect_size)
results$p_value = as.numeric(results$p_value)

# Adjust p-values for multiple testing using Benjamini-Hochberg
results$q_value = p.adjust(results$p_value, method = "BH")

# View results
print(results)

significance_threshold = 0.05
candidate_genes_with_significant_loss_of_function_mutation = results %>%
                                    filter(p_value < significance_threshold)
print(candidate_genes_with_significant_loss_of_function_mutation)
