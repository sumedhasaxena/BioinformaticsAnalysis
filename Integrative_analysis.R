# Integrative analysis

# 1.For each gene involved in a structural variant (SV), determine the expression level of the gene
# in the sample that harbours the SV
sv_data = read.table("sv_sclc_ucologne_2015.tsv", header = TRUE, sep = "\t")
transcriptomic_data = readRDS('expr_sclc_ucologne_2015.rds')

library(dplyr)

# Extract unique genes involved in SVs
unique_genes = unique(c(sv_data$site1_hugo_symbol, sv_data$site2_hugo_symbol))

# Initialize a data frame to store results
results = data.frame(sample_id = character(),
                     gene = character(),
                     expression_level = numeric(),
                     stringsAsFactors = FALSE)

for(sample in unique(sv_data$sample_id)){
  # Get the genes involved in SVs for the current sample
  sv_genes = sv_data %>%
    filter(sample_id == sample) %>%
    select(site1_hugo_symbol, site2_hugo_symbol) %>%
    unlist() %>%
    unique()
  
  # Retrieve expression levels for these genes in the current sample
  for(gene in sv_genes){
    if(gene %in% rownames(transcriptomic_data)){
      expression_level = transcriptomic_data[gene,sample]
      results = rbind(results, data.frame(sample_id = sample,
                                          gene = gene,
                                          expression_level = expression_level))
    }
  }
}

print(results)

# 2. Identify SVs that satisfy the following critiera:
# The involved pair of genes both have elevated expression levels in samples with the SV 
# compared to samples without the SV.

# Initialize a data frame to store results
results <- data.frame()

# Iterate through each unique SV
for (sample in unique(sv_data$sample_id)) {
  # Get the SV information for the current sample
  sv_data_sample = sv_data %>%
    filter(sample_id == sample)
  
  # Get the genes involved in the SV
  genes = unique(c(sv_data_sample$site1_hugo_symbol, sv_data_sample$site2_hugo_symbol))
  
  # upon testing, it was found that some genes in the 'genes' variable above are not present in original transcriptomic data.
  # hence, filtering the result to have only relevant genes
  
  present_genes = genes[genes %in% rownames(transcriptomic_data)]
  
  if (length(present_genes) > 0) {
    # Expression levels for the current sample
    expr_with_sv = transcriptomic_data[present_genes, sample]
    
    # Get samples without the SV
    samples_without_sv = setdiff(colnames(transcriptomic_data), sample)
    
    # Calculate average expression levels for samples without the SV
    avg_expr_without_sv = colMeans(transcriptomic_data[present_genes, samples_without_sv], na.rm = TRUE)
    
    # Check the criteria
    if (all(expr_with_sv > avg_expr_without_sv) &&
        any(sv_data_sample$site2_effect_on_frame == "in-frame")) {
      
      
        # Store the result for each gene pair
        for (i in seq_len(nrow(sv_data_sample))) {
          
          if (sv_data_sample$site1_hugo_symbol[i] %in% present_genes && 
              sv_data_sample$site2_hugo_symbol[i] %in% present_genes) {
            
              results = rbind(results, data.frame(
                sample_id = sample,
                gene1 = sv_data_sample$site1_hugo_symbol[i],
                gene2 = sv_data_sample$site2_hugo_symbol[i],
                expr_with_sv_gene1 = expr_with_sv[sv_data_sample$site1_hugo_symbol[i]],
                expr_with_sv_gene2 = expr_with_sv[sv_data_sample$site2_hugo_symbol[i]]))
          }
        }
      }
    }
  }


# View the resulting data frame with SVs meeting the criteria
print(results)

