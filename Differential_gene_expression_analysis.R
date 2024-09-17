clinical_data = read.table("pheno_sclc_ucologne_2015.tsv", header = TRUE, sep = "\t")

library(dplyr)

# 1. Define two groups of tumours as early stage (stages I-II) vs. advanced stage tumours (stages III-IV),
# while excluding samples missing stage information.

early_stage = c('I', 'Ia', 'Ib', 'IB', 'II', 'IIa', 'IIb')
advanced_stage = c('III', 'IIIa', 'IIIb', 'IV', 'IVa', 'IVb')


updated_clinical_data = clinical_data %>%
                 filter(!is.na(uicc_tumor_stage) & uicc_tumor_stage != '') %>%
                mutate(tumour_group = case_when(
                  uicc_tumor_stage %in%  early_stage ~ 'early',
                  uicc_tumor_stage %in%  advanced_stage ~ 'advanced'
                ))

#Load original transcriptomic data
transcriptomic_data = readRDS('expr_sclc_ucologne_2015.rds')

#2. Identify genes that differentially expressed in early vs. advanced stage tumours using an appropriate R package.

# For analyzing differential gene expression data derived from microarray experiments, the limma package in R is an appropriate choice.
# It is designed for the analysis of gene expression data and can handle continuous data such as floating-point numbers

BiocManager::install("limma")
library(limma)

# Subset the clinical data to include only relevant samples
# Since the clinical data and transciptomic data have a mismatch of samples,
# finding a common group of samples present in both data sets to proceed

relevant_clinical_data = updated_clinical_data[updated_clinical_data$patient_id %in% colnames(transcriptomic_data), ]
relevant_gene_data = transcriptomic_data[, colnames(transcriptomic_data) %in% relevant_clinical_data$patient_id]

# creating a design matrix for the linear model, including the tumour group
design = model.matrix(~ tumour_group, data = relevant_clinical_data)

# fit the linear model
fit = lmFit(relevant_gene_data, design)

# apply empirical Bayes moderation
fit = eBayes(fit)

results = topTable(fit, coef=2, number=Inf, adjust="fdr")

significant_genes_with_p_value = results[results$P.Value < 0.05, ]
significant_genes_with_adjusted_p_value = results[results$adj.P.Val < 0.05, ]

##Comment: No significant genes found with adjusted P value.
##Comment: 460 significant genes found with original P value.

print(significant_genes_with_p_value)
print(significant_genes_with_adjusted_p_value)

#visualize results with a volcano plot
library(ggplot2)

volcano_data = as.data.frame(results)
volcano_data$significant = volcano_data$P.Value < 0.05

ggplot(volcano_data, aes(x=logFC, y=-log10(P.Value), color=significant)) +
  geom_point(alpha=0.5) +
  theme_minimal() +
  labs(title="Volcano Plot of Differentially Expressed Genes",
       x="Log Fold Change",
       y="-Log10 P-Value") +
  scale_color_manual(values=c("black", "red"))