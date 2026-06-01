
install.packages("imsig")
install.packages(c("tidyverse", "ggplot2", "reshape2"))
install.packages("ggpubr")


# Load libraries
library(imsig)
library(tidyverse)
library(reshape2)
library(ggpubr)

# Load data
counts <- read.csv("imsig_input_counts_matrix.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("imsig_input_metadata.csv")

counts_g = counts[!counts$Gene == '',]

counts_nog = subset(counts_g, select = -c(Gene))

# Duplicate gene names we take the mean
counts_agg = aggregate(counts_nog,by=list(Gene=counts_g$Gene),data=counts_g,FUN=mean)
row.names(counts_agg) = counts_agg$Gene
counts_agg = subset(counts_agg, select = -c(Gene))

# Estimate immune cell abundance using imsig
abundance_scores <- imsig(counts_agg)
head(abundance_scores)

# Convert results to long format for easier plotting
imsig_long <- abundance_scores %>%
  rownames_to_column(var = "Code") %>%
  pivot_longer(-Code, names_to = "CellType", values_to = "Abundance")

head(imsig_long)

# Merge with metadata
imsig_annotated <- imsig_long %>%
  left_join(metadata, by = "Code")

head(imsig_annotated)

# Save results
write.csv(imsig_annotated, "imsig_abundance_results.csv")

# Remove categories that are not immune cells for the results
imsig_annotated_m = imsig_annotated[!imsig_annotated$CellType %in% c('Proliferation', 'Translation'),]

# Plot immune cell abundance grouped by "Baseline biopsy"
ggplot(imsig_annotated_m, aes(x = CellType, y = Abundance, fill = `Baseline.biopsy`)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 14) +
  labs(title = "Immune Cell Abundance per Sample Type",
       x = "Immune Cell Type",
       y = "Abundance Score",
       fill = "Baseline.biopsy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Kruskal–Wallis test per cell type
kruskal_results <- imsig_annotated_m %>%
  group_by(CellType) %>%
  summarise(
    p_value = kruskal.test(Abundance ~ `Baseline.biopsy`)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Show significant results
kruskal_results %>%
  filter(p_adj < 0.05)
kruskal_results

# post-hoc pairwise Wilcoxon tests for significant cell types
posthoc_results <- list()

signif_cells <- kruskal_results %>% filter(p_adj < 0.05) %>% pull(CellType)
signif_cells

for (cell in signif_cells) {
  df <- imsig_annotated_m %>% filter(CellType == cell)
  pw <- pairwise.wilcox.test(df$Abundance, df$`Baseline.biopsy`, p.adjust.method = "BH")
  posthoc_results[[cell]] <- pw$p.value
}
posthoc_results

# Save results
write.csv(kruskal_results, "imsig_kruskal_results.csv")
write.csv(posthoc_results, "imsig_pairwise_posthoc_results.csv")
saveRDS(posthoc_results, "imsig_pairwise_posthoc_results.rds")


