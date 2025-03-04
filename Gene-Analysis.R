# Install libraries 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(version='devel')
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}
if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db")
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}
if (!requireNamespace("M3C", quietly = TRUE)) {
  BiocManager::install("M3C")
}
if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}

# Load the libraries
library(gprofiler2)
library(DESeq2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(clusterProfiler)
library(ggplot2)
library(Rtsne)
library(umap)
library(magrittr)
library(M3C)
library(dplyr)

data_file <- file.path("./data/SRP075377.tsv")
metadata_file <- file.path("./data/metadata_SRP075377.tsv")
plotspath <- file.path("./plots")
resultspath <- file.path("./Results")

# Step 1: Load Files 
expression_data <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

metadata <- readr::read_tsv(metadata_file)

# Make the data in the order of the metadata
expression_df <- expression_data  %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

mapped_list <- mapIds(
  org.Hs.eg.db, 
  keys = expression_df$Gene,
  keytype = "ENSEMBL", 
  column = "SYMBOL", 
  multiVals = "list"
)

# Convert the list to a data frame 
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "HUGO") %>%
  tidyr::unnest(cols = HUGO)

# Add the HUGO symbols as a new column to the expression_df
expression_with_hugo_df <- expression_df %>%
  dplyr::mutate(HUGO = mapped_df$HUGO[match(Gene, mapped_df$Ensembl)])

# Reorder columns to make HUGO the second column
expression_with_hugo_df <- expression_with_hugo_df %>%
  dplyr::select(Gene, HUGO, everything())

# Now, `expression_with_hugo_df` has the HUGO column as the second column
head(expression_with_hugo_df)




#helps keep track of data
multi_mapped <- mapped_df %>%
  #count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "Hugo_count") %>%
  #Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(Hugo_count))
head(multi_mapped)


dim(expression_with_hugo_df) #check size of expression matrix

log_expression_data <- log2(expression_data + 2) #perform log analysis (didnt work with hugo expression?)
median_expression <- apply(log_expression_data, 1, median) #calculate median gene expr
log_expression_with_hugo_df <- expression_with_hugo_df %>%
  dplyr::mutate(across(-c(Gene, HUGO), ~ log2(. + 2)))

# View the log-scaled dataframe
head(log_expression_with_hugo_df)

expression_range <- apply(log_expression_data, 1, function(x) max(x) - min(x))
head(expression_range)
ggplot(data.frame(Expression_Range = expression_range), aes(x = Expression_Range)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(
    title = "Density Plot of Per-Gene Median Expression Ranges",
    subtitle = "Distribution of Expression Range Across All Genes",
    x = "Expression Range (log2 scale)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),      # Center title and adjust size
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"), # Center subtitle and adjust size
    axis.title = element_text(size = 12),                                 # Axis titles size
    axis.text = element_text(size = 10)                                   # Axis text size
  )

#---------------------------------------------------------------------------------------------------------------------
#Differential Analysis
#Set randomness for differential analysis
set.seed(12345)
# separate data into diabetic and non groups
metadata <- metadata %>%
  dplyr::mutate(Diabetic_status = dplyr::case_when(
    stringr::str_detect(refinebio_subject, "t2d") ~ "t2d",
    stringr::str_detect(refinebio_subject, "pancreatic islets|non t2d") ~ "non-diabetic"
  ))

# Check the new Diabetic_status column
table(metadata$Diabetic_status)
# List the sample and the status
dplyr::select(metadata, refinebio_title, Diabetic_status)

# Make mutation_status a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    Diabetic_status = factor(Diabetic_status, levels = c("non-diabetic", "t2d"))
  )
filtered_expression_df <- log_expression_data %>%
  dplyr::filter(rowSums(.) >= 4000)
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~Diabetic_status
)
deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
head(deseq_results)
# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))
head(deseq_df)
# Volcano Time!
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
volcano_plot

#-------------------------

top50_genes <- deseq_df %>%
  dplyr::filter(!is.na(padj) & padj < 0.05) %>% # Keep only significant results
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>% # Sort by absolute log2 fold change
  dplyr::slice_head(n = 50) # Take the top 50 genes
print(top50_genes)
# Define the path for saving the file
top_50_path <- file.path(resultspath, "top_50_differentially_expressed_genes.tsv")

# Write the top 50 genes to a TSV file
write.table(top50_genes, top_50_path, sep = "\t", row.names = FALSE, quote = FALSE)
# End Differential Analysis and id the top 50 genes



# Extract the Gene names from top50_genes and subset the expression matrix
top50_gene_names <- top50_genes$Gene  # Extract the list of top 50 gene names

# Subset log-transformed expression matrix to include only top 50 genes
heatmap_data <- expression_with_hugo_df %>%
  dplyr::filter(Gene %in% top50_gene_names) %>%
  dplyr::select(-Gene, -HUGO) %>%  # Remove the Gene and HUGO columns
  as.matrix()  # Convert to matrix

# Step 2: Prepare a color mapping for Diabetic status (groupings)
sample_grouping <- metadata$Diabetic_status
names(sample_grouping) <- metadata$refinebio_accession_code

# Create a color map for the groups
group_colors <- c("non-diabetic" = "blue", "t2d" = "red")

# Create a heatmap annotation for Diabetic_status
ha <- HeatmapAnnotation(
  df = data.frame(Diabetic_status = sample_grouping),
  col = list(Diabetic_status = group_colors)
)

# Step 3: Generate the heatmap with sample grouping sidebar
heatmap <- Heatmap(
  heatmap_data,
  name = "Expression",
  show_row_names = TRUE,  # Show gene names on rows
  show_column_names = FALSE,  # Hide sample names on columns
  cluster_rows = TRUE,  # Cluster genes (rows)
  cluster_columns = TRUE,  # Cluster samples (columns)
  top_annotation = ha,  # Add the sample grouping annotation
  row_names_gp = gpar(fontsize = 8),  # Adjust font size for row names
  column_title = "Samples",  # Title for the columns (samples)
  row_title = "Top 50 Genes"  # Title for the rows (genes)
)

# Draw the heatmap
draw(heatmap)

#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
# PCA
# Set rownames of metadata
rownames(metadata) <- metadata$refinebio_accession_code

# Check if rownames of metadata match colnames of countData
all(rownames(metadata) %in% colnames(expression_data))

# Ensure that metadata is a data frame
metadata <- as.data.frame(metadata)

metadata$Diabetic_status <- as.factor(metadata$Diabetic_status)
threshold <- 100000  # Minimum count threshold
keep <- rowSums(expression_data) > threshold
expression_filtered <- expression_data[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = round(expression_filtered), # Ensure counts are integers for DESeq2
  colData = metadata,
  design = ~ Diabetic_status
)


dds <- DESeq(dds)
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = "Diabetic_status", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
sum(is.na(metadata$Diabetic_status))

# Plot the PCA with ggplot2
ggplot(pcaData, aes(x = PC1, y = PC2, color = Diabetic_status)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Plot of Gene Expression",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("t2d" = "red", "non-diabetic" = "blue")) +  # Assign custom colors
  theme(legend.title = element_blank())
# END PCA
#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
# Plot the t-SNE 

# Make a matrix excluding the Gene and HUGO columns
expression_matrix <- as.matrix(log_expression_with_hugo_df %>% select(-Gene, -HUGO))
print(expression_matrix)

# Assign the diabetic status to labels to split up the data
labels <- as.factor(metadata$Diabetic_status)

# Make sure the labels are correct
print(levels(labels))

# Create the t-SNE then show the plot
tsne_result <- M3C::tsne(expression_matrix, labels)
tsne_result

#-----------
# Plot UMAP

# Split the data into the numbers and the categories
log_expression_numeric <- as.matrix(log_expression_with_hugo_df %>% select(-Gene, -HUGO))
log_expression_labels <- log_expression_with_hugo_df[, "HUGO"]

# Create the UMAP
log_expression_with_hugo_df_umap <- umap(log_expression_numeric)

# Take the data from the UMAP and turn it into data frames
umap_coordinates <- log_expression_with_hugo_df_umap$data
umap_coordinates_df <- as.data.frame(umap_coordinates)
# Match the accession codes to get the diabetic status of each
umap_coordinates_df$Diabetic_status <- metadata$Diabetic_status[match(rownames(umap_coordinates_df), metadata$refinebio_accession_code)]

# Creating the plot
ggplot(umap_coordinates_df, aes(x = X1, y = X2, color = Diabetic_status)) +
  geom_point(size = 3) +
  labs(title = "UMAP Plot of Gene Expression",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  scale_color_manual(values = c("t2d" = "red", "non-diabetic" = "blue")) +
  theme_minimal() +
  theme(legend.title = element_blank())

#End of t-SNE and UMAP
#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
# Wilcoxon Rank-Sum Test (Victoria)

# Grouping the data based on diabetic status
diabetic_s <- expression_with_hugo_df[, metadata$Diabetic_status == "t2d"]
nondiabetic_s <- expression_with_hugo_df[, metadata$Diabetic_status == "non-diabetic"]

# Run the test but exclude Gene and HUGO then run on the numeric data and extract the p-val
wilcoxon_results <- apply(expression_with_hugo_df[, -c(1, 2)], 1, function(gene_expression) {
  wilcox.test(gene_expression ~ metadata$Diabetic_status)$p.value
})

# Get the data frames of the Gene, HUGO, and p-val
wilcoxon_df <- data.frame(Gene = expression_with_hugo_df$Gene,
                          HUGO = expression_with_hugo_df$HUGO,
                          p_value = wilcoxon_results)

print(wilcoxon_df)

# Arrange the p-vals
wilcoxon_df <- wilcoxon_df %>%
  dplyr::arrange(p_value)
print(wilcoxon_df)

# Create the file with the arranged p-vals
wilcoxon_path <- file.path(resultspath, "wilcoxon_test_results.tsv")
write.table(wilcoxon_df, wilcoxon_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Get the significant p-vals
significant_genes <- wilcoxon_df %>%
  dplyr::filter(p_value < 0.05)

# Starting Gene Ontology--------------------------------------------------------

# Get just the HUGO
significant_hugo_genes <- significant_genes$HUGO
print(significant_genes$HUGO)

# Perform GO enrichment analysis using clusterProfiler
go_enrichment <- clusterProfiler::enrichGO(
  gene = significant_hugo_genes, 
  OrgDb = org.Hs.eg.db, 
  keyType = "SYMBOL", 
  ont = "BP", 
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05           
)

# Create data frames for the enrichment
go_df <- as.data.frame(go_enrichment)

print(go_df)

# Create the file
go_enrichment_path <- file.path(resultspath, "GO_enrichment_results.tsv")
write.table(as.data.frame(go_df), go_enrichment_path, sep = "\t", row.names = FALSE, quote = FALSE)

#End of Wilcoxon Rank-Sum Test and Gene Ontology
#---------------------------------------------------------------------------------------------------------------------
# Starting Disease Ontology Enrichment Analysis for Wilcoxon Rank-Sum Test (Dylan Fain)-------------------------------------

# Get just the HUGO for significant genes
significant_hugo_genes <- expression_with_hugo_df$HUGO
print(significant_hugo_genes)

#packages
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")


# Load required libraries
library(DOSE)  # Ensure this library is loaded since enrichDO is part of the DOSE package

# Assuming significant_hugo_genes contains the HUGO symbols of the significant genes
significant_hugo_genes <- expression_with_hugo_df$HUGO
print(significant_genes)
# Convert HUGO symbols to Entrez IDs if necessary
# You might need to use a package like biomaRt or org.Hs.eg.db to map HUGO to Entrez IDs
# For example, using org.Hs.eg.db:
library(org.Hs.eg.db)

# Get Entrez IDs for the HUGO symbols
entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = significant_hugo_genes, 
                     column = "ENTREZID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")
print(entrez_ids)
# Filter out any NA values (genes that couldn't be mapped)
entrez_ids <- na.omit(entrez_ids)

# Perform Disease Ontology enrichment analysis
do_enrichment <- enrichDO(
  gene = entrez_ids,               # Use the mapped Entrez IDs
  ont = "HDO",                      # Specify the ontology; "DO" for Disease Ontology
  organism = "hsa",               # Specify the organism; "hsa" for Homo sapiens
  pvalueCutoff = 0.05,            # P-value cutoff for significant results
  qvalueCutoff = 0.05,             # Q-value cutoff for FDR control
)
print(do_enrichment)
# Check the results of the enrichment analysis
summary(do_enrichment)
do_enrichment_df <- as.data.frame(do_enrichment@result)


# Write the results to a file
do_enrichment_path <- file.path(resultspath, "DO_enrichment_results.tsv")
write.table(do_enrichment_df, do_enrichment_path, sep = "\t", row.names = FALSE, quote = FALSE)
#End Section
#--------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
# gprofiler analysis and Gene Ontology (Kevin)
# Extract the gene symbols from the top50_genes dataframe
gene_list <- top50_genes$Gene

# Turn to HUGO
gene_list_hugo <- mapIds(
  org.Hs.eg.db, 
  keys = gene_list, 
  column = "SYMBOL", 
  keytype = "ENSEMBL", 
  multiVals = "first"
)

# Remove any NAs (genes that didn't map correctly)
gene_list_hugo <- na.omit(gene_list_hugo)

# Run gProfiler using HUGO symbols
gprofiler_results <- gost(
  query = gene_list_hugo,          
  organism = "hsapiens",      
  sources = c("GO:BP", "GO:CC", "GO:MF"),
  significant = FALSE
)

gsea_results_df <- gprofiler_results$result

# Combine the list so it can be made into a .tsv file
if ("parents" %in% colnames(gsea_results_df)) {
  gsea_results_df$parents <- sapply(gsea_results_df$parents, function(x) paste(x, collapse = ", "))
}

# Create the file
gsea_results_df_path <- file.path(resultspath, "gprofiler_results.tsv")
write.table(gsea_results_df, gsea_results_df_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Define a custom p-value threshold
custom_threshold <- 0.99

# Filter results based on the custom threshold
significant_terms <- gsea_results_df %>%
  filter(p_value < custom_threshold)

# View the filtered significant terms
print(significant_terms)

# Dotplot
# Create a dot plot for significantly enriched GO terms
ggplot(significant_terms, aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +
  geom_point(aes(size = -log10(p_value), color = source), alpha = 0.7) +
  scale_size_continuous(range = c(2, 8)) +  # Adjust size range of points
  labs(title = "Significantly Enriched Gene Ontology Terms",
       x = "GO Term",
       y = "-log10(P-value)",
       color = "Source") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#TopGO(Adam)

                                    

if (!requireNamespace("topGO", quietly = TRUE)) {
  BiocManager::install("topGO")
}
library(topGO)

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)



#TopGo(Adam)
#p_values extracted
p_value_of_genes <- deseq_df_with_hugo$pvalue
Hugo_Names <- deseq_df_with_hugo$HUGO


names(p_value_of_genes) <- Hugo_Names


#  gene list with p-values as a named vector
gene_list <- p_value_of_genes

#pvalues mapped to binary with 1 signifcating, 0 insignifcant 
gene_factor <- factor(as.integer(gene_list < 0.05))
names(gene_factor) <- names(gene_list)
#debug stuff
str(gene_list)
str(gene_factor)


#
GOdata <- new("topGOdata",
              ontology = "BP",  # You can change this to "MF" or "CC"
              allGenes = gene_factor,
              geneSel = function(x) (x == 1),
              annot = annFUN.org, 
              mapping = "org.Hs.eg.db", 
              ID = "symbol")

#fischer test 
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# top significant GO terms
topGOres <- GenTable(GOdata, 
                     classicFisher = resultFisher, 
                     orderBy = "classicFisher", 
                     topNodes = 10)  # top signifcant, can be adjusted to however many wanted

print(topGOres)

write.csv(topGOres, "GO_enrichment_results.csv", row.names = FALSE)


#CSV file
write.csv(topGOres, file = file.path(resultspath, "GO_enrichment_results.csv"), row.names = FALSE)


                                    
# End of gprofiler analysis and Gene Ontology
#---------------------------------------------------------------------------------------------------------------------
 # gprofiler analysis and Disease Ontology (Andrew)
# Extract the gene symbols from the top50_genes dataframe
gene_list <- top50_genes$Gene

# Turn to HUGO
gene_list_hugo <- mapIds(
  org.Hs.eg.db, 
  keys = gene_list, 
  column = "SYMBOL", 
  keytype = "ENSEMBL", 
  multiVals = "first"
)









# Remove any NAs (genes that didn't map correctly)
gene_list_hugo <- na.omit(gene_list_hugo)

gprofiler_results <- gost(
  query = gene_list_hugo,          
  organism = "hsapiens",      
  sources = c("HP"),  # Use "HPO" for Disease Ontology
  significant = FALSE  # Adjust as needed
)

# Convert the results into a data frame
gsea_results_df <- gprofiler_results$result


# Combine the list so it can be made into a .tsv file
if ("parents" %in% colnames(gsea_results_df)) {
  gsea_results_df$parents <- sapply(gsea_results_df$parents, function(x) paste(x, collapse = ", "))
}

# Create the file
gsea_results_df_path <- file.path(resultspath, "gprofiler_disease_results.tsv")
write.table(gsea_results_df, gsea_results_df_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Define a custom p-value threshold (example: 0.99)
custom_threshold <- 0.99





# Filter results based on the custom threshold
significant_terms <- gsea_results_df %>%
  filter(p_value < custom_threshold)

# View the filtered significant terms
print(significant_terms)




library(knitr)

# Select relevant columns for the table
significant_terms_table <- significant_terms %>%
  select(term_name, p_value, source, term_id)

# Display the static table
kable(significant_terms_table, caption = "Significantly Enriched Disease Ontology Terms")

significant_terms_path <- file.path(resultspath, "significant_terms_disease_ontology.tsv")
write.table(significant_terms_table, significant_terms_path, sep = "\t", row.names = FALSE, quote = FALSE)                                   

#---------------------------------------------------------------------------------------------------------------------
#(Dylan) - Unsupervised Analysis - Gaussian mixture models
# Install mclust package if not already installed
if (!require("mclust")) {
  install.packages("mclust")
}
library(mclust)

# Apply variance stabilizing transformation directly to the DESeq2 object
vsd <- varianceStabilizingTransformation(deseq_object)
# Calculate the variance for each gene from the transformed data
gene_variances_vst <- apply(assay(vsd), 1, var)

# Get the indices of the top 5000 most variable genes
top_5000_genes_vst <- order(gene_variances_vst, decreasing = TRUE)[1:5000]

# Subset the transformed data to include only the top 5000 genes
top_5000_gene_data_vst <- assay(vsd)[top_5000_genes_vst, ]

# Transpose the data for clustering
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <- top_5000_gene_data_vst_t[, 1:4591]



# Apply Gaussian Mixture Model clustering on the cleaned data
gmm_results_vst_clean <- Mclust(top_5000_gene_data_vst_t)
# View the summary of the model
summary(gmm_results_vst_clean)


#Vary the parameters to get different number of clusters and view the results
gmm_results_vst_clean <- Mclust(top_5000_gene_data_vst_t, G = 1:4)
summary(gmm_results_vst_clean)
gmm_results_vst_clean <- Mclust(top_5000_gene_data_vst_t, modelNames = "EII")
summary(gmm_results_vst_clean)

# Plot the classification results
plot(gmm_results_vst_clean, what = "classification")


# Define a function to run GMM on a specified number of top variable genes
run_gmm_on_genes <- function(num_genes) {
  # Select the top 'num_genes' from the vst-transformed data
  #selected_genes <- head(order(rowVars(top_5000_gene_data_vst), decreasing = TRUE), num_genes)
  top_gene_data <- top_5000_gene_data_vst[num_genes, ]
  
  # Transpose the data for clustering
  top_gene_data_t <- t(top_gene_data)
  
  # Fit GMM with varying parameters
  gmm_results <- Mclust(top_gene_data_t) # Adjust G based on expected clusters
  
  # Return the results
  return(gmm_results)
}

# Run GMM for different numbers of genes
results_10 <- run_gmm_on_genes(10)
results_100 <- run_gmm_on_genes(100)
results_1000 <- run_gmm_on_genes(1000)
results_4000 <- run_gmm_on_genes(4000)

# Store results for comparison
results_list <- list(
  genes_10 = results_10,
  genes_100 = results_100,
  genes_1000 = results_1000,
  genes_4000 = results_4000
)

# Print summaries of the results
for (name in names(results_list)) {
  cat("Results for", name, ":\n")
  print(summary(results_list[[name]]))
  cat("\n")
}

# Function to perform chi-squared test on clustering results, with checks for single clusters
chi_square_test_clusters <- function(clusters1, clusters2) {
  # Create a contingency table from the two cluster assignments
  contingency_table <- table(clusters1, clusters2)
  
  
  # Perform the chi-squared test
  chi_test <- chisq.test(contingency_table)
  
  # Return the test results
  return(chi_test)
}

# Extract the clustering assignments from each result
clusters_10 <- results_10$classification
clusters_100 <- results_100$classification
clusters_1000 <- results_1000$classification
clusters_4000 <- results_4000$classification

# Perform chi-squared tests for each pair of results
chi_10_100 <- chi_square_test_clusters(clusters_10, clusters_100)
chi_10_1000 <- chi_square_test_clusters(clusters_10, clusters_1000)
chi_10_4000 <- chi_square_test_clusters(clusters_10, clusters_4000)
chi_100_1000 <- chi_square_test_clusters(clusters_100, clusters_1000)
chi_100_4000 <- chi_square_test_clusters(clusters_100, clusters_4000)
chi_1000_4000 <- chi_square_test_clusters(clusters_1000, clusters_4000)

# Store the p-values from each test for the table
chi_results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 4000", 
                 "100 vs 1000", "100 vs 4000", "1000 vs 4000"),
  Chi_squared = c(chi_10_100$statistic, chi_10_1000$statistic, chi_10_4000$statistic, 
                  chi_100_1000$statistic, chi_100_4000$statistic, chi_1000_4000$statistic),
  p_value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_4000$p.value, 
              chi_100_1000$p.value, chi_100_4000$p.value, chi_1000_4000$p.value)
)

# Print the chi-squared test results
print(chi_results)

install.packages("ggalluvial")

# Load necessary libraries
library(ggplot2)
library(ggalluvial)

# Sample data (you will replace this with your actual clustering results)
data <- data.frame(
  sample = 1:1600,  # Assuming 1600 samples
  cluster_10_genes = clusters_10,    # Replace with actual cluster assignments
  cluster_100_genes = clusters_100,  # Replace with actual cluster assignments
  cluster_1000_genes = clusters_1000, # Replace with actual cluster assignments
  cluster_4000_genes = clusters_4000  # Replace with actual cluster assignments
)

# Melt the data into long format
library(reshape2)
long_data <- melt(data, id.vars = "sample")

# Create the alluvial plot
ggplot(long_data, aes(x = variable, stratum = value, alluvium = sample, fill = value, label = value)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray") +
  geom_stratum() +
  theme_minimal() +
  labs(title = "Alluvial Diagram of Cluster Memberships",
       x = "Clustering Setup", y = "Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#End of Dylan's Section
#-----------------------------------------------------------------------------------------

#-------------------------------------------------------------------
# (Andrew) Unsupervised analysis - PAM clustering
if (!require("cluster")) {
  install.packages("cluster")
}
library(cluster)

# Apply variance stabilizing transformation directly to the DESeq2 object
vsd <- varianceStabilizingTransformation(deseq_object)
gene_variances_vst <- apply(assay(vsd), 1, var)
top_5000_genes_vst <- order(gene_variances_vst, decreasing = TRUE)[1:5000]
top_5000_gene_data_vst <- assay(vsd)[top_5000_genes_vst, ]
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <- top_5000_gene_data_vst_t[, 1:4591]

# Function to run PAM clustering
run_pam_on_genes <- function(num_genes, k = 3) {
  top_gene_data <- top_5000_gene_data_vst[1:num_genes, ]
  top_gene_data_t <- t(top_gene_data)
  pam_results <- pam(top_gene_data_t, k = k)
  return(pam_results)
}

# Run PAM for different numbers of genes
results_10 <- run_pam_on_genes(10, k = 3)
results_100 <- run_pam_on_genes(100, k = 3)
results_1000 <- run_pam_on_genes(1000, k = 3)
results_4000 <- run_pam_on_genes(4000, k = 3)

# Store and print results
results_PAM <- list(
  genes_10 = results_10,
  genes_100 = results_100,
  genes_1000 = results_1000,
  genes_4000 = results_4000
)

for (name in names(results_PAM)) {
  cat("Results for", name, ":\n")
  print(results_PAM[[name]])
  cat("\n")
}

# Function to perform chi-squared test on clustering results
chi_square_test_clusters <- function(clusters1, clusters2) {
  contingency_table <- table(clusters1, clusters2)
  chi_test <- chisq.test(contingency_table)
  return(chi_test)
}

# Extract clustering assignments
clusters_10 <- results_10$clustering
clusters_100 <- results_100$clustering
clusters_1000 <- results_1000$clustering
clusters_4000 <- results_4000$clustering

# Perform chi-squared tests
chi_10_100 <- chi_square_test_clusters(clusters_10, clusters_100)
chi_10_1000 <- chi_square_test_clusters(clusters_10, clusters_1000)
chi_10_4000 <- chi_square_test_clusters(clusters_10, clusters_4000)
chi_100_1000 <- chi_square_test_clusters(clusters_100, clusters_1000)
chi_100_4000 <- chi_square_test_clusters(clusters_100, clusters_4000)
chi_1000_4000 <- chi_square_test_clusters(clusters_1000, clusters_4000)

# Store the chi-squared test results
chi_results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 4000", 
                 "100 vs 1000", "100 vs 4000", "1000 vs 4000"),
  Chi_squared = c(chi_10_100$statistic, chi_10_1000$statistic, chi_10_4000$statistic, 
                  chi_100_1000$statistic, chi_100_4000$statistic, chi_1000_4000$statistic),
  p_value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_4000$p.value, 
              chi_100_1000$p.value, chi_100_4000$p.value, chi_1000_4000$p.value)
)

# Print chi-squared test results
print(chi_results)

# Install and load ggalluvial for plotting
install.packages("ggalluvial")
library(ggplot2)
library(ggalluvial)
library(reshape2)

# Sample data for plotting
data <- data.frame(
  sample = 1:1600,
  cluster_10_genes = clusters_10,
  cluster_100_genes = clusters_100,
  cluster_1000_genes = clusters_1000,
  cluster_4000_genes = clusters_4000
)

# Melt the data into long format
long_data <- melt(data, id.vars = "sample")

# Create the alluvial plot
ggplot(long_data, aes(x = variable, stratum = value, alluvium = sample, fill = value, label = value)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray") +
  geom_stratum() +
  theme_minimal() +
  labs(title = "Alluvial Diagram of Cluster Memberships using PAM",
       x = "Clustering Setup", y = "Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#End of Andrew's Section
#-------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#(Kevin) - Unsupervised Analysis - Spectral Clustering

if (!require("cluster")) {
  install.packages("cluster")
}
if (!require("kernlab")) {
  install.packages("kernlab")
}
if (!require("ggalluvial")) {
  install.packages("ggalluvial")
}
if (!require("dplyr")) {
  install.packages("dplyr")
}
if (!require("tidyr")) {
  install.packages("tidyr")
}

library(dplyr)
library(tidyr)
library(ggalluvial)
library(cluster)
library(kernlab)

# Function to perform spectral clustering and PCA plot for a given number of genes/clusters
run_clustering_pca <- function(num_genes, data_matrix, num_clusters = 3) {
  # Select genes ranked by variance
  selected_genes <- order(apply(data_matrix, 1, var), decreasing = TRUE)[1:num_genes]
  gene_data <- data_matrix[selected_genes, ]
  
  # Transpose the data for clustering
  gene_data_t <- t(gene_data)
  
  # Spectral clustering
  spectral_clustering_result <- specc(gene_data_t, centers = num_clusters)
  cluster_labels <- spectral_clustering_result@.Data
  
  # PCA for 2d visualization
  pca_res <- prcomp(gene_data_t, scale. = TRUE)
  
  # Plot PCA with clusters
  pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Cluster = factor(cluster_labels))
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 3) +
    labs(title = paste("Spectral Clustering PCA with", num_genes, "Genes"), x = "PC1", y = "PC2") +
    theme_minimal()
  
  return(list(pca_plot = pca_plot, cluster_labels = cluster_labels))
}

# Clustering and plot PCA for different numbers of genes
results_10 <- run_clustering_pca(10, top_5000_gene_data_vst)
results_100 <- run_clustering_pca(100, top_5000_gene_data_vst)
results_1000 <- run_clustering_pca(1000, top_5000_gene_data_vst)
results_4500 <- run_clustering_pca(4500, top_5000_gene_data_vst)

# Store the previously computed results in a list
results_spectral <- list(
  results_10 = results_10,
  results_100 = results_100,
  results_1000 = results_1000,
  results_4500 = results_4500
)

# Show plots
print(results_10$pca_plot)
print(results_100$pca_plot)
print(results_1000$pca_plot)
print(results_4500$pca_plot)

# Chi Squared Test Function
perform_chi_squared_test <- function(cluster_labels1, cluster_labels2) {
  contingency_table <- table(cluster_labels1, cluster_labels2)
  
  # Perform chi-squared test
  chisq_result <- chisq.test(contingency_table)
  
  return(list(chisq_statistic = chisq_result$statistic, 
              p_value = chisq_result$p.value, 
              contingency_table = contingency_table))
}

chi_squared_results <- list()

# Perform chi-squared tests between each pair of clustering results
chi_squared_results[["10 vs 100"]] <- perform_chi_squared_test(results_10$cluster_labels, results_100$cluster_labels)
chi_squared_results[["10 vs 1000"]] <- perform_chi_squared_test(results_10$cluster_labels, results_1000$cluster_labels)
chi_squared_results[["10 vs 4500"]] <- perform_chi_squared_test(results_10$cluster_labels, results_4500$cluster_labels)
chi_squared_results[["100 vs 1000"]] <- perform_chi_squared_test(results_100$cluster_labels, results_1000$cluster_labels)
chi_squared_results[["100 vs 4500"]] <- perform_chi_squared_test(results_100$cluster_labels, results_4500$cluster_labels)
chi_squared_results[["1000 vs 4500"]] <- perform_chi_squared_test(results_1000$cluster_labels, results_4500$cluster_labels)

# Create a table of chi-squared test results
chi_squared_summary <- data.frame(
  Comparison = names(chi_squared_results),
  Chi_Squared_Statistic = sapply(chi_squared_results, function(x) x$chisq_statistic),
  P_Value = sapply(chi_squared_results, function(x) x$p_value)
)

# Store clustering results for alluvial diagram
alluvial_data <- data.frame(
  Sample = 1:length(results_10$cluster_labels),
  Clustering_10 = as.factor(results_10$cluster_labels),
  Clustering_100 = as.factor(results_100$cluster_labels),
  Clustering_1000 = as.factor(results_1000$cluster_labels),
  Clustering_4500 = as.factor(results_4500$cluster_labels)
)

# Data to long format for alluvial diagram
alluvial_long <- alluvial_data %>%
  pivot_longer(cols = -Sample, names_to = "Clustering_Setup", values_to = "Cluster")

# Create the alluvial diagram
ggplot(alluvial_long, aes(x = Clustering_Setup, stratum = Cluster, alluvium = Sample, fill = Cluster, label = Cluster)) +
  geom_flow(stat = "alluvium", aes.flow = "forward") +
  geom_stratum() +
  scale_x_discrete(limits = c("Clustering_10", "Clustering_100", "Clustering_1000", "Clustering_4500")) +
  theme_minimal() +
  labs(title = "Alluvial Diagram of Cluster Changes", x = "Clustering Setup", y = "Samples")

#End of Kevin's Section
#-----------------------------------------------------------------------------------------
#(Adam's section start) Hclust


vsd <- varianceStabilizingTransformation(deseq_object)
gene_variances_vst <- apply(assay(vsd), 1, var)
top_5000_genes_vst <- order(gene_variances_vst, decreasing = TRUE)[1:5000]
top_5000_gene_data_vst <- assay(vsd)[top_5000_genes_vst, ]
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <- top_5000_gene_data_vst_t[, 1:4591]
dist_matrix <- dist(top_5000_gene_data_vst_t, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

perform_hclust <- function(data, num_genes, num_clusters) {
  subset_data <- data[, 1:num_genes]
  dist_matrix <- dist(subset_data, method = "euclidean")
  hclust_result <- hclust(dist_matrix, method = "ward.D2")
  cluster_assignments <- cutree(hclust_result, k = num_clusters)
  return(cluster_assignments)
}

num_clusters <- 3
clusters_10_a <- perform_hclust(top_5000_gene_data_vst_t, 10, num_clusters)
clusters_100_a <- perform_hclust(top_5000_gene_data_vst_t, 100, num_clusters)
clusters_1000_a <- perform_hclust(top_5000_gene_data_vst_t, 1000, num_clusters)
clusters_4000_a <- perform_hclust(top_5000_gene_data_vst_t, 4000, num_clusters)

# Store results in a list
results_list_adam <- list(
  genes_10a = clusters_10_a,
  genes_100a = clusters_100_a,
  genes_1000a = clusters_1000_a,
  genes_4000a = clusters_4000_a
)

for (name in names(results_list_adam)) {
  cat("Results for", name, ":\n")
  print(summary(results_list_adam[[name]]))
  cat("\n")
}
                                

# End of Adam's section                                 
#-----------------------------------------------------------------------------------------
# (Victoria) Unsupervised analysis - Affinity Propagation Clustering
if (!require("apcluster")) {
  install.packages("apcluster")
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
}
if (!require("ggalluvial")) {
  install.packages("ggalluvial")
}
if (!require("reshape2")) {
  install.packages("reshape2")
}

library(apcluster)
library(ggplot2)
library(ggalluvial)
library(reshape2)

# Variance stabilizing transformation and transposing
vsd <- varianceStabilizingTransformation(deseq_object)
gene_variances_vst <- apply(assay(vsd), 1, var)
top_5000_genes_vst <- order(gene_variances_vst, decreasing = TRUE)[1:5000]
top_5000_gene_data_vst <- assay(vsd)[top_5000_genes_vst, ]
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <- top_5000_gene_data_vst_t[, 1:4591]

# NA value check
any(is.na(top_5000_gene_data_vst_t))

# Function for Affinity Propagation clustering
run_affinity_propagation <- function(num_genes) {
  top_gene_data <- top_5000_gene_data_vst[1:num_genes, ]
  top_gene_data_t <- t(top_gene_data)
  ap_result <- apcluster(negDistMat(r = 2), top_gene_data_t, maxits = 1000, convits = 100,lam = 0.7) 
  cluster_labels <- ap_result@clusters
  # Cluster labels
  formatted_clusters <- rep(NA, ncol(top_gene_data_t))
  for (i in seq_along(ap_result@clusters)) {
    formatted_clusters[ap_result@clusters[[i]]] <- i
  }
  return(list(ap_result = ap_result, cluster_labels = formatted_clusters))
}

# Affinity Propagation for different number of genes
results_10 <- run_affinity_propagation(10)
results_100 <- run_affinity_propagation(100)
results_1000 <- run_affinity_propagation(1000)
results_4500_Kmeans <- run_affinity_propagation(4500)

# Store results
results_list_Kmeans <- list(
  genes_10 = results_10,
  genes_100 = results_100,
  genes_1000 = results_1000,
  genes_4500 = results_4500_Kmeans
)
#print(results_4500_Kmeans$ap_result)

# Print results
for (name in names(results_list_Kmeans)) {
  cat("Results for", name, ":\n")
  print(results_list_Kmeans[[name]]$ap_result)
  cat("\n")
}

# Function for chi-squared test on results
chi_square_test_clusters <- function(clusters1, clusters2) {
  contingency_table <- table(clusters1, clusters2)
  chi_test <- chisq.test(contingency_table)
  return(chi_test)
}

# Extract
clusters_10 <- results_10$cluster_labels
clusters_100 <- results_100$cluster_labels
clusters_1000 <- results_1000$cluster_labels
clusters_4500 <- results_4500$cluster_labels

# Make sure they're the same size
clusters_10 <- clusters_10[1:1600]
clusters_100 <- clusters_100[1:1600]
clusters_1000 <- clusters_1000[1:1600]
clusters_4500 <- clusters_4500[1:1600] 

# Check lengths
length(clusters_10) 
length(clusters_100)
length(clusters_1000)
length(clusters_4500) 

# Chi-squared tests
chi_10_100 <- chi_square_test_clusters(clusters_10, clusters_100)
chi_10_1000 <- chi_square_test_clusters(clusters_10, clusters_1000)
chi_10_4500 <- chi_square_test_clusters(clusters_10, clusters_4500)
chi_100_1000 <- chi_square_test_clusters(clusters_100, clusters_1000)
chi_100_4500 <- chi_square_test_clusters(clusters_100, clusters_4500)
chi_1000_4500 <- chi_square_test_clusters(clusters_1000, clusters_4500)

# Store chi-squared test results
chi_results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 4500", 
                 "100 vs 1000", "100 vs 4500", "1000 vs 4500"),
  Chi_squared = c(chi_10_100$statistic, chi_10_1000$statistic, chi_10_4500$statistic, 
                  chi_100_1000$statistic, chi_100_4500$statistic, chi_1000_4500$statistic),
  p_value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_4500$p.value, 
              chi_100_1000$p.value, chi_100_4500$p.value, chi_1000_4500$p.value)
)

# Print
print(chi_results)

# Another check
table(clusters_10)
table(clusters_100)
table(clusters_1000)
table(clusters_4500)

# Use ggalluvial for plotting
library(ggplot2)
library(ggalluvial)
library(reshape2)

# Sample data for plotting
data <- data.frame(
  sample = 1:1600,
  cluster_10_genes = clusters_10,
  cluster_100_genes = clusters_100,
  cluster_1000_genes = clusters_1000,
  cluster_4500_genes = clusters_4500
)

# Melt data into long format
long_data <- melt(data, id.vars = "sample")

# Alluvial plot
ggplot(long_data, aes(x = variable, stratum = value, alluvium = sample, fill = value, label = value)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray") +
  geom_stratum() +
  theme_minimal() +
  labs(title = "Alluvial Diagram of Cluster Memberships using Affinity Propagation",
       x = "Clustering Setup", y = "Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# End of Victoria's Section
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# Begin Heatmap

install.packages("pheatmap")
# Load necessary libraries
library(pheatmap)

# Assuming 'data_matrix' is your gene expression matrix with rows as genes and columns as samples
# Placeholder variables for clustering results (to be provided by your groupmates)
# GMM is already done; K-means, hclust, and PAM will be added by your groupmates

# Placeholder cluster results (replace these with actual results later)
cluster_GMM <- results_list$genes_4000$classification  # GMM results already implemented by your group
cluster_AffinityPropagation <- results_4500_Kmeans$cluster_labels  # Placeholder for K-means clustering results
cluster_Hclust <- clusters_4000_a  # Placeholder for hierarchical clustering results
cluster_PAM <- results_PAM$genes_4000$clustering        # Placeholder for PAM clustering results
Cluster_spect <- results_4500$cluster_labels

# Placeholder for sample groups from Assignment 1
sample_groups <- sample_grouping  # Replace with actual sample group data

length(cluster_AffinityPropagation)
length(results_4500$cluster_labels)
length(cluster_Hclust)

# Create annotation data (each column corresponds to the cluster result from one method)
annotation_data <- data.frame(
  Cluster_GMM = cluster_GMM,
  cluster_AffinityPropagation = cluster_AffinityPropagation[1:1600],
  Cluster_Hclust = cluster_Hclust,
  Cluster_PAM = cluster_PAM,
  Cluster_spect = Cluster_spect,
  Sample_Group = sample_groups
)



# Define colors for annotations (this is optional, but improves readability)
annotation_colors <- list(
  Cluster_GMM = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3", "5" = "#FF6666", "6" = "#99CCFF"),
  Cluster_AffinityPropagation = rainbow(106),  # Use a rainbow palette for 106 clusters
  Cluster_PAM = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"),
  Cluster_spect = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"),
  Cluster_Hclust = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"),
  Sample_Group = c("non-diabetic" = "#FFD700", "t2d" = "#FF4500")  # Colors for sample group annotations
) 
heatmap_data <- expression_with_hugo_df %>%
  dplyr::filter(Gene %in% top50_gene_names) %>%
  dplyr::select(-Gene, -HUGO) %>%  # Remove the Gene and HUGO columns
  as.matrix()  # Convert to matrix

# Generate the heatmap with row/column dendrograms and annotations
pheatmap_result <-pheatmap(heatmap_data,
         annotation_col = annotation_data,  # Cluster annotations
         clustering_method = "ward.D2",     # Hierarchical clustering for rows and columns
         cluster_rows = T,                # Cluster rows (genes)
         cluster_cols = T,
         show_rownames = FALSE,             # Optionally hide row names if too many genes
         show_colnames = FALSE,             # Optionally hide column names if too many samples
         scale = "row",                     # Normalize expression levels across rows (genes)
         fontsize = 8,
         legend = TRUE,
         color = colorRampPalette(c("red", "white"))(100),  # Color gradient for heatmap
         annotation_colors = annotation_colors,  # Colors for cluster/sample group annotations
         main = "Heatmap of 5000 Genes with Clustering and Sample Group Annotations")

print(pheatmap_result)



#End Heatmap ---------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# Statistics
                                 
# Spectral

# Create a list to store the chi-squared test results comparing diabetic status with clustering results
chi_squared_diabetic_results <- list()

# Perform chi-squared tests for each clustering result
for (i in c("10", "100", "1000", "4500")) {
  cluster_labels <- results_spectral[[paste0("results_", i)]]$cluster_labels
  diabetic_subset <- diabetic_labels[1:length(cluster_labels)]
  contingency_table <- table(cluster_labels, diabetic_subset)
  chisq_test <- chisq.test(contingency_table)
  
  # Store results with updated names
  chi_squared_diabetic_results[[paste0("spectral_", i, " vs Diabetic Status")]] <- list(
    chisq_statistic = chisq_test$statistic,
    p_value = chisq_test$p.value,
    contingency_table = contingency_table
  )
}

# Create a summary table of chi-squared test results for diabetic status
chi_squared_diabetic_summary <- data.frame(
  Comparison = names(chi_squared_diabetic_results),
  Chi_Squared_Statistic = sapply(chi_squared_diabetic_results, function(x) x$chisq_statistic),
  P_Value = sapply(chi_squared_diabetic_results, function(x) x$p_value)
)

# Adjust p-values for multiple hypothesis testing
chi_squared_diabetic_summary$Adjusted_P_Value <- p.adjust(chi_squared_diabetic_summary$P_Value, method = "fdr")

# Print the final summary
print(chi_squared_diabetic_summary)

# Affinity

# Perform chi-squared tests for each clustering result
for (i in c("10", "100", "1000", "4500")) {
  cluster_labels <- results_list_Kmeans[[paste0("genes_", i)]]$cluster_labels
  diabetic_subset <- diabetic_labels[1:length(cluster_labels)]
  contingency_table <- table(cluster_labels, diabetic_subset)
  chisq_test <- chisq.test(contingency_table)
  
  # Store results with updated names
  chi_squared_diabetic_results[[paste0("affinity_", i, " vs Diabetic Status")]] <- list(
    chisq_statistic = chisq_test$statistic,
    p_value = chisq_test$p.value,
    contingency_table = contingency_table
  )
}


# Create a summary table of chi-squared test results for diabetic status
chi_squared_diabetic_summary <- data.frame(
  Comparison = names(chi_squared_diabetic_results),
  Chi_Squared_Statistic = sapply(chi_squared_diabetic_results, function(x) x$chisq_statistic),
  P_Value = sapply(chi_squared_diabetic_results, function(x) x$p_value)
)

# Adjust p-values for multiple hypothesis testing
chi_squared_diabetic_summary$Adjusted_P_Value <- p.adjust(chi_squared_diabetic_summary$P_Value, method = "fdr")

# Print the final summary
print(chi_squared_diabetic_summary)

# PAM

# Perform chi-squared tests for each clustering result against diabetic status
for (i in c("10", "100", "1000", "4000")) {
  cluster_labels1 <- results_PAM[[paste0("genes_", i)]]$clustering
  diabetic_subset <- diabetic_labels[1:length(cluster_labels1)]
  contingency_table <- table(cluster_labels1, diabetic_subset)
  chisq_test <- chisq.test(contingency_table)
  
  # Store results with updated names
  chi_squared_diabetic_results[[paste0("pam_", i, " vs Diabetic Status")]] <- list(
    chisq_statistic = chisq_test$statistic,
    p_value = chisq_test$p.value,
    contingency_table = contingency_table
  )
}

# Create a summary table of chi-squared test results for diabetic status
chi_squared_diabetic_summary <- data.frame(
  Comparison = names(chi_squared_diabetic_results),
  Chi_Squared_Statistic = sapply(chi_squared_diabetic_results, function(x) x$chisq_statistic),
  P_Value = sapply(chi_squared_diabetic_results, function(x) x$p_value)
)

# Adjust p-values for multiple hypothesis testing
chi_squared_diabetic_summary$Adjusted_P_Value <- p.adjust(chi_squared_diabetic_summary$P_Value, method = "fdr")

# Clean up the Comparison column
chi_squared_diabetic_summary$Comparison <- gsub("\\.X-squared$", "", chi_squared_diabetic_summary$Comparison) # If needed

# Print the final summary
print(chi_squared_diabetic_summary)

# Gaussian mixture models

# Perform chi-squared tests for each clustering result
for (i in c("10", "100", "1000", "4000")) {
  cluster_labels <- results_list[[paste0("genes_", i)]]$classification
  diabetic_subset <- diabetic_labels[1:length(cluster_labels)]
  contingency_table <- table(cluster_labels, diabetic_subset)
  chisq_test <- chisq.test(contingency_table)
  
  # Store results with updated names
  chi_squared_diabetic_results[[paste0("gmm_", i, " vs Diabetic Status")]] <- list(
    chisq_statistic = chisq_test$statistic,
    p_value = chisq_test$p.value,
    contingency_table = contingency_table
  )
}


# Create a summary table of chi-squared test results for diabetic status
chi_squared_diabetic_summary <- data.frame(
  Comparison = names(chi_squared_diabetic_results),
  Chi_Squared_Statistic = sapply(chi_squared_diabetic_results, function(x) x$chisq_statistic),
  P_Value = sapply(chi_squared_diabetic_results, function(x) x$p_value)
)

# Adjust p-values for multiple hypothesis testing
chi_squared_diabetic_summary$Adjusted_P_Value <- p.adjust(chi_squared_diabetic_summary$P_Value, method = "fdr")

# Print the final summary
print(chi_squared_diabetic_summary)

# hclust

# Perform chi-squared tests for each clustering result
for (i in c("10", "100", "1000", "4000")) {
  cluster_labels <- results_list_adam[[paste0("genes_",i,"a")]]
  diabetic_subset <- diabetic_labels[1:length(cluster_labels)]
  contingency_table <- table(cluster_labels, diabetic_subset)
  chisq_test <- chisq.test(contingency_table)
  
  # Store results with updated names
  chi_squared_diabetic_results[[paste0("hclust_", i, " vs Diabetic Status")]] <- list(
    chisq_statistic = chisq_test$statistic,
    p_value = chisq_test$p.value,
    contingency_table = contingency_table
  )
}


# Create a summary table of chi-squared test results for diabetic status
chi_squared_diabetic_summary <- data.frame(
  Comparison = names(chi_squared_diabetic_results),
  Chi_Squared_Statistic = sapply(chi_squared_diabetic_results, function(x) x$chisq_statistic),
  P_Value = sapply(chi_squared_diabetic_results, function(x) x$p_value)
)

# Adjust p-values for multiple hypothesis testing
chi_squared_diabetic_summary$Adjusted_P_Value <- p.adjust(chi_squared_diabetic_summary$P_Value, method = "fdr")

# Print the final summary
print(chi_squared_diabetic_summary)
#End Stats ---------------------------------------------------------------------------------------

#Supervised Analysis - Dylan ------------------------------------------------------

# Load necessary libraries
library(DESeq2)       # For working with DESeq2 objects
library(randomForest) # For Random Forest algorithm
library(caret)        # For confusionMatrix function and model evaluation
library(tidymodels)
library(dplyr)
library(ranger)

# Apply variance stabilizing transformation directly to the DESeq2 object
vsd <- varianceStabilizingTransformation(deseq_object)

# Calculate the variance for each gene from the transformed data
gene_variances_vst <- apply(assay(vsd), 1, var)

# Get the indices of the top 5000 most variable genes
top_5000_genes_vst <- order(gene_variances_vst, decreasing = TRUE)[1:5000]

# Subset the transformed data to include only the top 5000 genes
top_5000_gene_data_vst <- assay(vsd)[top_5000_genes_vst, ]

# Transpose the data for Random Forest
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <-top_5000_gene_data_vst_t[, 1:4591]


# Create a data frame with the top 5000 gene data and class labels
df <- data.frame(top_5000_gene_data_vst_t)
n_non_diabetic <- 547  # Replace with your actual number of non-diabetic samples
n_t2d <- 1053           # Replace with your actual number of type 2 diabetes samples
class_labels <- factor(c(rep("non-diabetic", n_non_diabetic), rep("t2d", n_t2d)))
df$class <- class_labels  # Add class labels to the data frame

# Set a seed for reproducibility
set.seed(123)

# Split the data into training and testing sets (70% train, 30% test)
train_index <- createDataPartition(df$class, p = 0.7, list = FALSE)
train_data <- df[train_index, ]
test_data <- df[-train_index, ]

# Train a Random Forest model
rf_model <- randomForest(class ~ ., data = train_data, importance = TRUE)

# Print the model details
print(rf_model)

# Make predictions on the test set
predictions <- predict(rf_model, newdata = test_data)

# Evaluate the model's performance
conf_matrix <- confusionMatrix(predictions, test_data$class)

# View the confusion matrix and overall accuracy
print(conf_matrix)


#Assign 1 portion ----------

# Define the recipe for preprocessing
recipe <- recipe(class ~ ., data = train_data) %>%
  step_normalize(all_predictors())  # Normalize features

# Split the data into training and testing sets (optional)
set.seed(123)
data_split <- initial_split(train_data, prop = 0.8)
train_set <- training(data_split)
test_set <- testing(data_split)
# Specify the random forest model
rf_model <- rand_forest() %>%
  set_engine("ranger") %>%
  set_mode("classification")

# Create a workflow
workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(rf_model)
# Fit the model
rf_fit <- fit(workflow, data = train_set)

# Make predictions on the test set
predictions <- predict(rf_fit, test_set) %>%
  bind_cols(test_set) %>%
  rename(.pred_class = .pred_class)  # Rename prediction column if needed

# View predictions
print(predictions)

# Evaluate the model's performance
conf_matrix <- conf_mat(predictions, truth = class, estimate = .pred_class)

# View the confusion matrix and overall accuracy
print(conf_matrix)

#Assign 3 portion------------

# Transpose the top_5000_gene_data_vst matrix to match the number of samples
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <- top_5000_gene_data_vst_t[, 1:4591]

# Add "features." prefix to each gene column name
colnames(top_5000_gene_data_vst_t) <- paste0("features.", colnames(top_5000_gene_data_vst_t))

# Extract cluster labels
cluster_labels <- as.factor(gmm_results_vst_clean$classification)  # Convert to factor for classification

# Combine features and labels to create the training data
train_data <- data.frame(top_5000_gene_data_vst_t, cluster_label = cluster_labels)

# Split the data into training and testing sets
set.seed(123)
data_split <- initial_split(train_data, prop = 0.8)
train_set <- training(data_split)
test_set <- testing(data_split)

# Ensure `test_set` contains `cluster_label`
str(test_set)

# Train the model
rf_fit <- fit(workflow, data = train_set)

# Make predictions on the test set
predictions <- predict(rf_fit, test_set) %>%
  bind_cols(test_set %>% select(cluster_label))

# Evaluate the model using a confusion matrix
conf_matrix <- conf_mat(predictions, truth = cluster_label, estimate = .pred_class)
print(conf_matrix)


#Part 4 ------

# Load necessary libraries
library(caTools)  # For train/test splitting
library(randomForest)  # For Random Forest model
library(pROC)  # For AUC calculation



# List of cluster variables
clusters_list <- list(
  clusters_10 = clusters_10,
  clusters_100 = clusters_100,
  clusters_1000 = clusters_1000,
  clusters_4000 = clusters_4000
)

# Function to train and evaluate a Random Forest model
train_and_evaluate_rf <- function(clusters, gene_data) {
  
  data_for_model <- data.frame(gene_data, cluster = as.factor(clusters))
  set.seed(123)
  
  # Split the data into training and test sets
  split <- sample.split(data_for_model$cluster, SplitRatio = 0.7)
  train_data <- subset(data_for_model, split == TRUE)
  test_data <- subset(data_for_model, split == FALSE)
  
  # Train Random Forest model
  rf_model <- randomForest(cluster ~ ., data = train_data, ntree = 100)
  
  # Predict probabilities for the test set
  pred_probs <- predict(rf_model, newdata = test_data, type = "prob")
  
  # Calculate AUC for each class (assuming binary or multiclass ROC)
  auc_values <- sapply(levels(test_data$cluster), function(label) {
    roc(test_data$cluster == label, pred_probs[, label])$auc
  })
  
  # Return the mean AUC across all classes
  mean_auc <- mean(auc_values)
  
  return(mean_auc)
}

# Initialize a list to store AUC results for each gene set
auc_results <- list()

# Iterate over the different clusters and calculate AUC
for (i in 1:length(clusters_list)) {
  clusters <- clusters_list[[i]]
  
  # Use the function for the given clusters and gene data
  auc <- train_and_evaluate_rf(clusters, top_5000_gene_data_vst_t)
  auc_results[[names(clusters_list)[i]]] <- auc
}

# Convert results to a data frame for better readability
auc_results_df <- data.frame(
  gene_count = c(10, 100, 1000, 4000),
  auc = unlist(auc_results)
)

# Print the AUC results
print(auc_results_df)

      
#End Dylan section -------------------------------------------------------------

#-------------------------------------------------------------------------------
# Supervised Analysis (Naïve Bayes) - Victoria

if (!require("tidymodels")) {
  install.packages("tidymodels")
}
if (!require("e1071")) {
  install.packages("e1071")
}
if (!require("caret")) {
  install.packages("caret")
}
if (!require("pROC")) {
  install.packages("pROC")
}
if (!require("caTools")) {
  install.packages("caTools")
}

library(tidymodels)
library(e1071)
library(caret)
library(pROC)
library(caTools)

#d. Predict the two groups from Assignment 1----------
sample_grouping <- metadata$Diabetic_status
names(sample_grouping) <- metadata$refinebio_accession_code
sample_grouping <- factor(sample_grouping, levels = c("non-diabetic", "t2d"))

table(sample_grouping)

labels <- sample_grouping

# Prepare the data
data_for_model <- data.frame(top_5000_gene_data_vst_t, label = as.factor(labels))
head(data_for_model)

# Split into training and testing sets
set.seed(123)
split <- sample.split(data_for_model$label, SplitRatio = 0.7)
train_data <- subset(data_for_model, split == "TRUE")
test_data <- subset(data_for_model, split == "FALSE")

# Scale the features
train_data_scaled <- scale(train_data[, -ncol(train_data)])
test_data_scaled <- scale(test_data[, -ncol(test_data)])

# Convert back to data frames
train_data_scaled <- data.frame(train_data_scaled, label = train_data$label)
test_data_scaled <- data.frame(test_data_scaled, label = test_data$label)

# Train 
set.seed(120)  
classifier_cl <- naiveBayes(label ~ ., data = train_data_scaled)

# Predict on test data
y_pred <- predict(classifier_cl, newdata = test_data_scaled)

# Confusion Matrix
cm <- table(test_data_scaled$label, y_pred)
confusionMatrix(cm)

#e. Retrain to predict clusters from Assignment 3----------

# Using 100 gene clusters
labels <- as.factor(clusters_100)
data_for_model <- data.frame(top_5000_gene_data_vst_t, cluster = labels)

# Split the data
set.seed(123)
split <- sample.split(data_for_model$cluster, SplitRatio = 0.7)
train_data <- subset(data_for_model, split == TRUE)
test_data <- subset(data_for_model, split == FALSE)

# Scale the features
train_data_scaled <- scale(train_data[, -ncol(train_data)])
test_data_scaled <- scale(test_data[, -ncol(test_data)])

# Convert back to data frames
train_data_scaled <- data.frame(train_data_scaled, cluster = train_data$cluster)
test_data_scaled <- data.frame(test_data_scaled, cluster = test_data$cluster)

# Train Naive Bayes classifier on the scaled data
set.seed(120)
classifier_cl <- naiveBayes(cluster ~ ., data = train_data_scaled)

# Predict on the test data
y_pred <- predict(classifier_cl, newdata = test_data_scaled)

cm <- table(test_data_scaled$cluster, y_pred)
confusionMatrix(cm)

#Part 4-----------------------------
# List of clusters for different gene subsets
clusters_list <- list(
  clusters_10 = results_10$classification,
  clusters_100 = results_100$classification,
  clusters_1000 = results_1000$classification,
  clusters_4000 = results_4000$classification
)

# Function to train and evaluate model with specific clusters
train_and_evaluate_auc <- function(clusters, gene_data) {
  
  # Combine the clusters with gene data
  data_for_model <- data.frame(gene_data, cluster = as.factor(clusters))
  
  # Split data 
  set.seed(123)
  split <- sample.split(data_for_model$cluster, SplitRatio = 0.7)
  train_data <- subset(data_for_model, split == TRUE)
  test_data <- subset(data_for_model, split == FALSE)
  
  # Standardize the training and test sets
  train_data_scaled <- scale(train_data[, -ncol(train_data)]) 
  test_data_scaled <- scale(test_data[, -ncol(test_data)])
  
  # Convert scaled data back to data frames and add the cluster column
  train_data_scaled <- data.frame(train_data_scaled, cluster = train_data$cluster)
  test_data_scaled <- data.frame(test_data_scaled, cluster = test_data$cluster)
  
  # Train Naive Bayes model
  set.seed(120)
  classifier_cl <- naiveBayes(cluster ~ ., data = train_data_scaled)
  
  # Predict 
  y_pred <- predict(classifier_cl, newdata = test_data_scaled, type = "raw")
  
  # Calculate AUC
  auc <- multiclass.roc(test_data_scaled$cluster, y_pred)
  
  # Return AUC
  return(list(
    auc = auc$auc
  ))
}

# Initialize AUC list
auc_results <- list()

# Iterate over each cluster subset
for (cluster_name in names(clusters_list)) {
  clusters <- clusters_list[[cluster_name]]
  result <- train_and_evaluate_auc(clusters, top_5000_gene_data_vst_t)
  auc_results[[cluster_name]] <- result$auc
}

# Convert AUC results to a data frame
auc_results_df <- data.frame(
  gene_count = c(10, 100, 1000, 4000),
  auc = unlist(auc_results)
)

print(auc_results_df)

#End Victoria's Section
#-------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Kevin - Supervised Analysis - Logistic Regression

if (!requireNamespace("tidymodels", quietly = TRUE)) {
  install.packages("tidymodels")
}
if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret")
}
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}
if (!requireNamespace("rsample", quietly = TRUE)) {
  install.packages("rsample")
}

# Load necessary libraries
library(DESeq2)
library(rsample)
library(tidymodels)
library(dplyr) 
library(kernlab)
library(glmnet)
library(caret)

# === T2D vs Non-T2D Classification === #

vsd <- varianceStabilizingTransformation(deseq_object)
gene_variances_vst <- apply(assay(vsd), 1, var)
top_5000_genes_vst <- order(gene_variances_vst, decreasing = TRUE)[1:5000]
top_5000_gene_data_vst <- assay(vsd)[top_5000_genes_vst, ]

# Transpose the data
top_5000_gene_data_vst_t <- t(top_5000_gene_data_vst)
top_5000_gene_data_vst_t <- top_5000_gene_data_vst_t[, 1:4591]

# Create data frame with T2D and Non-T2D labels
df <- as.data.frame(top_5000_gene_data_vst_t)
n_non_diabetic <- 547
n_t2d <- 1053
class_labels <- factor(c(rep("non-diabetic", n_non_diabetic), rep("t2d", n_t2d)))
df$class <- class_labels

# Split the data into training and testing sets
set.seed(123)
train_index <- createDataPartition(df$class, p = 0.7, list = FALSE)
train_data <- df[train_index, ]
test_data <- df[-train_index, ]


# Define logistic regression model for T2D vs Non-T2D classification
logistic_spec <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")

# Set up recipe and workflow
recipe <- recipe(class ~ ., data = train_data) %>%
  step_normalize(all_predictors())

logistic_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(logistic_spec)

# Train the model
logistic_fit <- fit(logistic_workflow, data = train_data)

# Make predictions on the test set
log_predictions <- predict(logistic_fit, test_data) %>%
  bind_cols(test_data)

# Evaluate the model's performance with a confusion matrix
log_conf_matrix <- conf_mat(log_predictions, truth = class, estimate = .pred_class)

# Print the confusion matrix
print("T2D vs Non-T2D Confusion Matrix")
print(log_conf_matrix)

                                 
# === Spectral Clustering-Based Classification === #

# Perform spectral clustering on the gene data
num_clusters <- 3  # Adjust the number of clusters as needed
spectral_clustering_result <- specc(top_5000_gene_data_vst_t, centers = num_clusters)
cluster_labels <- as.factor(spectral_clustering_result@.Data)

# Add cluster labels to the data frame
df$cluster_label <- cluster_labels

# Subset only numeric data for PCA
numeric_data <- df[, sapply(df, is.numeric)]  # Select only numeric columns

# Run PCA on numeric data (exclude the class and cluster label)
pca_res <- prcomp(numeric_data, scale. = TRUE)

# Get the top 50 principal components
pca_data <- as.data.frame(pca_res$x[, 1:50])
pca_data$cluster_label <- df$cluster_label  # Add cluster labels to PCA data

# Split the data into training and testing sets
set.seed(123)
data_split <- initial_split(pca_data, prop = 0.7, strata = cluster_label)
train_data <- training(data_split)
test_data <- testing(data_split)

# Define the logistic regression model with regularization (Lasso)
logistic_spec <- logistic_reg(penalty = 0.1, mixture = 0) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

train_one_vs_all <- function(train_data, cluster_id) {
  # Create binary outcome: Cluster cluster_id = 1, all other clusters = 0
  train_data$binary_cluster <- factor(ifelse(train_data$cluster_label == cluster_id, 1, 0))  # Convert to factor
  
  # Remove non-numeric columns for normalization (cluster_label and binary_cluster will be excluded later)
  numeric_train_data <- train_data %>% select(-cluster_label)
  
  recipe <- recipe(binary_cluster ~ ., data = numeric_train_data) %>%
    step_normalize(all_numeric_predictors())  # Normalize only numeric predictors
  
  # Create a workflow and add the logistic regression model
  logistic_workflow <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(logistic_spec)
  
  # Train the logistic regression model
  logistic_fit <- fit(logistic_workflow, data = train_data)
  
  return(logistic_fit)
}

# Train one-vs-all models for each cluster
models <- list()  # Initialize the models list

# Ensure that the factor levels for cluster_label are consistent across train and test datasets
train_data$cluster_label <- factor(train_data$cluster_label, levels = c("1", "2", "3"))
test_data$cluster_label <- factor(test_data$cluster_label, levels = c("1", "2", "3"))

for (cluster_id in levels(df$cluster_label)) {
  
  # Train the one-vs-all model for the current cluster
  model <- train_one_vs_all(train_data, cluster_id)

  models[[as.character(cluster_id)]] <- model
  
  # Check if model exists and proceed with prediction
  if (is.null(models[[as.character(cluster_id)]])) {
    print(paste("No model found for cluster", cluster_id))
    next  # Skip to the next cluster if no model is found
  }
  
  # Make predictions for the current model (cluster)
  predictions <- predict(models[[as.character(cluster_id)]], new_data = test_data) %>%
    bind_cols(test_data)
  
  # Ensure predicted_class has the same levels as cluster_label
  predictions$predicted_class <- factor(predictions$.pred_class, levels = c("1", "2", "3"))
  
  predictions$cluster_label <- factor(predictions$cluster_label, levels = levels(test_data$cluster_label))
  
  # Evaluate the model's performance with a confusion matrix
  conf_matrix <- conf_mat(predictions, truth = cluster_label, estimate = predicted_class)
  
  print(paste("Confusion Matrix for Cluster", cluster_id))
  print(conf_matrix)
}

# Select the top 10 genes by variance
top_10_genes <- order(gene_variances_vst, decreasing = TRUE)[1:10]
gene_data_subset_10 <- assay(vsd)[top_10_genes, ]

# Transpose the data for clustering
gene_data_subset_10_t <- t(gene_data_subset_10)

# Perform spectral clustering
num_clusters <- 3  # Adjust if needed
spectral_clustering_result_10 <- specc(gene_data_subset_10_t, centers = num_clusters)
cluster_labels_10 <- as.factor(spectral_clustering_result_10@.Data)

# Create data frame and add cluster labels
df_10 <- as.data.frame(gene_data_subset_10_t)
df_10$cluster_label <- cluster_labels_10

# Remove the cluster_label column and ensure only numeric columns are used for PCA
numeric_data_10 <- df_10[, sapply(df_10, is.numeric)]  # Only numeric columns
pca_res_10 <- prcomp(numeric_data_10, scale. = TRUE)  # Apply PCA on the numeric data

# Check the number of principal components available
num_pcs <- ncol(pca_res_10$x)  # Number of available principal components (10 in this case)

# Use all available principal components (10 in this case)
pca_data_10 <- as.data.frame(pca_res_10$x[, 1:num_pcs])  # Use the available principal components
pca_data_10$cluster_label <- df_10$cluster_label  # Add the cluster labels

# Split data into training and testing sets
set.seed(123)
data_split_10 <- initial_split(pca_data_10, prop = 0.7, strata = cluster_label)
train_data_10 <- training(data_split_10)
test_data_10 <- testing(data_split_10)

# Train one-vs-all logistic regression models for each cluster
models_10 <- list()
for (cluster_id in levels(train_data_10$cluster_label)) {
  train_data_10$binary_cluster <- factor(ifelse(train_data_10$cluster_label == cluster_id, 1, 0))
  
  # Set up recipe and model workflow
  recipe_10 <- recipe(binary_cluster ~ ., data = train_data_10 %>% select(-cluster_label)) %>%
    step_normalize(all_numeric_predictors())
  logistic_workflow_10 <- workflow() %>%
    add_recipe(recipe_10) %>%
    add_model(logistic_spec)
  
  # Train model and store
  models_10[[as.character(cluster_id)]] <- fit(logistic_workflow_10, data = train_data_10)
  
  # Predictions and confusion matrix for each cluster
  predictions_10 <- predict(models_10[[as.character(cluster_id)]], test_data_10) %>%
    bind_cols(test_data_10)
  predictions_10$predicted_class <- factor(predictions_10$.pred_class, levels = levels(test_data_10$cluster_label))
  
  conf_matrix_10 <- conf_mat(predictions_10, truth = cluster_label, estimate = predicted_class)
  print(paste("Confusion Matrix for Cluster", cluster_id, "using Top 10 Genes"))
  print(conf_matrix_10)
}

# Select the top 100 genes by variance
top_100_genes <- order(gene_variances_vst, decreasing = TRUE)[1:100]
gene_data_subset_100 <- assay(vsd)[top_100_genes, ]

# Transpose data
gene_data_subset_100_t <- t(gene_data_subset_100)

# Perform spectral clustering
spectral_clustering_result_100 <- specc(gene_data_subset_100_t, centers = num_clusters)
cluster_labels_100 <- as.factor(spectral_clustering_result_100@.Data)

# Create data frame and add cluster labels
df_100 <- as.data.frame(gene_data_subset_100_t)
df_100$cluster_label <- cluster_labels_100

# PCA and dimensionality reduction
pca_res_100 <- prcomp(df_100[, -ncol(df_100)], scale. = TRUE)
pca_data_100 <- as.data.frame(pca_res_100$x[, 1:50])
pca_data_100$cluster_label <- df_100$cluster_label

# Split data
data_split_100 <- initial_split(pca_data_100, prop = 0.7, strata = cluster_label)
train_data_100 <- training(data_split_100)
test_data_100 <- testing(data_split_100)

# One-vs-all logistic regression models
models_100 <- list()
for (cluster_id in levels(train_data_100$cluster_label)) {
  train_data_100$binary_cluster <- factor(ifelse(train_data_100$cluster_label == cluster_id, 1, 0))
  
  # Set up recipe and model
  recipe_100 <- recipe(binary_cluster ~ ., data = train_data_100 %>% select(-cluster_label)) %>%
    step_normalize(all_numeric_predictors())
  logistic_workflow_100 <- workflow() %>%
    add_recipe(recipe_100) %>%
    add_model(logistic_spec)
  
  # Train model
  models_100[[as.character(cluster_id)]] <- fit(logistic_workflow_100, data = train_data_100)
  
  # Predictions and evaluation
  predictions_100 <- predict(models_100[[as.character(cluster_id)]], test_data_100) %>%
    bind_cols(test_data_100)
  predictions_100$predicted_class <- factor(predictions_100$.pred_class, levels = levels(test_data_100$cluster_label))
  
  conf_matrix_100 <- conf_mat(predictions_100, truth = cluster_label, estimate = predicted_class)
  print(paste("Confusion Matrix for Cluster", cluster_id, "using Top 100 Genes"))
  print(conf_matrix_100)
}

# Select the top 1000 genes by variance
top_1000_genes <- order(gene_variances_vst, decreasing = TRUE)[1:1000]
gene_data_subset_1000 <- assay(vsd)[top_1000_genes, ]

# Transpose data
gene_data_subset_1000_t <- t(gene_data_subset_1000)

# Perform spectral clustering
spectral_clustering_result_1000 <- specc(gene_data_subset_1000_t, centers = num_clusters)
cluster_labels_1000 <- as.factor(spectral_clustering_result_1000@.Data)

# Create data frame and add cluster labels
df_1000 <- as.data.frame(gene_data_subset_1000_t)
df_1000$cluster_label <- cluster_labels_1000

# PCA and dimensionality reduction
pca_res_1000 <- prcomp(df_1000[, -ncol(df_1000)], scale. = TRUE)
pca_data_1000 <- as.data.frame(pca_res_1000$x[, 1:50])
pca_data_1000$cluster_label <- df_1000$cluster_label

# Split data
data_split_1000 <- initial_split(pca_data_1000, prop = 0.7, strata = cluster_label)
train_data_1000 <- training(data_split_1000)
test_data_1000 <- testing(data_split_1000)

# One-vs-all logistic regression models
models_1000 <- list()
for (cluster_id in levels(train_data_1000$cluster_label)) {
  train_data_1000$binary_cluster <- factor(ifelse(train_data_1000$cluster_label == cluster_id, 1, 0))
  
  # Set up recipe and model
  recipe_1000 <- recipe(binary_cluster ~ ., data = train_data_1000 %>% select(-cluster_label)) %>%
    step_normalize(all_numeric_predictors())
  logistic_workflow_1000 <- workflow() %>%
    add_recipe(recipe_1000) %>%
    add_model(logistic_spec)
  
  # Train model
  models_1000[[as.character(cluster_id)]] <- fit(logistic_workflow_1000, data = train_data_1000)
  
  # Predictions and evaluation
  predictions_1000 <- predict(models_1000[[as.character(cluster_id)]], test_data_1000) %>%
    bind_cols(test_data_1000)
  predictions_1000$predicted_class <- factor(predictions_1000$.pred_class, levels = levels(test_data_1000$cluster_label))
  
  conf_matrix_1000 <- conf_mat(predictions_1000, truth = cluster_label, estimate = predicted_class)
  print(paste("Confusion Matrix for Cluster", cluster_id, "using Top 1000 Genes"))
  print(conf_matrix_1000)
}

# ----- 4 --------
install.packages("pRoc")
# Load pROC library
library(pROC)

predictions_10_prob <- predict(models_10[[as.character(cluster_id)]], test_data_10, type = "prob")
# Match the row names of test_data_10 with the row names of df_10 to get the correct cluster_label
test_data_10$cluster_label <- df_10$cluster_label[match(rownames(test_data_10), rownames(df_10))]

# Check the distribution of cluster labels in the test set
table(test_data_10$cluster_label)

# Rename columns of predicted probabilities to match the levels of cluster_label
colnames(predictions_10_prob) <- c("1", "2", "3")

# Ensure that predictions_10_prob is a numeric matrix
predictions_10_prob_matrix <- as.matrix(predictions_10_prob)

# Recalculate the multiclass ROC
roc_10 <- multiclass.roc(test_data_10$cluster_label, predictions_10_prob_matrix)

# Calculate and print the AUC
auc_10 <- auc(roc_10)

predictions_100_prob <- predict(models_100[[as.character(cluster_id)]], test_data_100, type = "prob")
# Match the row names of test_data_10 with the row names of df_10 to get the correct cluster_label
test_data_100$cluster_label <- df_100$cluster_label[match(rownames(test_data_100), rownames(df_100))]

# Check the distribution of cluster labels in the test set
table(test_data_100$cluster_label)

# Rename columns of predicted probabilities to match the levels of cluster_label
colnames(predictions_100_prob) <- c("1", "2", "3")

# Ensure that predictions_10_prob is a numeric matrix
predictions_100_prob_matrix <- as.matrix(predictions_100_prob)

# Recalculate the multiclass ROC
roc_100 <- multiclass.roc(test_data_100$cluster_label, predictions_100_prob_matrix)

# Calculate and print the AUC
auc_100 <- auc(roc_100)

predictions_1000_prob <- predict(models_1000[[as.character(cluster_id)]], test_data_1000, type = "prob")
# Match the row names of test_data_10 with the row names of df_10 to get the correct cluster_label
test_data_1000$cluster_label <- df_1000$cluster_label[match(rownames(test_data_1000), rownames(df_1000))]

# Check the distribution of cluster labels in the test set
table(test_data_1000$cluster_label)

# Rename columns of predicted probabilities to match the levels of cluster_label
colnames(predictions_1000_prob) <- c("1", "2", "3")

# Ensure that predictions_10_prob is a numeric matrix
predictions_1000_prob_matrix <- as.matrix(predictions_1000_prob)

# Recalculate the multiclass ROC
roc_1000 <- multiclass.roc(test_data_1000$cluster_label, predictions_1000_prob_matrix)

# Calculate and print the AUC
auc_1000 <- auc(roc_1000)
print(paste("AUC for Top 10 Genes:", auc_10))
print(paste("AUC for Top 100 Genes:", auc_100))
print(paste("AUC for Top 1000 Genes:", auc_1000))
                                 
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
#Adam's section

if (!require("tidymodels")) {
  install.packages("tidymodels")
}
if (!require("e1071")) {
  install.packages("e1071")
}
if (!require("caret")) {
  install.packages("caret")
}
if (!require("pROC")) {
  install.packages("pROC")
}
if (!require("caTools")) {
  install.packages("caTools")
}

library(tidymodels)
library(e1071)
library(caret)
library(pROC)
library(caTools)

# PREDICTION ON DIABETIC VS NON DIABETIC (USING KNN)

sample_grouping_a <- metadata$Diabetic_status
names(sample_grouping_a) <- metadata$refinebio_accession_code
sample_grouping_a <- factor(sample_grouping_a, levels = c("non-diabetic", "t2d"))

table(sample_grouping_a)

labels_a <- sample_grouping_a

# Prepare the data
data_for_model_a <- data.frame(top_5000_gene_data_vst_t, label = as.factor(labels_a))
head(data_for_model_a)

# Split data
set.seed(123)
split_a <- sample.split(data_for_model_a$label, SplitRatio = 0.7)
train_data_a <- subset(data_for_model_a, split_a == "TRUE")
test_data_a <- subset(data_for_model_a, split_a == "FALSE")

# Scale data
train_data_scaled_a <- scale(train_data_a[, -ncol(train_data_a)])
test_data_scaled_a <- scale(test_data_a[, -ncol(test_data_a)])

train_data_scaled_a <- data.frame(train_data_scaled_a, label = train_data_a$label)
test_data_scaled_a <- data.frame(test_data_scaled_a, label = test_data_a$label)

library(caret)

set.seed(120)
ctrl_a <- trainControl(method = "cv", number = 10)
knn_model_a <- train(label ~ ., data = train_data_scaled_a, method = "knn",
                     tuneLength = 10, trControl = ctrl_a)

print(knn_model_a)

y_pred_knn_a <- predict(knn_model_a, newdata = test_data_scaled_a)
y_pred_knn_numeric <- ifelse(y_pred_knn_a == "non-diabetic", 0, 1)
cm_a <- confusionMatrix(test_data_scaled_a$label, y_pred_knn_a)
print(cm_a)

# PREDICTION ON CLUSTERS USING KN

#  clusters_100_a is the clustering result from Assignment 3

levels(cluster_labels) <- paste0("Class", levels(cluster_labels))
library(caret)
library(caTools)
library(pROC)

cluster_labels <- as.factor(clusters_100_a)

train_knn_with_genes_cluster <- function(num_genes, cluster_labels) {

  gene_data_subset <- top_5000_gene_data_vst_t[, 1:num_genes]
  
  data_for_model_clusters <- data.frame(gene_data_subset, cluster = cluster_labels)
  
  set.seed(123)
  split_clusters <- sample.split(data_for_model_clusters$cluster, SplitRatio = 0.7)
  train_data_clusters <- subset(data_for_model_clusters, split_clusters == TRUE)
  test_data_clusters <- subset(data_for_model_clusters, split_clusters == FALSE)
  
  train_data_clusters_scaled <- scale(train_data_clusters[, -ncol(train_data_clusters)])
  test_data_clusters_scaled <- scale(test_data_clusters[, -ncol(test_data_clusters)])

  train_data_clusters_scaled <- data.frame(train_data_clusters_scaled, cluster = train_data_clusters$cluster)
  test_data_clusters_scaled <- data.frame(test_data_clusters_scaled, cluster = test_data_clusters$cluster)
  
  set.seed(120)
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary)
  
  knn_model_clusters <- train(cluster ~ ., data = train_data_clusters_scaled, method = "knn",
                              tuneLength = 10, trControl = ctrl)
  
  y_pred_clusters <- predict(knn_model_clusters, newdata = test_data_clusters_scaled)
  cm <- confusionMatrix(test_data_clusters_scaled$cluster, y_pred_clusters)

  y_pred_probs <- predict(knn_model_clusters, newdata = test_data_clusters_scaled, type = "prob")
  multiclass_roc <- multiclass.roc(test_data_clusters_scaled$cluster, y_pred_probs)
  auc <- auc(multiclass_roc)

  print(cm)
  cat("AUC for", num_genes, "genes:", auc, "\n")
  
  return(list(accuracy = cm$overall["Accuracy"], auc = auc))
}

gene_subset_sizes <- c(10, 100, 1000, 4000)


results_list <- list()

for (num_genes in gene_subset_sizes) {
  cat("\nTraining with", num_genes, "genes to predict clusters_100_a:\n")
  result <- train_knn_with_genes_cluster(num_genes, cluster_labels)
  results_list[[paste0("genes_", num_genes)]] <- result
}

results_list

# KNN AUC RESULTS
train_knn_with_auc <- function(num_genes, cluster_labels) {
  gene_data_subset <- top_5000_gene_data_vst_t[, 1:num_genes]

  data_for_model_clusters <- data.frame(gene_data_subset, cluster = cluster_labels)
  

  set.seed(123)
  split_clusters <- sample.split(data_for_model_clusters$cluster, SplitRatio = 0.7)
  train_data_clusters <- subset(data_for_model_clusters, split_clusters == TRUE)
  test_data_clusters <- subset(data_for_model_clusters, split_clusters == FALSE)
  
  # scale features
  train_data_clusters_scaled <- scale(train_data_clusters[, -ncol(train_data_clusters)])
  test_data_clusters_scaled <- scale(test_data_clusters[, -ncol(test_data_clusters)])
  
  train_data_clusters_scaled <- data.frame(train_data_clusters_scaled, cluster = train_data_clusters$cluster)
  test_data_clusters_scaled <- data.frame(test_data_clusters_scaled, cluster = test_data_clusters$cluster)
  
  #  KNN Classifier
  set.seed(120)
  ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary)
  
  knn_model_clusters <- train(cluster ~ ., data = train_data_clusters_scaled, method = "knn",
                              tuneLength = 10, trControl = ctrl)
  
 
  y_pred_probs <- predict(knn_model_clusters, newdata = test_data_clusters_scaled, type = "prob")
  
  
  multiclass_roc <- multiclass.roc(test_data_clusters_scaled$cluster, y_pred_probs)
  auc <- auc(multiclass_roc)
  
  cat("AUC for", num_genes, "genes:", auc, "\n")
  return(list(auc = auc))
}

gene_subset_sizes <- c(10, 100, 1000, 4000)

auc_results <- list()

for (num_genes in gene_subset_sizes) {
  cat("\nTraining with", num_genes, "genes to predict clusters_100_a:\n")
  result <- train_knn_with_auc(num_genes, cluster_labels)
  auc_results[[paste0("genes_", num_genes)]] <- result
}

auc_results
                                 

#END OF ADAM'S SECTION: KNN MODEL


#-----------------------
#Begin Andrew's sections

# Load necessary packages (SVM implementation requires 'e1071' for svm function)
if (!require("tidymodels")) {
  install.packages("tidymodels")
}
if (!require("e1071")) {
  install.packages("e1071")
}
if (!require("caret")) {
  install.packages("caret")
}
if (!require("pROC")) {
  install.packages("pROC")
}
if (!require("caTools")) {
  install.packages("caTools")
}

library(tidymodels)
library(e1071)
library(caret)
library(pROC)
library(caTools)

# d. Predict the two groups from Assignment 1 ----------
# Prepare the grouping variable for diabetic status
sample_grouping  <- metadata$Diabetic_status
names(sample_grouping ) <- metadata$refinebio_accession_code
sample_grouping  <- factor(sample_grouping , levels = c("non-diabetic", "t2d"))

# Prepare the data
svm_data <- data.frame(top_5000_gene_data_vst_t, label = as.factor(sample_grouping ))

# Split into training and testing sets
set.seed(123)
data_split <- sample.split(svm_data$label, SplitRatio = 0.7)
train_set <- subset(svm_data, data_split == "TRUE")
test_set <- subset(svm_data, data_split == "FALSE")

# Scale the features
train_scaled <- scale(train_set[, -ncol(train_set)])
test_scaled <- scale(test_set[, -ncol(test_set)])

# Convert back to data frames
train_scaled <- data.frame(train_scaled, label = train_set$label)
test_scaled <- data.frame(test_scaled, label = test_set$label)

# Train SVM classifier
set.seed(120)
svm_classifier <- svm(label ~ ., data = train_scaled, kernel = "linear", cost = 1, scale = FALSE)

# Predict on test data
predictions <- predict(svm_classifier, newdata = test_scaled)

# Confusion Matrix
conf_matrix <- table(test_scaled$label, predictions)
confusionMatrix(conf_matrix)

# e. Retrain to predict clusters from Assignment 3 ----------
# Prepare labels based on clustering results from Assignment 3
cluster_labels <- as.factor(clusters_100)
cluster_data <- data.frame(top_5000_gene_data_vst_t, cluster = cluster_labels)

# Split the data
set.seed(123)
cluster_split <- sample.split(cluster_data$cluster, SplitRatio = 0.7)
train_cluster <- subset(cluster_data, cluster_split == TRUE)
test_cluster <- subset(cluster_data, cluster_split == FALSE)

# Scale the features
train_cluster_scaled <- scale(train_cluster[, -ncol(train_cluster)])
test_cluster_scaled <- scale(test_cluster[, -ncol(test_cluster)])

# Convert back to data frames
train_cluster_scaled <- data.frame(train_cluster_scaled, cluster = train_cluster$cluster)
test_cluster_scaled <- data.frame(test_cluster_scaled, cluster = test_cluster$cluster)

# Train SVM classifier on clusters
set.seed(120)
svm_cluster_classifier <- svm(cluster ~ ., data = train_cluster_scaled, kernel = "linear", cost = 1, scale = FALSE)

# Predict on test data
cluster_predictions <- predict(svm_cluster_classifier, newdata = test_cluster_scaled)

# Confusion Matrix for clusters
cluster_conf_matrix <- table(test_cluster_scaled$cluster, cluster_predictions)
confusionMatrix(cluster_conf_matrix)

# End of Andrew's SVM Section




#Begin P4 for SVM

calculate_auc <- function(gene_data, sample_grouping, top_genes = 10) {
  
  # Prepare the data
  svm_data <- data.frame(gene_data, label = as.factor(sample_grouping))
  
  # Split into training and testing sets
  set.seed(123)
  data_split <- sample.split(svm_data$label, SplitRatio = 0.7)
  train_set <- subset(svm_data, data_split == "TRUE")
  test_set <- subset(svm_data, data_split == "FALSE")
  
  # Scale the features
  train_scaled <- scale(train_set[, -ncol(train_set)])
  test_scaled <- scale(test_set[, -ncol(test_set)])
  
  # Convert back to data frames
  train_scaled <- data.frame(train_scaled, label = train_set$label)
  test_scaled <- data.frame(test_scaled, label = test_set$label)
  
  # Train SVM classifier
  set.seed(120)
  svm_classifier <- svm(label ~ ., data = train_scaled, kernel = "linear", cost = 1, scale = FALSE)
  
  # Predict on test data
  predictions <- predict(svm_classifier, newdata = test_scaled, type = "prob")
  
  # Prepare the test labels
  test_set$label <- as.factor(test_set$label)
  
  # Check the distribution of labels
  table(test_set$label)
  
  # Ensure that predictions are a numeric matrix
  predictions_matrix <- as.matrix(predictions)
  
  # Recalculate multiclass ROC
  roc_result <- multiclass.roc(test_set$label, predictions_matrix)
  
  # Calculate and return the AUC
  auc_value <- auc(roc_result)
  return(auc_value)
}


top_10_gene_data_vst_t <- top_5000_gene_data_vst_t[, order(-apply(top_5000_gene_data_vst_t, 2, var))[1:10]]
top_100_gene_data_vst_t <- top_5000_gene_data_vst_t[, order(-apply(top_5000_gene_data_vst_t, 2, var))[1:100]]
top_1000_gene_data_vst_t <- top_5000_gene_data_vst_t[, order(-apply(top_5000_gene_data_vst_t, 2, var))[1:1000]]


# Calculate AUC for Top 10, 100, and 1000 genes
auc_10 <- calculate_auc(top_10_gene_data_vst_t, metadata$Diabetic_status, top_genes = 10)
auc_100 <- calculate_auc(top_100_gene_data_vst_t, metadata$Diabetic_status, top_genes = 100)
auc_1000 <- calculate_auc(top_1000_gene_data_vst_t, metadata$Diabetic_status, top_genes = 1000)

# Print AUC values
print(paste("AUC for Top 10 Genes:", auc_10))
print(paste("AUC for Top 100 Genes:", auc_100))
print(paste("AUC for Top 1000 Genes:", auc_1000))
#End Andrew's SVM p4
#-----------------------------------



#-----------------------------------------------

#----------------------------------------------
#ROC curve portion
# Required libraries
if (!require("tidymodels")) install.packages("tidymodels")
if (!require("e1071")) install.packages("e1071")
if (!require("caret")) install.packages("caret")
if (!require("pROC")) install.packages("pROC")
if (!require("caTools")) install.packages("caTools")

library(tidymodels)
library(e1071)
library(caret)
library(pROC)
library(caTools)


log_predictions_a <- log_predictions[1:53, ]  # Rows 1 to 53 for all columns
predictions_a <- predictions[1:53, ]  # Rows 1 to 53 for all columns
y_pred_a <- y_pred[1:53]  # Rows 1 to 53 for all columns
y_pred_clusters_a_a <- y_pred_clusters_a[1:53]  # Rows 1 to 53 for all columns
cluster_predictions_a <- cluster_predictions[1:53]  # Rows 1 to 53 for all columns

# Now create the data frame
prediction_matrix <- data.frame(
  logistic = log_predictions_a,
  random_forest = predictions_a,
  naive_bayes = y_pred_a,
  knn = y_pred_clusters_a_a,
  svm = cluster_predictions_a
)

# Step 1: Extract the relevant columns from prediction_matrix based on the structure you provided

# Assuming we have the data frame 'prediction_matrix' loaded
# Extracting the first set of columns for logistic predictions (class + gene predictions)
logistic_predictions <- prediction_matrix[, grep("logistic\\.", colnames(prediction_matrix))]
logistic_class <- prediction_matrix$logistic.class

# Random forest: Class prediction and cluster label
random_forest_cluster <- prediction_matrix$random_forest.cluster_label

# Extracting principal component columns for random forest
random_forest_pcs <- prediction_matrix[, grep("random_forest.PC", colnames(prediction_matrix))]

# Naive Bayes, KNN, and SVM predictions (assuming the order is correct in prediction_matrix)
naive_bayes_pred <- prediction_matrix$naive_bayes
knn_pred <- prediction_matrix$knn
svm_pred <- prediction_matrix$svm

# Step 2: Create a consolidated prediction matrix for class labels (excluding gene names and PCs)
prediction_matrix_cleaned <- data.frame(
  logistic = logistic_class,
  random_forest_cluster = random_forest_cluster,
  naive_bayes = naive_bayes_pred,
  knn = knn_pred,
  svm = svm_pred
)





# Step 2a: Count Model Agreement on Each Class Label for Each Sample
class_counts <- prediction_matrix_cleaned %>%
  mutate(sample_id = row_number()) %>%  # Create sample_id based on row number
  pivot_longer(cols = c("logistic", "random_forest_cluster", "naive_bayes", "knn", "svm"),
               names_to = "model", values_to = "prediction") %>%
  group_by(sample_id, prediction) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = prediction, values_from = count, values_fill = list(count = 0))

# Ensure `count` column is numeric
class_counts <- class_counts %>%
  mutate(across(where(is.character), as.numeric, .names = "converted_{.col}"))

# Step 2b: Count Stability of Cluster Predictions (Models Predicting Same Cluster)
cluster_counts <- prediction_matrix_cleaned %>%
  mutate(sample_id = row_number()) %>%  # Create sample_id based on row number
  rowwise() %>%
  mutate(
    models_same_cluster = max(table(c(logistic, random_forest_cluster, naive_bayes, knn, svm))) / 5  # Proportion of agreement
  ) %>%
  ungroup() %>%
  select(sample_id, models_same_cluster)

# Step 3: Calculate Correlation Between Class Label Stability and Cluster Stability
# First, calculate class agreement as the proportion of models predicting the same class
class_agreement <- class_counts %>%
  mutate(sample_id = row_number()) %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(starts_with("non-diabetic"):starts_with("5"))),
         class_agreement = ifelse(total_count > 0, max(c_across(starts_with("non-diabetic"):starts_with("5"))) / total_count, 0)) %>%
  ungroup() %>%
  select(sample_id, class_agreement)
# Merge class agreement and cluster stability data
stability_df <- left_join(class_agreement, cluster_counts, by = "sample_id")

# Perform Spearman correlation test
cor_test <- cor.test(stability_df$class_agreement, stability_df$models_same_cluster, method = "spearman")

# Adjust p-values (example here uses Bonferroni correction)
adjusted_p_values <- p.adjust(cor_test$p.value, method = "bonferroni")

# Output results
print(class_counts)           # Number of models predicting each class label per sample
print(cluster_counts)         # Proportion of models agreeing on the same cluster per sample
print(cor_test)               # Correlation test result
print(adjusted_p_values)      # Adjusted p-values if multiple tests
#End ROC analysis
#----------------------------------------------------------------
#HEATMAP SECTION 5 
library(ComplexHeatmap)

combined_matrix <- cbind(
  KNN = y_pred_knn_a,             # Replace knn_predictions with your actual KNN data frame
  NaiveBayes = y_pred,       # Replace nb_predictions with actual Naive Bayes data
  LogReg = log_predictions,      # Replace log_reg_predictions with Logistic Regression data
  RandomForest = rf_predictions,     # Replace rf_predictions with Random Forest data
)

combined_matrix <- as.matrix(combined_matrix)


heatmap_colors <- colorRamp2(c(0, 1), c("blue", "red"))

Heatmap(
  combined_matrix,
  name = "Prediction",
  col = heatmap_colors,
  row_names_gp = gpar(fontsize = 8),   # Adjust font size for readability
  column_names_gp = gpar(fontsize = 8),
  show_row_names = TRUE,               # Show sample names or IDs
  show_column_names = TRUE,            # Show gene names or IDs
  cluster_rows = TRUE,                 # Cluster samples based on prediction similarity
  cluster_columns = TRUE,              # Cluster genes if needed
  row_title = "Samples",
  column_title = "Genes by Model",
  top_annotation = HeatmapAnnotation(
    Model = c(rep("KNN", ncol(y_pred_knn_a)),
              rep("Naive Bayes", ncol(y_pred)),
              rep("Logistic Regression", ncol(log_predictions)),
              rep("Random Forest", ncol(predictions)),
    col = list(Model = c("KNN" = "darkgreen", "Naive Bayes" = "purple", 
                         "Logistic Regression" = "orange", 
                         "Random Forest" = "blue"))
  )
)


                                 
