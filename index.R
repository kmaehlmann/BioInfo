# Load Libraries
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(org.Hs.eg.db)
library(magrittr)
library(dplyr)
library(tidyr)

# Load Datasets
metadata <- readr::read_tsv("metadata_SRP075377.tsv")
expression_df <- readr::read_tsv("SRP075377.tsv") %>%
  tibble::column_to_rownames("Gene")
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated HUGO symbols
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
