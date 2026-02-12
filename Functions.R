#library(hrbrthemes)
library(viridis)
library(umap)
library(ggplot2)
library(dplyr)
library(MASS)
library(rmarkdown)
library(httr)
library(jsonlite)
library(tidyr)
library(plotly)
search_gwasrapidd <- function(efo_id) { #perform GWAS search trait ->association
  library(gwasrapidd)
  library(dplyr)
  
  # Step 1: Retrieve associations for the given EFO ID
  my_associations <- get_associations(efo_id = efo_id)
  
  # Step 2: Check the number of associations retrieved
  cat("Total associations retrieved:", gwasrapidd::n(my_associations), "\n")
  
  # Step 3: Convert associations and risk genes slots to data frames
  associations_df <- as.data.frame(my_associations@associations)
  risk_genes_df <- as.data.frame(my_associations@ensembl_ids)
  
  # Step 4: Extract relevant columns from associations
  selected_columns <- c("association_id", "pvalue", "risk_frequency", "or_per_copy_number", "range","beta_number",'beta_direction')
  associations_df <- associations_df[, selected_columns, drop = FALSE]
  
  # Step 5: Merge datasets by association_id
  merged_data <- merge(risk_genes_df, associations_df, by = "association_id")
  
  # Drop identical columns in merged_data
  #risk_alleles_df <- as.data.frame(my_associations@risk_alleles)
  #merged_data <- merge(merged_data, risk_alleles_df, by = "association_id")
  
  # Drop rows with NA in ensembl_id after merging
  merged_data <- merged_data %>% filter(!is.na(ensembl_id))
  cat("No. of rowsafter dropping NAs in ensembl_id:", nrow(merged_data), "\n")
  
  # Step 6: Return the merged data frame
  return(merged_data)
}

process_association<- function(my_association) { #input an association object, merged desired column and return a dataset
  my_association<-my_association
  # Step 2: Check the number of associations retrieved
  cat("Total associations retrieved:", gwasrapidd::n(my_associations), "\n")
  
  # Step 3: Convert associations and risk genes slots to data frames
  associations_df <- as.data.frame(my_associations@associations)
  risk_genes_df <- as.data.frame(my_associations@ensembl_ids)
  
  # Step 4: Extract relevant columns from associations
  selected_columns <- c("association_id", "pvalue", "risk_frequency", "or_per_copy_number", "range","beta_number",'beta_direction')
  associations_df <- associations_df[, selected_columns, drop = FALSE]
  
  # Step 5: Merge datasets by association_id
  merged_data <- merge(risk_genes_df, associations_df, by = "association_id")
  
  # Drop identical columns in merged_data
  #risk_alleles_df <- as.data.frame(my_associations@risk_alleles)
  #merged_data <- merge(merged_data, risk_alleles_df, by = "association_id")
  
  # Drop rows with NA in ensembl_id after merging
  merged_data <- merged_data %>% filter(!is.na(ensembl_id))
  cat("No. of rowsafter dropping NAs in ensembl_id:", nrow(merged_data), "\n")
  
  # Step 6: Return the merged data frame
  return(merged_data)
}



# Define frequency levels function
assign_frequency_level <- function(data) {
  data %>% mutate(
    frequency_level = case_when(
      risk_frequency < 0.005 ~ "very rare",
      risk_frequency >= 0.005 & risk_frequency < 0.05 ~ "rare",
      risk_frequency >= 0.05 & risk_frequency < 0.1 ~ "uncommon",
      risk_frequency >= 0.1 ~ "common",
      TRUE ~ NA_character_ # Handle cases where risk_frequency is NA or missing
    )
  )
}

assign_beta_level <- function(data, column = "beta_number", new_column = "beta_level", quantiles = c(0, 0.25, 0.5, 0.75, 1)) {
  # Ensure the input column exists in the dataset
  if (!column %in% colnames(data)) {
    stop(paste("Column", column, "not found in the dataset."))
  }
  
  # Compute quantile thresholds
  thresholds <- quantile(data[[column]], probs = quantiles, na.rm = TRUE)
  
  # Assign beta_level based on quantiles
  data[[new_column]] <- cut(
    data[[column]],
    breaks = thresholds,
    include.lowest = TRUE,
    labels = paste("Q", seq_along(quantiles[-1]), sep = "")
  )
  
  return(data)
}


# Function to reorder clusters by numeric order
reorder_clusters <- function(data, column) {
  data %>%
    mutate(
      !!column := factor(
        .data[[column]],
        levels = unique(.data[[column]][order(as.numeric(gsub("Cluster (\\d+):.*", "\\1", .data[[column]])))])
      )
    )
}

cluster_id_to_name <- function(input_list, cluster = NULL) {
  if (is.null(cluster)) stop("Error: You haven't specified a cluster yet.")
  
  cluster_dict <- switch(cluster,
                         "Brain" = cluster_dict_brain,
                         "SingleCell" = cluster_dict_sc,
                         "Tissue" = cluster_dict_tissue,
                         stop("Error: Invalid cluster type. Choose 'Brain', 'SingleCell', or 'Tissue'."))
  
  cluster_dict$cluster_id <- as.character(cluster_dict$cluster_id)  # Ensure IDs are character for matching
  sapply(input_list, function(id) cluster_dict[[1]][match(as.character(id), cluster_dict$cluster_id)])
}


#Plotting cluster barplot
plot_cluster <- function(data, cluster_column, title=disease_name) {
  ggplot(data, aes_string(x = cluster_column, fill = "..count..")) +
    geom_bar() +
    labs(title = title, x = cluster_column, y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


plot_lollipop <- function(fisher_result, p_threshold = 0.05, title = "Fold Enrichment Across Clusters") {
  ggplot(data = fisher_result, 
         aes(x = reorder(cluster_name, fold_enrichment), 
             y = fold_enrichment)) +
    # Add lollipop sticks with dynamic color based on significance
    geom_segment(aes(xend = reorder(cluster_name, fold_enrichment), 
                     y = 0, yend = fold_enrichment, 
                     color = fisher_p < p_threshold), 
                 size = 1, show.legend = FALSE) +
    # Add dots with dynamic size and color based on significance
    geom_point(aes(color = fisher_p < p_threshold, 
                   size = fisher_p < p_threshold), 
               show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "orange", "FALSE" = "grey")) +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2)) +
    # Annotate p-values
    geom_text(aes(label = ifelse(fisher_p < p_threshold, paste0("p=", round(fisher_p, 4)), "")), 
              hjust = -0.3, vjust = 0, size = 3, color = "black") +
    # Add a reference line
    geom_hline(yintercept = 1, linetype = "dashed", color = "skyblue3") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = "Cluster Name", y = "Fold Enrichment", title = title) +
    theme_minimal()
}

#plot barplot using list column
plot_bar_list <- function(data, column_name) {
  ggplot(data.frame(Value = unlist(data[[column_name]])) %>%
           count(Value, sort = TRUE), aes(x = reorder(Value, -n), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_classic() +
    labs(x = column_name, y = "Frequency", title = paste("Barplot of", column_name)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

summarize_top_features <-function(data, column_name) {
  total <- nrow(data)
  data.frame(Value = unlist(data[[column_name]])) %>%
    count(Value, name = "Count", sort = TRUE) %>%
    arrange(desc(Count)) %>% # Order by count in descending order
    rename_with(~column_name, "Value") %>%
    mutate(Frequency =Count / total) # Add frequency column
}

plot_pie_list <- function(data, column_name) {
  data %>%
    mutate(!!column_name := sapply(.[[column_name]], toString)) %>%
    separate_rows(!!sym(column_name), sep = ",\\s*") %>%
    count(!!sym(column_name), name = "count") %>%
    mutate(label = ifelse(rank(-count) <= 5, !!sym(column_name), "")) %>%
    ggplot(aes(x = "", y = count, fill = !!sym(column_name))) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_viridis_d() +
    labs(title = "Protein Class Distribution (Top 5 Highlighted)", fill = column_name) +
    theme_void()
}

plot_boots_hist<- function(cluster.no, bootstrap_df=bootstrap_df, merged_data=merged_data) {
  ggplot(data.frame(Counts = bootstrap_df[, cluster.no]), aes(x = Counts)) +
    geom_histogram(fill = "steelblue", color = "black", binwidth = 1) +
    geom_vline(xintercept = merged_data$risk_Count[cluster.no], linetype = "dashed", color = "red") +
    labs(
      x = "Counts",
      y = "Frequency",
      title = paste("Histogram of Counts for Cluster", cluster.no)
    ) +
    theme_minimal()
}

run_fisher_test <- function(
    data = UMAP_data,
    gene_column = "gene",
    cluster_column,
    risk_genes = new_data_filt$gene,
    p_threshold = 0.05,
    cluster_type = NULL  # Specify cluster type for cluster_id_to_name
) {
  if (is.null(cluster_type)) stop("Error: You haven't specified a cluster type for translation.")
  
  # Perform Fisher test and calculate p-values
  result <- data %>%
    dplyr::select(all_of(c(gene_column, cluster_column))) %>%
    group_by(across(all_of(cluster_column))) %>%
    summarise(
      overlap = sum(.data[[gene_column]] %in% risk_genes),
      cluster_size = dplyr::n(),  # Explicitly qualify n() from dplyr
      fisher_p = {
        observed = sum(.data[[gene_column]] %in% risk_genes)
        total_genes = n_distinct(data[[gene_column]])
        matrix = matrix(c(
          observed,
          cluster_size - observed,
          length(risk_genes) - observed,
          total_genes - cluster_size - length(risk_genes) + observed
        ), nrow = 2)
        fisher.test(matrix, alternative = "greater", workspace = 2e7)$p.value
      },
      .groups = "drop"
    ) %>%
    filter(fisher_p < p_threshold)
  
  # Translate cluster IDs to names
  result <- result %>%
    mutate(cluster_name = cluster_id_to_name(.data[[cluster_column]], cluster = cluster_type))
  
  return(result)
}
