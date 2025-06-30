combine_tax_otu <- function(tax_file, otu_file, fasta_file) {
  
  # Step 1: Process the Tax table (taxonomy)
  Tax <- suppressWarnings({
    df <- read.table(tax_file, header = FALSE, sep = "\t")
    if (ncol(df) >= 4) df <- df[, -4] # Preliminary column removal if needed
    names(df) <- c("OTUs", "Taxonomy", "PID")[1:ncol(df)]
    
    df %>%
      filter(!grepl("Unassigned", Taxonomy, ignore.case = TRUE)) %>%
      separate(Taxonomy, 
               into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               sep = ";") %>%
      mutate(across(Kingdom:Species, ~ str_replace(., "^.*__", ""))) # Remove prefix patterns
  })
  
  # Step 2: Process the OTU table (abundance data)
  OTU <- fread(otu_file, sep = "\t", header = TRUE)
  colnames(OTU)[1] <- gsub("^#OTU ID", "OTUs", colnames(OTU)[1]) # Standardize OTU column name
  
  # Step 3: Import fasta file
  fasta_seq <- Biostrings::readDNAStringSet(fasta_file)
  raw_otus <- data.frame(
    OTUs = names(fasta_seq),
    Fasta = as.character(fasta_seq),
    stringsAsFactors = FALSE
  )
  
  # Step 4: Merge Tax + OTU tables
  combined_data <- merge(Tax, OTU, by = "OTUs")
  
  # Step 5: Return results as named list
  return(list(
    tax_otu_table = combined_data,
    raw_otus = raw_otus
  ))
}

parsing_taxa <- function(tax_otu_df, raw_otus_df, output_xlsx = "resultados_parciais_para_avaliação.xlsx") {
  # 0. Initial Merging
  raw_data_with_seqs <<- merge(
    tax_otu_df,
    raw_otus_df[, c("OTUs", "Fasta")],
    by = "OTUs",
    all.x = TRUE
  ) %>%
    relocate(Fasta, .after = OTUs) %>%
    mutate(Comp_Fasta = nchar(Fasta)) %>%
    relocate(Comp_Fasta, .after = Fasta)
  
  # 1. Initial data cleaning
filtered_data_with_flags <<- raw_data_with_seqs %>%
  # 1. Remove "sp." when Genus is empty/NA
  mutate(Species = ifelse(grepl("sp\\.", Species, ignore.case = TRUE) & (is.na(Genus) | Genus == ""), "", Species)) %>%
  
  # 2. Remove records without associated genus
  mutate(Species = ifelse(!is.na(Species) & Species != "" & (is.na(Genus) | Genus == ""), "", Species)) %>%
  
  # 3. Remove "sp." pattern followed only by numbers
  mutate(Species = ifelse(grepl("^sp\\.\\s*[0-9]+$", Species, ignore.case = TRUE), "", Species)) %>%
  
  # 4. Clean special characters when invalid patterns detected (uppercase/numbers)
  mutate(Species = ifelse(str_detect(Species, "[A-Z]{2}|[A-Z][0-9]|[0-9][A-Z]|[0-9]{2}"),
                         str_replace_all(Species, "[^a-záéíóúâêîôûãõç. ]", " "),
                         Species)) %>%
  
  # 5. Normalize multiple spaces
  mutate(Species = str_replace_all(Species, "\\s+", " ")) %>%
  
  # 6. Simplify invalid "sp." patterns to just "sp."
  mutate(Species = ifelse(str_detect(Species, "^sp\\.[:alnum:.]$") | 
                         str_detect(Species, "^sp\\..*[^[:alnum:][:space:]]"),
                       "sp.", Species)) %>%
  
  # 7. Remove all remaining "sp." or "sp"
  mutate(Species = str_replace_all(Species, "\\b(sp\\.?)(\\s|$)", "")) %>%
  
  # 8. Remove extra leading/trailing spaces
  mutate(Species = str_trim(Species)) %>%
  
  # 9. Remove records with only 1 character
  mutate(Species = ifelse(nchar(Species) == 1, "", Species)) %>%
  
  # 10. Create validation flag:
  mutate(Validation = case_when(
    str_count(Species, "\\S+") > 1 |
    str_detect(Species, "[^a-záéíóúâêîôûãõç]") ~ "Pending",
    Species != "" ~ "Validated",
    TRUE ~ "Validated"
  )) %>%
  relocate(Validation, .after = 1)

 ### NCBI VALIDATION SECTION ###
  
  # Select validated data
  validated_data <- filtered_data_with_flags %>%
    filter(
      Validation == "Validated",
      !is.na(Species),
      Species != ""
    ) %>%
    mutate(
      Full_Name = str_trim(paste(Genus, Species, sep = " "))
    ) %>%
    group_by(Full_Name) %>%
    ungroup()
  
  # NCBI validation
  db_download_ncbi()
  src <- src_ncbi()
  
  taxa_names <- unique(validated_data$Full_Name)
  taxa_names <- taxizedb::classification(taxa_names, db = "ncbi")

  # Helper function to check if object is a dataframe
  is_dataframe <- function(obj) {
    inherits(obj, "data.frame")
  }
  
  # Process each taxonomic name
  for (i in seq_along(taxa_names)) {
    if (!is_dataframe(taxa_names[[i]])) {
      taxa_names[[i]] <- data.frame(name = NA, rank = NA, id = NA)
    }
  }
  
  # Combine all results into one dataframe
  taxa_names_df <- do.call(rbind, lapply(names(taxa_names), function(name) {
    df <- taxa_names[[name]]
    df$taxid <- name
    return(df)
  }))
  
  # Clean and reorganize columns
  taxa_names_df <- taxa_names_df[, c("taxid", names(taxa_names_df)[1:(ncol(taxa_names_df) - 1)])]
  taxa_names_df <- taxa_names_df[, -ncol(taxa_names_df)]
  
  # Filter to only include standard taxonomic ranks
  taxa_names_df <- taxa_names_df %>%
    filter(rank %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species", NA)) %>%
    group_by(taxid) %>%
    distinct()
  
  # Reshape data to wide format
  identified_otus <- taxa_names_df %>%
    pivot_wider(names_from = rank, values_from = name)
  
  # Update validated data with NCBI results
  validated_data_updated <- validated_data %>%
    left_join(identified_otus %>% select(taxid, species), by = c("Full_Name" = "taxid")) %>%
    mutate(
      Validation = ifelse(
        is.na(species),
        "Validated - DB",
        Validation
      ),
      Species = ifelse(
        Validation == "Validated - DB",
        "",
        Species
      ),
      Validation = ifelse(
        Validation == "Validated - DB",
        "Validated",
        Validation
      )
    ) %>%
    select(-Full_Name, -species)
  
  # Update main dataset with validation results
  filtered_data_with_flags <<- filtered_data_with_flags %>%
    rows_update(
      validated_data_updated,
      by = "OTUs"
    )
  
  ### END OF NCBI VALIDATION SECTION ###

  # 2. Export results to XLSX file
  wb <- createWorkbook()
  
  # Create worksheets
  addWorksheet(wb, "raw_data_with_seqs")
  addWorksheet(wb, "filtered_data_with_flags")
  
  # Write data to worksheets
  writeData(wb, sheet = "raw_data_with_seqs", raw_data_with_seqs)
  writeData(wb, sheet = "filtered_data_with_flags", filtered_data_with_flags)
  
  # Save workbook

  dir.create(dirname(xlsx), recursive = TRUE, showWarnings = FALSE)

  saveWorkbook(wb, xlsx, overwrite = TRUE)
  
  # Print confirmation message
  message(paste0("File saved as: tax_assignments__to_validate.xlsx"))
}
