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

parsing_taxa <- function(tax_otu_df, raw_otus_df, output_xlsx = "tax_assignments__to_validate.xlsx") {
  
  # 0. Merge taxonomic table with OTU sequences
  raw_data_with_seqs <- merge(
    tax_otu_df,
    raw_otus_df,
    by = "OTUs",
    all.x = TRUE
  ) %>%
    relocate(Fasta, .after = OTUs) %>%
    mutate(Comp_Fasta = nchar(Fasta)) %>%
    relocate(Comp_Fasta, .after = Fasta)
  
  ### Preparing UNITE derived DB ###
  filtered_data_with_flags <<- raw_data_with_seqs %>%
  # Step 1: Replace underscores with spaces in both Species and Genus columns
  mutate(
    Species = str_replace_all(Species, "_", " "),
    Genus   = str_replace_all(Genus, "_", " ")
  ) %>%
  # Step 2: Remove any genus-derived substring from Species
  mutate(Species = mapply(
    function(sp, gn) {
      if (!is.na(sp) && !is.na(gn) && gn != "") {
        # split genus into words
        gn_words <- unlist(str_split(gn, "\\s+"))
        # sequentially remove all genus words from species
        for (w in gn_words) {
          if (nzchar(w)) {
            sp <- str_remove_all(sp, regex(w, ignore_case = TRUE))
          }
        }
        str_trim(str_squish(sp))
      } else {
        sp
      }
    },
    Species, Genus,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  ))
  
  # 1. Initial data cleaning
  filtered_data_with_flags <<- filtered_data_with_flags %>%
  mutate(across(c(Species),
                ~ ifelse(grepl("uncultured", ., ignore.case = TRUE), "", .))) %>%
  mutate(Species = ifelse(grepl("sp\\.", Species, ignore.case = TRUE) & (is.na(Genus) | Genus == ""), "", Species)) %>%
  mutate(Species = ifelse(!is.na(Species) & Species != "" & (is.na(Genus) | Genus == ""), "", Species)) %>%
  mutate(Species = ifelse(grepl("^sp\\.\\s*[0-9]+$", Species, ignore.case = TRUE), "", Species)) %>%
  mutate(Species = ifelse(str_detect(Species, "[A-Z]{2}|[A-Z][0-9]|[0-9][A-Z]|[0-9]{2}"),
                          str_replace_all(Species, "[^a-záéíóúâêîôûãõç. ]", " "),
                          Species)) %>%
  mutate(Species = str_replace_all(Species, "\\s+", " ")) %>%
  mutate(Species = ifelse(str_detect(Species, "^sp\\.[:alnum:.]$") | 
                            str_detect(Species, "^sp\\..*[^[:alnum:][:space:]]"),
                          "sp.", Species)) %>%
  mutate(Species = str_replace_all(Species, "\\b(sp\\.?)(\\s|$)", "")) %>%
  mutate(Species = str_trim(Species)) %>%
  mutate(Species = ifelse(nchar(Species) == 1, "", Species)) %>%
  mutate(Validation = case_when(
    str_count(Species, "\\S+") > 1 |
      str_detect(Species, "[^a-záéíóúâêîôûãõç]") ~ "Pending",
    Species != "" ~ "Validated",
    TRUE ~ "Validated"
  )) %>%
  relocate(Validation, .after = 1)
  
  ### NCBI VALIDATION SECTION ###
  validated_data <- filtered_data_with_flags %>%
    # Keep only validated records with non-empty Species
    filter(
      Validation == "Validated",
      !is.na(Species),
      Species != ""
    ) %>%
    # Build full name from Genus + Species
    mutate(
      Full_Name = str_trim(paste(Genus, Species, sep = " "))
    ) %>%
    group_by(Full_Name) %>%
    ungroup()
  
 if (nrow(validated_data) > 0) {
  # Download/update local NCBI taxonomy
  db_download_ncbi()
  src <- src_ncbi()
  
  # Query NCBI taxonomy classification
  taxa_names <- unique(validated_data$Full_Name)
  taxa_names <- taxizedb::classification(taxa_names, db = "ncbi")
  
  # Helper function: check if object is a data.frame
  is_dataframe <- function(obj) {
    inherits(obj, "data.frame")
  }
  
  # Ensure all taxa_names elements are data.frames
  for (i in seq_along(taxa_names)) {
    if (!is_dataframe(taxa_names[[i]])) {
      taxa_names[[i]] <- data.frame(name = NA, rank = NA, id = NA)
    }
  }
  
  # Combine results into a single dataframe
  taxa_names_df <- do.call(rbind, lapply(names(taxa_names), function(name) {
    df <- taxa_names[[name]]
    df$taxid <- name
    return(df)
  }))
  
  # Reorganize and filter ranks
  taxa_names_df <- taxa_names_df[, c("taxid", names(taxa_names_df)[1:(ncol(taxa_names_df) - 1)])]
  taxa_names_df <- taxa_names_df[, -ncol(taxa_names_df)]
  
  taxa_names_df <- taxa_names_df %>%
    filter(rank %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species", NA)) %>%
    group_by(taxid) %>%
    distinct()
  
  # Reshape and ensure 'species' column exists
  identified_otus <- taxa_names_df %>%
    pivot_wider(names_from = rank, values_from = name)
  
  # If 'species' column does not exist, create it
  if (!"species" %in% names(identified_otus)) {
    identified_otus$species <- NA_character_
  }
  
  # Rename 'species' to 'NCBI_species' for clarity
  names(identified_otus)[names(identified_otus) == "species"] <- "NCBI_species"
  
  # Merge NCBI results with validated data
  validated_data_updated <- validated_data %>%
    left_join(identified_otus %>% select(taxid, NCBI_species), by = c("Full_Name" = "taxid")) %>%
    mutate(
      Validation = ifelse(
        is.na(NCBI_species),
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
    select(-Full_Name, -NCBI_species)
  
    # Update main dataset with validation results
  filtered_data_with_flags <<- filtered_data_with_flags %>%
    rows_update(
      validated_data_updated,
      by = "OTUs"
    )
 }
  
    ### END OF NCBI VALIDATION SECTION ###
  
  # 2. Export results to XLSX file
  wb <- createWorkbook()
  addWorksheet(wb, "raw_data_with_seqs")
  addWorksheet(wb, "filtered_data_with_flags")
  writeData(wb, sheet = "raw_data_with_seqs", raw_data_with_seqs)
  writeData(wb, sheet = "filtered_data_with_flags", filtered_data_with_flags)
  
  dir.create(dirname(xlsx), recursive = TRUE, showWarnings = FALSE)
  saveWorkbook(wb, xlsx, overwrite = TRUE)
  
  message(paste0("File saved as: tax_assignments__to_validate.xlsx"))
  
  }

