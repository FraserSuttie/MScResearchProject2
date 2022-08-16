# create new data frame to count eQTL associations for each SNP
cred_eqtls <- data.frame(SNP = cred_snps$VARIANT_RS, CHR = str_sub(cred_snps$seqnames, 4), LOCUS = as.character(cred_snps$start))

# create an empty vector for snps_eqtls
snps_eqtls <- c()

# loop through each of the SNPs in "cred_eqtls"
for(i in 1:nrow(cred_eqtls)) {
  # paste the current SNP into the url to get the raw eQTL data
  raw_out <- httr::GET(url = paste("https://www.ebi.ac.uk/eqtl/api/associations?p_upper=0.05&size=1000&variant_id=",
                                   cred_eqtls[i, "SNP"], sep = ""))
  
  # if the request is not successful,
  # move to the next SNP
  if(raw_out$status_code != 200) {
    next
  }
  
  # store the json conent of "raw_out" as "raw_out_json"
  raw_out_json <- httr::content(raw_out, as = "parsed")
  
  # store the eQTL association data from the json as "eqtls"
  eqtls <- raw_out_json$"_embedded"$"associations"
  
  # if the number of eQTLs for the current SNP is greater than 0,
  # create a row for the SNP's eQTL associations
  if(length(eqtls) > 0) {
    # create a data frame of all the eQTL associations for the current SNP
    snp_eqtls <- bind_rows(raw_out_json$"_embedded"$"associations")
    
    # only keep the columns that we need
    snp_eqtls <- data.frame(SNP = snp_eqtls$rsid, GROUP = snp_eqtls$qtl_group, GENEID = snp_eqtls$gene_id, GENESYM = NA)
    
    # create a vector of all the unique gene IDs
    gene_ids <- unique(snp_eqtls$GENEID)
    
    # create a dictionary to link gene IDs with gene symbols
    gene_symbols <- ensembldb::select(EnsDb.Hsapiens.v86, keys = gene_ids, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
    
    # add the gene symbol to each eQTL association
    for(j in 1:nrow(snp_eqtls)) {
      gene_sym <- gene_symbols[which(gene_symbols$GENEID == snp_eqtls[j, "GENEID"]), "SYMBOL"]
      if(length(gene_sym) > 0) {
        snp_eqtls[j, "GENESYM"] <- gene_sym
      }
    }
    
    # remove gene IDs from the data frame
    snp_eqtls <- data.frame(SNP = snp_eqtls$SNP, GROUP = snp_eqtls$GROUP, GENE = snp_eqtls$GENESYM)
    
    # reorder "snp_eqtls" to categorise by "GROUP
    snp_eqtls <- snp_eqtls %>% dplyr::select(SNP, GROUP, GENE) %>% distinct() %>% group_by(SNP, GROUP) %>%
      summarise(GENE = paste(unique(GENE), collapse = ";"), .groups = "keep")
    
    # reorder "snp_eqtls" to have groups as column headings
    snp_eqtls <- snp_eqtls %>% spread(GROUP, GENE, fill = "")
    
    # if this is the first loop,
    # store "snp_eqtls" as "snps_eqtls",
    # otherwise, bind the rows together
    if(i == 1) {
      snps_eqtls <- snp_eqtls
    } else {
      snps_eqtls <- bind_rows(snps_eqtls, snp_eqtls)
    }
  }
  
  # update progress to console with a loading bar
  progress <- round((i/nrow(cred_eqtls)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("Credible SNPs: eQTLs ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# write "snps_eqtls" as a .csv file
write_csv(snps_eqtls, "cred_snps_data/cred_snps_38_eqtls.csv")

# remove rows from "cred_snps" that do not have any eQTL associations
cred_snps <- cred_snps[which(cred_snps$VARIANT_RS %in% snps_eqtls$SNP),]

#
if(sum(cred_snps$VARIANT_RS != snps_eqtls$SNP) == 0) {
  snps_eqtls <- snps_eqtls[, 2:ncol(snps_eqtls)]
  cred_snps <- cbind(cred_snps, snps_eqtls)
}

# remove the unnecessary objects from the environment
rm(cred_eqtls, eqtls, gene_symbols, raw_out, raw_out_json, snp_eqtls, gene_ids, gene_sym, i, j, load_bar, progress, snps_eqtls)