# store the URL for the Capture Hi-C zip file
url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/bin/mmc4.zip"

# download the zip file from the URL
download.file(url, basename(url))

# extract the files
unzip("mmc4.zip")
unzip("DATA_S1.zip")

# remove all the unncecessary files
unlink(x = c("mmc4.zip", "DATA_S1.zip", "ActivePromoterEnhancerLinks.tsv", "__MACOSX", "DataS1-README.pdf",
             "PCHiC_vs_rCHiC_peak_matrix.tsv"), recursive = T)

# store the Capture Hi-C .tsv file as "chic"
chic <- read_tsv("PCHiC_peak_matrix_cutoff5.tsv")

# add two new columns to "cred_snps" for gene symbols and cell types
cred_snps[, "CHi-C_SYMBOLS"] <- NA
cred_snps[, "CHi-C_CELLS"] <- NA

# loop through each of the SNPs
for(i in 1:nrow(cred_snps)) {
  # store the chromosome the SNP is on as "snp_chr"
  snp_chr <- str_sub(cred_snps[i, "seqnames"], 4)
  # store the SNP position as "snp_pos"
  snp_pos <- as.numeric(cred_snps[i, "start"])
  # store a subset of "chic" which matches the SNP chromosome as "temp_df"
  temp_df <- chic[which(chic$oeChr == snp_chr), ]
  # subset "temp_df" where the SNP position is between the start and end of the Capture Hi-C bounds
  temp_df <- temp_df[which(temp_df$oeStart < snp_pos & temp_df$oeEnd > snp_pos), ]
  # if there is no data in this subset,
  # move onto the next SNP
  if(nrow(temp_df) == 0) {
    next
  }
  
  # create empty vectors for genes and cells
  genes <- c()
  cells <- c()
  
  # loop through each row of "temp_df"
  for(j in 1:nrow(temp_df)) {
    # if one of the cell types has a score greater than 5,
    # append the gene symbols to "genes"
    if(sum(temp_df[j, 12:28] > 5) > 0) {
      genes <- append(genes, as.character(temp_df[j, "oeName"]))
    }
    
    # loop through each cell type
    for(k in 12:28) {
      # if the cell has a score greater than 5,
      # append the cell type to "cells"
      if(temp_df[j, k] > 5) {
        cells <- append(cells, colnames(temp_df)[k])
      }
    }
  }
  
  # collapse and add the gene symbols and cell types for each SNP to "cred_snps"
  cred_snps[i, "CHi-C_SYMBOLS"] <- paste(sort(unique(genes)), collapse = ";")
  cred_snps[i, "CHi-C_CELLS"] <- paste(sort(unique(cells)), collapse = ";")
  
  # update progress to console with a loading bar
  progress <- round((i/nrow(cred_snps)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("Credible SNPs: CHi-C ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# delete the CHi-C file from the directory
unlink("PCHiC_peak_matrix_cutoff5.tsv")

# organise the data frame into its final form
cred_snps_final <- data.frame(LEAD_SNP = cred_snps$LEAD_RS, NUMBER_IN_CREDIBLE_SET = cred_snps$CRED_SET,
                              SNP = cred_snps$VARIANT_RS, CHROMOSOME = str_sub(cred_snps$seqnames, 4), LOCUS = cred_snps$start,
                              ALLELES = cred_snps$ALLELES, PIP = cred_snps$PIP, VEP_DESCRIPTION = cred_snps$VEP_CONSEQUENCES,
                              VEP_GENES = cred_snps$VEP_SYMBOLS)

# bind the annotations together
cred_snps_final <- cbind(cred_snps_final, cred_snps[, 15:ncol(cred_snps)])

# convert all headings to upper case
colnames(cred_snps_final) <- toupper(colnames(cred_snps_final))

# write "cred_snps_final" to a .csv file
write.table(as.matrix(cred_snps_final), file = "cred_snps_data/cred_snps_annotation.csv",
            quote = F, sep = ",", row.names = F)

# remove unnecessary objects from the environment
rm(temp_df, cells, genes, i, j, k, load_bar, progress, snp_chr, snp_pos, url, chic, cred_snps_seqinfo, cred_snps, cred_snps_final)