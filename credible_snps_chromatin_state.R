# download the xref file for EpiMap and store it as "xref"
xref <- read_tsv(url("http://bartzabel.ls.manchester.ac.uk/worthingtonlab/functional_genomics_pub/epimap/epimap_xref.txt"))

# write "epimap_xref" to the current directory as a .tsv file
write_tsv(xref, "epimap_xref.tsv")

# read in the local "epimap_xref.tsv" file and add a column for the .bb download URL
xref <- read_tsv("epimap_xref.tsv") %>%
  mutate(url = paste0("https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/bigBed/",
                      BSSID, "_18_CALLS_segments.bb"))

# read in the FLS Sample Information .tsv file and add a column for the .bb download URL, store this in "sample_info"
sample_info <- read_tsv("http://bartzabel.ls.manchester.ac.uk/worthingtonlab/functional_genomics_pub/FLS/FLS_Sample_Info.tsv") %>%
  mutate(url = paste0("http://bartzabel.ls.manchester.ac.uk/worthingtonlab/functional_genomics_pub/FLS/ChromHMM/",
                      ID, "_18_segments.bb"))

# change the headings of "sample_info" to match "xref"
sample_info <- rename(sample_info, "BSSID" = "ID")

# bind the "sample_info" rows to the bottom of "xref"
xref <- rbind(xref, sample_info)

# creates an empty vector for results
results <- c()

# for each row of the "xref" data frame
for(i in 1:nrow(xref)) {
  
  # if the URL on the current row begins with "http",
  # get the URL header and store it in "status"
  if(startsWith(as.character(xref[i, "url"]), "http")) {
    status <- httr::HEAD(as.character(xref[i, "url"]))
  }
  
  # if the status code of the current URL is not 200,
  # move on to the next row
  if(status$status_code != 200) {
    next
  }
  
  # store the current row's URL as "xref_url"
  xref_url <- xref$url[i]
  
  # store the local destination and name of the .bb file as "bb_local"
  bb_local <- paste("xref_", as.character(i), ".bb",sep = "")
  
  # download the .bb file from the URL
  download.file(url = xref_url, destfile = bb_local, method = "curl")
  
  # update progress to console with a loading bar
  progress <- round((i/nrow(xref)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("Credible SNPs: EpiMap ", load_bar, " ", as.character(progress), "%", sep = ""))
  
  # import the file into the local environment as "bb"
  bb <- import(bb_local, format = "bigBed")
  
  # find where the credible SNPs and the .bb file overlap
  overlaps <- findOverlaps(cred_snps, bb)
  
  # store the number of overlaps as "len_over"
  len_over <- length(overlaps)
  
  # if there are no overlaps between the credible SNPs and the current .bb file,
  # move on to the next xref row
  if(len_over == 0) {
    # delete the downloaded .bb file from the working directory, as it takes up a lot of space
    unlink(bb_local)
    next
  }
  
  # for each instance of overlap between the credible SNPs and the .bb file
  for(j in 1:len_over) {
    # store the index for the query hits
    q_index <- queryHits(overlaps[j])
    # store the index for the subject hits
    s_index <- subjectHits(overlaps[j])
    # store the credible SNP associated with the query hit as "q"
    q <- cred_snps[q_index]
    # append this information to the results vector with the name of the SNP the xref information and the .bb subject name
    results <- rbind(results, c(q$VARIANT_RS, xref[i, c("Group", "Extended Info", "MNEMONIC")], bb[s_index]$name))
  }
  
  # delete the downloaded .bb file from the working directory, as it takes up a lot of space
  unlink(bb_local)
}

# store the results matrix in a data frame called "results_df"
results_df <- unique(as.data.frame(results))

# rename the columns in "results_df"
colnames(results_df) <- c("SNP","GROUP","NAME","TISSUE", "STATE")

# write "results_df" to a .csv file
write.table(as.matrix(results_df), file = "cred_snps_data/cred_snps_38_chrom_states.csv", quote = F, sep = ",",
            row.names = F)

# reorder "results_df" to categorise by tissue and store it as "results_by_tissue"
results_by_tissue <- results_df %>% dplyr::select(SNP, TISSUE, STATE) %>% distinct() %>% group_by(SNP, STATE) %>%
                     summarise(TISSUE = paste(unique(TISSUE), collapse = ";"), .groups = "keep")

# reorder "results_by_tissue" to have chromatin states as column headings and store as "epimap_results_by_tissue"
epimap_results_by_tissue <- results_by_tissue %>% spread(STATE, TISSUE, fill = "")

# write "epimap_results_by_tissue" to a .csv file
write.table(as.matrix(epimap_results_by_tissue), file = "cred_snps_data/cred_snps_38_chrom_states_by_tissue.csv",
            quote = F, sep = ",", row.names = F)

# reorder "results_df" to categorise by group and store it as "results_by_group"
results_by_group <- results_df %>% dplyr::select(SNP, GROUP, STATE) %>% distinct() %>% group_by(SNP, GROUP) %>%
                    summarise(STATE = paste(unique(STATE), collapse = ";"), .groups = "keep")

# reorder "results_by_group" to have groups as column headings and store as "epimap_results_by_group"
epimap_results_by_group <- results_by_group %>% spread(GROUP, STATE, fill = "")

# write "epimap_results_by_group" to a .csv file
write.table(as.matrix(epimap_results_by_group), file = "cred_snps_data/cred_snps_38_chrom_states_by_group.csv",
            quote = F, sep = ",", row.names = F)

# store "cred_snps" as a data frame
cred_snps <- as.data.frame(cred_snps)

# add the EpiMap annotations to "cred_snps"
cred_snps <- cbind(cred_snps, epimap_results_by_group[,2:ncol(epimap_results_by_group)])

# return SNP loci to original width
cred_snps$start <- cred_snps$start + 100
cred_snps$end <- cred_snps$end - 100
cred_snps$width <- cred_snps$end - cred_snps$start + 1

# remove the unnecessary objects from the environment
rm(bb, epimap_results_by_group, epimap_results_by_tissue, overlaps, q, results, results_by_group, results_by_tissue,
   results_df, sample_info, status, xref, bb_local, i, j, len_over, q_index, s_index, xref_url, load_bar, progress)