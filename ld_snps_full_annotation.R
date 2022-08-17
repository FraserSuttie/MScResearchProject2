# store lead SNP Rs ID's in "index_snps" vector
index_snps <- unique(ra_snps$LEAD_RS)

# creates an empty vector for LD SNPs
ld_snps <- c()

# loop to identify and annotate all of the LD SNPs for each of the lead SNPs
for(index_snp in index_snps) {
  # paste the URL to find SNPs in LD
  url_str <- paste("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?",
                   "window=500000&genome_build=grch38&token=e25744a9ef7d&pop=EUR&r2_d=r2&var=",
                   index_snp, sep = "")
  
  # get the raw data for the current SNP
  raw_out <- httr::GET(url = url_str)
  
  # if "raw_out" is not a successful request,
  # move onto the next lead SNP
  if(raw_out$status_code != 200) {
    next
  }
  
  # extract the data from the "raw_out" list
  data_out <- read.delim(textConnection(httr::content(raw_out, "text", encoding = "UTF-8")), header = T, sep = "\t")
  
  # if "data_out" contains the string "error",
  # we skip the whole LD SNPs analysis
  if((sum(grepl("error", data_out))) == 0) {
    # refine the number of SNPs in LD to only include those with R2 greater than or equal to 0.8
    data_out <- data_out %>% dplyr::filter(R2 >= 0.8)
    
    # for each LD SNP of the current index SNP
    for(i in 1:nrow(data_out)) {
      # select the i(th) row (SNP) of the refined data
      row <- data_out[i, ]
      
      # store these annotations in the GRanges class, which is designed to be a container for genomic annotations
      ld_snp <- as(as.character(row["Coord"]), "GRanges")
      
      # add the GRanges annotations to the "ld_snp" list for the current SNP
      values(ld_snp) <- DataFrame(name = row$RS_Number, row[c('Alleles', 'MAF', 'R2', 'RegulomeDB')],
                                  IndexSNP = index_snp)
      
      # paste the URL to find the variant effect predictor (VEP) information for the current LD SNP
      raw_out_json <- httr::GET(paste0("https://rest.ensembl.org/vep/human/id/", ld_snp$name,
                                       "?content-type=application/json"))
      
      # if "raw_out" is not a successful request,
      # move onto the next LD SNP
      if(raw_out_json$status_code != 200) {
        next
      }
      
      # store the parsed "raw_out_json" content
      data_out_json <- httr::content(raw_out_json, as = "parsed")
      
      # create an empty vector for consequences
      consequences <- c()
      
      # create an empty string for genes
      genes <- ""
      
      # if one of the components of the "data_out_json" list is "transcript_consequences",
      # store the data from this component in the "consequences" vector
      if("transcript_consequences" %in% names(data_out_json[[1]])) {
        temp_results <- bind_rows(data_out_json[[1]]$transcript_consequences)
        
        consequences <- unlist(unique(temp_results$consequence_terms))
        
        # if one of the components of the "temp_results" list is "gene_symbol",
        # paste together all the mentioned genes into the "genes" string
        if("gene_symbol" %in% names(temp_results)) {
          genes <- paste(unique(temp_results$gene_symbol[!is.na(temp_results$gene_symbol)]), collapse = ";")
        }
      }
      
      # if one of the components of the "data_out_json" list is "regulatory_feature_consequences",
      # store the data from this component in the "consequences" vector
      if("regulatory_feature_consequences" %in% names(data_out_json[[1]])) {
        temp_results <- bind_rows(data_out_json[[1]]$regulatory_feature_consequences)
        
        consequences <- c(consequences, unlist(unique(temp_results$consequence_terms)))
      }
      
      # paste all the consequences and separate them with a ";"
      consequences <- paste(consequences, collapse = ";")
      
      # add the VEP information to the "ld_snp" list for the current SNP
      values(ld_snp) <- DataFrame(values(ld_snp), VEP_CONSEQUENCES = consequences, VEP_SYMBOLS = genes)
      
      # append the current LD SNP annotation to the "ld_snps" list
      ld_snps <- append(ld_snps, ld_snp)
    }
  }
  # update progress to console with a loading bar
  progress <- round((match(index_snp, index_snps)/length(index_snps)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("LD SNPs: VEP ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# store the sequence information for "ld_snps"
ld_snps_seqinfo <- seqinfo(ld_snps)

# convert "ld_snps" from GRanges to a data frame
ld_snps <- as.data.frame(ld_snps) %>% mutate(chr_name = as.numeric(substr(as.character(seqnames), 4, length(seqnames))),
                                             chromosomal_region = paste(chr_name, start, end, sep = ":"))

# write ld_snps to a .csv file
write.table(as.matrix(ld_snps), file = "ld_snps_data/ld_snps_38_vep.csv",
            quote = F, sep = ",", row.names = F)

# remove the unnecessary objects from the environment
rm(data_out, data_out_json, ld_snp, raw_out, raw_out_json, row, temp_results, consequences, genes, i, index_snp,
   index_snps, url_str, load_bar, progress)

# store the number of LD SNP IDs that will need updating
num_ld_snps <- nrow(ld_snps)

# store the number of full blocks of 5 that can be updated and looped
full_repeats <- floor(num_ld_snps/5)

# create an empty vector for SNP Info
snp_info <- c()

# loop to update the LD SNP IDs in blocks of 5
for(i in 1:full_repeats) {
  # retrieve the SNP IDs for the current block
  snp_info_block <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
                          filters = "chromosomal_region", values = ld_snps$chromosomal_region[((5*i)-(4)):(5*i)],
                          mart = snp_mart)
  
  # add the current block to "snp_info"
  snp_info <- rbind(snp_info, snp_info_block)
  
  # update the console with progress
  message(paste("Block(", as.character(i), "/", as.character(full_repeats + 1), ")", sep = ""))
}

# store the number of remaining SNPs to be updated
snps_left <- num_ld_snps - ((full_repeats*5) + 1)

# if there are remaining SNPs to update,
# find this information and add it to "snp_info"
if(snps_left > 0) {
  # retrieve the SNP IDs for this demi-block
  snp_info_block <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
                          filters = "chromosomal_region",
                          values = ld_snps$chromosomal_region[(num_ld_snps - snps_left):(num_ld_snps)],
                          mart = snp_mart)
  
  # add it to "snp_info"
  snp_info <- rbind(snp_info, snp_info_block)
}

# update the console with progress
message(paste("Block(", as.character(full_repeats + 1), "/", as.character(full_repeats + 1), ")", sep = ""))

# create a duplicate of "ld_snps" to update specific Rs IDs called "ld_snps_id"
ld_snps_id <- data.frame(ld_snps, id_update = 0)

# store the single nucleotide results from "snp_info" as "snp_info_one"
snp_info_one <- snp_info[which(snp_info$chrom_start == snp_info$chrom_end), ]

# remove any SNPs that do not begin with "rs" from "snp_info_one"
snp_info_one <- snp_info_one[which(str_starts(snp_info_one$refsnp_id, "rs")), ]

# create a data frame for the frequency of each locus called "num_occur"
num_occur <- data.frame(table(snp_info_one$chrom_start))

# create a table of all the duplicated loci called "dup_table"
dup_table <- snp_info_one[snp_info_one$chrom_start %in% num_occur$Var1[num_occur$Freq > 1], ]

# add columns for matching the reference and variant alleles to "dup_table"
dup_table <- data.frame(dup_table, ref_match = 0, var_match = 0)

# for each SNP in "dup_table"
for(i in 1:nrow(dup_table)) {
  # get the raw information on the current SNP
  snpi <- httr::GET(paste0("https://rest.ensembl.org/vep/human/id/", dup_table$refsnp_id[i],
                           "?content-type=application/json"))
  
  # extract the data from this
  snpi_data <- httr::content(snpi, as = "parsed")
  
  # store the updated SNP alleles
  up_alls <- str_split(snpi_data[[1]]$allele_string, "/")
  
  # store the updated reference allele
  up_ref_all <- up_alls[[1]][1]
  
  # store the updated variant allele
  up_var_all <- up_alls[[1]][-1]
  
  # store the ld_snps alleles for current locus
  ld_snp_alls <- str_split(str_sub(ld_snps_id[which(ld_snps_id$start == dup_table$chrom_start[i]), "Alleles"], 2, -2), "/")
  
  # store the current reference allele
  ld_snp_ref_all <- ld_snp_alls[[1]][1]
  
  # store the current variant allele
  ld_snp_var_all <- ld_snp_alls[[1]][-1]
  
  # if the update and current tables have the same reference allele,
  # change "ref_match" from 0 to 1
  if(sum(up_ref_all %in% ld_snp_ref_all) > 0) {
    dup_table$ref_match[i] <- 1
  }
  
  # if the update and current tables have the same variant allele,
  # change "var_all" from 0 to 1
  if(sum(up_var_all %in% ld_snp_var_all) > 0) {
    dup_table$var_match[i] <- 1
  }
}

# create the column "allele_match" which is a sum of the other two
dup_table[, "allele_match"] <- dup_table$ref_match + dup_table$var_match

# store the duplicated loci
dup_loci <- unique(dup_table$chrom_start)

# create an empty vector called "rm_ids"
rm_ids <- c()

for(i in 1:length(dup_loci)) {
  # create a temporary table to compare the two duplicated loci
  temp_table <- dup_table[which(dup_table$chrom_start == dup_loci[i]), ]
  # store the least comparable SNP as "rm_id"
  rm_id <- temp_table[which(temp_table$allele_match == min(temp_table$allele_match)), "refsnp_id"]
  # add the SNP to remove to the vector "rm_ids"
  rm_ids <- append(rm_ids, rm_id)
}

# create a function to mean "not %in%"
"%!in%" <- function(x,y)!('%in%'(x,y))

# remove the SNPs in "rm_ids" from "snp_info_one"
snp_info_one <- snp_info_one[which(snp_info_one$refsnp_id %!in% rm_ids), ]

# loop through each row in "snp_info_one"
for(i in 1:nrow(snp_info_one)) {
  # store the index of the row in "ld_snps_id" that needs to be updated as "id_index"
  id_index <- match(x = snp_info_one$chrom_start[i], table = ld_snps_id$start)
  
  # change the name of the SNP in "ld_snps_id" to the update from "snp_info"
  ld_snps_id$name[id_index] <- snp_info_one$refsnp_id[i]
  
  # change the value in "id_update" from 0 to 1
  ld_snps_id$id_update[id_index] <- 1
}

# drop unnecessary columns from "ld_snps_id"
drops <- c("chr_name", "chromosomal_region")
ld_snps_id <- ld_snps_id[, !(names(ld_snps_id) %in% drops)]

# write "ld_snps_id" as a .csv file
write_csv(ld_snps_id, "ld_snps_data/ld_snps_full_id.csv")

# store a copy of "ld_snps_id" as "ld_snps_range"
ld_snps_range <- ld_snps_id

# decrease the start locus of all SNPs by 100
ld_snps_range$start <- ld_snps_range$start - 100

# increase the end locus of all SNPs by 100
ld_snps_range$end <- ld_snps_range$end + 100

# update the width of these new SNP ranges
ld_snps_range$width <- ld_snps_range$end - ld_snps_range$start + 1

# convert "ld_snps_range" into a GRanges object called "ld_snps
ld_snps <- makeGRangesFromDataFrame(ld_snps_range, keep.extra.columns = TRUE, seqinfo = ld_snps_seqinfo)

# remove the unnecessary objects from the environment
rm(dup_table, ld_snp_alls, ld_snps_id, ld_snps_range, num_occur, snp_info, snp_info_block, snp_info_one, snpi, snpi_data,
   temp_table, up_alls, drops, dup_loci, full_repeats, i, id_index, ld_snp_ref_all, ld_snp_var_all, num_ld_snps, rm_id,
   rm_ids, snps_left, up_ref_all, up_var_all, "%!in%")

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
  message(paste("LD SNPs: EpiMap ", load_bar, " ", as.character(progress), "%", sep = ""))
  
  # import the file into the local environment as "bb"
  bb <- import(bb_local, format = "bigBed")
  
  # find where the LD SNPs and the .bb file overlap
  overlaps <- findOverlaps(ld_snps, bb)
  
  # store the number of overlaps as "len_over"
  len_over <- length(overlaps)
  
  # if there are no overlaps between the LD SNPs and the current .bb file,
  # move on to the next xref row
  if(len_over == 0) {
    # delete the downloaded .bb file from the working directory, as it takes up a lot of space
    unlink(bb_local)
    next
  }
  
  # for each instance of overlap between the LD SNPs and the .bb file
  for(j in 1:len_over) {
    # store the index for the query hits
    q_index <- queryHits(overlaps[j])
    # store the index for the subject hits
    s_index <- subjectHits(overlaps[j])
    # store the LD SNP associated with the query hit as "q"
    q <- ld_snps[q_index]
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
write.table(as.matrix(results_df), file = "ld_snps_data/ld_snps_38_chrom_states.csv", quote = F, sep = ",",
            row.names = F)

# reorder "results_df" to categorise by tissue and store it as "results_by_tissue"
results_by_tissue <- results_df %>% dplyr::select(SNP, TISSUE, STATE) %>% distinct() %>% group_by(SNP, STATE) %>%
  summarise(TISSUE = paste(unique(TISSUE), collapse = ";"), .groups = "keep")

# reorder "results_by_tissue" to have chromatin states as column headings and store as "epimap_results_by_tissue"
epimap_results_by_tissue <- results_by_tissue %>% spread(STATE, TISSUE, fill = "")

# write "epimap_results_by_tissue" to a .csv file
write.table(as.matrix(epimap_results_by_tissue), file = "ld_snps_data/ld_snps_38_chrom_states_by_tissue.csv",
            quote = F, sep = ",", row.names = F)

# reorder "results_df" to categorise by group and store it as "results_by_group"
results_by_group <- results_df %>% dplyr::select(SNP, GROUP, STATE) %>% distinct() %>% group_by(SNP, GROUP) %>%
  summarise(STATE = paste(unique(STATE), collapse = ";"), .groups = "keep")

# reorder "results_by_group" to have groups as column headings and store as "epimap_results_by_group"
epimap_results_by_group <- results_by_group %>% spread(GROUP, STATE, fill = "")

# write "epimap_results_by_group" to a .csv file
write.table(as.matrix(epimap_results_by_group), file = "ld_snps_data/ld_snps_38_chrom_states_by_group.csv",
            quote = F, sep = ",", row.names = F)

# store "ld_snps" as a data frame
ld_snps <- as.data.frame(ld_snps)

# add the EpiMap annotations to "ld_snps"
ld_snps <- cbind(ld_snps, epimap_results_by_group[,2:ncol(epimap_results_by_group)])

# return SNP loci to original width
ld_snps$start <- ld_snps$start + 100
ld_snps$end <- ld_snps$end - 100
ld_snps$width <- ld_snps$end - ld_snps$start + 1

# remove the unnecessary objects from the environment
rm(bb, epimap_results_by_group, epimap_results_by_tissue, overlaps, q, results, results_by_group, results_by_tissue,
   results_df, sample_info, status, xref, bb_local, i, j, len_over, q_index, s_index, xref_url, load_bar, progress)

# create new data frame to count eQTL associations for each SNP
ld_eqtls <- data.frame(SNP = ld_snps$VARIANT_RS, CHR = str_sub(ld_snps$seqnames, 4), LOCUS = as.character(ld_snps$start))

# create an empty vector for snps_eqtls
snps_eqtls <- c()

# loop through each of the SNPs in "ld_eqtls"
for(i in 1:nrow(ld_eqtls)) {
  # paste the current SNP into the url to get the raw eQTL data
  raw_out <- httr::GET(url = paste("https://www.ebi.ac.uk/eqtl/api/associations?p_upper=0.05&size=1000&variant_id=",
                                   ld_eqtls[i, "SNP"], sep = ""))
  
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
  progress <- round((i/nrow(ld_eqtls)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("LD SNPs: eQTLs ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# write "snps_eqtls" as a .csv file
write_csv(snps_eqtls, "ld_snps_data/ld_snps_38_eqtls.csv")

# remove rows from "ld_snps" that do not have any eQTL associations
ld_snps <- ld_snps[which(ld_snps$VARIANT_RS %in% snps_eqtls$SNP),]

#
if(sum(ld_snps$VARIANT_RS != snps_eqtls$SNP) == 0) {
  snps_eqtls <- snps_eqtls[, 2:ncol(snps_eqtls)]
  ld_snps <- cbind(ld_snps, snps_eqtls)
}

# remove the unnecessary objects from the environment
rm(ld_eqtls, eqtls, gene_symbols, raw_out, raw_out_json, snp_eqtls, gene_ids, gene_sym, i, j, load_bar, progress, snps_eqtls)

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

# add two new columns to "ld_snps" for gene symbols and cell types
ld_snps[, "CHi-C_SYMBOLS"] <- NA
ld_snps[, "CHi-C_CELLS"] <- NA

# loop through each of the SNPs
for(i in 1:nrow(ld_snps)) {
  # store the chromosome the SNP is on as "snp_chr"
  snp_chr <- str_sub(ld_snps[i, "seqnames"], 4)
  # store the SNP position as "snp_pos"
  snp_pos <- as.numeric(ld_snps[i, "start"])
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
  
  # collapse and add the gene symbols and cell types for each SNP to "ld_snps"
  ld_snps[i, "CHi-C_SYMBOLS"] <- paste(sort(unique(genes)), collapse = ";")
  ld_snps[i, "CHi-C_CELLS"] <- paste(sort(unique(cells)), collapse = ";")
  
  # update progress to console with a loading bar
  progress <- round((i/nrow(ld_snps)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("LD SNPs: CHi-C ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# delete the CHi-C file from the directory
unlink("PCHiC_peak_matrix_cutoff5.tsv")

# organise the data frame into its final form
ld_snps_final <- data.frame(LEAD_SNP = ld_snps$LEAD_RS, NUMBER_IN_LD_SET = ld_snps$ld_SET,
                              SNP = ld_snps$VARIANT_RS, CHROMOSOME = str_sub(ld_snps$seqnames, 4), LOCUS = ld_snps$start,
                              ALLELES = ld_snps$ALLELES, PIP = ld_snps$PIP, VEP_DESCRIPTION = ld_snps$VEP_CONSEQUENCES,
                              VEP_GENES = ld_snps$VEP_SYMBOLS)

# bind the annotations together
ld_snps_final <- cbind(ld_snps_final, ld_snps[, 15:ncol(ld_snps)])

# convert all headings to upper case
colnames(ld_snps_final) <- toupper(colnames(ld_snps_final))

# write "ld_snps_final" to a .csv file
write.table(as.matrix(ld_snps_final), file = "ld_snps_data/ld_snps_annotation.csv",
            quote = F, sep = ",", row.names = F)

# remove unnecessary objects from the environment
rm(temp_df, cells, genes, i, j, k, load_bar, progress, snp_chr, snp_pos, url, chic, ld_snps_seqinfo, ld_snps, ld_snps_final)