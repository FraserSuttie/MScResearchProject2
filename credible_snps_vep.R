# make a copy of RA SNPs  to edit
cred_snps <- ra_snps

# create an empty vector for coord for all SNPs
coord_total <- c()

# loop through all cred SNPs to extract and store their details
for(i in 1:nrow(cred_snps)) {
  # store the chromosome of the current SNP as "snp_chr"
  snp_chr <- as.character(ra_snps[i, "seqnames"])
  # store the locus of the current SNP as "snp_loc"
  snp_loc <- as.character(ra_snps[i, "start"])
  # paste the SNP chromosome and locus together to make "coord"
  coord <- paste(snp_chr, ":", snp_loc, sep = "")
  # add the current "coord" to "coord_total"
  coord_total <- append(coord_total, coord)
}

# add the coord and allele data to "cred_snps"
cred_snps[, "COORD"] <- coord_total

# remove the GRanges data columns from "cred_snps"
cred_snps <- cred_snps[, 6:ncol(cred_snps)]

# remove the unnecessary objects from the environment
rm(coord, coord_total, i, snp_chr, snp_loc)

# store a copy of "cred_snps" as "data_out"
data_out <- cred_snps

# make "cred_snps" an empty vector
cred_snps <- c()

for(i in 1:nrow(data_out)) {
  # select the current row of "data_out"
  row <- data_out[i, ]
  # store these annotations in the GRanges class, which is designed to be a container for genomic annotations
  cred_snp <- as(as.character(row["COORD"]), "GRanges")
  # add the GRanges annotations to the "cred_snp" list for the current SNP
  values(cred_snp) <- DataFrame(VARIANT_ID = row$VARIANT_ID, VARIANT_RS = row$VARIANT_RS, CRED_SET = row$CRED_SET,
                                ALLELES = row$ALLELES, PIP = row$PIP, LEAD_ID = row$LEAD_ID, LEAD_RS = row$LEAD_RS)
  # paste the URL to find the variant effect predictor (VEP) information for the current LD SNP
  raw_out_json <- httr::GET(paste0("https://rest.ensembl.org/vep/human/id/", cred_snp$VARIANT_RS,
                                   "?content-type=application/json"))
  
  # if "raw_out_json" is not a successful request,
  # move on to annotating the next SNP
  if(raw_out_json$status_code != 200) {
    next
  }
  
  # store the parsed "raw_out_json" content
  data_out_json <- httr::content(raw_out_json,as = "parsed")
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
  
  # add the VEP information to the "cred_snp" list for the current SNP
  values(cred_snp) <- DataFrame(values(cred_snp), VEP_CONSEQUENCES = consequences, VEP_SYMBOLS = genes)
  
  # append the current LD SNP annotation to the "cred_snps" list
  cred_snps <- append(cred_snps, cred_snp)
  
  # update progress to console with a loading bar
  progress <- round((i/nrow(data_out)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("Credible SNPs: VEP ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# store the sequence information for "cred_snps"
cred_snps_seqinfo <- seqinfo(cred_snps)

# convert "cred_snps" from GRanges to a data frame
cred_snps <- as.data.frame(cred_snps)

# write cred_snps to a .csv file
write_csv(cred_snps, file = "cred_snps_data/cred_snps_38_vep.csv")

# decrease the start locus of all SNPs by 100
cred_snps$start <- cred_snps$start - 100
# increase the end locus of all SNPs by 100
cred_snps$end <- cred_snps$end + 100
# update the width of these new SNP ranges
cred_snps$width <- cred_snps$end - cred_snps$start + 1

# convert "cred_snps" into a GRanges object called "cred_snps
cred_snps <- makeGRangesFromDataFrame(cred_snps, keep.extra.columns = TRUE, seqinfo = cred_snps_seqinfo)

# remove the unnecessary objects from the environment
rm(cred_snp, data_out, data_out_json, raw_out_json, row, temp_results, consequences, genes, i, load_bar, progress)