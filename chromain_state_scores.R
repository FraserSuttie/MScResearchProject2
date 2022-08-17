# create a data frame to compare chromatin states for credible SNPs
chr_check <- data.frame(LEAD_SNP = cred_snps_annotation$LEAD_SNP, SET = cred_snps_annotation$NUMBER_IN_CREDIBLE_SET,
                        SNP = cred_snps_annotation$SNP, PIP = cred_snps_annotation$PIP,
                        BLOOD_STATE = cred_snps_annotation$`BLOOD & T-CELL`, STROMAL_STATE = cred_snps_annotation$STROMAL,
                        FLS_STATE = cred_snps_annotation$FLS, BLOOD_SCORE = NA, STROMAL_SCORE = NA, FLS_SCORE = NA,
                        RA_SCORE = NA)

# create a dictionary to assign scores to each chromatin state
chrom_states <- c("TssA" = 2,
                  "TssFlnk" = 1,
                  "TssFlnkU" = 1,
                  "TssFlnkD" = 1,
                  "Tx" = 2,
                  "TxWk" = 1,
                  "EnhG1" = 2,
                  "EnhG2" = 2,
                  "EnhA1" = 2,
                  "EnhA2" = 2,
                  "EnhWk" = 1,
                  "ZNF/Rpts" = 0,
                  "Het" = 0,
                  "TssBiv" = 0,
                  "EnhBiv" = 0,
                  "ReprPC" = 0,
                  "ReprPCWk" = 0,
                  "Quies" = 0)

# loop through each SNP in to calculate their blood state score
for(i in 1:nrow(chr_check)) {
  # extract the current chromatin states
  snp_states <- str_split(chr_check[i, "BLOOD_STATE"], ";")
  # take the average of the states as a value between 0-1
  state_score <- (mean(chrom_states[snp_states[[1]]]))/2
  # write the average score to "blood_chst"
  chr_check[i, "BLOOD_SCORE"] <- state_score
}

# loop through each SNP in to calculate their stromal state score
for(i in 1:nrow(chr_check)) {
  # extract the current chromatin states
  snp_states <- str_split(chr_check[i, "STROMAL_STATE"], ";")
  # take the average of the states as a value between 0-1
  state_score <- (mean(chrom_states[snp_states[[1]]]))/2
  # write the average score to "blood_chst"
  chr_check[i, "STROMAL_SCORE"] <- state_score
}

# loop through each SNP in to calculate their FLS state score
for(i in 1:nrow(chr_check)) {
  # extract the current chromatin states
  snp_states <- str_split(chr_check[i, "FLS_STATE"], ";")
  # take the average of the states as a value between 0-1
  state_score <- (mean(chrom_states[snp_states[[1]]]))/2
  # write the average score to "blood_chst"
  chr_check[i, "FLS_SCORE"] <- state_score
}

# calculate average score for each SNP
for(i in 1:nrow(chr_check)) {
  chr_check[i, "RA_SCORE"] <- mean(chr_check[i, "BLOOD_SCORE"], chr_check[i, "STROMAL_SCORE"], chr_check[i, "FLS_SCORE"])
}

# remove chromatin states, so just the scores are left
chr_check <- cbind(chr_check[, 1:4], chr_check[8:ncol(chr_check)])