# set the download timeout to 30 minutes
options(timeout = max(1800, getOption("timeout")))

# load the required R libraries
library(tidyverse)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(readr)
library(readxl)
library(gwascat)

# connect to a BioMart database hosted by Ensembl and store in "snp_mart"
snp_mart <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# read the RA SNPS.xlsx file as "ra_snps
ra_snps <- read_excel("ra_snps_data/RA SNPs.xlsx")

# remove the unnecessary rows and columns
ra_snps <- ra_snps[6:529, c(3, 4, 6, 7, 8, 13)]

# change the column names of "ra_snps"
colnames(ra_snps) <- c("LEAD_ID", "LEAD_RS", "CRED_SET", "VARIANT_ID", "RS_ID", "PIP")

# create empty vectors for coord and alleles for all SNPs
coord_total <- c()
alleles_total <- c()

# loop through all ra SNPs to extract and store their details
for(i in 1:nrow(ra_snps)) {
  # store the details of the current SNP as "snp_det"
  snp_det <- str_split(ra_snps$VARIANT_ID[i], ":")[[1]]
  # store the chromosome of the current SNP as "snp_chr"
  snp_chr <- snp_det[1]
  # store the locus of the current SNP as "snp_loc"
  snp_loc <- snp_det[2]
  # store the reference allele as "ref_all"
  ref_all <- snp_det[3]
  # store the alternative allele as "alt_all"
  alt_all <- snp_det[4]
  # paste the SNP chromosome and locus together to make "coord"
  coord <- paste("chr", snp_chr, ":", snp_loc, sep = "")
  # add the current "coord" to "coord_total"
  coord_total <- append(coord_total, coord)
  # paste the two refence and alternative alleles together to make "alleles"
  alleles <- paste("(", ref_all, "/", alt_all, ")", sep = "")
  # add the current "alleles" to "alleles_total"
  alleles_total <- append(alleles_total, alleles)
}

# add the coord and allele data to "ra_snps"
ra_snps[, "COORD"] <- coord_total
ra_snps[, "ALLELES"] <- alleles_total

# remove the unnecessary objects from the environment
rm(alleles, alleles_total, alt_all, coord, coord_total, i, ref_all, snp_chr, snp_det, snp_loc)

# create an empty vector for the RA SNPs GRanges object
ra_snps_granges <- c()

# loop through each row in "ra_snps"
for(i in 1:nrow(ra_snps)) {
  # store the current row of "ra_snps"
  row <- ra_snps[i, ]
  # convert the row to a GRanges object by the column "COORD"
  ra_snp <- as(as.character(row["COORD"]), "GRanges")
  # add the metadata from "ra_snps" to the GRanges object
  values(ra_snp) <- DataFrame(VARIANT_ID = row$VARIANT_ID, VARIANT_RS = row$RS_ID, CRED_SET = row$CRED_SET,
                              ALLELES = row$ALLELES, PIP = row$PIP, LEAD_ID = row$LEAD_ID, LEAD_RS = row$LEAD_RS)
  # add the GRanges object for the current row to the larger object for the whole dataset
  ra_snps_granges <- append(ra_snps_granges, ra_snp)
  
  # update progress to console with a loading bar
  progress <- round((i/nrow(ra_snps)) * 100, 0)
  load_bar <- paste(sample("■", progress, replace = T), collapse = "")
  cat("\014")
  message(paste("RA SNPs: GRanges ", load_bar, " ", as.character(progress), "%", sep = ""))
}

# store the full GRanges object as "ra_snps"
ra_snps <- ra_snps_granges

# remove the unnecessary objects from the environment
rm(i, load_bar, progress, ra_snp, ra_snps_granges, row)

# import chain file for converting human genome 19 to human genome 38
ch <- import.chain("hg19ToHg38.over.chain")

# rename the seqnames of "ra_snps" to be in the "UCSC" style
seqlevelsStyle(ra_snps) <- "UCSC"

# lift the intervals from build hg19 to hg38 and store as "ra_snps_38"
ra_snps_38 <- liftOver(ra_snps, ch)

# flatten the list
ra_snps_38 <- unlist(ra_snps_38)

# set the genome to be registered as "hg38"
genome(ra_snps_38) = "hg38"

# generate the object "ra_snps_38" from the class "gwasloc"
ra_snps_38 = new("gwaswloc", ra_snps_38)

# store the hg38 object as "ra_snps"
ra_snps <- ra_snps_38

# remove the unnecessary objects from the environment
rm(ra_snps_38, ch)

# convert the GRanges object to a data frame
ra_snps <- as.data.frame(ra_snps)

# remove the rows with a "LEAD_RS" that does not begin with "rs"
ra_snps <- ra_snps[which(str_starts(ra_snps$LEAD_RS, "rs")), ]
