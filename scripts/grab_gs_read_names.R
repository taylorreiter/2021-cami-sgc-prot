library(readr)
library(dplyr)

# gs = "/home/tereiter/github/2021-cami-annot/outputs/gs_read_annotations/CAMI_low_gs_read_annotations_with_eggnog.tsv",
# pfam_map = "inputs/Pfam-A.clans.tsv.gz"
# pfam_sgc = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.read_names.txt"

#pfam_map <- read_tsv("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz",
pfam_map <- read_tsv(snakemake@input[['pfam_map']],
                     col_names = c("pfam", "cl", "full_name", "name", "description"))

# get PFAM name for PFAM identifier
pfam_ident <- snakemake@wildcards[['pfam']]
#pfam_ident = "PF01297.19"
pfam_ident <- gsub("\\..*", "", pfam_ident)
pfam_name <- filter(pfam_map, pfam == pfam_ident)$name

# gs <- read_tsv("~/github/2021-cami-annot/outputs/osf_products/CAMI_low_gs_read_annotations_with_eggnog.tsv.gz")
gs <- read_tsv(snakemake@input[['gs']])

gs_single_pfam <- gs %>%
  filter(grepl(pfam_name, PFAMs))

write_tsv(gs_single_pfam, snakemake@output[['pfam_gs']])
