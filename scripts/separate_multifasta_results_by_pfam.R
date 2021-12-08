library(dplyr)
library(readr)

if (!dir.exists(snakemake@params[['outdir']])) {dir.create(snakemake@params[['outdir']])}

multifasta_results <- read_csv(snakemake@input[['cdbg_annot']])

multifasta_results %>%
  mutate(pfam = gsub("/home/tereiter/github/2021-pfam-shared-kmers/outputs/pfam_fastas/",
                     "", filename)) %>%
  mutate(pfam = gsub(".fa", "", pfam)) %>%
  select(pfam, cdbg_id) %>%
  group_by(pfam) %>%
  group_walk(~ write_tsv(.x, paste0(snakemake@params[['outdir']], .y$pfam, "_cdbg_nodes.tsv"), col_names = F))
