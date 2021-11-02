library(dplyr)
library(purrr)
library(readr)

files = Sys.glob("github/2021-cami-annot/outputs/eggnog_source_genomes/*.annotations")

eggnog <- files %>%
  set_names() %>%
  map_dfr(read_tsv, col_names = c('id', 'seed_ortholog', 'evalue', 'score', 
                                  'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 
                                  'Description', 'Preferred_name', 'GOs', 'EC',
                                  'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                                  'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
                                  'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'),
          col_types = "cccddccccccccccccccccc",
          comment = "#", na = "-", .id = "source_genome") %>%
  mutate(source_genome = gsub("outputs/eggnog_source_genomes/", "", source_genome)) %>%
  mutate(source_genome = gsub("\\.emapper\\.annotations", "", source_genome))

tmp <- eggnog %>%
  group_by(PFAMs) %>%
  tally()

# picked one of most represented, 
# middle, one that appeared once, and one that didn't appear at all.
