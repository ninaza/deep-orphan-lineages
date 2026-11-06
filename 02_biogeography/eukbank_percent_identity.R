library(tidyverse)

# Standard BLAST tabular output (-outfmt 6). If you used a custom format, change the names.
blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# One BLAST file per lineage (query = long-read DEEPOL sequences, subject = EukBank ASVs)
blast_files <- c(
  DEEPOL1 = "/Users/ninpo556/Mirror/PhD/Archaeplastids/PacBio_orphans/DEEPOL1/tree/cladeA_EukBank.out",
  DEEPOL2 = "/Users/ninpo556/Mirror/PhD/Archaeplastids/PacBio_orphans/DEEPOL2/tree/cladeB_EukBank.out",
  DEEPOL3 = "/Users/ninpo556/Mirror/PhD/Archaeplastids/PacBio_orphans/DEEPOL3/tree/cladeC_EukBank.out"
)

best_hits <- imap_dfr(blast_files, function(f, lin) {
  # ASVs that were retained after the phylogenetic validation
  kept <- read.table(file.path("data", "raw", paste0(lin, "_EukBank.list")),
                     header = FALSE, col.names = "amplicon")$amplicon
  
  read.delim(f, header = FALSE, col.names = blast_cols) %>%
    filter(sseqid %in% kept) %>%                  # only the final ASVs of this lineage
    group_by(sseqid) %>%                          # one row per ASV (it can hit several long reads)
    slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%   # keep the best hit
    ungroup() %>%
    mutate(lineage = lin)
})

pident_summary <- best_hits %>%
  group_by(lineage) %>%
  summarise(n_asv = n(),
            mean_pident   = mean(pident),
            sd_pident     = sd(pident),
            median_pident = median(pident),
            min_pident    = min(pident),
            max_pident    = max(pident),
            .groups = "drop")

print(pident_summary)
write.csv(pident_summary, file.path("results", "blast_pident_summary.csv"), row.names = FALSE)
