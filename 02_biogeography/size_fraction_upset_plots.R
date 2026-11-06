metadata_file <- here("data", "edited", "meta_with_presence_absence.csv")
message("[Data] Loading sample metadata:\n    ", metadata_file)
metadata_all <- read.csv(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE
)
cat("[OK] Loaded metadata: ", nrow(metadata_all), " rows\n")


# This block uses the prepare_lineage_data() helper to generate per-lineage filtered
# metadata and abundance (sample x ASV) matrices, for all downstream analyses.

filtered_data <- lapply(lineages, function(lin) {
  prepare_lineage_data(lin, metadata_all)
})
names(filtered_data) <- lineages

# Now for each lineage:
#    filtered_data[["DEEPOL1"]]$meta_f      # filtered meta for DEEPOL1
#    filtered_data[["DEEPOL1"]]$mat         # filtered PA matrix for DEEPOL1

# Print summary for each lineage
for(lin in lineages) {
  cat(sprintf("\n[%s] %d filtered samples, %d ASVs in matrix\n", lin,
              nrow(filtered_data[[lin]]$meta_f),
              ncol(filtered_data[[lin]]$mat)))
}

DEEPOL1 <- filtered_data[["DEEPOL1"]]
DEEPOL2 <- filtered_data[["DEEPOL2"]]
DEEPOL3 <- filtered_data[["DEEPOL3"]]

taxa_presence_by_size <- as.data.frame(DEEPOL3$mat) %>%
  mutate(Size = DEEPOL3$meta_f$size_fraction) %>%
  group_by(Size) %>%
  summarise(across(everything(), max)) %>% # 1 if present in at least one sample
  pivot_longer(-Size, names_to = "Taxon", values_to = "Present") %>%
  filter(Present > 0)

# 2. Convert to a list format for UpSetR
upset_list <- split(taxa_presence_by_size$Taxon, taxa_presence_by_size$Size)

# Create UpSet plot
png("/Users/ninpo556/Mirror/PhD/Archaeplastids/PacBio_orphans/EukBank/DEEPOL_Biogeography/results/DEEPOL3/DEEPOL3_size_upset.png", width = 2400, height = 1900, res = 300)
upset(fromList(upset_list), 
      order.by = "freq", 
      main.bar.color = "mediumseagreen", 
      sets.bar.color = "mediumseagreen",
      text.scale = c(3)) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
dev.off()
