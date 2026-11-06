library(dplyr)
library(ggplot2)

lineage_colours <- c(
  DEEPOL1 = "cornflowerblue",
  DEEPOL2 = "mediumpurple3",
  DEEPOL3 = "mediumseagreen"
)

nreads_all <- bind_rows(
  DEEPOL1_meta %>% select(nreads) %>% mutate(sample = "DEEPOL1"),
  DEEPOL2_meta %>% select(nreads) %>% mutate(sample = "DEEPOL2"),
  DEEPOL3_meta %>% select(nreads) %>% mutate(sample = "DEEPOL3")
)

nreads_box <- ggplot(nreads_all, aes(x = sample, y = nreads, colour = sample)) +
  geom_boxplot() +
  scale_colour_manual(values = lineage_colours) +
  theme_minimal() +
  theme(legend.position = "none")
nreads_box

ggsave("/Users/ninpo556/Mirror/PhD/Archaeplastids/PacBio_orphans/paper/supplementary files/nreads_box.pdf", plot = nreads_box, width = 8, height = 6)



library(dplyr)
library(ggplot2)

lineage_colours <- c(
  DEEPOL1 = "cornflowerblue",
  DEEPOL2 = "mediumpurple3",
  DEEPOL3 = "mediumseagreen"
)

# Subsets: samples in the big dataframe that are not in each detection dataframe
DEEPOL1_nodet <- anti_join(metadata_all, DEEPOL1_meta, by = "sample")
DEEPOL2_nodet <- anti_join(metadata_all, DEEPOL2_meta, by = "sample")
DEEPOL3_nodet <- anti_join(metadata_all, DEEPOL3_meta, by = "sample")

# Detection + no-detection combined, with a detection label for the pairing
nreads_all <- bind_rows(
  DEEPOL1_meta  %>% select(nreads) %>% mutate(sample = "DEEPOL1", detection = "Detected"),
  DEEPOL1_nodet %>% select(nreads) %>% mutate(sample = "DEEPOL1", detection = "Not detected"),
  DEEPOL2_meta  %>% select(nreads) %>% mutate(sample = "DEEPOL2", detection = "Detected"),
  DEEPOL2_nodet %>% select(nreads) %>% mutate(sample = "DEEPOL2", detection = "Not detected"),
  DEEPOL3_meta  %>% select(nreads) %>% mutate(sample = "DEEPOL3", detection = "Detected"),
  DEEPOL3_nodet %>% select(nreads) %>% mutate(sample = "DEEPOL3", detection = "Not detected")
) %>%
  mutate(detection = factor(detection, levels = c("Detected", "Not detected")))

nreads_box <- ggplot(nreads_all, aes(x = sample, y = nreads, colour = sample, fill = detection)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_colour_manual(values = lineage_colours) +
  scale_fill_manual(values = c("Detected" = "white", "Not detected" = "grey85")) +
  guides(colour = "none") +
  labs(x = NULL, y = "Number of reads", fill = NULL) +
  theme_minimal()
nreads_box

ggsave("/Users/ninpo556/Mirror/PhD/Archaeplastids/PacBio_orphans/paper/supplementary files/nreads_box.pdf",
       plot = nreads_box, width = 8, height = 6)

# polar sequencing depth
polar_oceans <- c("Arctic Ocean", "Southern Ocean")  # adjust to match your labels

nreads_polar <- metadata_all %>%
  filter(ocean %in% polar_oceans)

nreads_polar_box <- ggplot(nreads_polar, aes(x = "Polar oceans", y = nreads)) +
  geom_boxplot(fill = "grey85", colour = "grey20", width = 0.4) +
  labs(x = NULL, y = "Number of reads") +
  theme_minimal()
nreads_polar_box

mean(nreads_polar$nreads)
median(nreads_polar$nreads)
min(nreads_polar$nreads)
max(nreads_polar$nreads)
