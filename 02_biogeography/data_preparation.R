library(tidyverse)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(purrr)

# ============================================================
# CONFIGURATION: add / remove lineages here only
# (renamed from `lineages` to `lineage_cfg`, so it cannot clash with the character
#  vector `lineages` used by the coastal analysis script when both run in one session)
# ============================================================
lineage_cfg <- list(
  DEEPOL1 = list(
    asv_list = "data/raw/DEEPOL1_EukBank.list",
    out_dir  = "data/edited/DEEPOL1"
  ),
  DEEPOL2 = list(
    asv_list = "data/raw/DEEPOL2_EukBank.list",
    out_dir  = "data/edited/DEEPOL2"
  ),
  DEEPOL3 = list(
    asv_list = "data/raw/DEEPOL3_EukBank.list",
    out_dir  = "data/edited/DEEPOL3"
  )
)

# An ASV must be detected (nreads > 0) in at least this many samples to be retained.
# NOTE: the old code was `sum(nreads > 0) > 2`, i.e. at least 3 samples, while the comment said 2.
# The value is now 3 with `>=`, which reproduces the old numbers exactly. Set to 2 if you meant 2.
MIN_SAMPLES <- 3

POLAR_SOUTH_CUTOFF <- 60

ocean_polygon_file <- "data/World_Seas_IHO_v3/World_Seas_IHO_v3.shp"
snap_max_deg       <- 1           # points outside every polygon are snapped to the nearest one
# if within this many degrees (~100 km), otherwise left NA

# ============================================================
# METADATA PREPARATION (runs once, shared across all lineages)
# ============================================================
metadata_raw <- read.delim("data/raw/eukbank_18S_V4_samples.tsv",
                           header = TRUE, dec = ".")

# Per-sample sequencing depth = all reads in a sample. It is renamed to `sequencing_depth`
# everywhere, so it cannot be confused with `nreads` in the counts table (reads of ONE ASV in
# ONE sample). The old metadata column (`total_nreads` or `nreads`) is replaced by it.
depth_col <- intersect(c("total_nreads", "nreads"), names(metadata_raw))[1]
if (is.na(depth_col)) {
  stop("No per-sample depth column found. Columns containing 'read': ",
       paste(grep("read", names(metadata_raw), value = TRUE, ignore.case = TRUE), collapse = ", "))
}
cat("Per-sample depth taken from metadata column:", depth_col, "\n")
metadata_raw <- metadata_raw %>%
  mutate(sequencing_depth = as.numeric(.data[[depth_col]])) %>%
  select(-any_of(c("total_nreads", "nreads")))

marine_envplots <- c("marine_sediment", "marine_water", "marine_organism", "marine_ice")

# ============================================================
# OCEAN ASSIGNMENT
# ============================================================
# --- Legacy box/latitude method, kept ONLY to compare against the polygon method ---
assign_ocean_basin_legacy <- function(lat, lon) {
  lon <- ifelse(lon > 180, lon - 360, lon)
  case_when(
    lat >  66.5                                          ~ "Arctic Ocean",
    lat < -40                                            ~ "Southern Ocean",
    between(lat, 30, 46)  & between(lon, -6,  37)       ~ "Mediterranean Region",
    between(lat, 53, 66)  & between(lon, 12,  30)       ~ "Baltic Sea",
    between(lat, 36, 48)  & between(lon, 49,  55)       ~ NA_character_,
    between(lat, -40, 30) & between(lon,  20, 140)      ~ "Indian Ocean",
    lat >= 0  & (lon >= 100 | lon <= -70)               ~ "North Pacific Ocean",
    lat <  0  & (lon >= 140 | lon <= -70)               ~ "South Pacific Ocean",
    lat >= 0  & between(lon, -80, 20)                   ~ "North Atlantic Ocean",
    lat <  0  & between(lon, -70, 20)                   ~ "South Atlantic Ocean",
    TRUE                                                 ~ "Unassigned"
  )
}

# --- Polygon method ----------------------------------------------------------
get_ocean_polygons <- function() {
  if (!is.null(ocean_polygon_file) && file.exists(ocean_polygon_file)) {
    polys <- st_read(ocean_polygon_file, quiet = TRUE)
  } else {
    polys <- ne_download(scale = 10, type = "geography_marine_polys",
                         category = "physical", returnclass = "sf")
  }
  name_col <- intersect(c("NAME", "name", "name_en"), names(polys))[1]
  if (is.na(name_col)) stop("No name column found in ocean polygons: ", paste(names(polys), collapse = ", "))
  polys <- polys %>% st_transform(4326) %>% st_make_valid()
  polys$poly_name <- as.character(polys[[name_col]])
  polys[, "poly_name"]
}

# Regex rules mapping polygon names (seas, gulfs, oceans) to the basin labels used downstream.
# Edit these if the "unmapped polygon names" table printed below shows something you want to move.
rx_arctic   <- "^arctic ocean|barents|kara sea|laptev|east siberian|chukchi|beaufort|lincoln sea|greenland sea|baffin|hudson|foxe|white sea|amundsen gulf|lancaster|viscount melville|mcclure|nares|north ?west(ern)? passage"
rx_southern <- "southern ocean|weddell|ross sea|amundsen sea|bellingshausen"
rx_baltic   <- "baltic|bothni|gulf of finland|gulf of riga|kattegat|great belt|little belt|belt sea|oresund|øresund|archipelago sea|aland"
rx_med      <- "mediterranean|adriatic|ionian|aegean|tyrrhenian|ligurian|alboran|balearic|iberian sea|marmara|black sea|azov|levantine|libyan|sicil|sardinia|sea of crete|strait of gibraltar"
rx_indian   <- "indian ocean|arabian|bay of bengal|andaman|red sea|gulf of aden|persian|gulf of oman|laccadive|lakshadweep|mozambique|malacca|singapore|great australian bight|timor|bali sea|savu|gulf of mannar|palk|zanzibar"
rx_pacific  <- "pacific|south china|east(ern)? china|british columbia|alaska|yellow sea|japan|okhotsk|bering|philippine|celebes|sulu|banda|java sea|flores|molucca|halmahera|seram|arafura|coral sea|tasman|solomon|bismarck|thailand|tonkin|bohai|taiwan|korea|gulf of alaska|carpentaria|bass strait|salish|puget|california|panama|humboldt"
rx_atlantic <- "atlantic|north sea|norwegian sea|irish sea|english channel|biscay|celtic|labrador|davis strait|lawrence|fundy|gulf of maine|caribbean|gulf of mexico|sargasso|gulf of guinea|hebrides|inner seas|bristol channel|wadden|skagerrak|faroe|iceland|denmark strait|rio de la plata|gulf of cadiz|scotia sea|drake passage|falkland|patagon"

basin_from_name <- function(poly_name, lat) {
  n    <- tolower(poly_name)
  hemi <- ifelse(lat >= 0, "North", "South")   # IHO draws the Atlantic/Pacific split at the equator
  case_when(
    is.na(n)                  ~ NA_character_,
    grepl(rx_arctic,   n)     ~ "Arctic Ocean",
    grepl(rx_southern, n)     ~ "Southern Ocean",
    grepl(rx_baltic,   n)     ~ "Baltic Sea",
    grepl(rx_med,      n)     ~ "Mediterranean Region",
    grepl(rx_indian,   n)     ~ "Indian Ocean",
    grepl(rx_pacific,  n)     ~ paste(hemi, "Pacific Ocean"),
    grepl(rx_atlantic, n)     ~ paste(hemi, "Atlantic Ocean"),
    TRUE                      ~ "Unassigned"
  )
}

assign_ocean_polygon <- function(df) {
  old_s2 <- sf_use_s2()
  suppressMessages(sf_use_s2(FALSE))             # planar GEOS: robust for these polygons
  on.exit(suppressMessages(sf_use_s2(old_s2)), add = TRUE)
  
  polys <- get_ocean_polygons()
  poly_area <- as.numeric(st_area(st_set_crs(polys, NA)))   # planar, only used to rank overlaps
  
  coords <- df %>%
    distinct(longitude, latitude) %>%
    mutate(lon_n = ifelse(longitude > 180, longitude - 360, longitude))
  pts <- st_as_sf(coords, coords = c("lon_n", "latitude"), crs = 4326, remove = FALSE)
  
  # smallest polygon containing the point (most specific sea wins where polygons overlap)
  hits <- st_intersects(pts, polys)
  pick <- vapply(hits, function(ix) if (length(ix) == 0) NA_integer_ else ix[which.min(poly_area[ix])],
                 integer(1))
  coords$ocean_snapped <- FALSE
  
  # points outside every polygon (coastline rounding, estuaries): snap to the nearest if close
  miss <- which(is.na(pick))
  if (length(miss) > 0) {
    nn <- st_nearest_feature(pts[miss, ], polys)
    d  <- as.numeric(st_distance(st_set_crs(pts[miss, ], NA), st_set_crs(polys[nn, ], NA),
                                 by_element = TRUE))
    ok <- d <= snap_max_deg
    pick[miss[ok]] <- nn[ok]
    coords$ocean_snapped[miss[ok]] <- TRUE
  }
  
  coords$poly_name <- polys$poly_name[pick]
  coords$ocean     <- basin_from_name(coords$poly_name, coords$latitude)
  coords %>% select(longitude, latitude, poly_name, ocean, ocean_snapped)
}

marine_samples <- metadata_raw %>%
  filter(envplot %in% marine_envplots) %>%
  filter(!is.na(latitude), !is.na(longitude))

ocean_lookup <- assign_ocean_polygon(marine_samples)

marine_samples <- marine_samples %>%
  select(-any_of(c("ocean", "poly_name", "ocean_snapped"))) %>%   
  left_join(ocean_lookup, by = c("longitude", "latitude")) %>%
  mutate(ocean_legacy = assign_ocean_basin_legacy(latitude, longitude))

# ---- Ocean assignment QC ----
cat("\nPolygon names not mapped to a basin (edit the rx_* patterns above if needed):\n")
print(marine_samples %>% filter(ocean == "Unassigned") %>% count(poly_name, sort = TRUE))

cat("\nMarine samples outside all ocean polygons (ocean = NA):", sum(is.na(marine_samples$ocean)), "\n")
cat("Marine samples snapped to nearest polygon:", sum(marine_samples$ocean_snapped, na.rm = TRUE), "\n")

cat("\nOld (box/latitude) vs new (polygon) labels:\n")
print(table(old = marine_samples$ocean_legacy, new = marine_samples$ocean, useNA = "ifany"))

cat("\nSamples labelled 'Southern Ocean' per method, by latitude band:\n")
print(marine_samples %>%
        mutate(lat_band = cut(latitude, c(-91, -66.5, -60, -50, -40, 0, 91))) %>%
        group_by(lat_band) %>%
        summarise(old_southern = sum(ocean_legacy == "Southern Ocean", na.rm = TRUE),
                  new_southern = sum(ocean == "Southern Ocean", na.rm = TRUE),
                  n = n(), .groups = "drop"))

dir.create("data/edited", recursive = TRUE, showWarnings = FALSE)
write.csv(marine_samples %>%
            select(sample, latitude, longitude, envplot, poly_name, ocean, ocean_snapped, ocean_legacy),
          "data/edited/ocean_assignment_qc.csv", row.names = FALSE)

metadata_raw <- metadata_raw %>%
  select(-any_of(c("ocean", "sea_name", "polar_group"))) %>%   # so the join cannot create ocean.x / ocean.y
  left_join(marine_samples %>% select(sample, ocean, sea_name = poly_name), by = "sample") %>%
  mutate(polar_group = case_when(
    is.na(ocean) | is.na(latitude)                              ~ NA_character_,
    ocean == "Arctic Ocean" | latitude <= -POLAR_SOUTH_CUTOFF   ~ "polar",
    TRUE                                                         ~ "non-polar"
  ))
cat("\npolar_group (marine samples):\n"); print(table(metadata_raw$polar_group, useNA = "ifany"))

# ============================================================
# HABITAT & SIZE FRACTION PROCESSING
# ============================================================
metadata <- metadata_raw %>%
  mutate(habitat = case_when(
    envplot %in% marine_envplots                                         ~ "Marine",
    envplot %in% c("land_freshwater")                                    ~ "Freshwater",
    envplot %in% c("land_soil", "land_sediment",
                   "land_organism", "land_water")                        ~ "Terrestrial",
    envplot %in% c("none")                                               ~ "None"
  )) %>%
  mutate(
    size_fraction_lower_threshold = as.numeric(size_fraction_lower_threshold),
    size_fraction_upper_threshold = as.numeric(
      ifelse(is.na(size_fraction_upper_threshold), Inf, size_fraction_upper_threshold)
    ),
    size_fraction = case_when(
      size_fraction_lower_threshold <= 0.22 & size_fraction_upper_threshold == Inf ~ "Bulk_Total",
      size_fraction_lower_threshold >= 0.2  & size_fraction_upper_threshold <= 3   ~ "Pico",
      size_fraction_lower_threshold >= 3    & size_fraction_upper_threshold <= 20  ~ "Nano",
      # upper bound made inclusive (<= 200): the standard 20-200 um fraction was previously
      # excluded by `< 200` and ended up as Non_Standard
      size_fraction_lower_threshold >= 20   & size_fraction_upper_threshold <= 200 ~ "Micro",
      size_fraction_lower_threshold >= 200  & size_fraction_upper_threshold <  Inf ~ "Meso",
      TRUE ~ "Non_Standard"
    )
  ) %>%
  mutate(salinity_cat = case_when(
    is.na(salinity)                   ~ "unknown",
    salinity <= 0.5                   ~ "freshwater",
    salinity > 0.5  & salinity <= 30  ~ "brackish",
    salinity > 30   & salinity <= 50  ~ "saline",
    salinity > 50                     ~ "briny"
  )) %>%
  mutate(
    ocean_layer_fine = case_when(
      depth <= 200                ~ "photic",
      depth > 200 & depth <= 1000 ~ "twilight",
      depth > 1000                ~ "aphotic",
      TRUE                        ~ NA_character_
    ),
    ocean_layer = case_when(
      depth <= 200 ~ "sunlit",
      depth > 200  ~ "dark",
      TRUE         ~ NA_character_
    )
  )

# ---- Metadata QC (printed only, nothing is changed here) ----
cat("\n===== METADATA QC =====\n")
cat("Duplicated sample IDs:", sum(duplicated(metadata$sample)), "(must be 0, otherwise joins duplicate rows)\n")
cat("Latitude outside [-90,90] or longitude outside [-180,360]:",
    sum(abs(metadata$latitude) > 90 | metadata$longitude < -180 | metadata$longitude > 360, na.rm = TRUE), "\n")
cat("Samples at exactly 0,0 (likely missing coordinates):",
    sum(metadata$latitude == 0 & metadata$longitude == 0, na.rm = TRUE), "\n")
cat("Marine samples without coordinates (no ocean label possible):",
    sum(metadata$habitat == "Marine" & (is.na(metadata$latitude) | is.na(metadata$longitude)), na.rm = TRUE), "\n")
cat("Negative depths (would be called 'sunlit'):", sum(metadata$depth < 0, na.rm = TRUE), "\n")
cat("Marine samples by envplot (sediment depth is usually seafloor depth, not sampling depth):\n")
print(table(metadata$envplot[metadata$habitat == "Marine"], useNA = "ifany"))
cat("Marine samples by ocean_layer x envplot:\n")
print(table(metadata$envplot[metadata$habitat == "Marine"],
            metadata$ocean_layer[metadata$habitat == "Marine"], useNA = "ifany"))

cat("\nsize_fraction counts (marine):\n")
print(table(metadata$size_fraction[metadata$habitat == "Marine"], useNA = "ifany"))
cat("Most common lower/upper thresholds among marine 'Non_Standard' samples:\n")
print(metadata %>%
        filter(habitat == "Marine", size_fraction == "Non_Standard") %>%
        count(size_fraction_lower_threshold, size_fraction_upper_threshold, sort = TRUE) %>%
        head(15))

cat("\nSamples with NA or 0 sequencing_depth (relative abundance would be NA/Inf):",
    sum(is.na(metadata$sequencing_depth) | metadata$sequencing_depth <= 0), "\n")

# ============================================================
# LOAD FULL EUKBANK COUNTS
# ============================================================
eukbank_counts <- read.delim("data/raw/eukbank_18S_V4_counts.tsv",
                             header = TRUE, dec = ".")

# ============================================================
# PER-LINEAGE PROCESSING LOOP
# ============================================================
for (lineage_name in names(lineage_cfg)) {
  cfg <- lineage_cfg[[lineage_name]]
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  
  target_list <- read.table(cfg$asv_list, header = FALSE, col.names = "amplicon")
  
  lineage_counts <- eukbank_counts %>%
    filter(amplicon %in% target_list$amplicon)
  
  lineage_abundance <- lineage_counts %>%
    left_join(metadata %>% select(sample, sequencing_depth), by = "sample") %>%
    mutate(rel_abundance = nreads / sequencing_depth, rel_abundance_pct = rel_abundance * 100)
  
  cat(sprintf("\n%s: %d of %d listed ASVs found in counts; %d rows without sequencing_depth\n",
              lineage_name, n_distinct(lineage_counts$amplicon), nrow(target_list),
              sum(is.na(lineage_abundance$sequencing_depth))))
  
  out_abundance <- file.path(cfg$out_dir, paste0(lineage_name, "_abundance.csv"))
  write.csv(lineage_abundance, file = out_abundance, row.names = FALSE)
  
  lineage_filtered <- lineage_abundance %>%
    group_by(amplicon) %>%
    filter(sum(nreads > 0) >= MIN_SAMPLES) %>%
    ungroup()
  
  out_filtered <- file.path(cfg$out_dir, paste0(lineage_name, "_filtered_abundance.csv"))
  write.csv(lineage_filtered, file = out_filtered, row.names = FALSE)
  
  common_samples   <- intersect(lineage_filtered$sample, metadata$sample)
  metadata_lineage <- metadata[metadata$sample %in% common_samples, ]
  
  out_meta <- file.path(cfg$out_dir, paste0(lineage_name, "_samples_edited.csv"))
  write.csv(metadata_lineage, file = out_meta, row.names = FALSE)
}

# ============================================================================
# BUILD PRESENCE/ABSENCE TABLE
# (built from the UNFILTERED abundance files, so it includes ASVs seen in < MIN_SAMPLES samples;
#  the tree/abundance analyses use the filtered files. Make sure the analysis script uses the
#  same source for detection and for abundance.)
# ============================================================================
cat("\nBuilding presence/absence table...\n")

pa_list <- imap(lineage_cfg, function(cfg, name) {
  path <- file.path(cfg$out_dir, paste0(name, "_abundance.csv"))
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  
  read.csv(path, stringsAsFactors = FALSE) %>%
    filter(nreads > 0) %>%
    distinct(sample) %>%
    mutate(!!name := 1L)
})

pa_table <- metadata %>%
  select(sample) %>%
  reduce(pa_list, .init = ., function(acc, df) left_join(acc, df, by = "sample")) %>%
  mutate(across(all_of(names(lineage_cfg)), ~ replace_na(.x, 0L)))

cat(sprintf("  Presence/absence table: %d samples x %d lineages\n",
            nrow(pa_table), length(lineage_cfg)))

cat("  Detections per lineage:\n")
walk(names(lineage_cfg), function(nm) {
  cat(sprintf("    %s: %d samples (%.1f%%)\n", nm, sum(pa_table[[nm]]), 100 * mean(pa_table[[nm]])))
})

meta_pa <- metadata %>%
  left_join(pa_table, by = "sample")

write.csv(meta_pa,
          file      = "data/edited/meta_with_presence_absence.csv",
          row.names = FALSE)

cat("\nAll lineages processed and Presence/Absence table saved.\n")
