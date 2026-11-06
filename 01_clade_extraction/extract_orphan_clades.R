#!/usr/bin/env Rscript
# ==============================================================================
#  extract_orphan_clades.R
#  Identify candidate deep-branching orphan clades in a phylogeny and combine
#  them with the clique assignments from the sequence similarity network (SSN).
# ==============================================================================
#
# BACKGROUND
#   Orphan sequences are sequences that could not be assigned to an annotated
#   eukaryotic supergroup. In the SSN, "cliques" are groups of orphan sequences
#   that are all connected to each other by significant similarity. Here the
#   same orphan sequences are placed in a phylogeny containing annotated
#   reference sequences, and we ask which orphans form clades of their own.
#
# DEFINITION USED
#   An "orphan clade" is a monophyletic group (all descendants of one internal
#   node of the tree) in which
#     (1) EVERY tip is an orphan sequence, i.e. no annotated sequence is nested
#         inside it, and
#     (2) at least MIN_ORPHANS tips are present (default: 2).
#   By default only MAXIMAL orphan clades are reported, i.e. clades that are not
#   contained in a larger orphan-only clade (set KEEP_NESTED <- TRUE to report
#   every orphan-only node, including nested ones).
#
# INPUT FILES (tab-separated text, with a header line)
#   ORPHAN_IDS_FILE : one column "name" with the IDs of all orphan sequences.
#   CLIQUES_FILE    : columns "name" and "clade"; the clique each orphan sequence
#                     belongs to in the SSN (renamed to "clique" below).
#   TREE_FILE       : rooted phylogeny in Newick format; tip labels must match
#                     the IDs in "name".
#
# OUTPUT
#   OUTPUT_FILE     : one row per orphan sequence with columns
#                       name, clique, clade, clade_size
#                     "clade" is NA for orphans that do not form an orphan-only
#                     clade with at least MIN_ORPHANS members.
#
# REQUIREMENTS: R packages ape, dplyr (versions are printed at the end).
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
})

# ---- 1. Settings --------------------------------------------------------------
# Paths are relative to the repository root (edit if your layout differs).
ORPHAN_IDS_FILE <- "data/orphan_IDs.txt"
CLIQUES_FILE    <- "data/cover40_orphan_cliques.tsv"   # TODO: state what "cover40" means (e.g. coverage cut-off)
TREE_FILE       <- "data/RAxML_bestTree.all.concat.53.rooted_round.nw"
OUTPUT_FILE     <- "results/orphan_tree_SSN_results.tsv"

MIN_ORPHANS <- 2       # minimum number of orphan tips for a clade to be reported
KEEP_NESTED <- FALSE   # FALSE: only the largest (maximal) orphan-only clades

# ---- 2. Read and check the input ---------------------------------------------
orphan_ids <- read.delim(ORPHAN_IDS_FILE, header = TRUE, stringsAsFactors = FALSE)
cliques    <- read.delim(CLIQUES_FILE,    header = TRUE, stringsAsFactors = FALSE)
tree       <- read.tree(TREE_FILE)

stopifnot("`name` column missing in orphan ID file" = "name" %in% names(orphan_ids))
stopifnot("`name` and `clade` columns needed in clique file" = all(c("name", "clade") %in% names(cliques)))
stopifnot("The tree must be rooted (clades are defined relative to the root)" = is.rooted(tree))

orphan_tips <- as.character(orphan_ids$name)

missing_tips <- setdiff(orphan_tips, tree$tip.label)
if (length(missing_tips) > 0) {
  warning(length(missing_tips), " orphan IDs are not present in the tree, e.g.: ",
          paste(head(missing_tips, 3), collapse = ", "))
}

# ---- 3. Helper functions ------------------------------------------------------
# Tip labels of all tips below a given node of the tree (the node itself if it is a tip)
get_tips_below <- function(tree, node) {
  n_tips <- length(tree$tip.label)
  tips   <- integer(0)
  to_visit <- node
  while (length(to_visit) > 0) {
    current  <- to_visit[1]
    to_visit <- to_visit[-1]
    if (current <= n_tips) {
      tips <- c(tips, current)                                   # reached a tip
    } else {
      to_visit <- c(to_visit, tree$edge[tree$edge[, 1] == current, 2])   # add the children
    }
  }
  tree$tip.label[tips]
}

# Parent node of a node (NA for the root)
get_parent <- function(tree, node) tree$edge[match(node, tree$edge[, 2]), 1]

# Find orphan-only clades; returns one row per orphan tip in such a clade
find_orphan_clades <- function(tree, orphan_tips, min_orphans = 2, keep_nested = FALSE) {
  n_tips         <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)          # node IDs of all internal nodes

  # Step 1: internal nodes whose descendants are ALL orphans (and at least min_orphans of them)
  is_orphan_only <- vapply(internal_nodes, function(node) {
    tips_below <- get_tips_below(tree, node)
    length(tips_below) >= min_orphans && all(tips_below %in% orphan_tips)
  }, logical(1))
  orphan_nodes <- internal_nodes[is_orphan_only]

  # Step 2: optionally drop nodes nested inside a larger orphan-only clade
  if (!keep_nested) {
    parents      <- vapply(orphan_nodes, function(node) get_parent(tree, node), numeric(1))
    orphan_nodes <- orphan_nodes[is.na(parents) | !(parents %in% orphan_nodes)]
  }

  # Step 3: one row per tip, with a clade label and the clade size
  if (length(orphan_nodes) == 0) {
    return(data.frame(name = character(), clade = character(), clade_size = integer()))
  }
  bind_rows(lapply(seq_along(orphan_nodes), function(i) {
    data.frame(name  = get_tips_below(tree, orphan_nodes[i]),
               clade = paste0("clade", i),
               stringsAsFactors = FALSE)
  })) %>%
    group_by(clade) %>%
    mutate(clade_size = n()) %>%
    ungroup() %>%
    as.data.frame()
}

# ---- 4. Run -------------------------------------------------------------------
orphan_clades <- find_orphan_clades(tree, orphan_tips, MIN_ORPHANS, KEEP_NESTED)

cat("Orphan-only clades found:", n_distinct(orphan_clades$clade), "\n")
cat("Orphan sequences placed in these clades:", nrow(orphan_clades),
    "of", length(orphan_tips), "\n")
print(table(orphan_clades$clade_size))   # how many clades of each size

# ---- 5. Combine with the SSN cliques and write the result ---------------------
result <- orphan_ids %>%
  left_join(rename(cliques, clique = clade), by = "name") %>%   # SSN clique of each orphan
  left_join(orphan_clades, by = "name")                         # clade in the tree

if (nrow(result) != nrow(orphan_ids)) {
  warning("Some sequences appear more than once (overlapping cliques, or nested clades if KEEP_NESTED = TRUE).")
}

dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
write.table(result, file = OUTPUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Written to:", OUTPUT_FILE, "\n")

# Record the package versions used
cat("R", R.version.string, "| ape", as.character(packageVersion("ape")),
    "| dplyr", as.character(packageVersion("dplyr")), "\n")
