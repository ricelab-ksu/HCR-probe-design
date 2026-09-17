#!/usr/bin/env Rscript
# ============================================================================
# HCR PROBE DESIGNER (v41 - cross-species edition)
# ============================================================================
# Designs HCR v3.0 probe sets (amplifiers B1-B5 by default; v2/v3 selectable)
# from a gene-list CSV.
# Probe construction replicates the lab's "HCR PROBE MAKER Rice lab" sheet:
#   52-bp target = two 25-mers + 2-bp gap; each oligo = 18-nt initiator half
#   + 2-nt spacer + 25-nt antisense half (45 nt total, 90 bp per pair).
# The amplifier is NOT a user choice: the engine rotates B1-B5 per gene for
# the oligo naming, while the primary deliverable (HCR PROBE MAKER all-amp
# CSV) lists every amplifier in the SELECTED amplifier set for each 52-bp
# target. The lab's in-situ kit is HCR v2.0, so the set defaults to B1-B5;
# selecting "HCR v3.0" adds B7/B9/B10/B13-B15/B17. "Download all (ZIP)"
# bundles the complete output set (including HCR_design_log.txt).
#
# NEW IN v34: CROSS-SPECIES CHECK
#   To reuse the same probe set in other species, enter a comma-separated
#   species list in the sidebar ("Check probes in other species"). For each
#   designed probe pair, the app looks up the same gene in each species'
#   transcriptome (local rna.fna, or NCBI if the fallback is on), finds the
#   best ungapped match of each 25-nt binding half, and reports:
#     - identity of each half (e.g. 24/25), with mismatched bases highlighted
#     - whether the 2-bp gap between the two halves is preserved
#     - a verdict per pair: "exact" (both halves 25/25, gap = 2 bp),
#       "near" (<=2 mismatches per half), or "risky" (>2 mismatches).
#   When the other species' genomic.fna is also present, the panel additionally
#   draws a gene-model track plot of the ortholog transcripts (aligned to that
#   species' genome) with the designed probes' binding sites overlaid - the
#   same view as the design-species isoform plot.
#   NOTE: the initiator/spacer parts of the oligos are species-independent,
#   so only the two 25-nt antisense halves need to match. A "near" pair may
#   still work (HCR tolerates a few mismatches) - validate empirically.
#   Non-melanogaster RefSeq files index genes by protein name/LOC ID, so if
#   an ortholog is not found, try the full gene/protein name.
#
# SETUP
#   1. Install R (>= 4.1) and the shiny package: install.packages("shiny")
#   2. Keep this file where it is:
#         /Users/ricelab/Desktop/bioinformatics/HCR_probe_maker/HCR probe design.R
#      It automatically finds the transcript files one folder up:
#         /Users/ricelab/Desktop/bioinformatics/reference_files/
#            d_melanogaster/rna.fna
#            d_simulans/rna.fna
#            d_eugracilis/rna.fna
#      The app designs probes from the transcript sequences (rna.fna). If the
#      matching genomic.fna and genomic.gff/gtf are also present, each gene gets
#      an isoform-structure track plot, and the "shared exons" option can design
#      against exons common to all isoforms (shared CDS when annotated).
#      (The "species" column in your CSV should be e.g. "melanogaster",
#      "simulans", "eugracilis", "erecta", "ananassae", "teissieri".)
#
#      Nothing genome-sized ships with the app any more: each species' reference
#      set is downloaded from NCBI on first use into d_<species> (lazily, from
#      the "D. melanogaster reference" / "Add a species" panels or automatically
#      on the first run that needs it). A species that is already on disk is
#      detected automatically: any subfolder named d_<species> containing an
#      rna.fna(.gz). To add one by hand, download its RefSeq files into a new
#      subfolder, e.g.:
#        d_erecta:      GCF_003286155.1_DereRS2
#        d_ananassae:   GCF_017639315.1_ASM1763931v2
#        d_teissieri:   GCF_016746235.2_Prin_Dtei_1.1
#        d_virilis:     GCF_003285735.1_DvirRS2
#      from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/... (grab rna.fna.gz,
#      plus genomic.fna.gz and genomic.gff.gz for structure plots/shared CDS;
#      the app stores the NCBI annotation under its native genomic.gff.gz name
#      and reads both GTF and GFF3 attribute styles).
#      You can also DESIGN probes directly in any of these species: put the
#      species name in the CSV's species column (gene = protein name, LOC ID,
#      or accession) and uncheck all "Check probes in other species" boxes for
#      a single-species run with no cross-species check.
#   3a. Interactive app:  Rscript "HCR probe design.R"   (a browser window opens)
#   3b. Headless batch:   Rscript "HCR probe design.R" my_genes.csv
#       (flags: --cross to enable the cross-species check, --no-cross to
#       disable it; default is OFF, matching the interactive UI. Settings
#       come from default_settings(), overridable via the HCR_RELAX_FILTERS
#       and HCR_STRINGENCY environment variables.)
#      No browser/server needed - designs all genes in the CSV and writes an
#      HCR_batch_<timestamp>/ folder with the oligo order, all-amplifier table,
#      cross-species results, pool summary, run parameters, and plots PDF.
#      This is the recommended way to run BIG gene lists: it is immune to
#      browser/screen-lock issues. To keep your Mac awake and let it run
#      unattended (log output to a file, detach from the terminal):
#        caffeinate -i nohup Rscript "HCR probe design.R" my_genes.csv > run.log 2>&1 &
#      then check progress any time with:  tail -f run.log
#      (caffeinate prevents sleep; nohup ... & detaches it. Locking the
#      screen is then completely safe.)
#
# NOTE on freezes: all NCBI network calls have hard timeouts (30 s per
# request, 10 min for reference-genome downloads). A stalled connection now
# becomes a per-gene error message instead of hanging the whole run.
#
# Both FlyBase-style headers (name=dsx-RA; parent=FBgn...,dsx) and NCBI-style
# headers (... doublesex (dsx), transcript variant A, mRNA) are parsed
# automatically. Non-melanogaster RefSeq files annotate genes by protein
# name/LOC ID rather than symbol, so for those species use the full gene name
# (e.g. "doublesex", "retinin"), an accession in the transcript column, or a
# pasted sequence.
#
# If a gene is not found locally, the app can fall back to fetching the RefSeq
# transcript from NCBI over the internet (checkbox in the sidebar).
# ============================================================================

library(shiny)

PORT <- 7788
APP_VERSION <- "42.2.0"
APP_TITLE <- "HCR Probe Designer"

# Reference file directory (contains d_<species>/ subfolders with rna.fna).
# Resolution order: HCR_REFERENCE_DIR env var → original author default →
#   ~/Downloads/reference_files (common lab fallback). Set the env var to
#   override on any machine.
DATA_DIRS <- local({
  env <- Sys.getenv("HCR_REFERENCE_DIR", unset = "")
  if (nzchar(env)) return(env)
  candidates <- c("/Users/ricelab/Desktop/bioinformatics/reference_files",
                  file.path(path.expand("~"), "Downloads", "reference_files"))
  hit <- candidates[vapply(candidates, function(d) dir.exists(d), TRUE)]
  if (length(hit)) hit[1] else candidates[1]
})
if (!dir.exists(DATA_DIRS)) {
  message("WARNING: reference directory not found: ", DATA_DIRS,
          "\n         Set HCR_REFERENCE_DIR or install rna.fna files.")
} else {
  message("Using reference path: ", DATA_DIRS)
}

# Species available for the cross-species check: subfolders of DATA_DIRS named
# d_<species> that contain an rna.fna (or rna.fna.gz). Detected at startup so
# the UI can offer them as checkboxes; melanogaster is always listed (it is the
# design species) even if its folder is missing.
detect_species <- function() {
  spp <- character(0)
  if (dir.exists(DATA_DIRS)) {
    subs <- list.dirs(DATA_DIRS, recursive = FALSE, full.names = FALSE)
    for (s in subs) {
      sp <- sub("^d_", "", tolower(s))
      if (identical(sp, tolower(s))) next
      # accept both "rna.fna(.gz)" and NCBI-style "GCF_..._rna.fna.gz" names
      if (length(list.files(file.path(DATA_DIRS, s),
                            pattern = "(^|_)rna\\.fna(\\.gz)?$",
                            ignore.case = TRUE)))
        spp <- c(spp, sp)
    }
  }
  unique(c("melanogaster", sort(spp)))
}
SPECIES_CHOICES <- detect_species()

# ---- Probe design constants (from the lab spreadsheet) ----------------------

AMPLIFIERS <- list(
  B1 = list(half1 = "GAGGAGGGCAGCAAACGG", spacer1 = "AA",
            half2 = "GAAGAGTCTTCCTTTACG", spacer2 = "TA"),
  B2 = list(half1 = "CCTCGTAAATCCTCATCA", spacer1 = "AA",
            half2 = "ATCATCCAGTAAACCGCC", spacer2 = "AA"),
  B3 = list(half1 = "GTCCCTGCCTCTATATCT", spacer1 = "TT",
            half2 = "CCACTCAACTTTAACCCG", spacer2 = "TT"),
  B4 = list(half1 = "CCTCAACCTACCTCCAAC", spacer1 = "AA",
            half2 = "TCTCACCATATTCGCTTC", spacer2 = "AT"),
  B5 = list(half1 = "CTCACTCCCAATCTCTAT", spacer1 = "AA",
            half2 = "CTACCCTACAAATCCAAT", spacer2 = "AA"),
  # Additional amplifiers from Wang et al. 2020 (BioRxiv). Their published
  # spacers are "WW" (A/T); AA is used here.
  B7 = list(half1 = "CTTCAACCTCCACCTACC", spacer1 = "AA",
            half2 = "TCCAATCCCTACCCTCAC", spacer2 = "AA"),
  B9 = list(half1 = "CACGTATCTACTCCACTC", spacer1 = "AA",
            half2 = "TCAGCACACTCCCAACCC", spacer2 = "AA"),
  B10 = list(half1 = "CCTCAAGATACTCCTCTA", spacer1 = "AA",
             half2 = "CCTACTCGACTACCCTAG", spacer2 = "AA"),
  B13 = list(half1 = "AGGTAACGCCTTCCTGCT", spacer1 = "AA",
             half2 = "TTATGCTCAACATACAAC", spacer2 = "AA"),
  B14 = list(half1 = "AATGTCAATAGCGAGCGA", spacer1 = "AA",
             half2 = "CCCTATATTTCTGCACAG", spacer2 = "AA"),
  B15 = list(half1 = "CAGATTAACACACCACAA", spacer1 = "AA",
             half2 = "GGTATCTCGAACACTCTC", spacer2 = "AA"),
  B17 = list(half1 = "CGATTGTTTGTTGTGGAC", spacer1 = "AA",
             half2 = "GCATGCTAATCGGATGAG", spacer2 = "AA")
)
AMP_LIST <- names(AMPLIFIERS)   # every amplifier the engine knows about
# The lab's in-situ kit is HCR v2.0, so only B1-B5 are available by default.
# The v3 extras (B7, B9, B10, B13-B15, B17) are used only when the user
# explicitly selects "HCR v3.0" in the Settings block.
AMP_V2 <- c("B1", "B2", "B3", "B4", "B5")

# Active amplifier set for a run: B1-B5 (v2) unless the v3 set is requested.
# NULL-safe so batch/env-driven settings never error before amp_set is set.
active_amps <- function(settings) {
  amps <- if (!is.null(settings) && !is.null(settings$amp_set)) settings$amp_set else "v2"
  if (identical(amps, "v3")) AMP_LIST else AMP_V2
}
AMP_COLORS <- c(B1 = "#0ea5e9", B2 = "#8b5cf6", B3 = "#10b981",
                B4 = "#f59e0b", B5 = "#ef4444", B7 = "#ec4899",
                B9 = "#14b8a6", B10 = "#f97316", B13 = "#6366f1",
                B14 = "#84cc16", B15 = "#d946ef", B17 = "#78716c")
# Safe lookup: never errors on an unknown amplifier name (a bare [[ ]] lookup
# used to crash the plots with "subscript out of bounds").
amp_color <- function(a) {
  if (!is.null(a) && length(a) && a %in% names(AMP_COLORS)) AMP_COLORS[[a]]
  else "#64748b"
}

# ---- Sequence utilities -----------------------------------------------------

revcomp <- function(s) {
  s <- chartr("ACGTUacgtu", "TGCAAugcaa", s)
  chars <- strsplit(s, "")[[1]]
  paste0(rev(chars), collapse = "")
}

clean_sequence <- function(raw) {
  lines <- strsplit(raw, "\r?\n")[[1]]
  lines <- lines[!grepl("^\\s*>", lines)]
  s <- toupper(paste0(lines, collapse = ""))
  s <- gsub("[^ACGTUN]", "", s)
  gsub("U", "T", s)
}

gc_percent <- function(s) {
  if (!nchar(s)) return(0)
  100 * (nchar(gsub("[^GC]", "", s))) / nchar(s)
}

max_run <- function(s) {
  chars <- strsplit(s, "")[[1]]
  m <- 1; cur <- 1
  for (i in 2:length(chars)) {
    if (chars[i] == chars[i - 1]) { cur <- cur + 1; m <- max(m, cur) } else cur <- 1
  }
  m
}

# ---- FASTA loading / gene index ---------------------------------------------

parse_header <- function(h) {
  h <- sub("^>", "", h)
  id <- strsplit(h, "\\s+")[[1]][1]
  symbol <- NA_character_; variant <- NA_character_; gname <- NA_character_
  # NCBI style: ... gene name (symbol), transcript variant A, mRNA
  m <- regmatches(h, regexpr("\\([A-Za-z0-9_.-]+\\)", h))
  if (length(m) && nzchar(m)) symbol <- gsub("[()]", "", m[1])
  if (grepl("transcript variant [A-Za-z0-9]+", h))
    variant <- sub(".*transcript variant ([A-Za-z0-9]+).*", "\\1", h)
  # NCBI gene/protein name: text between "Drosophila <species> " and the
  # parenthetical, e.g. "protein doublesex (LOC122619035)" -> "doublesex".
  # (Non-melanogaster RefSeq files use LOC IDs as the symbol, so the name is
  # often the only way to find a gene.)
  nm2 <- sub("^\\S+\\s+", "", h)                       # drop accession
  nm2 <- sub("^(?:PREDICTED: )?Drosophila \\S+ ", "", nm2)
  nm2 <- sub("\\s*\\([^)]*\\).*$", "", nm2)
  nm2 <- trimws(sub(",\\s*(transcript|mRNA).*$", "", nm2))
  if (nzchar(nm2) && !grepl("^(uncharacterized|putative)", nm2, ignore.case = TRUE)) {
    gname <- sub("^protein ", "", nm2)
  }
  # FlyBase style: name=dsx-RA; ... parent=FBgn0000546,dsx
  if (grepl("(^|[ ;])name=[^;]+", h)) {
    nm <- sub(".*[ ;]name=([^;]+).*", "\\1", h)
    parts <- strsplit(nm, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 2 && grepl("^R[A-Z0-9]+$", parts[length(parts)])) {
      if (is.na(symbol)) symbol <- paste(parts[-length(parts)], collapse = "-")
      variant <- sub("^R", "", parts[length(parts)])
    }
  }
  if (is.na(symbol) && grepl("parent=[^;]+", h)) {
    par <- sub(".*parent=([^;]+).*", "\\1", h)
    genes <- strsplit(par, ",", fixed = TRUE)[[1]]
    if (length(genes) >= 2) symbol <- genes[2]
  }
  list(id = id, symbol = symbol, variant = variant, name = gname, title = h)
}

read_fasta <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  hdr <- grep("^>", lines)
  if (!length(hdr)) return(NULL)
  ends <- c(hdr[-1] - 1, length(lines))
  records <- lapply(seq_along(hdr), function(i) {
    info <- parse_header(lines[hdr[i]])
    info$seq <- toupper(gsub("[^A-Za-z]", "",
                             paste0(lines[(hdr[i] + 1):ends[i]], collapse = "")))
    info
  })
  records
}

# Locate the transcript FASTA for a species.
# Searches the hardcoded DATA_DIRS for subfolders matching the species name.
# Normalize species spellings to the folder names used in reference_files.
canonical_species <- function(sp) {
  s <- tolower(trimws(sp))
  if (s %in% c("tesseri", "tessieri", "teissier", "teissierii")) "teissieri"
  else s
}

find_species_file <- function(species) {
  key <- canonical_species(species)
  dir <- DATA_DIRS
  if (!dir.exists(dir)) {
    warning("Reference directory does not exist: ", dir)
    return(NULL)
  }
  # Look for subfolders like d_melanogaster, d_simulans, d_eugracilis
  subs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
  subs <- subs[grepl(key, tolower(basename(subs)), fixed = TRUE)]
  for (sd in subs) {
    # accept both "rna.fna(.gz)" and NCBI-style "GCF_..._rna.fna.gz" names
    hits <- list.files(sd, pattern = "(^|_)rna\\.fna(\\.gz)?$",
                       ignore.case = TRUE, full.names = TRUE)
    if (length(hits)) return(hits[1])
  }
  NULL
}

# Load (and cache) a species transcriptome; index transcripts by gene symbol.
load_species <- function(species, cache) {
  key <- canonical_species(species)
  if (!is.null(cache[[key]])) return(cache[[key]])
  path <- find_species_file(species)
  if (is.null(path)) return(NULL)
  recs <- read_fasta(path)
  idx <- new.env(hash = TRUE)
  byid <- new.env(hash = TRUE)
  for (r in recs) {
    byid[[tolower(sub("\\..*$", "", r$id))]] <- r   # accession, version-stripped
    keys <- character(0)
    if (!is.na(r$symbol)) keys <- c(keys, tolower(r$symbol))
    if (!is.na(r$name))   keys <- c(keys, tolower(r$name))
    for (k in unique(keys)) idx[[k]] <- c(idx[[k]], list(r))
  }
  idx[[".byid"]] <- byid
  cache[[key]] <- idx
  idx
}

# Choose one transcript from a gene's records: accession/isoform hint first,
# then NM_ accessions, then longest.
choose_transcript <- function(recs, hint = NULL) {
  warning_msg <- NULL
  chosen <- NULL
  if (!is.null(hint) && nzchar(hint)) {
    hint <- trimws(hint)
    # exact accession / record id
    for (r in recs) if (identical(r$id, hint)) chosen <- r
    # isoform hint: "RA", "RB", "A", "variant A"
    if (is.null(chosen)) {
      letter <- sub("^(?:R|variant\\s*)?([A-Za-z0-9]+)$", "\\1", hint, ignore.case = TRUE)
      for (r in recs)
        if (!is.na(r$variant) && toupper(r$variant) == toupper(letter)) chosen <- r
    }
    if (is.null(chosen))
      warning_msg <- paste0('Isoform hint "', hint, '" not matched; using default instead. ')
  }
  if (is.null(chosen)) {
    is_nm <- grepl("^NM_", vapply(recs, `[[`, "", "id"))
    ord <- order(!is_nm, -vapply(recs, function(r) nchar(r$seq), 0))
    chosen <- recs[[ord[1]]]
  }
  if (length(recs) > 1)
    warning_msg <- paste0(warning_msg, length(recs), " transcripts found; using ",
                          chosen$id, " (", nchar(chosen$seq),
                          " nt). Verify this is the isoform you want.")
  # Build isoform summary table
  isoform_df <- do.call(rbind, lapply(recs, function(r) {
    data.frame(
      ID = r$id,
      Symbol = if (!is.na(r$symbol)) r$symbol else "",
      Name = if (!is.na(r$name)) r$name else "",
      Variant = if (!is.na(r$variant)) r$variant else "",
      Length = nchar(r$seq),
      Selected = identical(r$id, chosen$id),
      stringsAsFactors = FALSE
    )
  }))
  list(rec = chosen, warning = warning_msg, all_isoforms = isoform_df)
}

# ---- JBrowse / FlyBase link helpers ----------------------------------------

jbrowse_url <- function(rec) {
  # FlyBase gene report if we have a symbol
  if (!is.na(rec$symbol) && nzchar(rec$symbol)) {
    return(paste0("https://flybase.org/reports/", rec$symbol, ".html"))
  }
  # NCBI nuccore if we have an accession
  if (grepl("^[NX]M_[0-9]+", rec$id)) {
    return(paste0("https://www.ncbi.nlm.nih.gov/nuccore/", rec$id))
  }
  NULL
}

# ---- NCBI online fallback ---------------------------------------------------

# Map local file base names to NCBI organism names (extend as needed).
SPECIES_ORGANISMS <- c(
  melanogaster = "Drosophila melanogaster",
  teissieri    = "Drosophila teissieri",
  tesseri      = "Drosophila teissieri",   # common misspelling
  ananassae    = "Drosophila ananassae",
  erecta       = "Drosophila erecta",
  eugracilis   = "Drosophila eugracilis",
  simulans     = "Drosophila simulans",
  virilis      = "Drosophila virilis"
)
organism_for <- function(species) {
  s <- tolower(trimws(species))
  if (s %in% names(SPECIES_ORGANISMS)) return(unname(SPECIES_ORGANISMS[[s]]))
  if (grepl(" ", s)) return(species)          # already a full organism name
  paste("Drosophila", species)                # sensible default for this lab
}

# Reference assemblies for the one-click species downloader (RefSeq).
SPECIES_GCF <- c(
  melanogaster = "GCF_000001215.4_Release_6_plus_ISO1_MT",
  simulans     = "GCF_016746395.2_Prin_Dsim_3.1",
  eugracilis   = "GCF_018153835.1_ASM1815383v1",
  erecta       = "GCF_003286155.1_DereRS2",
  ananassae    = "GCF_017639315.1_ASM1763931v2",
  teissieri    = "GCF_016746235.2_Prin_Dtei_1.1",
  virilis      = "GCF_003285735.1_DvirRS2"
)

# NCBI compressed download size estimate per species (rna + genome + annotation,
# MB); used only for user-facing copy.
SPECIES_DL_MB <- c(
  melanogaster =  73, simulans = 61, eugracilis = 66, erecta = 62,
  ananassae = 86, teissieri = 61, virilis = 75
)

# ---- on-demand species references (installed on first use, not bundled) ------

# Directory where a species' downloaded reference files land.
species_ref_dir <- function(species) {
  file.path(DATA_DIRS, paste0("d_", canonical_species(species)))
}

# Has an actual rna.fna[.gz] on disk for the species? (The ground-truth check;
# detect_species() always lists melanogaster for the UI, so it can't be used
# here.)
species_reference_installed <- function(species) {
  suppressWarnings(!is.null(find_species_file(canonical_species(species))))
}

# Species (of the ones the app can auto-install) currently on disk.
installed_species <- function() {
  names(SPECIES_GCF)[vapply(names(SPECIES_GCF),
                            species_reference_installed, TRUE)]
}

# Checkbox labels for the cross-species section: installed species show their
# plain name; anything else advertises the on-demand NCBI download.
species_checkbox_labels <- function() {
  installed <- installed_species()
  setNames(names(SPECIES_GCF), vapply(names(SPECIES_GCF), function(sp) {
    if (sp %in% installed) return(sp)
    paste0(sp, " (on demand ~", SPECIES_DL_MB[[sp]], " MB unless you tick it)")
  }, character(1)))
}

# Download a species' reference set from NCBI into DATA_DIRS/d_<species>.
# Files stay compressed (.gz); the app's FASTA/GTF readers handle that natively.
# The NCBI GFF3 annotation is saved as genomic.gff.gz (its own format) and is
# picked up by find_gtf_file(): parse_gtf_exons()/gtf_locus() already read
# GFF3 attribute syntax (gene=, transcript_id=, Parent=), so no conversion is
# needed. Existing components are kept; only missing ones are fetched.
# Returns a status message.
download_species_ref <- function(species, with_genome = TRUE,
                                 progress = message) {
  # large files (up to ~150 MB) get 15 min; still finite so a stalled
  # connection errors out instead of freezing the app forever
  old_to <- getOption("timeout")
  on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = 900)
  sp <- canonical_species(species)
  gcf <- SPECIES_GCF[[sp]]
  if (is.null(gcf)) return(paste0("No known assembly for '", species, "'."))
  acc <- sub("^(GCF_[0-9]+)\\..*$", "\\1", gcf)   # GCF_016746235
  nums <- sub("^GCF_", "", acc)                   # 016746235
  segs <- substring(nums, c(1, 4, 7), c(3, 6, 9)) # 016 / 746 / 235
  base <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/", segs[1],
                 "/", segs[2], "/", segs[3], "/", gcf, "/")
  dest <- species_ref_dir(sp)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  components <- list(
    rna        = paste0(gcf, "_rna.fna.gz"),
    genome     = paste0(gcf, "_genomic.fna.gz"),
    annotation = paste0(gcf, "_genomic.gff.gz")
  )
  if (!with_genome) components <- components["rna"]
  saved <- character(0)
  for (nm in names(components)) {
    f <- switch(nm, rna = "rna.fna.gz", genome = "genomic.fna.gz",
                annotation = "genomic.gff.gz")
    out <- file.path(dest, f)
    if (file.exists(out) && file.size(out) > 0) { saved <- c(saved, f); next }
    part <- paste0(out, ".part")
    on.exit(if (file.exists(part)) unlink(part), add = TRUE)
    url <- paste0(base, components[[nm]])
    progress(paste0("Downloading ", sp, " ", components[[nm]], " ..."))
    ok <- tryCatch({
      utils::download.file(url, part, mode = "wb", quiet = TRUE)
      file.exists(part) && file.size(part) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) {
      if (file.exists(part)) unlink(part)
      if (nm == "rna")
        return(paste0("Download failed for ", sp, " (rna.fna.gz). Check your ",
                      "internet connection, or download manually from ", base))
      progress(paste0("  (", components[[nm]], " unavailable for ", sp, "; continuing)"))
      next
    }
    file.rename(part, out)        # atomic: never leave half-written files
    saved <- c(saved, f)
  }
  paste0(sp, " reference saved to ", dest, " (", paste(saved, collapse = ", "), ")")
}

# Fetch a URL as text lines with a HARD timeout. Raw readLines(url) has no
# timeout: if the network stalls mid-request (VPN drop, Mac sleeping, NCBI
# hiccup) R waits forever and the app appears frozen. download.file honours
# options("timeout"), so we route every small HTTP read through it.
http_lines <- function(url, timeout = 30) {
  old <- getOption("timeout")
  on.exit(options(timeout = old), add = TRUE)
  options(timeout = timeout)
  tmp <- tempfile()
  ok <- tryCatch({
    utils::download.file(url, tmp, quiet = TRUE, mode = "wt")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(tmp)) return(NULL)
  lines <- tryCatch(readLines(tmp, warn = FALSE), error = function(e) NULL)
  unlink(tmp)
  lines
}

read_fasta_from_accession <- function(acc) {
  fetch <- function(a) {
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",
                  utils::URLencode(a), "&rettype=fasta&retmode=text")
    r <- http_lines(url)
    if (is.null(r)) character(0) else r
  }
  fa <- fetch(acc)
  # NCBI 400s on unknown version suffixes (e.g. NM_169202.3) — retry without it
  if (!length(fa) || !any(grepl("^>", fa))) fa <- fetch(sub("\\..*$", "", acc))
  if (!length(fa) || !any(grepl("^>", fa))) stop(paste("Could not fetch accession", acc))
  tmp <- tempfile(fileext = ".fasta")
  writeLines(fa, tmp)
  read_fasta(tmp)
}

ncbi_fetch_gene <- function(gene, organism) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  if (grepl("^[A-Z]{2,3}_[0-9]+(\\.[0-9]+)?$", gene))
    return(read_fasta_from_accession(gene))  # direct accession

  esearch <- function(term) {
    url <- paste0(base, "esearch.fcgi?db=nuccore&retmax=20&term=",
                  utils::URLencode(term, reserved = TRUE))
    xml <- http_lines(url)
    if (is.null(xml)) stop("Could not reach NCBI (timed out).")
    Sys.sleep(0.4)  # stay under NCBI's 3 req/s unauthenticated limit
    ids <- regmatches(xml, regexpr("(?<=<Id>)\\d+(?=</Id>)", xml, perl = TRUE))
    ids[nzchar(ids)]
  }
  efetch_records <- function(ids) {
    fa <- http_lines(paste0(base, "efetch.fcgi?db=nuccore&id=",
                            paste(ids, collapse = ","),
                            "&rettype=fasta&retmode=text"))
    if (is.null(fa)) return(NULL)
    Sys.sleep(0.4)
    tmp <- tempfile(fileext = ".fasta")
    writeLines(fa, tmp)
    read_fasta(tmp)
  }

  filt <- 'srcdb_refseq[PROP] AND biomol_mrna[PROP]'
  ids <- esearch(sprintf('%s[Gene Name] AND "%s"[Organism] AND %s', gene, organism, filt))
  if (length(ids)) return(efetch_records(ids))

  # alias fallback: the gene may be indexed under a different symbol
  # (e.g. "Drop" is NCBI symbol "Dr") - search all fields, then keep only
  # records whose parsed symbol or gene/protein name matches the query.
  ids <- esearch(sprintf('%s[All Fields] AND "%s"[Organism] AND %s', gene, organism, filt))
  if (!length(ids)) return(NULL)
  recs <- efetch_records(ids)
  if (is.null(recs)) return(NULL)
  g <- tolower(gene)
  pat <- paste0("\\b", gsub("([^A-Za-z0-9])", "\\\\\\1", g), "\\b")
  keep <- Filter(function(r)
    (!is.na(r$symbol) && tolower(r$symbol) == g) ||
    (!is.na(r$name) && grepl(pat, tolower(r$name))), recs)
  if (length(keep)) keep else NULL
}

# Look a gene up in NCBI's Gene database for a given organism. This is how
# you get the symbol/LOC ID a gene is filed under in another species when
# the local RefSeq file doesn't know the melanogaster symbol (e.g. "Abd-B"
# is "LOC6504858" in D. ananassae). Returns list(symbol, description, uid)
# or NULL. One esearch + one esummary call.
ncbi_gene_lookup <- function(gene, organism) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  term <- sprintf('%s[Gene Name] AND "%s"[Organism]', gene, organism)
  url <- paste0(base, "esearch.fcgi?db=gene&retmax=5&term=",
                utils::URLencode(term, reserved = TRUE))
  xml <- http_lines(url)
  if (is.null(xml)) return(NULL)
  Sys.sleep(0.4)   # stay under NCBI's 3 req/s unauthenticated limit
  uid <- regmatches(xml, regexpr("(?<=<Id>)\\d+(?=</Id>)", xml, perl = TRUE))
  if (!length(uid)) return(NULL)
  uid <- uid[1]
  sxml <- http_lines(paste0(base, "esummary.fcgi?db=gene&id=", uid))
  Sys.sleep(0.4)
  if (is.null(sxml)) return(list(symbol = NULL, description = NULL, uid = uid))
  txt <- paste(sxml, collapse = " ")
  grab <- function(tag) {
    m <- regmatches(txt, regexpr(paste0("(?<=<", tag, ">)[^<]+(?=</", tag, ">)"),
                                 txt, perl = TRUE))
    if (length(m) && nzchar(m)) m else NULL
  }
  list(symbol = grab("Name"), description = grab("Description"), uid = uid)
}

# ---- GTF / exon utilities -----------------------------------------------

# Read (and cache) the raw lines of a GTF file. The GTF is reused by every
# gene in a run, so we only pay the file-read cost once per species.
get_gtf_lines <- function(gtf_path, cache = NULL) {
  key <- paste0(".gtf:", gtf_path)
  if (!is.null(cache) && !is.null(cache[[key]])) return(cache[[key]])
  if (!file.exists(gtf_path)) return(NULL)
  con <- if (grepl("\\.gz$", gtf_path)) gzfile(gtf_path, "rt") else file(gtf_path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  if (!is.null(cache)) cache[[key]] <- lines
  lines
}

# Parse GTF file and extract exon info for a gene
parse_gtf_exons <- function(gtf_path, gene_symbol, cache = NULL) {
  lines <- get_gtf_lines(gtf_path, cache)
  if (is.null(lines)) return(NULL)
  # Filter for exon features matching the gene
  g <- tolower(gene_symbol)
  # GTF is tab-delimited; look for lines with "exon" in the 3rd column
  exon_lines <- lines[grepl("^[^#].*	exon	", lines, perl = TRUE)]
  # Match by gene symbol in attributes (column 9). GTF style is
  # gene_name "x" / gene_id "x" / gene "x"; NCBI GFF3 uses ;gene=x; style.
  pat <- paste0('gene_name "', g, '"|gene_id "', g, '"|gene "', g, '"',
                '|(^|;)gene=', g, '([;,]|$)')
  matched <- exon_lines[grepl(pat, exon_lines, ignore.case = TRUE)]
  if (!length(matched)) return(NULL)

  # Parse each tab-delimited line
  exons <- lapply(matched, function(line) {
    parts <- strsplit(line, "	")[[1]]
    if (length(parts) < 9) return(NULL)
    attr <- parts[9]
    # Extract transcript_id from attributes (GTF "..." or GFF3 = style)
    tid <- sub('.*transcript_id "([^"]+)".*', "\\1", attr)
    if (identical(tid, attr)) tid <- sub(".*transcript_id=([^;]+).*", "\\1", attr)
    if (identical(tid, attr)) tid <- sub(".*Parent=([^;,]+).*", "\\1", attr)
    list(chr = parts[1], start = as.integer(parts[4]), end = as.integer(parts[5]),
         strand = parts[7], transcript_id = tid)
  })
  exons <- exons[!vapply(exons, is.null, TRUE)]
  if (!length(exons)) return(NULL)

  # Group by transcript
  by_tx <- split(exons, vapply(exons, `[[`, "", "transcript_id"))
  by_tx
}

# Find shared exon regions across all isoforms (intersection on genome)
find_shared_exons <- function(by_tx) {
  if (length(by_tx) < 2) {
    # Only one isoform - return all its exons
    out <- do.call(rbind, lapply(by_tx[[1]], function(e) 
      data.frame(chr = e$chr, start = e$start, end = e$end, 
                 strand = e$strand, stringsAsFactors = FALSE)))
    # Deduplicate in case the single isoform has overlapping exons
    if (!is.null(out) && nrow(out)) {
      out <- out[!duplicated(paste(out$chr, out$start, out$end, out$strand)), ]
    }
    return(out)
  }

  # Get all exon regions from first transcript
  shared <- data.frame(
    chr = vapply(by_tx[[1]], `[[`, "", "chr"),
    start = vapply(by_tx[[1]], `[[`, 0, "start"),
    end = vapply(by_tx[[1]], `[[`, 0, "end"),
    strand = vapply(by_tx[[1]], `[[`, "", "strand"),
    stringsAsFactors = FALSE
  )

  # Intersect with each subsequent transcript's exons
  for (i in 2:length(by_tx)) {
    tx_exons <- data.frame(
      chr = vapply(by_tx[[i]], `[[`, "", "chr"),
      start = vapply(by_tx[[i]], `[[`, 0, "start"),
      end = vapply(by_tx[[i]], `[[`, 0, "end"),
      stringsAsFactors = FALSE
    )

    # Find overlaps between shared and current transcript exons
    new_shared <- data.frame(chr = character(0), start = integer(0), 
                              end = integer(0), strand = character(0),
                              stringsAsFactors = FALSE)
    for (j in seq_len(nrow(shared))) {
      for (k in seq_len(nrow(tx_exons))) {
        # Check if on same chromosome
        if (shared$chr[j] != tx_exons$chr[k]) next
        # Find intersection
        o_start <- max(shared$start[j], tx_exons$start[k])
        o_end <- min(shared$end[j], tx_exons$end[k])
        if (o_start <= o_end) {
          new_shared <- rbind(new_shared, data.frame(
            chr = shared$chr[j], start = o_start, end = o_end,
            strand = shared$strand[j], stringsAsFactors = FALSE
          ))
        }
      }
    }
    shared <- new_shared
    if (!nrow(shared)) break
  }
  # Deduplicate: exact same coordinates from overlapping isoform intersections
  if (!is.null(shared) && nrow(shared)) {
    shared <- shared[!duplicated(paste(shared$chr, shared$start, shared$end, shared$strand)), ]
  }
  shared
}

# Read (and cache) the genome as a chromosome -> sequence index.
get_genome_index <- function(genome_path, cache = NULL) {
  key <- paste0(".genome:", genome_path)
  if (!is.null(cache) && !is.null(cache[[key]])) return(cache[[key]])
  if (!file.exists(genome_path)) return(NULL)
  fa <- read_fasta(genome_path)
  if (is.null(fa)) return(NULL)
  chr_idx <- new.env(hash = TRUE)
  for (r in fa) chr_idx[[r$id]] <- r$seq
  if (!is.null(cache)) cache[[key]] <- chr_idx
  chr_idx
}

# Extract sequence from genomic FASTA for given regions
extract_regions <- function(genome_path, regions, cache = NULL) {
  if (!nrow(regions)) return("")
  chr_idx <- get_genome_index(genome_path, cache)
  if (is.null(chr_idx)) return("")

  seqs <- character(nrow(regions))
  for (i in seq_len(nrow(regions))) {
    chr <- regions$chr[i]
    if (!exists(chr, chr_idx)) next
    seq <- chr_idx[[chr]]
    s <- regions$start[i]
    e <- regions$end[i]
    if (s > 0 && e <= nchar(seq)) {
      frag <- substr(seq, s, e)
      if (regions$strand[i] == "-") frag <- revcomp(frag)
      seqs[i] <- frag
    }
  }
  paste0(seqs, collapse = "")
}

# Find genomic FASTA for a species (parallel to find_species_file)
find_genome_file <- function(species) {
  key <- canonical_species(species)
  dir <- DATA_DIRS
  if (!dir.exists(dir)) return(NULL)
  subs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
  subs <- subs[grepl(key, tolower(basename(subs)), fixed = TRUE)]
  for (sd in subs) {
    # Accept "genomic.fna(.gz)" and NCBI-style "GCF_..._genomic.fna(.gz)",
    # but never "cds_from_genomic.fna" (that's CDS sequence, not chromosomes).
    hits <- list.files(sd, pattern = "(^|_)genomic\\.fna(\\.gz)?$",
                       ignore.case = TRUE, full.names = TRUE)
    hits <- hits[!grepl("cds_from", basename(hits), ignore.case = TRUE)]
    if (length(hits)) return(hits[1])
  }
  NULL
}

# Find GTF file for a species
find_gtf_file <- function(species) {
  key <- canonical_species(species)
  dir <- DATA_DIRS
  if (!dir.exists(dir)) return(NULL)
  subs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
  subs <- subs[grepl(key, tolower(basename(subs)), fixed = TRUE)]
  for (sd in subs) {
    hits <- list.files(sd, pattern = "(^|_)genomic\\.(gtf|gff)(\\.gz)?$",
                       ignore.case = TRUE, full.names = TRUE)
    if (length(hits)) return(hits[1])
  }
  NULL
}

# ---- Sequence-based exon mapping ----------------------------------------------

# Extract a padded locus sequence on the gene's strand for alignment.
get_locus_seq <- function(genome_path, locus, pad = 25000, cache = NULL) {
  chr_idx <- get_genome_index(genome_path, cache)
  if (is.null(chr_idx) || !exists(locus$chr, envir = chr_idx)) return(NULL)
  g <- chr_idx[[locus$chr]]
  s <- max(1, locus$start - pad)
  e <- min(nchar(g), locus$end + pad)
  sq <- substr(g, s, e)
  if (locus$strand == "-") sq <- revcomp(sq)
  list(seq = sq, chr = locus$chr, gstart = s, gend = e, strand = locus$strand)
}

# Greedy exact-match chaining of a transcript onto a locus (mini-BLAT for
# same-assembly sequences). Returns exon blocks in GENOMIC coordinates in
# transcript order, plus the fraction of the transcript covered. Splice
# junctions are found implicitly: the transcript simply has no intron bases,
# so each new seed lands at the next exon start.
chain_tx_to_locus <- function(tx_seq, locus) {
  t <- tx_seq; g <- locus$seq
  nt <- nchar(t); ng <- nchar(g)
  if (!nt || !ng) return(NULL)
  pos_t <- 1; pos_g <- 1; fails <- 0
  rows <- NULL; covered <- 0; guard <- 0
  while (pos_t <= nt && guard < 20000) {
    guard <- guard + 1
    k <- min(30, nt - pos_t + 1)
    if (k < 12) break
    hit <- regexpr(substr(t, pos_t, pos_t + k - 1),
                   substr(g, pos_g, ng), fixed = TRUE)[1]
    if (hit < 0 && k > 12) {                       # retry with a shorter seed
      k <- 12
      hit <- regexpr(substr(t, pos_t, pos_t + k - 1),
                     substr(g, pos_g, ng), fixed = TRUE)[1]
    }
    if (hit < 0) {              # SNP/indel base, tiny exon, or unaligned tail
      fails <- fails + 1
      pos_t <- pos_t + if (fails > 100) 50 else 1
      next
    }
    fails <- 0
    g0 <- pos_g + hit - 1
    # extend the exact match to the right, chunk by chunk
    mlen <- k
    repeat {
      rem_t <- nt - (pos_t + mlen) + 1
      rem_g <- ng - (g0 + mlen) + 1
      if (rem_t <= 0 || rem_g <= 0) break
      step <- min(4000, rem_t, rem_g)
      ct <- substr(t, pos_t + mlen, pos_t + mlen + step - 1)
      cg <- substr(g, g0 + mlen, g0 + mlen + step - 1)
      if (identical(ct, cg)) { mlen <- mlen + step; next }
      neq <- which(strsplit(ct, "", fixed = TRUE)[[1]] !=
                   strsplit(cg, "", fixed = TRUE)[[1]])
      mlen <- mlen + if (length(neq)) neq[1] - 1 else step
      break
    }
    g1 <- g0 + mlen - 1
    if (locus$strand == "-") {
      gs <- locus$gstart + ng - g1; ge <- locus$gstart + ng - g0
    } else {
      gs <- locus$gstart + g0 - 1; ge <- locus$gstart + g1 - 1
    }
    rows <- rbind(rows, data.frame(tstart = pos_t, tend = pos_t + mlen - 1,
                                   chr = locus$chr, start = gs, end = ge,
                                   strand = locus$strand,
                                   stringsAsFactors = FALSE))
    covered <- covered + mlen
    pos_t <- pos_t + mlen
    pos_g <- g1 + 1
  }
  if (is.null(rows)) return(NULL)
  list(exons = rows, coverage = covered / nt)
}

# Merge blocks separated by tiny gaps (<=3 bp on both the transcript and the
# genome). These come from SNPs/reference edits, not real introns.
merge_near_blocks <- function(ex) {
  if (nrow(ex) < 2) return(ex)
  out <- ex[1, ]; j <- 1
  for (r in 2:nrow(ex)) {
    tgap <- ex$tstart[r] - out$tend[j] - 1
    ggap <- if (ex$strand[1] == "-") out$start[j] - ex$end[r] - 1
            else ex$start[r] - out$end[j] - 1
    if (tgap <= 3 && ggap <= 3) {
      out$tend[j] <- ex$tend[r]
      if (ex$strand[1] == "-") out$start[j] <- ex$start[r] else out$end[j] <- ex$end[r]
    } else {
      out <- rbind(out, ex[r, ]); j <- j + 1
    }
  }
  out
}

# Align all transcript records of a gene to its locus. Returns by_tx-style
# exon lists keyed by transcript id, plus per-transcript coverage.
map_isoforms_seq <- function(recs, locus, min_cov = 0.7) {
  by_tx <- list(); cov <- numeric(0)
  for (r in recs) {
    m <- chain_tx_to_locus(r$seq, locus)
    if (is.null(m) || !nrow(m$exons)) next
    cov[r$id] <- m$coverage
    if (m$coverage < min_cov) next
    ex <- merge_near_blocks(m$exons)
    by_tx[[r$id]] <- lapply(seq_len(nrow(ex)), function(j)
      list(chr = ex$chr[j], start = ex$start[j], end = ex$end[j],
           strand = ex$strand[j], transcript_id = r$id))
  }
  list(by_tx = by_tx, coverage = cov)
}

# Gene locus (chr/start/end/strand) from GTF exon lines; coordinates only,
# transcript grouping not needed.
gtf_locus <- function(gtf_path, gene_symbol, cache = NULL) {
  lines <- get_gtf_lines(gtf_path, cache)
  if (is.null(lines)) return(NULL)
  g <- tolower(gene_symbol)
  exon_lines <- lines[grepl("^[^#].*\texon\t", lines, perl = TRUE)]
  # GTF: gene_name "x"/gene_id "x"/gene "x"; NCBI GFF3: ;gene=x; style.
  pat <- paste0('gene_name "', g, '"|gene_id "', g, '"|gene "', g, '"',
                '|(^|;)gene=', g, '([;,]|$)')
  matched <- exon_lines[grepl(pat, exon_lines, ignore.case = TRUE)]
  if (!length(matched)) return(NULL)
  starts <- integer(0); ends <- integer(0); chrs <- character(0); strands <- character(0)
  for (line in matched) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 8) next
    s <- suppressWarnings(as.integer(parts[4]))
    e <- suppressWarnings(as.integer(parts[5]))
    if (is.na(s) || is.na(e)) next
    starts <- c(starts, s); ends <- c(ends, e)
    chrs <- c(chrs, parts[1]); strands <- c(strands, parts[7])
  }
  if (!length(starts)) return(NULL)
  chr <- names(which.max(table(chrs)))
  keep <- chrs == chr
  strand <- names(which.max(table(strands[keep])))
  list(chr = chr, start = min(starts[keep]), end = max(ends[keep]), strand = strand)
}

# Fallback locus finder: exact 60-mer anchors from the transcript searched
# across all chromosomes. Used when the GTF has no entry for the gene.
find_locus_by_sequence <- function(genome_path, tx_seq, cache = NULL) {
  chr_idx <- get_genome_index(genome_path, cache)
  if (is.null(chr_idx)) return(NULL)
  nt <- nchar(tx_seq)
  if (nt < 80) return(NULL)
  k <- 60
  spos <- unique(round(seq(1, nt - k + 1, length.out = min(5, nt - k + 1))))
  seeds <- vapply(spos, function(p) substr(tx_seq, p, p + k - 1), "")
  best <- NULL
  for (chr in ls(chr_idx)) {
    g <- chr_idx[[chr]]
    if (nchar(g) < 10000) next                    # skip tiny scaffolds
    for (st in c("+", "-")) {
      anchor_g <- integer(0); anchor_t <- integer(0)
      for (si in seq_along(seeds)) {
        s <- if (st == "+") seeds[si] else revcomp(seeds[si])
        pos <- regexpr(s, g, fixed = TRUE)[1]
        if (pos > 0) { anchor_g <- c(anchor_g, pos); anchor_t <- c(anchor_t, spos[si]) }
      }
      if (length(anchor_g) >= 2 && (is.null(best) || length(anchor_g) > best$n))
        best <- list(chr = chr, strand = st, n = length(anchor_g),
                     start = min(anchor_g - anchor_t),
                     end   = max(anchor_g + (nt - anchor_t)))
    }
  }
  if (!is.null(best)) best$start <- max(1, best$start)
  best
}

# ---- CDS annotation -----------------------------------------------------------

# Parse CDS features for a gene from the GTF, grouped by transcript_id.
gtf_cds_by_tx <- function(gtf_path, gene_symbol, cache = NULL) {
  lines <- get_gtf_lines(gtf_path, cache)
  if (is.null(lines)) return(NULL)
  g <- tolower(gene_symbol)
  cds_lines <- lines[grepl("^[^#].*\tCDS\t", lines, perl = TRUE)]
  # GTF: gene_name "x"/gene_id "x"/gene "x"; NCBI GFF3: ;gene=x; style.
  pat <- paste0('gene_name "', g, '"|gene_id "', g, '"|gene "', g, '"',
                '|(^|;)gene=', g, '([;,]|$)')
  matched <- cds_lines[grepl(pat, cds_lines, ignore.case = TRUE)]
  if (!length(matched)) return(NULL)
  out <- list()
  for (line in matched) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 9) next
    attr <- parts[9]
    tid <- sub('.*transcript_id "([^"]+)".*', "\\1", attr)
    if (identical(tid, attr)) tid <- sub(".*transcript_id=([^;]+).*", "\\1", attr)
    if (identical(tid, attr)) tid <- sub(".*Parent=([^;,]+).*", "\\1", attr)
    s <- suppressWarnings(as.integer(parts[4]))
    e <- suppressWarnings(as.integer(parts[5]))
    if (is.na(s) || is.na(e)) next
    out[[tid]] <- c(out[[tid]], list(list(chr = parts[1], start = s, end = e,
                                          strand = parts[7], transcript_id = tid)))
  }
  if (!length(out)) return(NULL)
  out
}

# Flatten a list of segment lists into a data.frame.
segs_to_df <- function(segs) {
  do.call(rbind, lapply(segs, function(x)
    data.frame(chr = x$chr, start = x$start, end = x$end, strand = x$strand,
               stringsAsFactors = FALSE)))
}

# Genomic overlap of two segment data.frames (same chr required).
intersect_two_dfs <- function(a, b) {
  out <- NULL
  for (j in seq_len(nrow(a))) for (k in seq_len(nrow(b))) {
    if (a$chr[j] != b$chr[k]) next
    s <- max(a$start[j], b$start[k]); e <- min(a$end[j], b$end[k])
    if (s <= e) out <- rbind(out, data.frame(chr = a$chr[j], start = s, end = e,
                                             strand = a$strand[j],
                                             stringsAsFactors = FALSE))
  }
  out
}

# Merge overlapping or nearly adjacent (gap <= 3 bp) rows. Adjacent fragments
# arise from alignment boundary wobble and would otherwise duplicate or drop
# sequence when extracted.
merge_overlaps <- function(df) {
  if (is.null(df) || nrow(df) < 2) return(df)
  df <- df[order(df$chr, df$start), ]
  out <- df[1, ]; j <- 1
  for (r in 2:nrow(df)) {
    if (df$chr[r] == out$chr[j] && df$start[r] <= out$end[j] + 4) {
      out$end[j] <- max(out$end[j], df$end[r])
    } else {
      out <- rbind(out, df[r, ]); j <- j + 1
    }
  }
  rownames(out) <- NULL
  out
}

# Regions that are CODING in every isoform: intersect per-isoform CDS segments,
# where each isoform's CDS = its exons clipped to its GTF CDS (fallback: the
# gene-level CDS union, for isoforms missing from the GTF).
shared_cds_regions <- function(by_tx, cds_by_tx) {
  cds_union <- merge_overlaps(segs_to_df(unlist(cds_by_tx, recursive = FALSE)))
  per_tx <- lapply(names(by_tx), function(id) {
    m <- match(sub("\\..*$", "", id), sub("\\..*$", "", names(cds_by_tx)))
    cds <- if (!is.na(m)) segs_to_df(cds_by_tx[[m]]) else cds_union
    intersect_two_dfs(segs_to_df(by_tx[[id]]), cds)
  })
  per_tx <- Filter(function(d) !is.null(d) && nrow(d) > 0, per_tx)
  if (!length(per_tx)) return(NULL)
  out <- per_tx[[1]]
  if (length(per_tx) > 1)
    for (i in 2:length(per_tx)) {
      out <- intersect_two_dfs(out, per_tx[[i]])
      if (is.null(out) || !nrow(out)) return(NULL)
    }
  # Deduplicate: exact same coordinates from overlapping isoform CDS
  if (!is.null(out) && nrow(out)) {
    out <- out[!duplicated(paste(out$chr, out$start, out$end, out$strand)), ]
  }
  out
}

# ---- Isoform track plot ------------------------------------------------------

# Defensive cleanup: drop transcripts with no exons and exons lacking
# coordinates. Malformed entries used to crash the plots with
# "subscript out of bounds".
sanitize_by_tx <- function(by_tx) {
  if (is.null(by_tx) || !length(by_tx)) return(NULL)
  if (is.null(names(by_tx)) || any(!nzchar(names(by_tx))))
    names(by_tx) <- paste0("tx", seq_along(by_tx))
  by_tx <- lapply(by_tx, function(es)
    Filter(function(x) !is.null(x$start) && !is.null(x$end) &&
             length(x$start) > 0 && length(x$end) > 0 &&
             is.finite(x$start) && is.finite(x$end), es))
  by_tx <- Filter(function(es) length(es) > 0, by_tx)
  if (!length(by_tx)) NULL else by_tx
}

# Flip a set of start/end coordinates around [lo, hi] (for drawing
# minus-strand loci 5'->3' left to right).
flip_se <- function(s, e, lo, hi) list(s = lo + hi - e, e = lo + hi - s)

# Gene-model style diagram: one row per transcript, exons as blocks on a gray
# intron line, genomic coordinates on the x-axis (like IGV/Ensembl tracks).
#   by_tx       list of transcripts, each a list of exons (chr/start/end/strand)
#   selected_id accession of the transcript chosen for design (highlighted)
#   shared      optional data.frame of shared exon regions (extra green track)
#   labels      optional named character vector of row labels (names = tx ids)
#   sel_col     color for the selected transcript (amplifier color)
plot_isoform_tracks <- function(by_tx, selected_id = NULL, shared = NULL,
                                labels = NULL, sel_col = "#0ea5e9",
                                probe_frags = NULL, probe_col = "#0f172a",
                                cds_by_tx = NULL, shared_label = "shared") {
  by_tx <- sanitize_by_tx(by_tx)
  if (is.null(by_tx)) return(invisible(NULL))
  tx_ids <- names(by_tx)
  # order transcripts by genomic start so similar isoforms line up
  tx_start <- vapply(by_tx, function(es) min(vapply(es, `[[`, 0, "start")), 0)
  tx_end   <- vapply(by_tx, function(es) max(vapply(es, `[[`, 0, "end")), 0)
  ord <- order(tx_start, tx_end)
  by_tx <- by_tx[ord]; tx_ids <- tx_ids[ord]

  sel_base <- if (!is.null(selected_id)) sub("\\..*$", "", selected_id) else ""
  is_sel <- vapply(tx_ids, function(id)
    identical(sub("\\..*$", "", id), sel_base), TRUE)

  has_shared <- !is.null(shared) && nrow(shared) > 0
  n_rows <- length(by_tx) + if (has_shared) 1 else 0

  x_min <- min(vapply(by_tx, function(es) min(vapply(es, `[[`, 0, "start")), 0))
  x_max <- max(vapply(by_tx, function(es) max(vapply(es, `[[`, 0, "end")), 0))
  if (has_shared) {
    x_min <- min(x_min, min(shared$start)); x_max <- max(x_max, max(shared$end))
  }
  pad <- 0.02 * (x_max - x_min)

  # Draw minus-strand loci 5'->3' (left to right), so the view reads in the
  # direction of transcription and matches flipped loci in other species.
  locus_strand <- {
    ss <- vapply(by_tx, function(es) {
      st <- es[[1]]$strand; if (is.null(st) || !length(st)) "" else st
    }, "")
    pick <- ss[ss %in% c("+", "-")]
    if (!is.null(selected_id)) {
      hit <- which(sub("\\..*$", "", names(by_tx)) ==
                     sub("\\..*$", "", selected_id))
      if (length(hit) && ss[hit[1]] %in% c("+", "-")) unname(ss[hit[1]])
      else if (length(pick)) unname(pick[1]) else "+"
    } else if (length(pick)) unname(pick[1]) else "+"
  }
  if (identical(locus_strand, "-")) {
    lo <- x_min; hi <- x_max
    by_tx <- lapply(by_tx, function(es) lapply(es, function(x) {
      f <- flip_se(x$start, x$end, lo, hi)
      x$start <- f$s; x$end <- f$e; x$strand <- "+"; x
    }))
    if (has_shared) {
      f <- flip_se(shared$start, shared$end, lo, hi)
      shared$start <- f$s; shared$end <- f$e; shared$strand <- "+"
    }
    if (!is.null(cds_by_tx))
      cds_by_tx <- lapply(cds_by_tx, function(es) lapply(es, function(x) {
        f <- flip_se(x$start, x$end, lo, hi)
        x$start <- f$s; x$end <- f$e; x$strand <- "+"; x
      }))
    if (!is.null(probe_frags) && length(probe_frags))
      probe_frags <- lapply(probe_frags, function(fr)
        data.frame(gstart = lo + hi - fr$gend, gend = lo + hi - fr$gstart))
    tx_ids <- names(by_tx)
    is_sel <- vapply(tx_ids, function(id)
      identical(sub("\\..*$", "", id), sel_base), TRUE)
  }

  op <- par(mar = c(3, 14, 1, 1)); on.exit(par(op))
  plot(0, type = "n", xlim = c(x_min - pad, x_max + pad),
       ylim = c(0.5, n_rows + 0.5), axes = FALSE, xlab = "", ylab = "")
  axis(1, col = "gray40", col.axis = "gray40", cex.axis = 0.8)
  mtext(if (identical(locus_strand, "-"))
          "genomic position (bp) - minus strand, shown 5' to 3'"
        else "genomic position (bp)",
        side = 1, line = 1.9, cex = 0.7, col = "gray40")

  # gene-level CDS union, used for isoforms missing from the GTF
  cds_union_df <- NULL
  if (!is.null(cds_by_tx))
    cds_union_df <- merge_overlaps(segs_to_df(unlist(cds_by_tx, recursive = FALSE)))

  # exons: UTR portions as thin gray boxes, CDS portions as thick colored
  # boxes (JBrowse-style); falls back to uniform thick boxes without CDS data
  draw_track <- function(y, exons, col, strand = NULL, cds = NULL) {
    s <- vapply(exons, `[[`, 0, "start"); e <- vapply(exons, `[[`, 0, "end")
    segments(min(s), y, max(e), y, col = "gray70", lwd = 2, lend = 2)
    # small direction arrows along the intron line
    if (!is.null(strand) && strand %in% c("+", "-")) {
      xs <- seq(min(s), max(e), length.out = 7)[2:6]
      text(xs, rep(y, length(xs)), labels = if (strand == "+") ">" else "<",
           col = "gray55", cex = 0.55)
    }
    if (!is.null(cds) && nrow(cds)) {
      rect(s, y - 0.09, e, y + 0.09, col = "#cbd5e1", border = NA)
      rect(cds$start, y - 0.18, cds$end, y + 0.18, col = col, border = NA)
    } else {
      rect(s, y - 0.18, e, y + 0.18, col = col, border = NA)
    }
  }

  y_sel <- NULL
  for (i in seq_along(by_tx)) {
    y <- n_rows - i + 1                     # first transcript on top
    if (is_sel[i]) y_sel <- y
    col <- if (is_sel[i]) sel_col else "#7dd3fc"
    # this transcript's CDS, clipped to its (aligned) exons
    cds_df <- NULL
    if (!is.null(cds_by_tx)) {
      m <- match(sub("\\..*$", "", tx_ids[i]), sub("\\..*$", "", names(cds_by_tx)))
      cds_df <- if (!is.na(m)) segs_to_df(cds_by_tx[[m]]) else cds_union_df
      if (!is.null(cds_df))
        cds_df <- intersect_two_dfs(segs_to_df(by_tx[[i]]), cds_df)
    }
    tx_strand <- by_tx[[i]][[1]]$strand
    if (is.null(tx_strand) || !length(tx_strand)) tx_strand <- "+"
    draw_track(y, by_tx[[i]], col, strand = tx_strand, cds = cds_df)
    lab <- if (!is.null(labels) && tx_ids[i] %in% names(labels)) labels[[tx_ids[i]]]
           else tx_ids[i]
    text(par("usr")[1], y, labels = lab, adj = c(1, 0.5), xpd = NA, cex = 0.65,
         font = if (is_sel[i]) 2 else 1,
         col  = if (is_sel[i]) sel_col else "gray30")
  }

  if (has_shared) {
    sh <- lapply(seq_len(nrow(shared)), function(j)
      list(start = shared$start[j], end = shared$end[j]))
    sh_strand <- if ("strand" %in% names(shared)) shared$strand[1] else "+"
    if (is.null(sh_strand) || !length(sh_strand)) sh_strand <- "+"
    draw_track(1, sh, "#10b981", strand = sh_strand)
    text(par("usr")[1], 1, labels = shared_label, adj = c(1, 0.5), xpd = NA,
         cex = 0.65, font = 2, col = "#10b981")
  }

  # overlay designed probe targets on the selected transcript's track
  if (!is.null(probe_frags) && length(probe_frags) && !is.null(y_sel)) {
    for (fr in probe_frags) {
      if (is.null(fr) || is.null(fr$gstart) || is.null(fr$gend) ||
          !length(fr$gstart)) next
      rect(fr$gstart, y_sel - 0.11, fr$gend, y_sel + 0.11,
           col = adjustcolor(probe_col, alpha.f = 0.85), border = NA)
    }
    mtext("dark blocks = designed probe targets (selected isoform)",
          side = 3, line = 0, cex = 0.6, col = "gray40", adj = 0)
  }
}

# Combined cross-species track plot: the design species' isoforms on top, then
# each checked species' transcripts, all drawn in LOCUS-RELATIVE coordinates
# (each species' locus shifted to start at 0) so exon structures line up for
# comparison. Track labels carry the species name; the transcript used for the
# probe match is saturated, other transcripts faded, and the designed probes'
# binding sites are overlaid as dark blocks on each species' selected track.
plot_combined_species_tracks <- function(res) {
  groups <- list()
  ref_by_tx <- sanitize_by_tx(res$isoform_exons)
  if (!is.null(ref_by_tx))
    groups[[length(groups) + 1]] <- list(
      species = res$species, by_tx = ref_by_tx,
      selected = res$accession, frags = res$probe_frags,
      shared = res$shared_exons,
      shared_label = if (isTRUE(res$shared_is_cds)) "shared CDS" else "shared exons")
  if (!is.null(res$cross_species))
    for (e in res$cross_species) {
      sp_by_tx <- sanitize_by_tx(e$isoform_exons)
      if (!is.null(sp_by_tx))
        groups[[length(groups) + 1]] <- list(
          species = e$species, by_tx = sp_by_tx,
          selected = e$accession, frags = e$probe_frags)
    }
  if (!length(groups)) return(invisible(NULL))
  
  # Flip minus-strand loci 5'->3' BEFORE any alignment math, so the probe-site
  # anchoring below is computed in the same orientation for every species.
  # (Flipping after the offset calculation misaligns the groups, because each
  # group would flip around its own extent.)
  flip_group <- function(g) {
    strds <- vapply(g$by_tx, function(es) {
      st <- es[[1]]$strand; if (is.null(st) || !length(st)) "" else st
    }, "")
    strd <- unname(strds[strds %in% c("+", "-")][1])
    if (is.na(strd) || !identical(strd, "-")) return(g)
    lo <- min(vapply(g$by_tx, function(es) min(vapply(es, `[[`, 0, "start")), 0))
    hi <- max(vapply(g$by_tx, function(es) max(vapply(es, `[[`, 0, "end")), 0))
    g$by_tx <- lapply(g$by_tx, function(es) lapply(es, function(x) {
      f <- flip_se(x$start, x$end, lo, hi)
      x$start <- f$s; x$end <- f$e; x$strand <- "+"; x
    }))
    if (!is.null(g$frags) && length(g$frags))
      g$frags <- lapply(g$frags, function(fr) {
        if (is.null(fr) || is.null(fr$gstart) || !length(fr$gstart)) return(fr)
        data.frame(gstart = lo + hi - fr$gend, gend = lo + hi - fr$gstart)
      })
    if (!is.null(g$shared) && nrow(g$shared)) {
      f <- flip_se(g$shared$start, g$shared$end, lo, hi)
      g$shared$start <- f$s; g$shared$end <- f$e
      if ("strand" %in% names(g$shared)) g$shared$strand <- "+"
    }
    g
  }
  groups <- lapply(groups, flip_group)
  
  pal <- c("#0ea5e9", "#8b5cf6", "#10b981", "#f59e0b", "#ef4444")
  
  # Reference = design species (group 1). Its locus is shifted to start at 0.
  # Other species are then shifted so their MATCHED PROBE SITES line up with
  # the design species' probe sites - this makes the shared-CDS guide lines
  # meaningful across species (they mark the same biological regions).
  # midpoint of each probe fragment; robust to empty/malformed fragments
  frag_mids <- function(frags) {
    if (is.null(frags) || !length(frags)) return(NULL)
    mids <- vapply(frags, function(d) {
      if (is.null(d) || is.null(d$gstart) || is.null(d$gend) ||
          !length(d$gstart) || !length(d$gend)) NA_real_
      else mean(c(d$gstart, d$gend))
    }, 0)
    mids <- mids[is.finite(mids)]
    if (!length(mids)) NULL else mids
  }
  
  off_ref <- min(vapply(groups[[1]]$by_tx,
                        function(es) min(vapply(es, `[[`, 0, "start")), 0))
  ref_anchors <- frag_mids(groups[[1]]$frags)
  if (!is.null(ref_anchors)) ref_anchors <- ref_anchors - off_ref
  
  # Flatten to a track list in display order (top first)
  tracks <- list()
  for (gi in seq_along(groups)) {
    g <- groups[[gi]]
    off <- min(vapply(g$by_tx, function(es) min(vapply(es, `[[`, 0, "start")), 0))
    if (gi == 1) off <- off_ref
    # anchor ortholog loci to the design species via the probe match sites
    anchored <- FALSE
    if (gi > 1 && !is.null(ref_anchors)) {
      sp_anchors <- frag_mids(g$frags)
      if (!is.null(sp_anchors)) {
        # pair anchors by pair NUMBER (frag names), so pairs that were not
        # found in this species don't shift the correspondence
        common <- intersect(names(ref_anchors), names(sp_anchors))
        if (length(common)) {
          off <- stats::median(sp_anchors[common] - ref_anchors[common])
          anchored <- TRUE
        } else {
          n <- min(length(sp_anchors), length(ref_anchors))
          off <- stats::median(sp_anchors[seq_len(n)] - ref_anchors[seq_len(n)])
          anchored <- TRUE
        }
      }
    }
    tx_ids <- names(g$by_tx)
    tx_ids <- tx_ids[order(vapply(g$by_tx, function(es)
      min(vapply(es, `[[`, 0, "start")), 0))]
    sel_base <- if (is.null(g$selected)) "" else sub("\\..*$", "", g$selected)
    for (k in seq_along(tx_ids)) {
      tid <- tx_ids[k]
      es <- g$by_tx[[tid]]
      is_sel <- identical(sub("\\..*$", "", tid), sel_base)
      fr <- NULL
      if (is_sel && !is.null(g$frags) && length(g$frags))
        fr <- lapply(g$frags, function(d) {
          if (is.null(d) || is.null(d$gstart) || is.null(d$gend) ||
              !length(d$gstart)) return(NULL)
          data.frame(s = d$gstart - off, e = d$gend - off)
        })
      if (!is.null(fr)) fr <- Filter(Negate(is.null), fr)
      if (!length(fr)) fr <- NULL
      st <- es[[1]]$strand
      if (is.null(st) || !length(st) || !st %in% c("+", "-")) st <- "+"
      tracks[[length(tracks) + 1]] <- list(
        species = g$species, id = tid,
        s = vapply(es, `[[`, 0, "start") - off,
        e = vapply(es, `[[`, 0, "end") - off,
        strand = st, sel = is_sel,
        col = pal[((gi - 1) %% length(pal)) + 1],
        frags = fr, group_first = (k == 1), gi = gi)
    }
    # extra track: the shared regions the probes were designed against
    if (!is.null(g$shared) && nrow(g$shared)) {
      sh_st <- if ("strand" %in% names(g$shared)) g$shared$strand[1] else "+"
      if (is.null(sh_st) || !length(sh_st) || !sh_st %in% c("+", "-"))
        sh_st <- "+"
      tracks[[length(tracks) + 1]] <- list(
        species = g$species, id = g$shared_label,
        s = g$shared$start - off, e = g$shared$end - off,
        strand = sh_st, sel = TRUE,
        col = "#16a34a", frags = NULL, group_first = FALSE, gi = gi)
    }
  }
  
  n_tracks <- length(tracks)
  n_groups <- length(groups)
  gap <- 0.9
  total_h <- n_tracks + gap * (n_groups - 1)
  x_max <- max(vapply(tracks, function(t) max(t$e), 0))
  x_min <- min(0, vapply(tracks, function(t) min(t$s), 0))
  pad <- 0.02 * max(x_max - x_min, 1)
  
  op <- par(mar = c(3, 24, 1, 1)); on.exit(par(op))
  plot(0, type = "n", xlim = c(x_min - pad, x_max + pad),
       ylim = c(0.5, total_h + 0.5), axes = FALSE, xlab = "", ylab = "")
  axis(1, col = "gray40", col.axis = "gray40", cex.axis = 0.8)
  mtext("aligned position (bp) - all loci shown 5' to 3'; ortholog loci anchored at probe match sites",
        side = 1, line = 1.9, cex = 0.7, col = "gray40")
  
  # Faint vertical guides at the design species' shared-region boundaries, so
  # you can see whether ortholog matches fall inside the conserved blocks.
  sh <- Filter(function(t) t$gi == 1 && t$id %in% c("shared CDS", "shared exons"),
               tracks)
  if (length(sh)) {
    for (b in unique(c(sh[[1]]$s, sh[[1]]$e)))
      segments(b, 0.5, b, total_h + 0.5, col = adjustcolor("#16a34a", 0.18),
               lty = 3, lwd = 1)
  }
  
  # ---- PASS 1: compute y positions for every track ----
  cursor <- total_h + 1
  prev_gi <- 0
  for (i in seq_along(tracks)) {
    t <- tracks[[i]]
    if (t$gi != prev_gi) {
      cursor <- cursor - if (prev_gi == 0) 1 else gap + 1
      prev_gi <- t$gi
    } else {
      cursor <- cursor - 1
    }
    tracks[[i]]$y <- cursor
  }
  
  # ---- SYNTENY RIBBONS: connect matched probe fragments across species ----
  # Quadratic Bezier through start / control / end points
  qbezier <- function(p0, p1, p2, n = 40) {
    t <- seq(0, 1, length.out = n)
    x <- (1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1]
    y <- (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2]
    list(x = x, y = y)
  }
  
  for (gi in 2:length(groups)) {
    sel_prev <- Filter(function(t) t$sel && t$gi == gi - 1, tracks)
    sel_curr <- Filter(function(t) t$sel && t$gi == gi, tracks)
    if (!length(sel_prev) || !length(sel_curr)) next
    t1 <- sel_prev[[1]]; t2 <- sel_curr[[1]]
    if (is.null(t1$frags) || is.null(t2$frags)) next
    
    common <- intersect(names(t1$frags), names(t2$frags))
    if (!length(common)) next
    
    for (fn in common) {
      f1 <- t1$frags[[fn]]; f2 <- t2$frags[[fn]]
      x1 <- mean(c(f1$s, f1$e))
      x2 <- mean(c(f2$s, f2$e))
      y1 <- t1$y; y2 <- t2$y
      
      # Control point sits midway in x, midway in y (gives a smooth S-curve)
      xc <- (x1 + x2) / 2
      yc <- (y1 + y2) / 2
      bz <- qbezier(c(x1, y1), c(xc, yc), c(x2, y2))
      # Use the upstream group's colour, very faint so exons stay readable
      lines(bz$x, bz$y, col = adjustcolor(t1$col, alpha.f = 0.22), lwd = 2.2)
    }
  }
  
  # ---- PASS 2: draw tracks ----
  for (t in tracks) {
    y <- t$y
    # group separator line
    if (t$group_first && t$gi > 1)
      abline(h = y + 0.5 + gap / 2, col = "gray88", lty = 2)
    # intron line + direction arrows
    segments(min(t$s), y, max(t$e), y, col = "gray70", lwd = 2, lend = 2)
    if (t$strand %in% c("+", "-")) {
      xs <- seq(min(t$s), max(t$e), length.out = 7)[2:6]
      text(xs, rep(y, length(xs)), labels = if (t$strand == "+") ">" else "<",
           col = "gray55", cex = 0.55)
    }
    col <- if (t$sel) t$col else adjustcolor(t$col, alpha.f = 0.4)
    rect(t$s, y - 0.18, t$e, y + 0.18, col = col, border = NA)
    # probe binding sites on the selected track
    if (!is.null(t$frags))
      for (frl in t$frags)
        rect(frl$s, y - 0.11, frl$e, y + 0.11,
             col = adjustcolor("#0f172a", alpha.f = 0.9), border = NA)
    text(par("usr")[1], y, labels = paste0(t$species, " | ", t$id),
         adj = c(1, 0.5), xpd = NA, cex = 0.6,
         font = if (t$sel) 2 else 1, col = if (t$sel) t$col else "gray30")
  }
  mtext("dark blocks = designed probe binding sites (selected transcript per species); faint curves = synteny between matched probe sites",
        side = 3, line = 0, cex = 0.6, col = "gray40", adj = 0)
}


# Reorder a transcript's exons into 5'->3' (transcript) order as a data.frame.
# Minus-strand transcripts run from high to low genomic coordinates.
exons_to_tx_segments <- function(exons) {
  df <- data.frame(start  = vapply(exons, `[[`, 0, "start"),
                   end    = vapply(exons, `[[`, 0, "end"),
                   strand = vapply(exons, `[[`, "", "strand"),
                   stringsAsFactors = FALSE)
  df <- if (df$strand[1] == "-") df[order(-df$start), ] else df[order(df$start), ]
  rownames(df) <- NULL
  df
}

# Map probe pairs from design-sequence coordinates onto genomic coordinates.
# segs must be in design-sequence (5'->3') order; a probe spanning a splice
# junction yields one fragment per exon it overlaps.
map_pairs_to_genome <- function(pairs, segs) {
  lens  <- segs$end - segs$start + 1
  cum0  <- c(0, cumsum(lens))
  total <- sum(lens)
  frags <- list()
  # name fragments by pair index so cross-species anchors pair up correctly
  # even when some pairs were not found in the ortholog
  nms <- names(pairs)
  if (is.null(nms) || any(!nzchar(nms))) nms <- as.character(seq_along(pairs))
  for (j in seq_along(pairs)) {
    p <- pairs[[j]]
    ps <- p$start; pe <- p$end
    if (ps < 1 || pe > total) next
    rows <- NULL
    for (k in seq_len(nrow(segs))) {
      ov_s <- max(ps, cum0[k] + 1); ov_e <- min(pe, cum0[k + 1])
      if (ov_s > ov_e) next
      if (segs$strand[k] == "-") {
        g_s <- segs$end[k] - (ov_e - cum0[k] - 1)
        g_e <- segs$end[k] - (ov_s - cum0[k] - 1)
      } else {
        g_s <- segs$start[k] + (ov_s - cum0[k] - 1)
        g_e <- segs$start[k] + (ov_e - cum0[k] - 1)
      }
      rows <- rbind(rows, data.frame(gstart = min(g_s, g_e), gend = max(g_s, g_e)))
    }
    if (!is.null(rows)) frags[[nms[j]]] <- rows
  }
  frags
}

# ---- Thermodynamics (SantaLucia 1998 RNA/DNA nearest-neighbor) --------------
# Same parameters and formulas as the Jeff Lee Lab HCRv3 designer.

NN_dH <- c(AA=-7.8, AC=-5.9, AG=-9.1, AT=-8.3, CA=-9.0, CC=-9.3, CG=-16.3,
           CT=-7.0, GA=-5.5, GC=-8.0, GG=-12.8, GT=-7.8, TA=-7.8, TC=-8.6,
           TG=-10.4, TT=-11.5)
NN_dS <- c(AA=-21.9, AC=-12.3, AG=-23.5, AT=-23.9, CA=-26.1, CC=-23.2,
           CG=-47.1, CT=-19.7, GA=-13.5, GC=-17.1, GG=-31.9, GT=-21.6,
           TA=-23.2, TC=-22.9, TG=-28.4, TT=-36.4)

# Returns list(Tm, dG, dH, dS) for a target RNA/DNA sequence.
calc_thermo <- function(seq, temp_c = 37, na_m = 0.3, oligo_conc = 5e-5) {
  n <- nchar(seq)
  if (n < 2) return(list(Tm = NA_real_, dG = NA_real_, dH = NA_real_, dS = NA_real_))
  din <- substring(seq, 1:(n - 1), 2:n)
  dH <- sum(NN_dH[din], na.rm = TRUE) + 1.9     # kcal/mol, + initiation
  dS <- sum(NN_dS[din], na.rm = TRUE) - 3.9     # cal/(mol K)
  dG <- (dH * 1000 - (temp_c + 273.15) * dS) / 1000
  Tm <- (dH * 1000 / (dS + 1.9872 * log(oligo_conc / 4))) - 273.15 +
        16.6 * log10(na_m)
  list(Tm = Tm, dG = dG, dH = dH, dS = dS)
}

# ---- Composition rules (Choi et al. / Jeff Lee Lab) -------------------------

comp_a_ok  <- function(s) nchar(gsub("[^A]", "", s)) / nchar(s) < 0.28
comp_c_ok  <- function(s) { f <- nchar(gsub("[^C]", "", s)) / nchar(s)
                            f > 0.22 & f < 0.28 }
no_astack  <- function(s) !grepl("AAAA", s)
no_cstack  <- function(s) !grepl("CCCC", s)
# C-stacking: no 6-nt window among the first 12 nt with >50% C
c_spec_ok  <- function(s) {
  for (i in 1:7) {
    w <- substr(s, i, i + 5)
    if (nchar(w) == 6 && nchar(gsub("[^C]", "", w)) / 6 > 0.5) return(FALSE)
  }
  TRUE
}

# Annotate candidates with thermodynamics + composition, then score.
# score = -(deviation from target dG); higher is better.
annotate_candidates <- function(cands, settings) {
  # Fall back to defaults if a setting is missing: assigning numeric(0) to a
  # list element would silently DELETE $score and crash dp_select downstream.
  tdg  <- if (is.null(settings$target_dg)) -80 else settings$target_dg
  tdgh <- if (is.null(settings$target_dg_half)) -38 else settings$target_dg_half
  for (i in seq_along(cands)) {
    w  <- cands[[i]]$target52
    h1 <- substr(w, 1, 25)
    h2 <- substr(w, 28, 52)
    tf <- calc_thermo(w,  settings$temp_c, settings$na_m, settings$oligo_conc)
    t1 <- calc_thermo(h1, settings$temp_c, settings$na_m, settings$oligo_conc)
    t2 <- calc_thermo(h2, settings$temp_c, settings$na_m, settings$oligo_conc)
    cands[[i]]$dg  <- tf$dG; cands[[i]]$tm <- tf$Tm
    cands[[i]]$dg1 <- t1$dG; cands[[i]]$dg2 <- t2$dG
    cands[[i]]$comp_ok <-
      (!settings$comp_a  || (comp_a_ok(h1)  && comp_a_ok(h2)))  &&
      (!settings$comp_c  || (comp_c_ok(h1)  && comp_c_ok(h2)))  &&
      (!settings$comp_stack || (no_astack(h1) && no_astack(h2) &&
                                no_cstack(h1)  && no_cstack(h2))) &&
      (!settings$comp_cstack || (c_spec_ok(h1) && c_spec_ok(h2)))
    cands[[i]]$score <- -(abs(tf$dG - tdg) +
                          abs(t1$dG - tdgh) +
                          abs(t2$dG - tdgh))
  }
  cands
}

# Keep candidates inside the dG/Tm ranges (only when filter_thermo is on)
# and passing the individually-enabled composition rules.
filter_scored <- function(cands, settings) {
  if (!length(cands)) return(cands)
  keep <- vapply(cands, function(c) {
    (!isTRUE(settings$filter_thermo) ||
       ((is.na(c$dg) || (c$dg >= settings$dg_min && c$dg <= settings$dg_max)) &&
        (is.na(c$tm) || (c$tm >= settings$tm_min && c$tm <= settings$tm_max)))) &&
    isTRUE(c$comp_ok)
  }, TRUE)
  cands[keep]
}

# ---- Optimal probe-set selection (dynamic programming) ----------------------
# Weighted interval scheduling with a count cap: choose up to n_want
# non-overlapping candidates (min `spacing` bp between windows) maximizing
# total score, i.e. best-matched thermodynamics at even coverage.

dp_select <- function(cands, n_want, spacing) {
  m <- length(cands)
  if (m <= 1) return(cands)
  starts <- vapply(cands, `[[`, 0, "start")
  ends   <- vapply(cands, `[[`, 0, "end")
  scores <- vapply(cands, function(c)
    if (is.null(c$score) || !length(c$score) || !is.finite(c$score)) 0
    else c$score, 0)
  ord <- order(starts)
  starts <- starts[ord]; ends <- ends[ord]; scores <- scores[ord]
  cands <- cands[ord]
  # large positive count bonus: always prefer more probes, tie-break on quality
  scores <- 1e4 + scores
  # p[i] = largest j with ends[j] <= starts[i] - spacing (0 if none)
  p <- findInterval(starts - spacing, ends)
  p <- pmin(p, seq_len(m) - 1L)
  K <- min(n_want, m)
  neg <- -1e15
  dp <- matrix(neg, nrow = m + 1, ncol = K + 1)
  dp[, 1] <- 0
  take <- matrix(FALSE, nrow = m, ncol = K)
  for (i in seq_len(m)) {
    for (k in 1:K) {
      with_i <- if (p[i] >= 0) dp[p[i] + 1, k] + scores[i] else neg
      if (with_i > dp[i, k + 1]) {   # dp[i, k+1] currently = best of first i-1
        dp[i + 1, k + 1] <- with_i
        take[i, k] <- TRUE
      } else {
        dp[i + 1, k + 1] <- dp[i, k + 1]
      }
    }
  }
  # backtrack from the best feasible count (may be < K when windows overlap)
  sel <- integer(0)
  i <- m; k <- which.max(dp[m + 1, ]) - 1L
  while (i >= 1 && k >= 1) {
    if (take[i, k]) { sel <- c(i, sel); i <- p[i]; k <- k - 1 }
    else i <- i - 1
    if (i == 0) break
  }
  cands[sel]
}

# ---- Off-target screening (pure-R k-mer seed scan) --------------------------
# Approximates the Jeff Lee Lab BLAST screen without needing BLAST+: a
# transcript is flagged when it contains an exact k-mer from BOTH 25-nt
# halves of a probe (split probes need both halves to initiate HCR, so a
# single-half match is only a weak risk and is reported but not dropped).

get_ot_index <- function(species, cache) {
  key <- paste0(".otindex.", tolower(species))
  if (exists(key, envir = cache)) return(get(key, envir = cache))
  idx <- load_species(species, cache)
  ot <- NULL
  if (!is.null(idx)) {
    keys <- setdiff(ls(idx), ".byid")
    seen <- new.env(hash = TRUE)
    tx_id <- character(0); gene <- character(0); seqs <- character(0)
    for (gk in keys) {
      for (r in idx[[gk]]) {
        if (exists(r$id, envir = seen)) next
        assign(r$id, TRUE, envir = seen)
        tx_id <- c(tx_id, r$id); gene <- c(gene, gk); seqs <- c(seqs, r$seq)
      }
    }
    if (length(seqs)) {
      big <- paste(seqs, collapse = "\f")
      tx_start <- c(1, cumsum(nchar(seqs) + 1)[-length(seqs)] + 1)
      ot <- list(big = big, tx_start = tx_start, gene = gene, tx_id = tx_id)
    }
  }
  assign(key, list(ot = ot), envir = cache)   # cache even NULL result
  list(ot = ot)
}

# Distinct genes whose transcripts contain any k-mer of seq (exact match).
# exclude_tx: accession bases (e.g. "NM_001234") of the target gene's own
# transcripts, which would otherwise flag under protein-name index keys.
ot_genes_hit <- function(ot, seq, k, exclude_tx = character(0)) {
  n <- nchar(seq)
  seeds <- unique(substring(seq, 1:(n - k + 1), k:n))
  hits <- gregexpr(paste(seeds, collapse = "|"), ot$big)[[1]]
  if (hits[1] == -1) return(character(0))
  tx <- findInterval(hits, ot$tx_start)
  keep <- !(sub("\\..*$", "", ot$tx_id[tx]) %in% exclude_tx)
  unique(ot$gene[tx[keep]])
}

# Annotate one candidate with off-target genes (both halves) and weak hits
# (one half only). exclude: lowercase gene keys to ignore (the target gene).
ot_screen_candidate <- function(cand, ot, k, exclude = character(0),
                                exclude_tx = character(0)) {
  h1 <- substr(cand$target52, 1, 25)
  h2 <- substr(cand$target52, 28, 52)
  g1 <- setdiff(ot_genes_hit(ot, h1, k, exclude_tx), exclude)
  g2 <- setdiff(ot_genes_hit(ot, h2, k, exclude_tx), exclude)
  cand$ot_genes <- intersect(g1, g2)
  cand$ot_weak  <- length(union(g1, g2)) - length(cand$ot_genes)
  cand$ot_n     <- length(cand$ot_genes)
  cand
}

# ---- Probe design engine ----------------------------------------------------

find_candidates <- function(seq, spacing, gc_min, gc_max, max_hp, exon_bounds = NULL,
                            stride = NULL) {
  n <- nchar(seq)
  if (is.null(stride)) stride <- 52 + spacing
  stride <- max(1, stride)
  starts <- seq(1, n - 51, by = stride)
  out <- list()
  for (st in starts) {
    w <- substr(seq, st, st + 51)
    if (grepl("N", w)) next

    # --- Splice junction filter: 52-bp window must lie entirely within one exon ---
    if (!is.null(exon_bounds) && nrow(exon_bounds) > 0) {
      # Check if this window [st, st+51] is fully contained in any single exon
      en <- st + 51
      inside_one_exon <- any(
        exon_bounds$start <= st & exon_bounds$end >= en
      )
      if (!inside_one_exon) next  # window spans a junction; skip it
    }

    half_a <- substr(w, 1, 25)
    half_b <- substr(w, 28, 52)
    ga <- gc_percent(half_a); gb <- gc_percent(half_b)
    if (ga < gc_min || ga > gc_max || gb < gc_min || gb > gc_max) next
    if (max_run(half_a) > max_hp || max_run(half_b) > max_hp) next
    out[[length(out) + 1]] <- list(start = st, end = st + 51, target52 = w,
                                   gc1 = ga, gc2 = gb)
  }
  out
}

pick_evenly <- function(cands, n) {
  if (length(cands) <= n) return(cands)
  idx <- unique(round(seq(1, length(cands), length.out = n)))
  cands[idx]
}

build_pair <- function(gene, amp, pair_index, cand) {
  a <- AMPLIFIERS[[amp]]
  rc <- revcomp(cand$target52)

  # --- Intermediate components for visual display ---
  # The 52-bp target: two 25-mers + 2-bp gap
  half_a <- substr(cand$target52, 1, 25)      # 5' half of target (sense strand)
  half_b <- substr(cand$target52, 28, 52)     # 3' half of target (sense strand)
  gap    <- substr(cand$target52, 26, 27)      # 2-bp gap

  # Reverse complement of each half (these are what the probes bind to)
  rc_half_a <- revcomp(half_a)   # 25-nt antisense for odd probe
  rc_half_b <- revcomp(half_b)    # 25-nt antisense for even probe

  # Odd probe = initiator half (18 nt) + spacer (2 nt) + antisense half_a (25 nt)
  odd_antisense <- substr(rc, nchar(rc) - 24, nchar(rc))  # = rc_half_a
  odd  <- paste0(a$half1, a$spacer1, odd_antisense)

  # Even probe = antisense half_b (25 nt) + spacer (2 nt) + initiator half (18 nt)
  even_antisense <- substr(rc, 1, 25)  # = rc_half_b
  even <- paste0(even_antisense, a$spacer2, a$half2)

  pad <- sprintf("%02d", pair_index)
  gv <- function(x) if (is.null(x)) NA_real_ else x
  list(pair = pair_index, start = cand$start, end = cand$end,
       target52 = cand$target52, gc1 = cand$gc1, gc2 = cand$gc2,
       dg = gv(cand$dg), tm = gv(cand$tm), dg1 = gv(cand$dg1), dg2 = gv(cand$dg2),
       ot_n = if (is.null(cand$ot_n)) NA_integer_ else cand$ot_n,
       ot_genes = if (is.null(cand$ot_genes) || !length(cand$ot_genes)) ""
                  else paste(cand$ot_genes, collapse = "; "),
       ot_ok = if (is.null(cand$ot_ok)) NA else cand$ot_ok,
       # --- Raw components ---
       half_a = half_a, half_b = half_b, gap = gap,
       rc_half_a = rc_half_a, rc_half_b = rc_half_b,
       # --- Amplifier components ---
       amp = amp,
       odd_initiator = a$half1, odd_spacer = a$spacer1,
       even_initiator = a$half2, even_spacer = a$spacer2,
       # --- Final probes ---
       odd_name = paste0(gene, "-", amp, "-", pad, "a"), odd = odd,
       even_name = paste0(gene, "-", amp, "-", pad, "b"), even = even)
}

# Resolve sequence for one gene job, then design its probe pairs.
design_gene <- function(job, settings, index, cache) {
  amps <- active_amps(settings)
  # Explicit CSV "amplifier" column wins, but only if that amplifier is in the
  # selected set (v2 mode = B1-B5 only). A v3 amplifier named in the CSV while
  # v2 mode is active falls back to auto B1-B5 rotation.
  amp <- if (!is.null(job$amplifier) && job$amplifier %in% amps) job$amplifier
         else if (settings$default_amp == "auto") amps[(index %% length(amps)) + 1]
         else settings$default_amp
  res <- list(gene = job$gene, species = job$species, amplifier = amp, status = "ok",
              warnings = character(0), accession = NA, seq_len = 0,
              n_candidates = 0, pairs = list(), error = NULL,
              all_isoforms = NULL, shared_exons = NULL, used_shared = FALSE,
              isoform_exons = NULL, probe_frags = NULL, design_seq_len = 0,
              struct_source = NULL, map_coverage = NULL,
              settings_used = settings,
              protein_name = NULL,
              cds_by_tx = NULL, shared_is_cds = FALSE,
              avoid_splice_junctions = settings$avoid_splice_junctions)
  tryCatch({
    rec <- NULL
    recs <- NULL
    if (!is.null(job$sequence) && nchar(clean_sequence(job$sequence)) >= 60) {
      rec <- list(id = "user-provided", seq = clean_sequence(job$sequence),
                  title = "pasted sequence", variant = NA, symbol = job$gene)
    } else if (!is.null(job$transcript) &&
               grepl("^[A-Z]{2,3}_[0-9]+(\\.[0-9]+)?$", trimws(job$transcript))) {
      # explicit accession: exact match in local files, else fetch from NCBI
      acc <- trimws(job$transcript)
      acc_base <- sub("\\..*$", "", acc)   # ignore version suffix when matching
      idx <- load_species(job$species, cache)
      if (!is.null(idx) && !is.null(idx[[".byid"]]))
        rec <- idx[[".byid"]][[tolower(acc_base)]]
      if (is.null(rec)) {
        if (!settings$ncbi_fallback)
          stop(paste0("Accession ", acc, " not in local files and NCBI fallback is off."))
        fa <- tryCatch(read_fasta_from_accession(acc), error = function(e) NULL)
        if (is.null(fa) || !length(fa)) stop(paste("Could not fetch accession", acc))
        rec <- fa[[1]]
      }
    } else {
      idx <- load_species(job$species, cache)
      recs <- if (is.null(idx)) NULL else idx[[tolower(job$gene)]]
      if (is.null(recs) && !is.null(idx)) {
        # fallback: word-boundary search over gene/protein names
        # (e.g. "ovo" matches "transcriptional regulator ovo")
        pat <- paste0("\\b", gsub("([^A-Za-z0-9])", "\\\\\\1", tolower(job$gene)), "\\b")
        keys <- setdiff(ls(idx), ".byid")
        hits <- keys[vapply(keys, function(k) grepl(pat, k), TRUE)]
        if (length(hits))
          recs <- unique(unlist(lapply(hits, function(k) idx[[k]]), recursive = FALSE))
      }
      if (is.null(recs)) {
        if (!settings$ncbi_fallback)
          stop(paste0('Gene "', job$gene, '" not found in local "', job$species,
                      '" file, and NCBI fallback is off. Note: non-melanogaster ',
                      "RefSeq files are annotated by protein name/LOC ID, not gene ",
                      "symbol - try the full gene name, put an XM_/NM_ accession in ",
                      "the transcript column, or paste the sequence. ",
                      "(Searched in: ", DATA_DIRS, ")"))
        recs <- ncbi_fetch_gene(job$gene, organism_for(job$species))
        if (is.null(recs))
          stop(paste0('No RefSeq mRNA found for "', job$gene, '" in ',
                      organism_for(job$species),
                      ". Paste the transcript sequence into the CSV instead."))
      }
      picked <- choose_transcript(recs, job$transcript)
      rec <- picked$rec
      if (!is.null(picked$warning)) res$warnings <- c(res$warnings, picked$warning)
      res$all_isoforms <- picked$all_isoforms
    }

    res$accession <- rec$id
    res$seq_len <- nchar(rec$seq)
    res$seq <- rec$seq          # needed by the FASTA export (full transcript)
    # Protein/gene-product name of the design transcript (e.g. "homeotic
    # protein abdominal-B") - used to find the ortholog in other species'
    # RefSeq files, which are indexed by protein name rather than symbol.
    if (!is.null(rec$name) && !is.na(rec$name) && nzchar(rec$name))
      res$protein_name <- rec$name

    # --- Exon structures -----------------------------------------------------
    # Preferred source: sequence alignment of every isoform onto the genome
    # (verifiable and independent of annotation quirks). Fallback: GTF
    # coordinates. All of this is best-effort and must never kill probe design.
    genome_path <- find_genome_file(job$species)
    gtf_path <- find_gtf_file(job$species)
    by_tx <- NULL
    all_recs <- if (!is.null(recs) && length(recs)) recs else list(rec)
    if (!is.null(genome_path)) {
      locus <- NULL
      if (!is.null(gtf_path))
        locus <- tryCatch(gtf_locus(gtf_path, job$gene, cache),
                          error = function(e) NULL)
      if (is.null(locus))
        locus <- tryCatch(find_locus_by_sequence(genome_path, rec$seq, cache),
                          error = function(e) NULL)
      if (!is.null(locus)) {
        loc <- tryCatch(get_locus_seq(genome_path, locus, cache = cache),
                        error = function(e) NULL)
        mapped <- if (!is.null(loc))
          tryCatch(map_isoforms_seq(all_recs, loc), error = function(e) NULL)
          else NULL
        if (!is.null(mapped) && length(mapped$by_tx)) {
          by_tx <- mapped$by_tx
          res$struct_source <- "sequence alignment"
          res$map_coverage <- mapped$coverage
          if (length(mapped$by_tx) < length(all_recs))
            res$warnings <- c(res$warnings, paste0(
              length(all_recs) - length(mapped$by_tx), " of ", length(all_recs),
              " transcripts could not be aligned to the genome (>30% unmappable) ",
              "and were excluded from the isoform plot and shared-exon set."))
        }
      }
    }
    if (is.null(by_tx) && !is.null(gtf_path)) {
      by_tx <- tryCatch(parse_gtf_exons(gtf_path, job$gene, cache = cache),
                        error = function(e) NULL)
      if (!is.null(by_tx)) res$struct_source <- "GTF"
    }
    if (!is.null(by_tx)) res$isoform_exons <- by_tx

    # CDS annotation for this gene (for the plot's UTR/CDS rendering and for
    # shared-CDS design); NULL when the GTF lacks this gene.
    cds_by_tx <- NULL
    if (!is.null(by_tx) && !is.null(gtf_path))
      cds_by_tx <- tryCatch(gtf_cds_by_tx(gtf_path, job$gene, cache),
                            error = function(e) NULL)
    res$cds_by_tx <- cds_by_tx

    # If requested, use shared exons across all isoforms
    seq_for_design <- rec$seq
    if (isTRUE(settings$use_shared_exons) && !is.null(res$all_isoforms) && nrow(res$all_isoforms) >= 1) {
      genome_path <- find_genome_file(job$species)
      if (!is.null(by_tx) && !is.null(genome_path)) {
        shared <- find_shared_exons(by_tx)
        if (nrow(shared)) {
          shared <- merge_overlaps(shared)
          # Restrict to regions CODING in all isoforms when CDS annotation
          # exists; otherwise keep the full shared exons (UTRs included).
          if (!is.null(cds_by_tx)) {
            shared_cds <- tryCatch(shared_cds_regions(by_tx, cds_by_tx),
                                   error = function(e) NULL)
            if (!is.null(shared_cds) && nrow(shared_cds)) {
              shared <- merge_overlaps(shared_cds)
              res$shared_is_cds <- TRUE
            } else {
              res$warnings <- c(res$warnings, paste0(
                "No coding region is common to all ", length(by_tx),
                " isoforms; using shared exons (UTRs included) instead."))
            }
          } else {
            res$warnings <- c(res$warnings,
              "No CDS annotation for this gene in the GTF; shared regions include UTRs.")
          }
          # order regions 5'->3' for the gene's strand, so the concatenated
          # sequence matches the mRNA (minus strand runs high->low coords)
          if (shared$strand[1] == "-") shared <- shared[order(-shared$start), ]
          shared_seq <- extract_regions(genome_path, shared, cache = cache)
          if (nchar(shared_seq) >= 60) {
            seq_for_design <- shared_seq
            res$used_shared <- TRUE
            res$shared_exons <- shared
            # --- Rebuild by_tx from shared_exons for consistent plotting ---
            # When there's only 1 isoform, by_tx from sequence alignment may
            # have slightly different coordinates than shared_exons (due to
            # alignment wobble / merge_near_blocks). Rebuild by_tx so the
            # isoform plot and probe mapping use identical exon boundaries.
            if (!is.null(by_tx) && length(by_tx) == 1) {
              tx_id <- names(by_tx)[1]
              by_tx[[tx_id]] <- lapply(seq_len(nrow(shared)), function(i) {
                list(chr = shared$chr[i], start = shared$start[i],
                     end = shared$end[i], strand = shared$strand[i],
                     transcript_id = tx_id)
              })
              res$isoform_exons <- by_tx
            }
            res$warnings <- c(res$warnings, paste0(
              "Using shared ", if (isTRUE(res$shared_is_cds)) "CDS" else "exon",
              " sequence (", nrow(shared), " regions, ",
              nchar(shared_seq), " nt) across ", nrow(res$all_isoforms), " isoforms."))
          } else {
            res$warnings <- c(res$warnings, paste0(
              "Shared exons only ", nchar(shared_seq), " nt; using full transcript instead."))
          }
        } else {
          res$warnings <- c(res$warnings, "No shared exons found across isoforms; using full transcript.")
        }
      } else {
        missing <- c(if (is.null(genome_path)) "genomic.fna" else NULL,
                     if (is.null(gtf_path)) "genomic.gtf"
                     else if (is.null(by_tx)) paste0('exons for "', job$gene, '" in genomic.gtf')
                     else NULL)
        res$warnings <- c(res$warnings, paste0(
          "Cannot compute shared exons - missing: ", paste(missing, collapse = " and "),
          ". (Looked in ", file.path(DATA_DIRS, paste0("d_", tolower(job$species))), ")"))
      }
    }

    # --- Build exon boundary map for splice-junction avoidance ---
    exon_bounds <- NULL
    if (isTRUE(settings$avoid_splice_junctions)) {
      if (isTRUE(res$used_shared) && !is.null(res$shared_exons) && nrow(res$shared_exons) > 0) {
        # Shared exons: each region is an "exon" in the concatenated design sequence.
        # Cumulative positions in seq_for_design.
        bounds <- data.frame(start = integer(0), end = integer(0), stringsAsFactors = FALSE)
        pos <- 1
        for (i in seq_len(nrow(res$shared_exons))) {
          reg_len <- res$shared_exons$end[i] - res$shared_exons$start[i] + 1
          bounds <- rbind(bounds, data.frame(start = pos, end = pos + reg_len - 1,
                                             stringsAsFactors = FALSE))
          pos <- pos + reg_len  # no spacer; design sequence is pure concatenation
        }
        exon_bounds <- bounds
      } else if (!is.null(by_tx)) {
        # Full transcript: need transcript-relative exon coordinates.
        # Get the selected transcript's exons in 5'->3' order.
        ti <- match(sub("\\..*$", "", rec$id), sub("\\..*$", "", names(by_tx)))
        if (!is.na(ti)) {
          segs <- exons_to_tx_segments(by_tx[[ti]])
          # These are genomic coords; convert to transcript-relative coords
          # by cumulative sum of exon lengths.
          bounds <- data.frame(start = integer(0), end = integer(0), stringsAsFactors = FALSE)
          pos <- 1
          for (j in seq_len(nrow(segs))) {
            ex_len <- segs$end[j] - segs$start[j] + 1
            bounds <- rbind(bounds, data.frame(start = pos, end = pos + ex_len - 1,
                                               stringsAsFactors = FALSE))
            pos <- pos + ex_len
          }
          exon_bounds <- bounds
        }
      }
    }

    # Denser stride (13 bp) when quality scoring is on so the DP selector has
    # more windows to choose from; original coarse stride otherwise.
    stride <- if (isTRUE(settings$use_thermo)) 13 else NULL
    cands <- find_candidates(seq_for_design, settings$spacing, settings$gc_min,
                             settings$gc_max, settings$max_hp, exon_bounds,
                             stride = stride)
    # short/AT-rich transcripts may yield few windows - relax one step if allowed
    if (isTRUE(settings$relax_filters) && length(cands) < settings$pairs_per_gene) {
      gc_min2 <- max(settings$gc_min - 5, 0)
      gc_max2 <- min(settings$gc_max + 5, 100)
      hp2 <- settings$max_hp + 1
      cands2 <- find_candidates(seq_for_design, settings$spacing, gc_min2, gc_max2, hp2, exon_bounds,
                                stride = stride)
      if (length(cands2) > length(cands)) {
        res$warnings <- c(res$warnings, paste0(
          "Filters auto-relaxed to GC ", gc_min2, "\u2013", gc_max2,
          "%, homopolymer \u2264 ", hp2, " to find more windows (",
          length(cands), " \u2192 ", length(cands2), ")."))
        cands <- cands2
      }
    }

    # --- Thermodynamic scoring + composition rules ---------------------------
    if (isTRUE(settings$use_thermo) && length(cands)) {
      cands_all <- annotate_candidates(cands, settings)
      n_before <- length(cands_all)
      cands <- filter_scored(cands_all, settings)
      if (!length(cands) && isTRUE(settings$relax_filters) && n_before > 0) {
        # relax: widen dG/Tm ranges and drop composition rules
        s2 <- settings
        s2$dg_min <- settings$dg_min - 15; s2$dg_max <- settings$dg_max + 15
        s2$tm_min <- settings$tm_min - 10; s2$tm_max <- settings$tm_max + 10
        s2$comp_a <- s2$comp_c <- s2$comp_stack <- s2$comp_cstack <- FALSE
        cands <- filter_scored(cands_all, s2)
        if (length(cands)) {
          res$settings_used <- s2          # record the relaxed criteria actually applied
          res$warnings <- c(res$warnings,
            "dG/Tm ranges widened and composition rules dropped to find enough windows.")
        }
      }
      if (length(cands) && length(cands) < n_before)
        res$warnings <- c(res$warnings, paste0(
          n_before - length(cands), " windows removed by dG/Tm/composition filters."))
    } else if (length(cands)) {
      # still attach thermo values for reporting, but don't filter
      cands <- annotate_candidates(cands, settings)
    }
    res$n_candidates <- length(cands)
    if (!length(cands)) {
      if (isTRUE(settings$avoid_splice_junctions) && !is.null(exon_bounds)) {
        stop("No valid 52-bp windows passed the filters while avoiding splice junctions. Try unchecking 'Avoid splice junctions', relaxing the dG/Tm ranges, or reducing pairs per gene.")
      } else {
        stop("No valid 52-bp windows passed the filters.")
      }
    }

    # --- Optimal set selection (DP) with off-target screen -------------------
    picked <- if (isTRUE(settings$use_thermo))
      dp_select(cands, settings$pairs_per_gene, settings$spacing)
    else
      pick_evenly(cands, settings$pairs_per_gene)

    if (isTRUE(settings$offtarget_screen) && length(picked)) {
      ot <- tryCatch(get_ot_index(job$species, cache)$ot, error = function(e) NULL)
      if (!is.null(ot)) {
        exclude <- tolower(job$gene)
        # Exclude every transcript of the target gene, whichever index key
        # (symbol vs full gene/protein name) it was filed under.
        exclude_tx <- sub("\\..*$", "", rec$id)
        if (!is.null(recs) && length(recs))
          exclude_tx <- c(exclude_tx, sub("\\..*$", "", vapply(recs, `[[`, "", "id")))
        if (!is.null(res$all_isoforms) && nrow(res$all_isoforms)) {
          ai <- res$all_isoforms
          acc_col <- which(tolower(names(ai)) %in% c("accession", "id"))
          if (length(acc_col))
            exclude_tx <- c(exclude_tx, sub("\\..*$", "", ai[[acc_col[1]]]))
        }
        exclude_tx <- unique(exclude_tx)
        for (round in 1:3) {
          picked <- lapply(picked, ot_screen_candidate, ot = ot,
                           k = settings$ot_seed, exclude = exclude,
                           exclude_tx = exclude_tx)
          bad <- vapply(picked, function(p) p$ot_n > settings$max_offtarget, TRUE)
          if (!any(bad)) break
          if (round == 3) {
            # Per-pair explicitness: mark offenders so downstream can see which
            # pairs exceed the cap even though cleaner windows were unavailable.
            for (i in which(bad)) picked[[i]]$ot_ok <- FALSE
            for (i in which(!bad)) picked[[i]]$ot_ok <- TRUE
            res$warnings <- c(res$warnings, paste0(
              sum(bad), " pair(s) have off-target genes (", 
              paste(unique(unlist(lapply(picked[bad], `[[`, "ot_genes"))), collapse = "; "),
              ") but no cleaner windows were available: pair(s) ",
              paste(vapply(picked[bad], function(p) p$start, 0), collapse = ", "),
              " (design-sequence start positions)."))
            break
          }
          # drop offenders and re-select
          drop_ids <- vapply(picked[bad], `[[`, 0, "start")
          cands <- Filter(function(c) !(c$start %in% drop_ids), cands)
          if (!length(cands)) break
          picked <- if (isTRUE(settings$use_thermo))
            dp_select(cands, settings$pairs_per_gene, settings$spacing)
          else
            pick_evenly(cands, settings$pairs_per_gene)
        }
      } else {
        res$warnings <- c(res$warnings,
          "Off-target screen skipped: no local transcriptome index for this species.")
      }
    }
    if (length(picked) < settings$pairs_per_gene)
      res$warnings <- c(res$warnings, paste0(
        "Only ", length(picked), " non-overlapping pairs could be placed (< ",
        settings$pairs_per_gene, " requested; ", length(cands),
        " valid windows available). All were used."))
    res$pairs <- lapply(seq_along(picked),
                        function(i) build_pair(job$gene, amp, i, picked[[i]]))
    res$design_seq_len <- nchar(seq_for_design)

    # Map the designed probes from design-sequence to genomic coordinates so
    # they can be drawn on the isoform track plot. Any mismatch between the
    # GTF exons and the transcript just skips the overlay.
    res$probe_frags <- tryCatch({
      segs <- NULL
      if (isTRUE(res$used_shared) && !is.null(res$shared_exons)) {
        segs <- res$shared_exons            # row order = design sequence order
      } else if (!is.null(by_tx)) {
        ti <- match(sub("\\..*$", "", rec$id), sub("\\..*$", "", names(by_tx)))
        if (!is.na(ti)) segs <- exons_to_tx_segments(by_tx[[ti]])
      }
      if (is.null(segs) || !length(res$pairs)) NULL
      else if (identical(res$struct_source, "GTF") &&
               abs(sum(segs$end - segs$start + 1) - res$design_seq_len) >
               0.20 * res$design_seq_len) NULL   # GTF/sequence version mismatch
      else map_pairs_to_genome(res$pairs, segs)
    }, error = function(e) NULL)
  }, error = function(e) {
    res$status <<- "error"
    res$error <<- conditionMessage(e)
  })
  res
}

# oPool packing: 90 bp per probe pair, keep each pool under the bp limit.
build_pools <- function(results, bp_limit) {
  ok <- Filter(function(r) r$status == "ok", results)
  pools <- list()
  cur <- list(genes = character(0), bp = 0)
  for (r in ok) {
    bp <- length(r$pairs) * 90
    if (bp > bp_limit)
      warning("Gene '", r$gene, "' requires ", bp, " bp of probes, exceeding the pool ",
              "limit of ", bp_limit, " bp; its pairs all land in the same pool.")
    if (length(cur$genes) && cur$bp + bp > bp_limit) {
      pools <- c(pools, list(cur)); cur <- list(genes = character(0), bp = 0)
    }
    cur$genes <- c(cur$genes, r$gene)
    cur$bp <- cur$bp + bp
  }
  if (length(cur$genes)) pools <- c(pools, list(cur))
  lapply(seq_along(pools), function(i)
    list(index = i, genes = pools[[i]]$genes, bp = pools[[i]]$bp))
}

# ---- Cross-species probe conservation check ----------------------------------
#
# The initiator and spacer parts of an HCR oligo are species-independent
# (they come from the amplifier, not the target). To reuse a probe set in
# another species, only the two 25-nt antisense halves need to match the
# orthologous transcript. These helpers score that.

# Best ungapped match of a short pattern against a (sense) transcript.
# Returns list(pos, mismatches, mm_pos, matched) or NULL if sequence too short.
best_ungapped_match <- function(pat, seq) {
  np <- nchar(pat); ns <- nchar(seq)
  if (is.na(np) || is.na(ns) || ns < np || np < 8) return(NULL)
  starts <- seq_len(ns - np + 1)
  subs <- substring(seq, starts, starts + np - 1)
  ok <- !grepl("N", subs, fixed = TRUE)
  if (!any(ok)) return(NULL)
  subs <- subs[ok]; starts <- starts[ok]
  m <- do.call(rbind, strsplit(subs, "", fixed = TRUE))
  p <- strsplit(pat, "", fixed = TRUE)[[1]]
  mm <- rowSums(m != matrix(rep(p, each = nrow(m)), nrow = nrow(m), ncol = np))
  i <- which.min(mm)[1]
  list(pos = starts[i], mismatches = mm[i],
       mm_pos = which(m[i, ] != p), matched = subs[i])
}

# Check every probe pair of a result against one ortholog transcript sequence.
# half_a / half_b are the 5' and 3' 25-mers of the 52-bp target (sense strand),
# i.e. exactly what the odd/even probe antisense halves must base-pair with.
check_pairs_in_seq <- function(pairs, ortho_seq) {
  lapply(pairs, function(p) {
    a <- best_ungapped_match(p$half_a, ortho_seq)
    b <- best_ungapped_match(p$half_b, ortho_seq)
    if (is.null(a) || is.null(b))
      return(list(pair = p$pair, found = FALSE, verdict = "no match",
                  a = a, b = b, gap = NA_integer_))
    gap <- b$pos - (a$pos + 25)          # expected: 2 (as in the design)
    verdict <- if (a$mismatches == 0 && b$mismatches == 0 && gap == 2) "exact"
               else if (a$mismatches <= 2 && b$mismatches <= 2) "near"
               else "risky"
    list(pair = p$pair, found = TRUE, verdict = verdict, a = a, b = b, gap = gap)
  })
}

# Match a gene-product name ("abdominal B") against a species index whose
# keys may use a longer/differently-worded annotation ("homeobox protein
# abdominal-b"). Splits both sides into alphanumeric tokens, drops generic
# annotation words, and keeps keys containing ALL of the query's informative
# tokens as whole words. Among matching keys the most specific one (fewest
# extra tokens) wins. Returns the key's records or NULL.
match_by_protein_tokens <- function(protein_name, idx) {
  if (is.null(idx) || is.null(protein_name) || !nzchar(protein_name)) return(NULL)
  stop_words <- c("protein", "homeobox", "transcription", "factor", "putative",
                  "uncharacterized", "like", "isoform", "family", "member",
                  "drosophila", "melanogaster", "mrna", "gene", "partial",
                  "homolog", "of", "and", "the", "a")
  toks_of <- function(s) {
    t <- strsplit(tolower(s), "[^a-z0-9]+")[[1]]
    t[nzchar(t)]
  }
  qt <- setdiff(toks_of(protein_name), stop_words)
  if (!length(qt)) qt <- toks_of(protein_name)
  if (!length(qt)) return(NULL)
  keys <- setdiff(ls(idx), ".byid")
  best <- NULL; best_extra <- Inf
  for (k in keys) {
    kt <- toks_of(k)
    # every query token must appear as a whole token on the candidate side
    if (!all(qt %in% kt)) next
    extra <- length(setdiff(kt, qt))
    if (extra < best_extra) { best <- k; best_extra <- extra }
  }
  if (is.null(best)) return(NULL)
  recs <- idx[[best]]
  if (length(recs)) { attr(recs, "how") <- "protein-name"; recs } else NULL
}

# Find the same gene's transcripts in another species: local transcriptome
# first (symbol, then word-boundary name search), NCBI fallback optional.
# The "how" attribute on the returned list records the match path:
# "symbol"/"accession" = exact index hit; "name-search" = loose word-boundary
# match on protein names (CAN grab the wrong gene, e.g. cad -> "CAD protein";
# flagged in the UI); "ncbi" = fetched online.
resolve_gene_records <- function(gene, species, cache, ncbi_fallback = FALSE,
                                 allow_name_fallback = TRUE) {
  idx <- load_species(species, cache)
  how <- NULL
  recs <- if (is.null(idx)) NULL else idx[[tolower(gene)]]
  if (!is.null(recs)) how <- "symbol"
  # direct accession: exact match in the local file's accession index
  if (is.null(recs) && !is.null(idx) && !is.null(idx[[".byid"]]) &&
      grepl("^[A-Z]{2,3}_[0-9]+(\\.[0-9]+)?$", trimws(gene))) {
    r <- idx[[".byid"]][[tolower(sub("\\..*$", "", trimws(gene)))]]
    if (!is.null(r)) { recs <- list(r); how <- "accession" }
  }
  if (is.null(recs) && !is.null(idx) && isTRUE(allow_name_fallback)) {
    pat <- paste0("\\b", gsub("([^A-Za-z0-9])", "\\\\\\1", tolower(gene)), "\\b")
    keys <- setdiff(ls(idx), ".byid")
    hits <- keys[vapply(keys, function(k) grepl(pat, k), TRUE)]
    if (length(hits)) {
      recs <- unique(unlist(lapply(hits, function(k) idx[[k]]), recursive = FALSE))
      how <- "name-search"
    }
  }
  if (is.null(recs) && isTRUE(ncbi_fallback)) {
    recs <- tryCatch(ncbi_fetch_gene(gene, organism_for(species)),
                     error = function(e) NULL)
    if (!is.null(recs)) how <- "ncbi"
  }
  if (!is.null(recs)) attr(recs, "how") <- how
  recs
}

# Align a gene's ortholog transcripts to the other species' genome (same
# chaining engine as the main design path) and map the probe-binding positions
# - the best matches of the two 25-nt halves - onto genomic coordinates, so the
# cross-species panel can draw the same gene-model track plot as the design
# species. Best-effort: returns NULL when no local genome exists for the
# species (the identity table still works without it).
ortholog_structure <- function(recs, rec, checked_pairs, species, cache) {
  genome_path <- find_genome_file(species)
  if (is.null(genome_path)) return(NULL)
  locus <- tryCatch(find_locus_by_sequence(genome_path, rec$seq, cache),
                    error = function(e) NULL)
  if (is.null(locus)) return(NULL)
  loc <- tryCatch(get_locus_seq(genome_path, locus, cache = cache),
                  error = function(e) NULL)
  if (is.null(loc)) return(NULL)
  mapped <- tryCatch(map_isoforms_seq(recs, loc), error = function(e) NULL)
  if (is.null(mapped) || !length(mapped$by_tx)) return(NULL)

  # Probe fragments: span from the half-a match start to the half-b match end
  # in the selected ortholog transcript, mapped through its exon segments.
  frags <- NULL
  ti <- match(sub("\\..*$", "", rec$id), sub("\\..*$", "", names(mapped$by_tx)))
  if (!is.na(ti)) {
    segs <- exons_to_tx_segments(mapped$by_tx[[ti]])
    pseudo <- list()
    for (cp in checked_pairs) {
      if (!isTRUE(cp$found)) next
      # keep the pair number as the name so anchors line up with the design
      # species' fragments even when some pairs are not found here
      pseudo[[as.character(cp$pair)]] <-
        list(start = min(cp$a$pos, cp$b$pos), end = max(cp$a$pos, cp$b$pos) + 24)
    }
    if (length(pseudo))
      frags <- tryCatch(map_pairs_to_genome(pseudo, segs), error = function(e) NULL)
  }
  list(by_tx = mapped$by_tx, frags = frags, coverage = mapped$coverage,
       src = "sequence alignment")
}

# Run the conservation check for one gene result across a species list.
cross_species_check <- function(res, species_list, cache, ncbi_fallback = FALSE) {
  out <- list()
  for (sp in species_list) {
    sp <- trimws(sp)
    if (!nzchar(sp) || tolower(sp) == tolower(res$species)) next
    entry <- list(species = sp, accession = NA_character_, tx_len = NA_integer_,
                  n_isoforms = NA_integer_, pairs = NULL, error = NULL,
                  isoform_exons = NULL, probe_frags = NULL,
                  struct_source = NULL, map_coverage = NULL)
    tryCatch({
      # Ortholog search terms: the optional CSV "ortholog" column may contain
      # plain terms (apply to every species) and/or species-tagged terms like
      # "simulans:LOC6727147; eugracilis:LOC108115147" - tagged terms are used
      # only for the matching species and take priority over the gene symbol.
      orth_terms <- character(0)
      if (!is.null(res$ortholog) && nzchar(res$ortholog)) {
        raw <- trimws(strsplit(res$ortholog, ";", fixed = TRUE)[[1]])
        for (rt in raw[nzchar(raw)]) {
          if (grepl("^[A-Za-z]+\\s*:", rt)) {
            tag <- tolower(trimws(sub(":.*$", "", rt)))
            val <- trimws(sub("^[A-Za-z]+\\s*:\\s*", "", rt))
            if (grepl(tag, tolower(sp), fixed = TRUE) && nzchar(val))
              orth_terms <- c(orth_terms, val)
          } else orth_terms <- c(orth_terms, rt)
        }
      }
      terms <- unique(c(orth_terms, res$gene))
      terms <- terms[!is.na(terms) & nzchar(terms)]
      recs <- NULL; matched_term <- NULL; match_how <- NULL
      for (tm in terms) {
        # When the user gave explicit ortholog terms, disable the loose
        # protein-name fallback: a wrong-gene match (e.g. cad -> "CAD protein")
        # is worse than an honest "not found".
        r <- resolve_gene_records(tm, sp, cache, ncbi_fallback,
                                  allow_name_fallback = !length(orth_terms))
        if (!is.null(r) && length(r)) {
          recs <- r; matched_term <- tm; match_how <- attr(r, "how"); break
        }
      }
      # Auto-resolution layer 1 (offline): the design transcript's protein
      # name ("abdominal B") usually reappears in the other species' RefSeq
      # annotation, possibly worded differently ("homeobox protein
      # abdominal-b") - exact key first, then all-tokens match.
      if (is.null(recs) && !is.null(res$protein_name) &&
          nzchar(res$protein_name)) {
        idx_sp <- load_species(sp, cache)
        r <- resolve_gene_records(res$protein_name, sp, cache, FALSE,
                                  allow_name_fallback = FALSE)
        if (is.null(r) && !is.null(idx_sp))
          r <- match_by_protein_tokens(res$protein_name, idx_sp)
        if (!is.null(r) && length(r)) {
          recs <- r; matched_term <- res$protein_name
          match_how <- "protein-name"
        }
      }
      # Auto-resolution layer 2 (online): ask NCBI's Gene database what
      # symbol/LOC ID this gene has in the check species, then use that
      # against the local file (or fetch the transcripts from NCBI).
      if (is.null(recs) && isTRUE(ncbi_fallback)) {
        gl <- tryCatch(ncbi_gene_lookup(res$gene, organism_for(sp)),
                       error = function(e) NULL)
        if (!is.null(gl)) {
          for (tm2 in unique(c(gl$symbol, gl$description))) {
            if (is.null(tm2) || !nzchar(tm2)) next
            r <- resolve_gene_records(tm2, sp, cache, FALSE,
                                      allow_name_fallback = FALSE)
            if (!is.null(r) && length(r)) {
              recs <- r; matched_term <- tm2; match_how <- "ncbi-gene"; break
            }
          }
          if (is.null(recs) && !is.null(gl$symbol) && nzchar(gl$symbol)) {
            r <- tryCatch(ncbi_fetch_gene(gl$symbol, organism_for(sp)),
                          error = function(e) NULL)
            if (!is.null(r) && length(r)) {
              recs <- r; matched_term <- gl$symbol; match_how <- "ncbi-gene"
            }
          }
        }
      }
      if (is.null(recs))
        stop(paste0('No transcript found for "', res$gene, '" in ', sp,
                    " (tried: ", paste(terms, collapse = ", "),
                    if (!is.null(res$protein_name))
                      paste0(", protein name '", res$protein_name, "'") else "",
                    if (isTRUE(ncbi_fallback)) ", NCBI Gene lookup" else "", "). ",
                    "Non-melanogaster RefSeq files index genes by protein ",
                    "name/LOC ID - put the protein name, LOC ID, or an XM_/NM_ ",
                    "accession for this species in the CSV's ortholog column",
                    if (!isTRUE(ncbi_fallback))
                      ", or enable 'NCBI fallback' to auto-look-up the LOC ID"
                    else "", "."))
      if (!identical(tolower(matched_term), tolower(res$gene)))
        entry$ortholog_term <- matched_term
      entry$match_how <- match_how
      entry$loose <- identical(match_how, "name-search")
      picked <- choose_transcript(recs, NULL)
      entry$accession <- picked$rec$id
      entry$tx_len <- nchar(picked$rec$seq)
      entry$n_isoforms <- length(recs)
      entry$pairs <- check_pairs_in_seq(res$pairs, picked$rec$seq)
      # Gene-model track plot data for this species (needs its genomic.fna)
      st <- ortholog_structure(recs, picked$rec, entry$pairs, sp, cache)
      if (!is.null(st)) {
        entry$isoform_exons <- st$by_tx
        entry$probe_frags <- st$frags
        entry$struct_source <- st$src
        entry$map_coverage <- st$coverage
      }
    }, error = function(e) entry$error <<- conditionMessage(e))
    out[[length(out) + 1]] <- entry
  }
  if (length(out)) out else NULL
}

# Verdict color (shared by the UI table and the CSV export logic).
VERDICT_COLORS <- c(exact = "#16a34a", near = "#d97706", risky = "#dc2626",
                    "no match" = "#64748b")

# Flatten the cross-species check into one row per pair x species (for CSV).
cross_species_rows <- function(results) {
  rows <- list()
  for (r in results) {
    if (is.null(r$cross_species)) next
    for (e in r$cross_species) {
      if (!is.null(e$error)) {
        rows[[length(rows) + 1]] <- data.frame(
          Gene = r$gene, Check_species = e$species, Accession = "",
          Match_method = "", Pair = NA,
          HalfA_identity = "", HalfA_mismatches = NA, HalfA_mm_pos = "",
          HalfB_identity = "", HalfB_mismatches = NA, HalfB_mm_pos = "",
          Gap_bp = NA, Verdict = paste0("ERROR: ", e$error),
          stringsAsFactors = FALSE)
        next
      }
      for (cp in e$pairs) {
        fmt <- function(x) if (is.null(x)) list(id = "", mm = NA, pos = "") else
          list(id = paste0(25 - x$mismatches, "/25"), mm = x$mismatches,
               pos = paste(x$mm_pos, collapse = ";"))
        fa <- fmt(cp$a); fb <- fmt(cp$b)
        rows[[length(rows) + 1]] <- data.frame(
          Gene = r$gene, Check_species = e$species, Accession = e$accession,
          Match_method = if (is.null(e$match_how)) "" else e$match_how,
          Pair = cp$pair,
          HalfA_identity = fa$id, HalfA_mismatches = fa$mm, HalfA_mm_pos = fa$pos,
          HalfB_identity = fb$id, HalfB_mismatches = fb$mm, HalfB_mm_pos = fb$pos,
          Gap_bp = if (is.null(cp$gap)) NA else cp$gap, Verdict = cp$verdict,
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ---- CSV parsing ------------------------------------------------------------

# Minimal RFC-4180-ish parser: handles quoted cells with commas/newlines and
# escaped "" quotes. (read.csv misaligns empty fields before multiline quoted
# cells, so we don't use it.)
parse_csv <- function(text) {
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  n <- length(chars)
  rows <- list(); row <- character(0); cell <- ""
  in_q <- FALSE; i <- 1
  while (i <= n) {
    ch <- chars[i]
    if (in_q) {
      if (ch == "\"") {
        if (i < n && chars[i + 1] == "\"") { cell <- paste0(cell, "\""); i <- i + 1 }
        else in_q <- FALSE
      } else cell <- paste0(cell, ch)
    } else if (ch == "\"") {
      in_q <- TRUE
    } else if (ch == ",") {
      row <- c(row, cell); cell <- ""
    } else if (ch == "\n" || ch == "\r") {
      if (ch == "\r" && i < n && chars[i + 1] == "\n") i <- i + 1
      row <- c(row, cell); cell <- ""
      if (any(nzchar(trimws(row)))) rows <- c(rows, list(row))
      row <- character(0)
    } else cell <- paste0(cell, ch)
    i <- i + 1
  }
  row <- c(row, cell)
  if (any(nzchar(trimws(row)))) rows <- c(rows, list(row))
  rows
}

parse_jobs <- function(text) {
  rows <- parse_csv(text)
  if (!length(rows)) stop("CSV is empty.")
  header <- tolower(trimws(rows[[1]]))
  if (!("gene" %in% header)) stop('CSV must have a "gene" column header.')
  known <- c("gene", "species", "amplifier", "transcript", "sequence", "ortholog")
  extra <- setdiff(header, known)
  if (length(extra))
    warning("Ignoring CSV column(s) not used by the designer: ",
            paste(sprintf('"%s"', extra), collapse = ", "),
            ". Recognized columns: ", paste(known, collapse = ", "), ".")
  jobs <- list()
  for (r in rows[-1]) {
    getv <- function(nm) {
      j <- match(nm, header)
      if (is.na(j) || j > length(r)) "" else trimws(r[j])
    }
    gene <- getv("gene")
    if (!nzchar(gene)) next
    amp <- toupper(getv("amplifier"))
    sp  <- getv("species")
    tr  <- getv("transcript")
    sq  <- getv("sequence")
    orth <- getv("ortholog")
    jobs <- c(jobs, list(list(
      gene = gene,
      species = if (nzchar(sp)) canonical_species(sp) else "melanogaster",
      amplifier = if (amp %in% AMP_LIST) amp else NULL,
      transcript = if (nzchar(tr)) tr else NULL,
      sequence = if (nzchar(sq)) sq else NULL,
      ortholog = if (nzchar(orth)) orth else NULL)))
  }
  jobs
}

# ---- Cross-species results panel ---------------------------------------------

# Monospace sequence with mismatched bases highlighted in red.
highlight_mm <- function(seq25, mm_pos, base_col = "#0f172a") {
  chars <- strsplit(seq25, "", fixed = TRUE)[[1]]
  spans <- lapply(seq_along(chars), function(i)
    span(style = if (i %in% mm_pos)
           paste0("color:#dc2626; font-weight:bold; text-decoration:underline;")
         else paste0("color:", base_col, ";"), chars[i]))
  tagList(spans)
}

# Compact pair x species matrix plus per-pair mismatch detail.
cross_species_ui <- function(cs) {
  if (is.null(cs) || !length(cs)) return(NULL)

  # Union of pair indices across species, in order
  pair_ids <- sort(unique(unlist(lapply(cs, function(e)
    vapply(e$pairs, function(x) x$pair, 0)))))
  # Cell lookup: [[species]][[pair]] -> check result
  cell_of <- function(e, pid) {
    if (!is.null(e$error) || is.null(e$pairs)) return(NULL)
    hit <- Filter(function(x) x$pair == pid, e$pairs)
    if (length(hit)) hit[[1]] else NULL
  }
  cell_tag <- function(e, pid) {
    if (!is.null(e$error))
      return(tags$td(style = "color:#64748b; font-size:10px;", "not found"))
    cp <- cell_of(e, pid)
    if (is.null(cp)) return(tags$td("-"))
    if (!isTRUE(cp$found))
      return(tags$td(style = "color:#64748b;", "no match"))
    col <- VERDICT_COLORS[[cp$verdict]]
    tags$td(style = paste0("color:", col, "; font-weight:bold; font-size:11px;"),
            paste0(25 - cp$a$mismatches, "/25 + ", 25 - cp$b$mismatches, "/25"),
            tags$br(),
            span(style = "font-weight:normal; font-size:10px;",
                 paste0(cp$verdict, if (!is.na(cp$gap) && cp$gap != 2)
                   paste0(", gap ", cp$gap, " bp!") else "")))
  }

  detail_blocks <- lapply(cs, function(e) {
    if (!is.null(e$error))
      return(p(style = "color:#dc2626; font-size:12px;",
               paste0(e$species, ": ", e$error)))
    interesting <- Filter(function(x) !isTRUE(x$found) || x$verdict != "exact",
                          e$pairs)
    if (!length(interesting))
      return(p(style = "color:#16a34a; font-size:12px;",
               paste0(e$species, " (", e$accession, "): all pairs exact matches. ",
                      "This probe set should work as-is.")))
    tagList(
      if (isTRUE(e$loose))
        p(style = "color:#d97706; font-size:12px; font-weight:bold;",
          "⚠ Ortholog matched by loose protein-name search only - this may be ",
          "the WRONG gene (e.g. 'cad' can match 'CAD protein'). Verify the ",
          "accession before trusting these verdicts; for a definitive lookup, ",
          "put a LOC ID or accession in the CSV ortholog column."),
      p(style = "color:#475569; font-size:12px; margin-bottom:2px;",
        strong(e$species), " (", code(e$accession),
        if (!is.null(e$ortholog_term))
          paste0(", found via '", e$ortholog_term, "'") else "", "): ",
        length(e$pairs) - length(interesting), " of ", length(e$pairs),
        " pairs are exact; mismatched pairs below. ",
        "Ortholog sequence shown, mismatched bases in red:"),
      lapply(interesting, function(cp) {
        if (!isTRUE(cp$found))
          return(div(style = "font-family:'Courier New',monospace; font-size:12px; color:#64748b;",
                     paste0("pair ", cp$pair, ": no match found in this transcript")))
        div(style = "font-family:'Courier New',monospace; font-size:12px; background:#f8fafc; padding:6px 8px; border-radius:4px; margin-bottom:4px;",
          div(span(style = "color:#64748b; font-size:10px;",
                   paste0("pair ", cp$pair, " · half-a (odd probe binds) @", cp$a$pos,
                          " · ", 25 - cp$a$mismatches, "/25:  ")),
              highlight_mm(cp$a$matched, cp$a$mm_pos)),
          div(span(style = "color:#64748b; font-size:10px;",
                   paste0("pair ", cp$pair, " · half-b (even probe binds) @", cp$b$pos,
                          " · ", 25 - cp$b$mismatches, "/25:  ")),
              highlight_mm(cp$b$matched, cp$b$mm_pos)),
          if (!is.na(cp$gap) && cp$gap != 2)
            div(style = "color:#dc2626; font-size:10px;",
                paste0("gap between halves is ", cp$gap,
                       " bp in this species (2 bp in the design) - probes will still bind, but spacing differs"))
        )
      })
    )
  })

  tagList(
    tags$table(class = "table table-condensed",
      style = "margin-bottom:8px; width:auto;",
      tags$thead(tags$tr(
        tags$th("pair"),
        lapply(cs, function(e) tags$th(e$species)))),
      tags$tbody(lapply(pair_ids, function(pid)
        tags$tr(tags$td(paste0("#", pid)),
                lapply(cs, function(e) cell_tag(e, pid)))))),
    detail_blocks
  )
}

# ---- Shiny UI ---------------------------------------------------------------

SAMPLE_CSV <- paste(
  "gene,species,amplifier,transcript,sequence",
  "Abd-B,melanogaster,,,",
  "cad,melanogaster,,,",
  sep = "\n")

ui <- fluidPage(
  titlePanel("HCR Probe Designer"),
  p("High-throughput HCR v3.0 probe-set design (B1\u2013B5) from a gene-list CSV. ",
    "Sequences come from the local transcript files (", code("reference_files/"),
    "), or from NCBI if the fallback is enabled."),
  sidebarLayout(
    sidebarPanel(
      h4("1. Gene list (CSV)"),
      fileInput("csv_file", NULL, accept = c(".csv", "text/csv")),
      textAreaInput("csv_text", NULL, rows = 7,
                    placeholder = "gene,species,amplifier,transcript,sequence\nAbd-B,melanogaster,,,"),
      actionButton("load_sample", "Load sample", class = "btn-default btn-sm"),
      downloadButton("dl_template", "Download template", class = "btn-default btn-sm"),
      helpText("Columns: gene (required - symbol for melanogaster, symbol or full ",
               "protein name for other species), species (any detected: ",
               paste(SPECIES_CHOICES, collapse = ", "),
               "), amplifier (B1-B17, optional), transcript (isoform ",
               "hint like RA, or an accession; optional), sequence (optional ",
               "DNA/FASTA override), ortholog (optional - protein name, LOC ID, or ",
               "accession to use when looking up this gene in OTHER species for the ",
               "cross-species check)."),
      hr(),
      h4("2. Settings"),
      radioButtons("stringency", "Filtering stringency",
                   choices = c("Strict (recommended)" = "strict",
                               "Loose (mirrors auto-relax: GC +/-5, homopolymer +1, dG +/-15, Tm +/-10, compositions off)" = "loose"),
                   selected = "strict", inline = TRUE),
      radioButtons("amp_set", "Amplifier set",
                   choices = c("HCR v2.0 \u2014 B1\u2013B5 only" = "v2",
                               "HCR v3.0 \u2014 all B1\u2013B17"  = "v3"),
                   selected = "v2", inline = TRUE),
      fluidRow(
        column(6, numericInput("pairs_per_gene", "Pairs per gene", 6, min = 1, max = 40)),
        column(6, numericInput("spacing", "Spacing between targets (bp)", 5, min = 0))
      ),
      fluidRow(
        column(6, numericInput("gc_min", "GC min (%)", 40)),
        column(6, numericInput("gc_max", "GC max (%)", 60))
      ),
      fluidRow(
        column(6, numericInput("max_hp", "Max homopolymer run", 4, min = 3)),
        column(6, numericInput("pool_bp", "Pool bp limit", 3300))
      ),
      fluidRow(
        column(6, numericInput("tm_min", "Tm min (\u00b0C)", 60)),
        column(6, numericInput("tm_max", "Tm max (\u00b0C)", 90))
      ),
      fluidRow(
        column(6, numericInput("dg_min", "dG min (kcal/mol)", -100)),
        column(6, numericInput("dg_max", "dG max (kcal/mol)", -60))
      ),
      fluidRow(
        column(6, numericInput("target_dg", "Target dG (pair)", -80)),
        column(6, numericInput("target_dg_half", "Target dG (half)", -38))
      ),
      fluidRow(
        column(6, numericInput("temp_c", "Temperature (\u00b0C)", 37, min = 0, max = 100)),
        column(6, numericInput("na_m", "Na+ (M)", 0.3, min = 0.01, max = 5, step = 0.1))
      ),
      checkboxInput("use_thermo",
                    "Score + rank windows by dG/Tm (nearest-neighbor)", TRUE),
      checkboxInput("filter_thermo",
                    "Reject windows outside the dG/Tm ranges above", FALSE),
      checkboxInput("comp_a", "Composition: A < 28% in each half", TRUE),
      checkboxInput("comp_c", "Composition: C 22-28% in each half (strict)", FALSE),
      checkboxInput("comp_stack", "Composition: no AAAA/CCCC runs", TRUE),
      checkboxInput("comp_cstack", "Composition: C-stacking limit (first 12 nt)", TRUE),
      checkboxInput("offtarget_screen",
                    "Screen probes for off-targets in the transcriptome (slower)", TRUE),
      fluidRow(
        column(6, numericInput("ot_seed", "Off-target seed length (nt)", 18, min = 12, max = 25)),
        column(6, numericInput("max_offtarget", "Max off-target genes", 0, min = 0))
      ),
      checkboxInput("ncbi_fallback",
                    "Fetch from NCBI if not found locally (needs internet)", FALSE),
      helpText(style = "font-size:11px; color:#64748b; margin-top:-6px;",
               "Orthologs are auto-resolved by protein name from the local ",
               "files (offline); with NCBI fallback on, the NCBI Gene ",
               "database is also queried to find the LOC ID when the local ",
               "annotation doesn't mention the gene."),
      checkboxInput("relax_filters",
                    "Auto-relax GC/homopolymer filters when few windows found", FALSE),
      checkboxInput("use_shared_exons",
                    "Design against shared regions only (CDS when annotated; requires GTF + genome)", TRUE),
      checkboxInput("avoid_splice_junctions",
                    "Avoid splice junctions (52-bp window must lie within one exon)", TRUE),
      hr(),
      h4("Cross-species check (optional)"),
      checkboxGroupInput("check_species", "Check probes in other species",
                         choices = species_checkbox_labels(),
                         selected = character(0)),
      textOutput("ref_status", inline = TRUE),
      # One-click download of a species' RefSeq files into reference_files/
      div(style = "background:#f1f5f9; border-radius:6px; padding:8px 10px; margin-bottom:8px;",
        tags$b(style = "font-size:12px;", "Add a species"),
        selectInput("dl_species", NULL,
                    choices = setdiff(names(SPECIES_GCF), detect_species()),
                    width = "100%"),
        checkboxInput("dl_genome_too",
                      "Also download genome + annotation GFF (large, ~50-90 MB; needed for structure plots)",
                      TRUE),
        actionButton("dl_species_btn", "Download from NCBI",
                     class = "btn-default btn-sm"),
        textOutput("dl_species_status", inline = FALSE)
      ),
      helpText("OFF by default: tick any species to run the cross-species check. ",
               "Each probe's two 25-nt binding halves are then matched against ",
               "the same gene's transcript in the ticked species, so you can see ",
               "which pairs are reusable across species. Species that are not on ",
               "disk yet are fetched from NCBI automatically on the first run that ",
               "uses them (unless the NCBI fallback above is on). Leave everything ",
               "unticked to design for the CSV's species only, with no ",
               "cross-species check."),
      # First install: only melanogaster is required, and it is fetched on demand.
      div(style = "background:#fff7ed; border:1px solid #fdba74; border-radius:6px; padding:8px 10px; margin-bottom:8px;",
        tags$b(style = "font-size:12px;", "D. melanogaster reference"),
        p(style = "font-size:11px; margin:4px 0;", "Downloaded from NCBI on first ",
          "use into ", code("d_melanogaster/"), " inside the reference folder. ",
          "The ~2 GB bundled genome set is no longer shipped with the app."),
        actionButton("dl_mel_btn", "Download D. melanogaster (~73 MB) now",
                     class = "btn-primary btn-sm", width = "100%"),
        textOutput("dl_mel_status", inline = TRUE)
      ),
      actionButton("run", "3. Design probes", class = "btn-primary btn-lg",
                   width = "100%")
    ),
    mainPanel(
      uiOutput("summary"),
      uiOutput("pools_ui"),
      uiOutput("results_ui")
    )
  )
)

# ---- Shiny server -----------------------------------------------------------

# ---- FASTA export helper ----------------------------------------------------

# Build a FASTA file for a single gene result:
#   1. Full transcript sequence (header = gene|accession|full)
#   2. Shared CDS regions concatenated with NNNNNN spacers (header = gene|shared_CDS)
#   3. Each 52-bp target window (header = gene|pair_N|start-end)
build_gene_fasta <- function(r) {
  if (r$status != "ok" || !length(r$pairs)) return(NULL)

  lines <- character(0)

  # 1. Full transcript
  lines <- c(lines, paste0(">", r$gene, "|", r$accession, "|full_transcript"))
  lines <- c(lines, r$seq)

  # 2. Shared CDS regions (if available)
  if (!is.null(r$shared_exons) && nrow(r$shared_exons) > 0) {
    genome_path <- find_genome_file(r$gene)
    if (!is.null(genome_path)) {
      # Build concatenated sequence with NNNNNN spacers between regions
      region_seqs <- character(nrow(r$shared_exons))
      for (i in seq_len(nrow(r$shared_exons))) {
        region_seqs[i] <- extract_regions(
          genome_path,
          r$shared_exons[i, , drop = FALSE]
        )
      }
      # Filter out empty sequences
      region_seqs <- region_seqs[nzchar(region_seqs)]
      if (length(region_seqs)) {
        concat_seq <- paste(region_seqs, collapse = "NNNNNN")
        lines <- c(lines, paste0(">", r$gene, "|shared_", 
                    if (isTRUE(r$shared_is_cds)) "CDS" else "exons"))
        lines <- c(lines, concat_seq)
      }
    }
  }

  # 3. Each 52-bp target window
  for (p in r$pairs) {
    lines <- c(lines, paste0(">", r$gene, "|pair_", p$pair, "|", p$start, "-", p$end))
    lines <- c(lines, p$target52)
  }

  paste(lines, collapse = "\n")
}

# Build a combined FASTA for all genes
build_combined_fasta <- function(results) {
  ok <- Filter(function(r) r$status == "ok" && length(r$pairs) > 0, results)
  if (!length(ok)) return("")
  fastas <- vapply(ok, build_gene_fasta, "")
  fastas <- fastas[nzchar(fastas)]
  if (!length(fastas)) return("")
  paste(fastas, collapse = "\n\n")
}

# ---- All-amplifier probe table helper ----------------------------------------

# Build a row for the Excel-style all-amplifier table.
# For each 52-bp target, output probe sequences for ALL amplifiers (B1-B5).
# Columns: Probe name, 52bp_Reverse_Comp, 52bp_Gene_Orientation,
#          B1-1, B1-2, B2-1, B2-2, ... for every amplifier in `amps`.
build_all_amp_row <- function(gene, amp, pair_index, target52, accession = "", start = 0, end = 0, gc1 = 0, gc2 = 0,
                              dg = NA_real_, tm = NA_real_, dg1 = NA_real_, dg2 = NA_real_,
                              ot_n = NA_integer_, ot_genes = "", amps = AMP_V2) {
  rc52 <- revcomp(target52)
  pad <- sprintf("%02d", pair_index)
  probe_name <- paste0(gene, " #", pair_index)

  # --- Antisense halves (no initiator) for Benchling primer visualization ---
  odd_antisense  <- substr(rc52, nchar(rc52) - 24, nchar(rc52))   # 25 nt
  even_antisense <- substr(rc52, 1, 25)                            # 25 nt

  # Build probe sequences for every amplifier in the selected set
  probes <- character(2 * length(amps))
  names(probes) <- as.vector(rbind(paste0(amps, "-1"), paste0(amps, "-2")))

  for (a_name in amps) {
    a <- AMPLIFIERS[[a_name]]
    # Odd probe: initiator_half1 + spacer1 + last 25 nt of RC(52bp)
    odd <- paste0(a$half1, a$spacer1, odd_antisense)
    # Even probe: first 25 nt of RC(52bp) + spacer2 + initiator_half2
    even <- paste0(even_antisense, a$spacer2, a$half2)
    probes[[paste0(a_name, "-1")]] <- odd
    probes[[paste0(a_name, "-2")]] <- even
  }

  c(list(
    Gene = gene,
    Accession = accession,
    Pair = pair_index,
    Start = start,
    End = end,
    Target_52bp = target52,
    GC_half1 = round(gc1, 1),
    GC_half2 = round(gc2, 1),
    dG_full = round(dg, 1),
    Tm_full = round(tm, 1),
    dG_half1 = round(dg1, 1),
    dG_half2 = round(dg2, 1),
    Offtarget_n = ot_n,
    Offtarget_genes = ot_genes,
    Probe_name = probe_name,
    Target_52bp_RC = rc52,
    Target_52bp_sense = target52,
    odd_antisense = odd_antisense,
    even_antisense = even_antisense
  ), as.list(probes))
}

# Build the full all-amplifier table from results (one row per target, with a
# probe-sequences column pair for every amplifier in `amps`).
build_all_amp_table <- function(results, amps = AMP_V2) {
  ok <- Filter(function(r) r$status == "ok" && length(r$pairs) > 0, results)
  if (!length(ok)) return(NULL)

  rows <- do.call(rbind, lapply(ok, function(r) {
    do.call(rbind, lapply(r$pairs, function(p) {
      as.data.frame(build_all_amp_row(r$gene, r$amplifier, p$pair, p$target52,
                    r$accession, p$start, p$end, p$gc1, p$gc2,
                    p$dg, p$tm, p$dg1, p$dg2, p$ot_n, p$ot_genes,
                    amps = amps),
                    stringsAsFactors = FALSE)
    }))
  }))
  rows
}

server <- function(input, output, session) {
  results <- reactiveVal(list())
  pools <- reactiveVal(list())
  last_settings <- reactiveVal(NULL)      # for the run-params download
  refs_used <- reactiveVal(list())        # per-species reference provenance (log)
  cross_requested <- reactiveVal(character(0))  # cross-species species list (log)
  jobs_run <- reactiveVal(NULL)           # parsed gene list of the last run (log)
  species_cache <- new.env(hash = TRUE)   # session-level transcriptome cache

  # ---- on-demand reference management (melanogaster on first use, the rest on
  # demand; nothing genome-sized is shipped with the app anymore) ----
  ref_state <- reactiveVal(0L)      # bumped whenever a download completes
  output$ref_status <- renderText({
    ref_state()
    local <- installed_species()
    remote <- setdiff(names(SPECIES_GCF), local)
    if (!length(local))
      return(paste0("On demand: ", paste(remote, collapse = ", ")))
    paste0("Local: ", paste(local, collapse = ", "),
           if (length(remote))
             paste0("  \u00b7  On demand: ", paste(remote, collapse = ", ")))
  })
  # Refresh checkbox labels + species picker whenever install state changes.
  observe({
    ref_state()
    updateCheckboxGroupInput(session, "check_species",
                             choices = species_checkbox_labels(),
                             selected = isolate(input$check_species))
    updateSelectInput(session, "dl_species",
                      choices = setdiff(names(SPECIES_GCF), detect_species()))
  })
  # Shared downloader used by the sidebar button, the first-install modal and
  # the run-guard. Idempotent: only missing components are fetched.
  dl_install_species <- function(sp, with_genome) {
    msg <- tryCatch(
      withProgress(message = paste("Downloading", sp, "from NCBI"), value = 0.1, {
        download_species_ref(sp, with_genome,
                             progress = function(m) setProgress(0.5, detail = m))
      }),
      error = function(e) paste("Download failed:", conditionMessage(e)))
    if (grepl("^Download failed", msg))
      showNotification(msg, type = "error", duration = 15)
    else
      showNotification(msg, type = "message", duration = 8)
    ref_state(ref_state() + 1L)
    invisible(msg)
  }
  output$dl_mel_status <- renderText({
    ref_state()                                  # refresh after a download
    if (species_reference_installed("melanogaster")) "melanogaster installed"
    else "not installed yet \u2013 needed for melanogaster designs"
  })
  observeEvent(input$dl_mel_btn, dl_install_species("melanogaster", TRUE))
  # First launch: offer the melanogaster download right away.
  if (!species_reference_installed("melanogaster")) {
    showModal(modalDialog(
      title = "Install the D. melanogaster genome?",
      p("This build no longer ships the ~2 GB genome bundle \u2013 each species' ",
        "reference files are downloaded from NCBI on first use."),
      p("The D. melanogaster reference (~73 MB, compressed) is needed before ",
        "you can design melanogaster probes or run cross-species checks."),
      footer = tagList(
        actionButton("mel_dl_go", "Download D. melanogaster (~73 MB)",
                     class = "btn-primary"),
        modalButton("Skip for now")
      ),
      size = "m", easyClose = FALSE)
    )
  }
  observeEvent(input$mel_dl_go, {
    removeModal()
    dl_install_species("melanogaster", TRUE)
  })

  observeEvent(input$load_sample,
               updateTextAreaInput(session, "csv_text", value = SAMPLE_CSV))
  observeEvent(input$csv_file, {
    req(input$csv_file)
    updateTextAreaInput(session, "csv_text",
                        value = paste(readLines(input$csv_file$datapath, warn = FALSE),
                                      collapse = "\n"))
  })

  output$dl_template <- downloadHandler(
    filename = function() "HCR_template.csv",
    content = function(file) writeLines(SAMPLE_CSV, file)
  )

  # ---- one-click species download ----
  dl_species_status <- reactiveVal("")
  output$dl_species_status <- renderText(dl_species_status())
  observeEvent(input$dl_species_btn, {
    sp <- input$dl_species
    req(sp)
    msg <- dl_install_species(sp, isTRUE(input$dl_genome_too))
    dl_species_status(msg)
  })

  observeEvent(input$run, {
    jobs <- tryCatch(parse_jobs(input$csv_text), error = function(e) e)
    if (inherits(jobs, "error")) {
      showNotification(conditionMessage(jobs), type = "error"); return()
    }
    if (!length(jobs)) {
      showNotification("No genes found. Paste a CSV or load the sample.", type = "error")
      return()
    }
    jobs_run(jobs)                          # kept for HCR_design_log.txt
    # Settings: single source of truth is default_settings() (batch & env vars).
    # UI inputs override the defaults at runtime so the GUI is always live.
    settings <- default_settings()
    settings$pairs_per_gene <- input$pairs_per_gene
    settings$spacing        <- input$spacing
    settings$gc_min         <- input$gc_min
    settings$gc_max         <- input$gc_max
    settings$max_hp         <- input$max_hp
    # Amplifiers are NOT user-facing: the engine rotates B1-B5 per gene for the
    # oligo naming, and the primary deliverable (HCR PROBE MAKER all-amp CSV)
    # carries every combination in the selected amplifier set so the choice can
    # be made at the bench.
    settings$default_amp    <- "auto"
    # Which amplifiers are available. Default v2 (B1-B5); v3 adds B7/B9/B10/
    # B13-B15/B17. NULL-guarded for older saved UI state.
    settings$amp_set        <- if (is.null(input$amp_set)) "v2" else input$amp_set
    settings$ncbi_fallback  <- input$ncbi_fallback
    settings$relax_filters  <- input$relax_filters
    settings$use_shared_exons <- input$use_shared_exons
    settings$avoid_splice_junctions <- input$avoid_splice_junctions
    settings$use_thermo     <- input$use_thermo
    settings$filter_thermo  <- input$filter_thermo
    settings$temp_c         <- input$temp_c
    settings$na_m           <- input$na_m
    settings$oligo_conc     <- 5e-5             # fixed reagent concentration
    settings$tm_min         <- input$tm_min
    settings$tm_max         <- input$tm_max
    settings$dg_min         <- input$dg_min
    settings$dg_max         <- input$dg_max
    settings$target_dg      <- input$target_dg
    settings$target_dg_half <- input$target_dg_half
    settings$comp_a         <- input$comp_a
    settings$comp_c         <- input$comp_c
    settings$comp_stack     <- input$comp_stack
    settings$comp_cstack    <- input$comp_cstack
    settings$offtarget_screen <- input$offtarget_screen
    settings$ot_seed        <- input$ot_seed
    settings$max_offtarget  <- input$max_offtarget
    # Stringency preset: "loose" mirrors the auto-relax rescues and is applied
    # to every gene for the whole run. Defaults to strict when the input is
    # missing (e.g. older saved UI state).
    settings$stringency <- if (is.null(input$stringency)) "strict"
                           else input$stringency
    settings <- apply_stringency(settings)
    last_settings(settings)
    check_sps <- input$check_species
    if (is.null(check_sps)) check_sps <- character(0)
    # Make sure every species this run touches is on disk. melanogaster is
    # downloaded on first use; ticked cross-species species are fetched on
    # demand. Skipped when the NCBI fallback is on (that path fetches the
    # transcripts one by one instead of pulling whole compressed genomes).
    needed <- unique(c(check_sps,
                       vapply(jobs, function(j) canonical_species(j$species), "")))
    needed <- intersect(needed, names(SPECIES_GCF))
    cross_requested(check_sps)                    # recorded for the design log
    pre_installed <- needed[vapply(needed, species_reference_installed, TRUE)]
    missing <- needed[!vapply(needed, species_reference_installed, TRUE)]
    if (length(missing) && !isTRUE(settings$ncbi_fallback)) {
      for (sp in missing)
        if (grepl("^Download failed", dl_install_species(sp, TRUE))) return()
    }
    # Reference provenance for HCR_design_log.txt: local vs freshly downloaded,
    # per species, with the NCBI assembly accession (SPECIES_GCF).
    refs_used(setNames(lapply(needed, function(sp) {
      if (sp %in% pre_installed) "local (already installed)"
      else if (species_reference_installed(sp)) "freshly downloaded from NCBI"
      else "not on disk (NCBI fallback enabled: transcripts fetched per gene)"
    }), needed))
    out <- vector("list", length(jobs))
    withProgress(message = "Designing probes", value = 0, {
      for (i in seq_along(jobs)) {
        setProgress(value = (i - 1) / length(jobs), detail = jobs[[i]]$gene)
        out[[i]] <- design_gene(jobs[[i]], settings, i - 1, species_cache)
        out[[i]]$ortholog <- jobs[[i]]$ortholog  # needed by cross_species_check
        if (length(check_sps) && out[[i]]$status == "ok" && length(out[[i]]$pairs)) {
          setProgress(detail = paste0(jobs[[i]]$gene, ": cross-species check"))
          out[[i]]$cross_species <- tryCatch(
            cross_species_check(out[[i]], check_sps, species_cache,
                                settings$ncbi_fallback),
            error = function(e) NULL)
        }
      }
    })
    results(out)
    pools(build_pools(out, input$pool_bp))
  })

  total_oligos <- reactive(sum(vapply(results(), function(r) length(r$pairs) * 2, 0)))

  output$summary <- renderUI({
    req(length(results()) > 0)
    wellPanel(
      style = "display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap;",
      span(strong(total_oligos()), " oligos \u00b7 ",
           strong(length(pools())), " oPool(s)"),
      span(
        downloadButton("dl_zip", "Download all (ZIP)", class = "btn-primary btn-sm"),
        downloadButton("dl_all_amp", "HCR PROBE MAKER all-amp CSV (selected amplifier set)", class = "btn-primary btn-sm"),
        downloadButton("dl_order", "Oligo order CSV", class = "btn-default btn-sm"),
        downloadButton("dl_antisense_fasta", "Antisense halves CSV file (For Benchling primer upload)", class = "btn-default btn-sm"),
        downloadButton("dl_pools", "Pool summary", class = "btn-default btn-sm"),
        downloadButton("dl_params", "Run parameters (params.txt)", class = "btn-default btn-sm"),
        downloadButton("dl_cross_species", "Cross-species conservation CSV", class = "btn-default btn-sm"),
        downloadButton("dl_shared_cds_fasta", "Shared CDS sequences (FASTA)", class = "btn-default btn-sm"),
        downloadButton("dl_fasta", "FASTA (52 bp target sequences shared across all isoforms)", class = "btn-default btn-sm")
      ),
      helpText(style = "font-size:11px; color:#64748b; width:100%; margin-top:6px;",
               "The HCR PROBE MAKER CSV lists every amplifier in the selected set (B1\u2013B5 for ",
               "HCR v2.0 by default, or all B1\u2013B17 for v3.0) as odd/even oligos for each 52-bp ",
               "target, so the amplifier can be chosen at the bench; the oligo order CSV keeps the ",
               "auto-assigned set. 'Download all (ZIP)' bundles every output file of this run (plus ",
               "HCR_design_log.txt) for lab record keeping.")
    )
  })

  output$pools_ui <- renderUI({
    req(length(pools()) > 0)
    wellPanel(
      h4("oPool packing"),
      tags$table(class = "table table-condensed",
        tags$tbody(lapply(pools(), function(p)
          tags$tr(
            tags$td(strong(paste("Pool", p$index)), ": ", paste(p$genes, collapse = ", ")),
            tags$td(style = "text-align:right;", paste0(p$bp, " bp"))
          ))))
    )
  })

  output$results_ui <- renderUI({
    res <- results()
    if (!length(res))
      return(wellPanel(style = "border-style:dashed; text-align:center; color:#94a3b8; padding:40px;",
                       "Results will appear here. Load the sample CSV and hit \u201cDesign probes\u201d to try it out."))

    # --- Color legend for probe construction ---
    legend <- wellPanel(
      style = "background:#f8fafc; border:1px solid #e2e8f0; padding:10px 14px; margin-bottom:16px;",
      div(style = "font-size:12px; color:#475569; font-weight:bold; margin-bottom:6px;", 
          "Probe construction color key:"),
      div(style = "display:flex; flex-wrap:wrap; gap:12px; font-size:11px;",
        div(span(style = "color:#0f172a; font-weight:bold;", "██"), " 52-bp target"),
        div(span(style = "color:#f59e0b; font-weight:bold;", "██"), " 2-nt gap"),
        div(span(style = "color:#8b5cf6; font-weight:bold;", "██"), " Antisense (RC)"),
        div(span(style = "color:#0ea5e9; font-weight:bold;", "██"), " Initiator half"),
        div(span(style = "color:#94a3b8; font-weight:bold;", "██"), " Spacer"),
        div(span(style = "color:#ef4444; font-weight:bold;", "██"), " Odd probe"),
        div(span(style = "color:#10b981; font-weight:bold;", "██"), " Even probe")
      )
    )
    tagList(
      legend,
      lapply(seq_along(res), function(i) {
      r <- res[[i]]
      plot_id <- paste0("map_", i)
      table_id <- paste0("tbl_", i)
      iso_plot_id <- paste0("isomap_", i)
      local({
        rr <- r
        output[[plot_id]] <- renderPlot({
          if (rr$status != "ok" || !length(rr$pairs)) return(NULL)
          par(mar = c(2, 1, 1, 1))
          plot(0, type = "n", xlim = c(1, max(rr$seq_len, 1)), ylim = c(0, 1),
               axes = FALSE, xlab = "", ylab = "")
          axis(1)
          segments(1, 0.5, rr$seq_len, 0.5, col = "#94a3b8", lwd = 5, lend = 2)
          col <- amp_color(rr$amplifier)
          for (p in rr$pairs)
            rect(p$start, 0.25, p$end, 0.75, col = col, border = NA)
        }, height = 70)
        output[[table_id]] <- renderUI({
          if (rr$status != "ok" || !length(rr$pairs)) return(NULL)

          # Color scheme for visual display
          col_target <- "#0f172a"    # dark slate - target sequence
          col_gap    <- "#f59e0b"    # amber - gap
          col_rc     <- "#8b5cf6"    # purple - reverse complement (antisense)
          col_init   <- "#0ea5e9"    # sky blue - initiator half
          col_spacer <- "#94a3b8"    # slate gray - spacer
          col_odd    <- "#ef4444"    # red - odd probe
          col_even   <- "#10b981"    # green - even probe

          tagList(
            lapply(rr$pairs, function(p) {
              tagList(
                hr(style = "border-top: 1px solid #e2e8f0; margin: 12px 0;"),

                # --- STEP 1: The 52-bp target ---
                h5(style = paste0("color:", col_target, "; margin-top:8px;"),
                   paste0("Pair ", p$pair, " · 52-bp target (positions ", p$start, "–", p$end, ")")),
                div(style = "color:#64748b; font-size:11px; margin-bottom:2px;",
                    paste0(
                      if (!is.null(p$dg) && !is.na(p$dg))
                        paste0("\u0394G ", round(p$dg, 1), " kcal/mol · Tm ", round(p$tm, 1),
                               " \u00b0C · halves \u0394G ", round(p$dg1, 1), " / ",
                               round(p$dg2, 1), " · ") else "",
                      if (!is.null(p$ot_n) && !is.na(p$ot_n))
                        paste0("off-targets: ", if (p$ot_n == 0) "none"
                               else paste0(p$ot_n, " gene(s): ", p$ot_genes))
                      else "off-target screen not run")),
                div(style = "font-family: 'Courier New', monospace; font-size: 13px; line-height: 1.6; background:#f8fafc; padding:10px; border-radius:6px;",
                  # Target sense strand with gap highlighted
                  div(style = "margin-bottom:4px;",
                    span(style = "color:#64748b; font-size:11px;", "Target (sense 5'→3'): "),
                    span(style = paste0("color:", col_target, "; font-weight:bold;"), p$half_a),
                    span(style = paste0("color:", col_gap, "; font-weight:bold;"), p$gap),
                    span(style = paste0("color:", col_target, "; font-weight:bold;"), p$half_b)
                  ),
                  # Labels under the target
                  div(style = "margin-bottom:8px; padding-left:0px;",
                    span(style = "color:#64748b; font-size:10px;", "                25-nt half-a   2-nt gap   25-nt half-b")
                  ),

                  # --- STEP 2: Reverse complement each half ---
                  div(style = "margin:8px 0; border-left:3px solid #cbd5e1; padding-left:8px;",
                    div(style = "color:#64748b; font-size:11px; margin-bottom:4px;", "Step 1: Reverse-complement each half (antisense strands for probe binding):"),
                    div(
                      span(style = "color:#64748b; font-size:11px;", "  odd probe binds → "),
                      span(style = paste0("color:", col_rc, "; font-weight:bold;"), p$rc_half_a)
                    ),
                    div(
                      span(style = "color:#64748b; font-size:11px;", "  even probe binds → "),
                      span(style = paste0("color:", col_rc, "; font-weight:bold;"), p$rc_half_b)
                    )
                  ),

                  # --- STEP 3: Add amplifier initiator sequences ---
                  div(style = "margin:8px 0; border-left:3px solid #cbd5e1; padding-left:8px;",
                    div(style = "color:#64748b; font-size:11px; margin-bottom:4px;", 
                        paste0("Step 2: Add ", p$amp, " amplifier initiator halves + spacers:")),

                    # Odd probe construction
                    div(style = "margin-bottom:6px;",
                      span(style = "color:#64748b; font-size:11px;", "  Odd probe (5'→3'): "),
                      span(style = paste0("color:", col_init, "; font-weight:bold;"), p$odd_initiator),
                      span(style = paste0("color:", col_spacer, ";"), p$odd_spacer),
                      span(style = paste0("color:", col_rc, "; font-weight:bold;"), p$rc_half_a),
                      br(),
                      span(style = "color:#64748b; font-size:10px; padding-left:18px;", 
                           paste0("  └─ ", nchar(p$odd_initiator), "-nt initiator  +  ", nchar(p$odd_spacer), "-nt spacer  +  25-nt antisense"))
                    ),

                    # Even probe construction
                    div(style = "margin-bottom:6px;",
                      span(style = "color:#64748b; font-size:11px;", "  Even probe (5'→3'): "),
                      span(style = paste0("color:", col_rc, "; font-weight:bold;"), p$rc_half_b),
                      span(style = paste0("color:", col_spacer, ";"), p$even_spacer),
                      span(style = paste0("color:", col_init, "; font-weight:bold;"), p$even_initiator),
                      br(),
                      span(style = "color:#64748b; font-size:10px; padding-left:18px;", 
                           paste0("  └─ 25-nt antisense  +  ", nchar(p$even_spacer), "-nt spacer  +  ", nchar(p$even_initiator), "-nt initiator"))
                    )
                  ),

                  # --- STEP 4: Final sequences ---
                  div(style = "margin-top:8px; padding:8px; background:#f1f5f9; border-radius:4px;",
                    div(style = "color:#64748b; font-size:11px; margin-bottom:6px; font-weight:bold;", "Final oligo sequences (45 nt each, 90 bp per pair):"),
                    div(style = "margin-bottom:4px;",
                      span(style = paste0("color:", col_odd, "; font-weight:bold;"), p$odd_name), 
                      span(style = "color:#64748b; font-size:11px;", " ("), 
                      span(style = paste0("color:", col_odd, ";"), paste0(nchar(p$odd), " nt")),
                      span(style = "color:#64748b; font-size:11px;", "): "),
                      span(style = paste0("color:", col_odd, "; font-family:monospace;"), p$odd)
                    ),
                    div(
                      span(style = paste0("color:", col_even, "; font-weight:bold;"), p$even_name), 
                      span(style = "color:#64748b; font-size:11px;", " ("), 
                      span(style = paste0("color:", col_even, ";"), paste0(nchar(p$even), " nt")),
                      span(style = "color:#64748b; font-size:11px;", "): "),
                      span(style = paste0("color:", col_even, "; font-family:monospace;"), p$even)
                    )
                  )
                )
              )
            })
          )
        })
        output[[iso_plot_id]] <- renderPlot({
          tryCatch({
          if (is.null(rr$isoform_exons) || !length(rr$isoform_exons)) return(NULL)
          # label rows with accession + variant letter + alignment coverage
          ids <- names(rr$isoform_exons)
          labs <- ids
          if (!is.null(rr$all_isoforms)) {
            m <- match(sub("\\..*$", "", ids), sub("\\..*$", "", rr$all_isoforms$ID))
            vr <- rr$all_isoforms$Variant[m]
            labs <- ifelse(!is.na(vr) & nzchar(vr), paste0(labs, " (", vr, ")"), labs)
          }
          if (!is.null(rr$map_coverage)) {
            pct <- round(100 * as.numeric(rr$map_coverage[ids]))
            labs <- ifelse(!is.na(pct), paste0(labs, " ", pct, "%"), labs)
          }
          labs <- stats::setNames(labs, ids)
          plot_isoform_tracks(rr$isoform_exons, selected_id = rr$accession,
                              shared = rr$shared_exons, labels = labs,
                              sel_col = amp_color(rr$amplifier),
                              probe_frags = rr$probe_frags,
                              cds_by_tx = rr$cds_by_tx,
                              shared_label = if (isTRUE(rr$shared_is_cds)) "shared CDS"
                                             else "shared")
          }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot skipped:", conditionMessage(e)),
                 col = "#dc2626", cex = 0.8)
          })
        })
        # Combined cross-species gene-model plot: all species in one figure,
        # with the designed probes' binding sites overlaid on each ortholog.
        if (!is.null(rr$cross_species)) {
          output[[paste0("csplot_all_", i)]] <- renderPlot({
            tryCatch(
              plot_combined_species_tracks(rr),
              error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Plot skipped:", conditionMessage(e)),
                     col = "#dc2626", cex = 0.8)
              })
          })
        }
      })
      badge <- if (r$status == "ok")
        span(class = "label label-primary", paste(length(r$pairs), "pairs"))
      else span(class = "label label-danger", "failed")
      shared_badge <- if (isTRUE(r$used_shared))
        span(class = "label label-success",
             if (isTRUE(r$shared_is_cds)) "shared CDS" else "shared exons")
      else NULL
      junction_badge <- if (isTRUE(r$avoid_splice_junctions))
        span(class = "label label-info", "intron-free")
      else NULL
      wellPanel(
        h4(r$gene, " ", badge, " ", shared_badge, " ", junction_badge, " ",
           if (!is.na(r$accession)) code(r$accession),
           if (r$status == "ok")
             tags$small(style = "color:#64748b;",
                        paste0(" \u00b7 ", format(r$seq_len, big.mark = ","), " nt \u00b7 ",
                               r$n_candidates, " valid windows"))),
                if (!is.null(r$error)) p(style = "color:#dc2626;", r$error),
        lapply(r$warnings, function(w) p(style = "color:#d97706;", paste0("⚠ ", w))),
        if (!is.null(r$isoform_exons) && length(r$isoform_exons)) tagList(
          tags$details(open = "open",
            tags$summary(style = "color:#0284c7; cursor:pointer;",
                         paste0("Isoform structures (", length(r$isoform_exons),
                                " transcripts, genomic view",
                                if (!is.null(r$struct_source))
                                  paste0(" · ", r$struct_source) else "", ")")),
            plotOutput(iso_plot_id,
                       height = paste0(70 + 26 * (length(r$isoform_exons) +
                                         if (!is.null(r$shared_exons)) 1 else 0), "px")))
        ),
        if (!is.null(r$all_isoforms) && nrow(r$all_isoforms) > 1) tagList(
          tags$details(
            tags$summary(style = "color:#64748b; cursor:pointer;",
                         paste0("Show ", nrow(r$all_isoforms), " isoforms")),
            tags$table(class = "table table-condensed table-striped",
              tags$thead(tags$tr(
                tags$th("ID"), tags$th("Symbol"), tags$th("Name"),
                tags$th("Variant"), tags$th("Length"), tags$th("Selected"), tags$th("Link")
              )),
              tags$tbody(lapply(seq_len(nrow(r$all_isoforms)), function(j) {
                iso <- r$all_isoforms[j, ]
                url <- jbrowse_url(list(id = iso$ID, symbol = iso$Symbol))
                tags$tr(
                  tags$td(code(iso$ID)),
                  tags$td(iso$Symbol),
                  tags$td(iso$Name),
                  tags$td(iso$Variant),
                  tags$td(iso$Length),
                  tags$td(if(iso$Selected) "✓" else ""),
                  tags$td(if (!is.null(url)) a(href = url, target = "_blank", "View") else "")
                )
              }))
            )
          )
        ),
        if (!is.null(r$shared_exons) && nrow(r$shared_exons) > 0) tagList(
          tags$details(
            tags$summary(style = "color:#16a34a; cursor:pointer;",
                         paste0("Show ", nrow(r$shared_exons),
                                if (isTRUE(r$shared_is_cds)) " shared CDS regions"
                                else " shared exon regions")),
            tags$table(class = "table table-condensed table-striped",
              tags$thead(tags$tr(
                tags$th("Chr"), tags$th("Start"), tags$th("End"), tags$th("Length"), tags$th("Strand")
              )),
              tags$tbody(lapply(seq_len(nrow(r$shared_exons)), function(j) {
                ex <- r$shared_exons[j, ]
                tags$tr(
                  tags$td(ex$chr),
                  tags$td(format(ex$start, big.mark = ",")),
                  tags$td(format(ex$end, big.mark = ",")),
                  tags$td(ex$end - ex$start + 1),
                  tags$td(ex$strand)
                )
              }))
            )
          )
        ),
        if (!is.null(r$cross_species) && length(r$cross_species)) tagList(
          tags$details(open = "open",
            tags$summary(style = "color:#7c3aed; cursor:pointer; font-weight:bold;",
                         paste0("Cross-species conservation (",
                                paste(vapply(r$cross_species, function(e) e$species, ""),
                                      collapse = ", "), ")")),
            cross_species_ui(r$cross_species),
            {
              n_rows <- if (!is.null(r$isoform_exons)) length(r$isoform_exons) else 0
              if (!is.null(r$shared_exons) && nrow(r$shared_exons))
                n_rows <- n_rows + 1   # shared-region track
              n_grps <- as.integer(n_rows > 0)
              for (e in r$cross_species)
                if (!is.null(e$isoform_exons) && length(e$isoform_exons)) {
                  n_rows <- n_rows + length(e$isoform_exons); n_grps <- n_grps + 1
                }
              if (n_rows > 0) tagList(
                p(style = "color:#475569; font-size:12px; margin:10px 0 2px;",
                  strong("Transcript structures across species"),
                  " (relative coordinates; saturated track = transcript used for the probe match):"),
                plotOutput(paste0("csplot_all_", i),
                           height = paste0(90 + 26 * n_rows + 14 * n_grps, "px"))
              )
            })
        ),
        if (r$status == "ok" && length(r$pairs)) tagList(
          plotOutput(plot_id, height = "80px"),
          tags$details(tags$summary(style = "color:#64748b; cursor:pointer;",
                                    paste("Show", length(r$pairs) * 2, "oligo sequences — with design steps")),
                       uiOutput(table_id))
        )
      )
    })
  )
  })

  # ---- downloads ----
  # Provenance line for the design log: where the gene list came from.
  csv_source_label <- reactive({
    if (!is.null(input$csv_file) && length(input$csv_file$name) &&
        nzchar(input$csv_file$name %||% ""))
      paste0("uploaded CSV: ", input$csv_file$name)
    else
      paste0("pasted CSV text")
  })

  order_rows <- reactive(build_order_table(results(), pools()))

  output$dl_order <- downloadHandler(
    filename = function() "HCR_oligo_order.csv",
    content = function(file) utils::write.csv(order_rows(), file, row.names = FALSE,
                                              quote = FALSE)
  )

  # ---- Download all (ZIP): every output file of this run, for lab records ----
  output$dl_zip <- downloadHandler(
    filename = function() paste0("HCR_probes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      tmp <- file.path(tempdir(),
        paste0("hcr_zip_", format(Sys.time(), "%H%M%S"), "_",
               paste(sample(letters, 6, replace = FALSE), collapse = "")))
      dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
      write_all_run_files(tmp, results(), last_settings(), pools(),
                          jobs_run(), cross_requested(), refs_used(),
                          length(cross_requested()) > 0, run_params_txt = FALSE,
                          source_label = csv_source_label())
      msg <- zip_run_files(file, tmp)
      unlink(tmp, recursive = TRUE)
      if (!grepl("^ZIP created", msg)) {
        # No zip utility on this machine: leave a small explanation in place of
        # the archive instead of failing the download outright.
        writeLines(msg, file)
        showNotification(msg, type = "error", duration = 15)
      }
    }
  )

  output$dl_pools <- downloadHandler(
    filename = function() "HCR_pool_summary.csv",
    content = function(file)
      utils::write.csv(build_pool_table(pools(), total_oligos() * 45), file,
                       row.names = FALSE)
  )

  output$dl_params <- downloadHandler(
    filename = function() "HCR_run_params.txt",
    content = function(file) {
      s <- last_settings()
      if (is.null(s)) {
        writeLines("No run parameters yet.", file)
        return()
      }
      writeLines(build_params_text(s, paste0("HCR Probe Designer run: ", Sys.time()),
                                   paste0("csv_source = ", csv_source_label())),
                 file)
    }
  )

  # ---- CSV download: cross-species conservation per pair ----
  output$dl_cross_species <- downloadHandler(
    filename = function() "HCR_cross_species_conservation.csv",
    content = function(file) {
      rows <- build_cross_species_table(results())
      if (is.null(rows))
        rows <- data.frame(Message = paste("No cross-species checks run.",
          "Enter species in 'Check probes in other species' and re-run."))
      utils::write.csv(rows, file, row.names = FALSE)
    }
  )

  # ---- FASTA download: shared CDS sequences per gene ----
  output$dl_shared_cds_fasta <- downloadHandler(
    filename = function() "HCR_shared_CDS_sequences.fasta",
    content = function(file) writeLines(build_shared_cds_fasta_text(results()), file)
  )

  # ---- FASTA download: full transcript + shared CDS + 52-bp targets per gene ----
  output$dl_fasta <- downloadHandler(
    filename = function() "HCR_probe_targets.fasta",
    content = function(file) {
      fa <- build_combined_fasta(results())
      if (!nzchar(fa)) {
        fa <- ">no_data\nNo valid probe designs found."
      }
      writeLines(fa, file)
    }
  )

  # ---- PRIMARY deliverable: HCR PROBE MAKER all-amplifier CSV -----------------
  # One row per 52-bp target with probe sequences for every amplifier in the
  # selected amplifier set (B1-B5 by default; all B1-B17 when v3.0 chosen).
  output$dl_all_amp <- downloadHandler(
    filename = function() "HCR_probe_maker_all_amplifiers.csv",
    content = function(file) {
      tbl <- build_all_amp_table(results(), active_amps(last_settings()))
      if (is.null(tbl)) {
        tbl <- data.frame(Message = "No valid probe designs found.")
      }
      utils::write.csv(tbl, file, row.names = FALSE)
    }
  )

  # ---- CSV download: antisense halves for Benchling "Import DNA/RNA oligos" ----
  output$dl_antisense_fasta <- downloadHandler(
    filename = function() "HCR_antisense_halves_for_Benchling.csv",
    content = function(file) {
      rows <- build_antisense_table(results())
      if (is.null(rows)) {
        writeLines("Name,Sequence\nno_data,No valid probe designs found.", file)
        return()
      }
      utils::write.csv(rows, file, row.names = FALSE, quote = FALSE)
    }
  )
}

# ============================================================================
# Shared run-output builders (used by BOTH the Shiny download handlers and
# run_batch(), so the interactive ZIP download and the batch folder always
# contain identical files). No chemistry here - pure re-formatting of the
# design results.
# ============================================================================

# --- Oligo order CSV rows (specific amplifier assignments) -------------------
build_order_table <- function(results, pools) {
  pool_of <- new.env(hash = TRUE)
  for (p in pools) for (g in p$genes)
    pool_of[[strsplit(g, " ")[[1]][1]]] <- p$index
  do.call(rbind, unlist(lapply(results, function(r)
    lapply(r$pairs, function(p) data.frame(
      Pool = { v <- pool_of[[r$gene]]; if (is.null(v)) "" else v },
      Name = c(p$odd_name, p$even_name),
      Sequence = c(p$odd, p$even), stringsAsFactors = FALSE))),
    recursive = FALSE))
}

# --- Pool summary rows (with the TOTAL row) ----------------------------------
build_pool_table <- function(pools, total_oligos_bp) {
  rows <- if (length(pools)) do.call(rbind, lapply(pools, function(p) data.frame(
    Pool = p$index, Genes = paste(p$genes, collapse = "; "),
    Total_bp = p$bp, stringsAsFactors = FALSE))) else NULL
  rbind(rows, data.frame(Pool = "TOTAL", Genes = "",
                         Total_bp = total_oligos_bp))
}

# --- Antisense halves CSV rows (for Benchling primer upload) ------------------
build_antisense_table <- function(results) {
  ok <- Filter(function(r) r$status == "ok" && length(r$pairs) > 0, results)
  if (!length(ok)) return(NULL)
  do.call(rbind, lapply(ok, function(r) {
    do.call(rbind, lapply(r$pairs, function(p) {
      rc52 <- revcomp(p$target52)
      odd_as  <- substr(rc52, nchar(rc52) - 24, nchar(rc52))
      even_as <- substr(rc52, 1, 25)
      data.frame(
        Name = c(paste0(r$gene, "_pair", p$pair, "_odd_antisense"),
                 paste0(r$gene, "_pair", p$pair, "_even_antisense")),
        Sequence = c(odd_as, even_as),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# --- Shared-CDS/shared-exon FASTA text ----------------------------------------
build_shared_cds_fasta_text <- function(results) {
  ok <- Filter(function(r) r$status == "ok" && !is.null(r$shared_exons) && nrow(r$shared_exons) > 0, results)
  if (!length(ok)) return(paste0(">no_data\nNo shared CDS data available."))
  lines <- character(0)
  for (r in ok) {
    region_seqs <- character(nrow(r$shared_exons))
    genome_path <- find_genome_file(r$species)
    if (is.null(genome_path)) next
    for (i in seq_len(nrow(r$shared_exons))) {
      region_seqs[i] <- extract_regions(
        genome_path,
        r$shared_exons[i, , drop = FALSE]
      )
    }
    region_seqs <- region_seqs[nzchar(region_seqs)]
    if (!length(region_seqs)) next
    concat_seq <- paste(region_seqs, collapse = "NNNNNN")
    label <- if (isTRUE(r$shared_is_cds)) "shared_CDS" else "shared_exons"
    lines <- c(lines, paste0(">", r$gene, "|", label, "|", nrow(r$shared_exons), "regions|", nchar(concat_seq), "nt"))
    lines <- c(lines, concat_seq)
  }
  if (!length(lines)) lines <- ">no_data\nNo shared CDS sequences extracted."
  paste(lines, collapse = "\n")
}

# --- Cross-species conservation table (NULL when the check was not run) -------
build_cross_species_table <- function(results) {
  rows <- cross_species_rows(results)
  if (!is.null(rows) && nrow(rows)) rows else NULL
}

# --- Plot pages (one per gene: isoform track + cross-species track) -----------
# Shared by run_batch() and the ZIP download. Draw to the active pdf() device.
draw_plots_pdf_pages <- function(results) {
  for (r in results) {
    if (r$status != "ok") next
    if (!is.null(r$isoform_exons) && length(r$isoform_exons)) {
      ids <- names(r$isoform_exons)
      labs <- stats::setNames(ids, ids)
      tryCatch(
        plot_isoform_tracks(r$isoform_exons, selected_id = r$accession,
                            shared = r$shared_exons, labels = labs,
                            sel_col = amp_color(r$amplifier),
                            probe_frags = r$probe_frags, cds_by_tx = r$cds_by_tx,
                            shared_label = if (isTRUE(r$shared_is_cds)) "shared CDS"
                                          else "shared"),
        error = function(e) NULL)
      title(main = paste0(r$gene, " - isoforms (", r$accession, ")"),
            cex.main = 0.9, line = 0.5)
    }
    if (!is.null(r$cross_species))
      tryCatch({ plot_combined_species_tracks(r)
                 title(main = paste0(r$gene, " - cross-species comparison"),
                       cex.main = 0.9, line = 0.5) },
               error = function(e) NULL)
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Run parameters text (shared by HCR_run_params.txt and the design log) ----
build_params_text <- function(settings, header = NULL, extra = character(0)) {
  c(header,
    extra,
    paste0(names(settings), " = ",
           vapply(settings, function(x) paste(x, collapse = ","), "")))
}

# --- Species reference provenance for the log ---------------------------------
# Records, per species used by the run, whether the reference set was already
# on disk (local) or freshly downloaded from NCBI during the run, plus the
# NCBI assembly accession. refs_used is a named list of source strings.
build_reference_lines <- function(refs_used, jobs, check_sps) {
  used <- unique(c(vapply(jobs, function(j) canonical_species(j$species), ""),
                   vapply(check_sps, canonical_species, "")))
  used <- intersect(used, names(SPECIES_GCF))
  if (!length(used)) return(character(0))
  header <- c("", "SPECIES REFERENCES (RefSeq, NCBI)", "---------------------------------")
  c(header, unlist(lapply(used, function(sp) {
    src <- refs_used[[sp]]
    if (is.null(src))
      src <- if (species_reference_installed(sp)) "local (already installed)"
             else "not installed (no local reference)"
    gcf <- SPECIES_GCF[[sp]]
    sprintf("  %-14s %-40s %s", sp, src, if (is.null(gcf)) "unknown assembly" else gcf)
  })))
}

# ============================================================================
# HCR_design_log.txt - the definitive, self-contained record of one run.
# Anyone holding this file (plus the app version it names) can reproduce the
# exact run: every setting, the stringency preset and its effective values,
# which species' references were used and how they got there, per-gene design
# summaries with all warnings, and the cross-species verdicts.
# ============================================================================
build_design_log <- function(results, settings, jobs, check_sps, refs_used,
                             cross_requested,
                             source_label = "interactive (Shiny app)") {
  raw <- default_settings()          # pre-preset defaults, for the delta table
  lines <- c(
    paste0(APP_TITLE, " design log"),
    "=============================",
    paste0("Run timestamp      : ", format(Sys.time())),
    paste0("App version        : ", APP_VERSION),
    paste0("R version          : ", R.version.string),
    paste0("Platform           : ", R.version$platform),
    paste0("Run source         : ", source_label),
    paste0("Stringency mode    : ", settings$stringency),
    if (identical(settings$stringency, "loose"))
      "  (loose preset applied: GC +/-5, homopolymer +1, dG +/-15, Tm +/-10, composition filters off)"
    else
      "  (strict: thresholds as entered, no preset adjustments)",
    "",
    "RUN PARAMETERS (effective values used for this run)",
    "---------------------------------------------------",
    build_params_text(settings),
    "",
    "STRINGENCY PRESET - strict defaults vs effective values",
    "-------------------------------------------------------",
    sprintf("  %-14s %10s -> %10s", "parameter", "strict", "effective"),
    sprintf("  %-14s %10s -> %10s", "gc_min", raw$gc_min, settings$gc_min),
    sprintf("  %-14s %10s -> %10s", "gc_max", raw$gc_max, settings$gc_max),
    sprintf("  %-14s %10s -> %10s", "max_hp", raw$max_hp, settings$max_hp),
    sprintf("  %-14s %10s -> %10s", "tm_min", raw$tm_min, settings$tm_min),
    sprintf("  %-14s %10s -> %10s", "tm_max", raw$tm_max, settings$tm_max),
    sprintf("  %-14s %10s -> %10s", "dg_min", raw$dg_min, settings$dg_min),
    sprintf("  %-14s %10s -> %10s", "dg_max", raw$dg_max, settings$dg_max),
    paste0("  composition filters (comp_a/comp_c/comp_stack/comp_cstack): ",
           if (identical(settings$stringency, "loose")) "OFF (loose)" else "ON (strict)"),
    build_reference_lines(refs_used, jobs, check_sps),
    "",
    "PER-GENE DESIGN SUMMARY",
    "-----------------------")
  for (r in results) {
    amp_txt <- if (is.null(r$amplifier)) "-" else r$amplifier
    lines <- c(lines,
      paste0("  Gene ", r$gene, " (", r$species, "): status = ", r$status),
      paste0("    amplifier assigned : ", amp_txt,
             if (identical(settings$default_amp, "auto"))
               paste0(" (auto-rotation over the active amplifier set: ",
                      paste(active_amps(settings), collapse = "/"), ", by gene order)")
             else ""),
      paste0("    pairs designed     : ", length(r$pairs)),
      paste0("    candidates found   : ", r$n_candidates),
      paste0("    accession          : ", if (is.na(r$accession)) "-" else r$accession),
      paste0("    design seq length  : ", r$design_seq_len, " nt"),
      paste0("    regions used       : ",
             if (isTRUE(r$used_shared) && !is.null(r$shared_exons) && nrow(r$shared_exons))
               paste0(nrow(r$shared_exons), " ",
                      if (isTRUE(r$shared_is_cds)) "shared CDS" else "shared exon",
                      " region(s)")
             else "full transcript (no shared regions)"))
    if (!is.null(r$error))
      lines <- c(lines, paste0("    error              : ", r$error))
    if (length(r$warnings)) {
      lines <- c(lines, paste0("    warnings (", length(r$warnings), "):"),
                 paste0("      ! ", r$warnings))
    } else lines <- c(lines, "    warnings           : none")
    lines <- c(lines, "")
  }
  lines <- c(lines, "CROSS-SPECIES CHECK", "-------------------")
  if (!cross_requested || !length(check_sps)) {
    lines <- c(lines, "  OFF for this run (no species ticked / --no-cross).")
  } else {
    lines <- c(lines, paste0("  Species checked: ", paste(check_sps, collapse = ", ")))
    for (r in results) {
      if (is.null(r$cross_species)) next
      for (e in r$cross_species) {
        lines <- c(lines, paste0("  ", r$gene, " -> ", e$species, ":"))
        if (!is.null(e$error)) {
          lines <- c(lines, paste0("    ERROR: ", e$error))
          next
        }
        id <- function(x) if (is.null(x)) "-" else paste0(25 - x$mismatches, "/25")
        lines <- c(lines,
          paste0("    ortholog matched   : ", e$accession,
                 if (!is.null(e$ortholog_term)) paste0(" (term: ", e$ortholog_term, ")") else ""),
          paste0("    match method       : ",
                 if (is.null(e$match_how)) "-" else e$match_how),
          paste0("    isoforms in target : ", e$n_isoforms),
          paste0("    transcript length  : ", e$tx_len, " nt"))
        if (length(e$pairs)) {
          lines <- c(lines, "    per-pair verdicts:")
          for (cp in e$pairs)
            lines <- c(lines, sprintf(
              "      pair %s: halfA %s, halfB %s, gap %s bp -> %s",
              cp$pair, id(cp$a), id(cp$b),
              if (is.null(cp$gap)) "-" else cp$gap, cp$verdict))
        }
        lines <- c(lines, "")
      }
    }
  }
  lines
}



# --- Zip helper ---------------------------------------------------------------
# Locates a usable zip binary without installing anything (base R only).
find_zip_binary <- function() {
  hit <- Sys.which("zip")
  if (nzchar(hit)) return(hit)
  candidates <- c("/usr/bin/zip",                                   # macOS
                  file.path(R.home("bin"), "zip.exe"),              # bundled R for Windows
                  "C:/Rtools/bin/zip.exe",                          # Rtools
                  "C:/Program Files/Git/usr/bin/zip.exe",           # Git for Windows
                  file.path(Sys.getenv("RTOOLS40_HOME", ""), "usr", "bin", "zip.exe"))
  hit <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(hit)) hit[1] else NULL
}

# Zip the contents of dir (flat archive, no leading folder) into zipfile.
# Returns a status message; never errors out the caller.
zip_run_files <- function(zipfile, dir) {
  zb <- find_zip_binary()
  if (is.null(zb))
    return(paste0("Could not create the ZIP: no 'zip' binary found on this system. ",
                  "The individual files are still available via the other buttons."))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  files <- list.files(dir, all.files = FALSE, no.. = TRUE)
  res <- utils::zip(zipfile = zipfile, files = files, zip = zb, flags = "-j")
  if (is.null(res) || !file.exists(zipfile))
    paste("ZIP creation failed:", res)
  else
    paste0("ZIP created: ", basename(zipfile), " (", length(files), " files)")
}

# --- Write the complete, canonical run-output file set into out_dir -----------
# Single source of truth for the Download-all (ZIP) button and run_batch().
write_all_run_files <- function(out_dir, results, settings, pools, jobs,
                                check_sps = character(0), refs_used = list(),
                                cross_requested = FALSE, run_params_txt = FALSE,
                                source_label = "interactive (Shiny app)",
                                csv_path = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. PRIMARY deliverable: HCR PROBE MAKER all-amp CSV (selected amp set).
  amp_tbl <- build_all_amp_table(results, active_amps(settings))
  if (!is.null(amp_tbl))
    utils::write.csv(amp_tbl,
                     file.path(out_dir, "HCR_probe_maker_all_amplifiers.csv"),
                     row.names = FALSE)
  # 2. Oligo order CSV (specific amplifier assignments, if needed).
  order <- build_order_table(results, pools)
  if (!is.null(order))
    utils::write.csv(order, file.path(out_dir, "HCR_oligo_order.csv"),
                     row.names = FALSE, quote = FALSE)
  # 3. Antisense halves CSV (Benchling primer upload).
  as_rows <- build_antisense_table(results)
  if (!is.null(as_rows))
    utils::write.csv(as_rows, file.path(out_dir, "HCR_antisense_halves.csv"),
                     row.names = FALSE, quote = FALSE)
  # 4. Pool summary (with TOTAL row).
  utils::write.csv(
    build_pool_table(pools,
                     sum(vapply(results, function(r) length(r$pairs) * 2, 0)) * 45),
    file.path(out_dir, "HCR_pool_summary.csv"), row.names = FALSE)
  # 5. Cross-species conservation (only when the check was run).
  if (cross_requested) {
    cs <- build_cross_species_table(results)
    if (!is.null(cs))
      utils::write.csv(cs,
                       file.path(out_dir, "HCR_cross_species_conservation.csv"),
                       row.names = FALSE)
  }
  # 6. The definitive run log.
  writeLines(build_design_log(results, settings, jobs, check_sps, refs_used,
                              cross_requested, source_label),
             file.path(out_dir, "HCR_design_log.txt"))
  # 7. Shared-CDS FASTA.
  writeLines(build_shared_cds_fasta_text(results),
             file.path(out_dir, "HCR_shared_CDS_sequences.fasta"))
  # 8. Probe-target FASTA.
  fa <- build_combined_fasta(results)
  if (!nzchar(fa)) fa <- ">no_data\nNo valid probe designs found."
  writeLines(fa, file.path(out_dir, "HCR_probe_targets.fasta"))
  # 9. Plots PDF (one page per gene: isoform track + cross-species track).
  pdf(file.path(out_dir, "HCR_plots.pdf"), width = 11, height = 8.5)
  tryCatch(draw_plots_pdf_pages(results), error = function(e) NULL)
  dev.off()
  # Batch mode also keeps the classic HCR_run_params.txt.
  if (run_params_txt)
    writeLines(build_params_text(settings,
                                 paste0("HCR probe design batch run: ", Sys.time()),
                                 c(paste0("mode = batch"),
                                   if (!is.null(csv_path)) paste0("input_csv = ", csv_path))),
               file.path(out_dir, "HCR_run_params.txt"))
  invisible(file.path(out_dir, "HCR_probe_maker_all_amplifiers.csv"))
}



# Default settings, matching the Shiny UI defaults. Used by batch mode;
# edit values here to change batch behavior.
# Batch-mode env-var overrides (documented, for automation):
#   HCR_RELAX_FILTERS=1/true  -> force relax_filters=TRUE (e.g. to reproduce
#                                pre-v41 batch default or for AT-rich inputs)
default_settings <- function() {
  s <- list(
    pairs_per_gene = 6, spacing = 5, gc_min = 40, gc_max = 60, max_hp = 4,
    default_amp = "auto", pool_bp = 3300, ncbi_fallback = FALSE,
    amp_set = "v2",
    relax_filters = FALSE, use_shared_exons = TRUE, avoid_splice_junctions = TRUE,
    use_thermo = TRUE, filter_thermo = FALSE,
    temp_c = 37, na_m = 0.3, oligo_conc = 5e-5,
    tm_min = 60, tm_max = 90, dg_min = -100, dg_max = -60,
    target_dg = -80, target_dg_half = -38,
    comp_a = TRUE, comp_c = FALSE, comp_stack = TRUE, comp_cstack = TRUE,
    offtarget_screen = TRUE, ot_seed = 18, max_offtarget = 0,
    stringency = "strict")
  rf <- Sys.getenv("HCR_RELAX_FILTERS", unset = "")
  if (nzchar(rf))
    s$relax_filters <- tolower(rf) %in% c("1", "true", "yes", "on")
  sy <- Sys.getenv("HCR_STRINGENCY", unset = "")
  if (nzchar(sy)) {
    s$stringency <- tolower(sy)
    if (!s$stringency %in% c("strict", "loose"))
      stop("HCR_STRINGENCY must be 'strict' or 'loose' (got '", sy, "')")
  }
  as_ <- Sys.getenv("HCR_AMP_SET", unset = "")
  if (nzchar(as_)) {
    s$amp_set <- tolower(as_)
    if (!s$amp_set %in% c("v2", "v3"))
      stop("HCR_AMP_SET must be 'v2' or 'v3' (got '", as_, "')")
  }
  s
}

# Apply the "loose" filtering preset, mirroring the auto-relax rescues
# (GC +/-5, homopolymer +1, dG +/-15, Tm +/-10, drop ALL composition rules).
# Strict mode leaves the caller's thresholds untouched.
apply_stringency <- function(settings) {
  if (!identical(settings$stringency, "loose")) return(settings)
  settings$gc_min  <- max(settings$gc_min - 5, 0)
  settings$gc_max  <- min(settings$gc_max + 5, 100)
  settings$max_hp  <- settings$max_hp + 1
  settings$tm_min  <- settings$tm_min - 10
  settings$tm_max  <- settings$tm_max + 10
  settings$dg_min  <- settings$dg_min - 15
  settings$dg_max  <- settings$dg_max + 15
  settings$comp_a      <- FALSE
  settings$comp_c      <- FALSE
  settings$comp_stack  <- FALSE
  settings$comp_cstack <- FALSE
  settings
}

# Headless batch mode: design probes for a CSV of genes and write all output
# files to a folder - no server, no browser. Usage:
#   Rscript HCR_probe_design_v41.R my_genes.csv            (cross OFF, like the UI)
#   Rscript HCR_probe_design_v41.R my_genes.csv --cross    (enable the check)
run_batch <- function(csv_path, cross = FALSE) {
  jobs <- parse_jobs(paste(readLines(csv_path, warn = FALSE), collapse = "\n"))
  settings <- apply_stringency(default_settings())
  cache <- new.env(hash = TRUE)
  out_dir <- paste0("HCR_batch_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  dir.create(out_dir)
  cat("Batch mode: ", length(jobs), " genes from ", csv_path, "\n", sep = "")
  cat("Stringency: ", settings$stringency,
      if (identical(settings$stringency, "loose"))
        " (GC +/-5, homopolymer +1, dG +/-15, Tm +/-10, compositions off)\n"
      else " (strict defaults)\n", sep = "")
  if (cross)
    cat("Cross-species check: all detected species except each gene's own (",
        paste(SPECIES_CHOICES, collapse = ", "), ")\n", sep = "")
  else cat("Cross-species check: OFF (--no-cross)\n", sep = "")

  # Lazy reference ensure (same policy as the GUI run-guard): the design
  # species' reference set is downloaded from NCBI on first use unless the
  # NCBI fallback is on. Records provenance for HCR_design_log.txt.
  needed <- intersect(unique(vapply(jobs, function(j) canonical_species(j$species), "")),
                      names(SPECIES_GCF))
  pre_installed <- needed[vapply(needed, species_reference_installed, TRUE)]
  missing <- needed[!vapply(needed, species_reference_installed, TRUE)]
  if (length(missing) && !isTRUE(settings$ncbi_fallback)) {
    for (sp in missing) {
      cat("Downloading reference set for ", sp, " from NCBI ...\n", sep = "")
      tryCatch(download_species_ref(sp, TRUE),
               error = function(e) cat("  download failed: ",
                                       conditionMessage(e), "\n", sep = ""))
    }
  }
  refs_used <- setNames(lapply(needed, function(sp) {
    if (sp %in% pre_installed) "local (already installed)"
    else if (species_reference_installed(sp)) "freshly downloaded from NCBI"
    else "not on disk (NCBI fallback enabled: transcripts fetched per gene)"
  }), needed)

  results <- vector("list", length(jobs))
  for (i in seq_along(jobs)) {
    cat(sprintf("[%d/%d] %s ... ", i, length(jobs), jobs[[i]]$gene))
    r <- design_gene(jobs[[i]], settings, i - 1, cache)
    r$ortholog <- jobs[[i]]$ortholog
    check_sps <- if (cross)
      setdiff(SPECIES_CHOICES, tolower(jobs[[i]]$species)) else character(0)
    if (r$status == "ok" && length(r$pairs) && length(check_sps))
      r$cross_species <- tryCatch(
        cross_species_check(r, check_sps, cache, settings$ncbi_fallback),
        error = function(e) NULL)
    results[[i]] <- r
    cat(r$status, "-", length(r$pairs), "pairs\n")
    for (w in r$warnings) cat("    !", w, "\n")
  }

  pools <- build_pools(results, settings$pool_bp)

  # --- The complete, canonical output file set (identical to the ZIP) -------
  write_all_run_files(out_dir, results, settings, pools, jobs,
                      check_sps = if (cross) SPECIES_CHOICES else character(0),
                      refs_used = refs_used,
                      cross_requested = cross, run_params_txt = TRUE,
                      source_label = paste0("batch (input_csv = ", csv_path, ")"),
                      csv_path = csv_path)
  cat("\nDone. Output folder: ", normalizePath(out_dir), "\n", sep = "")
}

# Run when executed directly (Rscript app.R); skip when source()d (e.g. tests)
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  csv_arg <- args[!grepl("^--", args)][1]
  if (!is.na(csv_arg) && file.exists(csv_arg)) {
    # Cross-species check defaults to OFF (matching the interactive UI);
    # --cross enables it, --no-cross is accepted for backwards compatibility.
    cross <- "--cross" %in% args && !("--no-cross" %in% args)
    run_batch(csv_path = csv_arg, cross = cross)
    quit(save = "no")
  }
  cat("\n============================================================\n")
  cat("  HCR Probe Designer\n")
  cat("  Open your browser to: http://127.0.0.1:", PORT, "\n", sep = "")
  cat("  (The browser should open automatically)\n")
  cat("============================================================\n\n")
  cat("  Reference folder:", DATA_DIRS, "\n")
  if (dir.exists(DATA_DIRS)) {
    cat("  Reference folder found: YES\n")
    subs <- list.dirs(DATA_DIRS, recursive = FALSE, full.names = TRUE)
    if (length(subs)) {
      cat("  Species subfolders detected:\n")
      for (d in subs) {
        has_rna <- length(list.files(d, pattern = "rna\\.fna(\\.gz)?$", ignore.case = TRUE)) > 0
        has_gen <- length(list.files(d, pattern = "(^|_)genomic\\.fna(\\.gz)?$", ignore.case = TRUE)) > 0
        has_gtf <- length(list.files(d, pattern = "(^|_)genomic\\.(gtf|gff)(\\.gz)?$", ignore.case = TRUE)) > 0
        cat("   -", basename(d), ":",
            if (has_rna) "rna.fna OK" else "NO rna.fna!",
            "|", if (has_gen) "genomic.fna OK" else "no genomic.fna",
            "|", if (has_gtf) "genomic.gtf/gff OK" else "no genomic.gtf/gff", "\n")
      }
    } else {
      cat("  WARNING: no species subfolders found in reference directory\n")
    }
  } else {
    cat("  WARNING: reference folder not found - NCBI fallback will be used.\n")
    cat("  Expected:", DATA_DIRS, "\n")
  }
  cat("\n")
  shinyApp(ui, server, options = list(port = PORT, launch.browser = TRUE))
}
