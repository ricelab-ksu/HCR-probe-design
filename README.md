# HCR Probe Designer

High-throughput **HCR v2.0/v3.0 probe-set design** (amplifiers B1–B5 by
default; B1–B17 optional) from a gene-list CSV, with cross-species probe
conservation checks and a standalone Electron desktop app (bundled R runtime —
no R install required for users).

## What it does

For each gene in your CSV it designs 52-bp probe targets (two 25-mers + a
2-bp gap) and generates the full probe set. Construction replicates the lab's
"HCR PROBE MAKER Rice lab" sheet:

- **52-bp target** = two 25-mers + a 2-bp gap
- **Each oligo** = 18-nt initiator half + 2-nt spacer + 25-nt antisense half
  (45 nt total, 90 bp per pair)
- **Amplifier is *not* a user choice** — the engine rotates B1–B5 internally
  for oligo naming, and the **primary deliverable** (HCR PROBE MAKER all-amp
  CSV) lists every amplifier in the selected set for each 52-bp target so the
  amplifier can be chosen at the bench. The lab's in-situ kit is HCR v2.0, so
  **B1–B5 is the default**; HCR v3.0 adds B7/B9/B10/B13–B15/B17.

## Downloads

Pre-built macOS apps (signed ad-hoc; Gatekeeper will ask to confirm on first
open) are published on the **[Releases page][releases]**:

- `HCR Probe Designer-<ver>-arm64.dmg` — Apple Silicon
- `HCR Probe Designer-<ver>.dmg`      — Intel (x64)

## Features

- **CSV input** (`gene`, optional `species`, `amplifier`, `transcript`,
  `sequence`, `ortholog` columns) — see `HCR_template.csv`.
- **Amplifier set toggle** — the lab's in-situ kit is **HCR v2.0**, so probes
  use **B1–B5 by default**; select **HCR v3.0** in Settings to also load
  B7/B9/B10/B13–B15/B17. The HCR PROBE MAKER all-amp CSV lists every
  amplifier in the chosen set.
- **Lazy NCBI reference download** — no genomes ship with the app; each
  species' RefSeq set (`rna.fna.gz` + `genomic.fna.gz` + `genomic.gff.gz`) is
  fetched on first use. Both GTF and GFF3 attribute styles are parsed.
- **Filtering stringency** (Strict / Loose) plus full thermodynamic filters
  (GC%, homopolymer, dG, Tm, composition).
- **Cross-species conservation check** — reuse a probe set in other species;
  per-pair verdicts (`exact` / `near` / `risky`) with ortholog accessions.
- **Download all (ZIP)** — every output file of a run, for lab record keeping:
  - `HCR_probe_maker_all_amplifiers.csv` (primary; every amplifier in the set)
  - `HCR_oligo_order.csv`
  - `HCR_antisense_halves.csv` (Benchling primer upload)
  - `HCR_pool_summary.csv`
  - `HCR_cross_species_conservation.csv`
  - `HCR_design_log.txt` (definitive reproducibility record)
  - `HCR_shared_CDS_sequences.fasta`
  - `HCR_probe_targets.fasta`
  - `HCR_plots.pdf` (one page per gene: isoform track + cross-species track)

## Repository layout

```
HCR_probe_design_v41.R     # the design engine (shiny UI/server + batch)
HCR_template.csv           # input CSV template
desktop/                   # Electron shell + build scripts
  main.js, preload.js      # Electron main/preload
  package.json             # version + scripts
  electron-builder.yml     # packaging config
  resources/app/           # synced payload (engine + runner.R)
  scripts/                 # R bundling / payload sync / icon tooling
```

## Running from source (no Electron)

Requirements: R (≥ 4.1) with the `shiny` package.

```sh
# Interactive app (opens a browser):
Rscript HCR_probe_design_v41.R

# Headless batch mode (cross-species check OFF by default, like the UI):
Rscript HCR_probe_design_v41.R my_genes.csv
Rscript HCR_probe_design_v41.R my_genes.csv --cross   # enable the check

# Batch settings overrides:
HCR_STRINGENCY=loose Rscript HCR_probe_design_v41.R my_genes.csv
HCR_RELAX_FILTERS=1  Rscript HCR_probe_design_v41.R my_genes.csv
```

Reference files are read from the directory named by `HCR_REFERENCE_DIR`, or
`~/Downloads/reference_files`, or the lab's reference store. Each species lives
in a `d_<species>/` subfolder with `rna.fna(.gz)` (+ `genomic.fna/.gff` for
isoform plots and shared-CDS design).

## Build the desktop app

```sh
cd desktop
npm install
npm run bundle:mac        # (re)download + bundle the R runtime
npm run dist:mac          # produces the .dmg in desktop/dist/
```

## Probe chemistry guardrails

Initiator/spacer sequences, the 18+2+25 / 52-bp / 2-bp-gap construction rules,
and the input CSV contract are fixed by design and must not be changed casually.

[releases]: https://github.com/ricelab-ksu/HCR-probe-design/releases
