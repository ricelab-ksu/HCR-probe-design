# HCR Probe Designer — Desktop Edition

A standalone, one-click desktop application wrapping the **HCR Probe Designer v41**
R/Shiny design engine (cross-species edition) with Electron and an embedded R runtime
+ reference data, so end users don't need R, Node, or any data files to design probes.

```
┌─────────────────────────────┐   ┌──────────────────────────────┐
│  Electron (this repo)        │   │  bundled R.framework (mac)  │
│  main.js spawns R            │──▶│  exec/R --vanilla runner.R  │
│  reads HCR_PORT=NNNN         │   │  sources HCR_probe_design_  │
│  BrowserWindow → Shiny       │◀──│  v41.R, starts Shiny on 127 │
│  (127.0.0.1:NNNN)            │   │  .0.0.1:<random port>       │
└───────────────┬─────────────┘   └──────────────┬───────────────┘
                └─────────────────────────────────┘
        httpuv websocket = interactive design UI (unchanged)
```

## Why this works on macOS (the hard part)

CRAN's R framework ships with dylibs stamped with absolute
`/Library/Frameworks/...` install names. You can't just copy it into an app:

- `install_name_tool` rewrites the load commands, **but R is code-signed**, so the
  edit invalidates the signature and macOS refuses to run the binaries.
- Fix (implemented in `scripts/fix_r_framework_mac.sh`): rewrite every R dylib /
  shared-object load command to be **`@loader_path`-relative** (computed per file,
  so nested `library/*/libs/*.so` resolve back to `../lib` correctly), then
  **re-sign everything ad hoc** (`codesign -s -`, which works without a Developer
  account and clears the invalidated signature on both Intel and Apple Silicon).

Windows needs no such trick: `Rscript.exe` derives `R_HOME` from its own location
(the same mechanism R-Portable uses), so we just extract the official installer.

## Repository layout

```
desktop/
├── package.json            # electron + electron-builder
├── main.js                 # Electron main (spawn R, port handshake, window, ipc)
├── preload.js              # minimal contextBridge surface
├── electron-builder.yml    # mac DMG / win NSIS targets + extraResources
├── assets/icon.png         # generated icon
├── resources/
│   ├── app/                # runner.R + HCR_probe_design_v41.R + HCR_template.csv
│   ├── r/                  # bundled R runtime (produced by bundle scripts)
│   └── reference_files/    # 2 GB species FASTA store
├── scripts/
│   ├── bundle_r_mac.sh     # download CRAN R → extract → relocate → install pkgs
│   ├── bundle_r_win.ps1    # download R win installer → silent extract → pkgs
│   ├── fix_r_framework_mac.sh   # the install_name_tool + code-sign trick
│   ├── copy_reference_files.js  # stage reference store into resources/
│   └── make_icon.js        # zero-dependency icon generator
└── build/                  # transient (cache, temp scripts) — gitignored
```

## Quick start (development)

Requires: Node ≥ 20, macOS R ≥ 4.0 installed (or `HCR_RSCRIPT` override), the
reference store on the lab volume.

```bash
npm install
npm run dev        # runs Electron against the SYSTEM R + live reference store
```

Dev mode intentionally skips the bundled runtime and the reference copy: it
symlinks `resources/reference_files` to the live 2 GB store and calls your
`Rscript`.

## Building distributable bundles

### macOS (per arch)

```bash
npm install
./scripts/bundle_r_mac.sh 4.6.0 arm64   # or x86_64 for Intel
npm run dist:mac
```

`bundle_r_mac.sh` downloads the official CRAN `R-4.6.0-<arch>.pkg` (~170 MB),
extracts it with `pkgutil --expand-full` (no admin needed), relocates +
re-signs it, then installs `shiny` + deps into the bundled library from CRAN
macOS binaries. The result is a fully self-contained `.app` / `.dmg`.

### Windows (from a Windows machine / CI)

```powershell
npm install
powershell -ExecutionPolicy Bypass -File scripts\bundle_r_win.ps1
npm run dist:win
```

## End-user data location

Reference data ships inside the app, but on first launch Electron copies it to
`~/Library/Application Support/HCR Probe Designer/reference_files` (or
`%APPDATA%\HCR Probe Designer`) so that the app's built-in *one-click download
of species RefSeq files* can add new `d_<species>/rna.fna` folders **without
rewriting the (signed) app bundle**.

Probe-design CSVs download to the system **Downloads** folder, as in the
original lab workflow.

## Notes & limits

- **Unsigned builds**: out of the box the mac build is ad-hoc signed
  (`identity: null`). Gatekeeper will ask to confirm first-run for internal
  builds; for wider distribution add your Developer ID + notarization.
- **First launch** stages ~2 GB of reference data (instant on APFS thanks to
  clonefile); subsequent launches skip it.
- The HD port (`HCR_PORT=`) handshake makes the R server port collision-proof;
  `runner.R` picks a free port in `[30000, 60000]` via `httpuv::randomPort`.