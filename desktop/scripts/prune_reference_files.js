#!/usr/bin/env node
/**
 * Prune the bundled reference payload so the installer ships no genomes.
 *
 * The reference_files folder inside the .app used to carry ~2.1 GB of
 * uncompressed FASTA/GTF plus marker CSVs. None of that is needed anymore:
 * every species' transcriptome/genome/annotation is downloaded from NCBI on
 * first use (~60-90 MB compressed) into the app's user-data directory.
 *
 * This script replaces resources/reference_files with a tiny stub. It refuses
 * to delete a symlink unless --force is given, because dev setups commonly
 * symlink the lab's live reference store and we must never delete that.
 */
'use strict';
const fs = require('fs');
const path = require('path');

const DEST = path.join(__dirname, '..', 'resources', 'reference_files');

function main() {
  if (fs.existsSync(DEST)) {
    const st = fs.lstatSync(DEST);
    if (st.isSymbolicLink()) {
      if (!process.argv.includes('--force')) {
        console.error('[prune] refusing to prune symlink at', DEST);
        console.error('[prune]   (dev setups symlink the lab reference store).');
        console.error('[prune] pass --force to replace it with the stub.');
        process.exit(1);
      }
      fs.unlinkSync(DEST);
    } else {
      fs.rmSync(DEST, { recursive: true, force: true });
    }
  }
  fs.mkdirSync(DEST, { recursive: true });
  const readme = [
    'Reference files are NOT bundled with this app anymore.',
    '',
    "Each species' reference set (rna.fna.gz + genomic.fna.gz + genomic.gff.gz)",
    'is downloaded from NCBI on first use and stored in the app user-data',
    'directory under reference_files/ (one d_<species>/ folder each).',
    '',
    'Nothing else is required: the marker CSVs, XLSX/RTF exports, fbgn_map and',
    'dmel_gene_reference that used to ship here are not read by the app.',
    '',
  ].join('\n');
  fs.writeFileSync(path.join(DEST, 'README.txt'), readme + '\n');
  console.log('[prune] resources/reference_files is now a stub (README.txt only).');
}

main();