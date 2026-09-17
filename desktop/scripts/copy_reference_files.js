#!/usr/bin/env node
/**
 * copy_reference_files.js — synchronize the reference-file store into
 * resources/reference_files for packaging.
 *
 * The dev workspace keeps `resources/reference_files` as a symlink to the live
 * reference dir. electron-builder dereferences symlinks on package, but for
 * clean, explicit builds this script performs a real copy (cloning moves fast
 * on APFS). Safe to run repeatedly.
 *
 * Usage:
 *   node scripts/copy_reference_files.js [--skip-if-present] [--no-symlink]
 */
'use strict';

const fs = require('fs');
const path = require('path');

const ROOT = path.join(__dirname, '..');
const DEST = path.join(ROOT, 'resources', 'reference_files');
const SOURCES = [
  process.env.HCR_REFERENCE_SOURCE,
  '/Volumes/Backup_Plus/HCR_probe_designer/reference_files', // canonical lab volume
  path.join(process.env.HOME || '.', 'Downloads', 'reference_files'),
].filter(Boolean);

const args = process.argv.slice(2);
const skipIfPresent = args.includes('--skip-if-present');
const disallowRelink = args.includes('--no-symlink');

function humanSize(n) {
  const u = ['B', 'KB', 'MB', 'GB'];
  let i = 0;
  while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
  return `${n.toFixed(1)} ${u[i]}`;
}

function findSource() {
  for (const s of SOURCES) if (fs.existsSync(s)) return s;
  return null;
}

(async () => {
  const src = findSource();
  if (!src) {
    console.error('copy_reference_files: no source reference_files found. ' +
      'Set HCR_REFERENCE_SOURCE or place files at one of:');
    SOURCES.forEach((s) => console.error('  ' + s));
    process.exit(1);
  }

  const stat = fs.lstatSync(DEST, { throwIfNoEntry: false });
  if (skipIfPresent && (stat || fs.existsSync(DEST))) {
    console.log(`copy_reference_files: ${DEST} already present, skipping.`);
    return;
  }

  // Dev mode: use a symlink so redesigns don't re-copy 2 GB on every run.
  if (!disallowRelink && process.platform !== 'win32') {
    fs.rmSync(DEST, { recursive: true, force: true });
    fs.symlinkSync(src, DEST, 'dir');
    console.log(`symlinked ${DEST}\n        -> ${src}`);
    return;
  }

  console.log(`copying ${src} -> ${DEST}`);
  fs.rmSync(DEST, { recursive: true, force: true });
  fs.mkdirSync(path.dirname(DEST), { recursive: true });
  fs.cpSync(src, DEST, { recursive: true });
  const bytes = fs.readdirSync(DEST).reduce(
    (acc, e) => acc + fs.statSync(path.join(DEST, e)).size, 0);
  console.log(`copied ${humanSize(bytes)} of reference data.`);
})();