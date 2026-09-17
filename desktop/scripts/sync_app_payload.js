#!/usr/bin/env node
/**
 * Copy the app payload into resources/app before packaging.
 *
 * The design engine (../HCR_probe_design_v41.R) and the CSV template live at
 * the repo root and are edited there, while Electron packages whatever sits in
 * resources/app. Without this step the .app can silently ship an older engine
 * than the one in the working tree (the two files are otherwise unrelated).
 */
'use strict';
const fs = require('fs');
const path = require('path');

const DESKTOP = path.join(__dirname, '..');
const ROOT = path.join(DESKTOP, '..');
const DEST = path.join(DESKTOP, 'resources', 'app');

// [source path, file name inside resources/app]
const PAYLOAD = [
  [path.join(ROOT, 'HCR_probe_design_v41.R'), 'HCR_probe_design_v41.R'],
  [path.join(ROOT, 'HCR_template.csv'), 'HCR_template.csv'],
  [path.join(DESKTOP, 'runner.R'), 'runner.R'],
];

fs.mkdirSync(DEST, { recursive: true });
for (const [src, name] of PAYLOAD) {
  const dst = path.join(DEST, name);
  if (!fs.existsSync(src)) {
    if (fs.existsSync(dst)) continue;      // nothing to sync for this entry
    console.error('[sync] missing payload source:', src);
    process.exit(1);
  }
  if (fs.existsSync(dst) && fs.readFileSync(dst).equals(fs.readFileSync(src))) {
    console.log('[sync] up to date:', name);
    continue;
  }
  fs.copyFileSync(src, dst);
  console.log('[sync] copied ->', dst);
}
console.log('[sync] resources/app payload ready.');
