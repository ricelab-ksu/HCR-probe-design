#!/usr/bin/env node
/**
 * make_icon.js — generate a 1024×1024 PNG app icon with zero dependencies.
 *
 * Draws a rounded-rect "gene/DNA chip" motif: dark indigo background, a
 * double-helix gradient, and the marker "HCR" band. Encodes a PNG by hand
 * (zlib + crc32) so there is no asset pipeline to install.
 *
 * Usage: node scripts/make_icon.js [output.png]
 */
'use strict';

const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

const SIZE = 1024;
const outFile = process.argv[2] || path.join(__dirname, '..', 'assets', 'icon.png');

// ---- minimal PNG encoder (RGB8 + deflate, no filters, no alpha) ----
const CRC_TABLE = new Int32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  CRC_TABLE[n] = c;
}
function crc32(buf) {
  let c = 0xffffffff;
  for (let i = 0; i < buf.length; i++) c = CRC_TABLE[(c ^ buf[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}
function chunk(type, data) {
  const len = Buffer.alloc(4);
  len.writeUInt32BE(data.length, 0);
  const typeBuf = Buffer.from(type, 'ascii');
  const crcBuf = Buffer.alloc(4);
  crcBuf.writeUInt32BE(crc32(Buffer.concat([typeBuf, data])), 0);
  return Buffer.concat([len, typeBuf, data, crcBuf]);
}
function encodePNG(width, height, rgba) {
  const raw = Buffer.alloc(height * (1 + width * 4));
  let off = 0;
  for (let y = 0; y < height; y++) {
    raw[off++] = 0; // filter type 0
    for (let x = 0; x < width; x++) {
      const i = (y * width + x) * 4;
      raw[off++] = rgba[i];
      raw[off++] = rgba[i + 1];
      raw[off++] = rgba[i + 2];
      raw[off++] = rgba[i + 3];
    }
  }
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(width, 0);
  ihdr.writeUInt32BE(height, 4);
  ihdr[8] = 8;  // bit depth
  ihdr[9] = 6;  // color type RGBA
  return Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk('IHDR', ihdr),
    chunk('IDAT', zlib.deflateSync(raw, { level: 9 })),
    chunk('IEND', Buffer.alloc(0)),
  ]);
}

// ---- drawing helpers on an RGBA canvas ----
const px = new Uint8Array(SIZE * SIZE * 4);

function setPx(x, y, r, g, b, a = 255) {
  const i = (y * SIZE + x) * 4;
  px[i] = r; px[i + 1] = g; px[i + 2] = b; px[i + 3] = a;
}
function blend(x, y, r, g, b, a) {
  const i = (y * SIZE + x) * 4;
  const sa = a / 255;
  const da = px[i + 3] / 255;
  const oa = sa + da * (1 - sa);
  if (oa <= 0) return;
  px[i] = Math.round((r * sa + px[i] * da * (1 - sa)) / oa);
  px[i + 1] = Math.round((g * sa + px[i + 1] * da * (1 - sa)) / oa);
  px[i + 2] = Math.round((b * sa + px[i + 2] * da * (1 - sa)) / oa);
  px[i + 3] = Math.round(oa * 255);
}
const lerp = (a, b, t) => a + (b - a) * t;
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

function inRoundedRect(x, y, x0, y0, x1, y1, rad) {
  if (x < x0 || x > x1 || y < y0 || y > y1) return false;
  const cx = Math.max(x0 + rad, Math.min(x, x1 - rad));
  const cy = Math.max(y0 + rad, Math.min(y, y1 - rad));
  const dx = x - cx, dy = y - cy;
  return dx * dx + dy * dy <= rad * rad;
}

// ---- compose the icon ----
const M = 64;   // margin
const x0 = M, y0 = M, x1 = SIZE - M, y1 = SIZE - M;
const rad = 210;

// Background: vertical indigo gradient
for (let y = 0; y < SIZE; y++) {
  const t = y / SIZE;
  const r = Math.round(lerp(30, 17, t));
  const g = Math.round(lerp(34, 22, t));
  const b = Math.round(lerp(59, 36, t));
  for (let x = 0; x < SIZE; x++) setPx(x, y, r, g, b);
}

// Rounded-square clip mask shape (the "chip")
for (let y = 0; y < SIZE; y++) {
  for (let x = 0; x < SIZE; x++) {
    if (!inRoundedRect(x, y, x0, y0, x1, y1, rad)) {
      for (let k = 0; k < 4; k++) px[(y * SIZE + x) * 4 + k] = 0; // transparent corner
    }
  }
}

// DNA double-helix motif: two sine strands with crossing circles
const cy = SIZE / 2;
const amp = 210;
const period = 240;
const strandW = 10;
for (let y = M + 70; y < SIZE - M - 70; y++) {
  const t = (y - (M + 70)) / (SIZE - 2 * M - 140);
  const s = Math.sin((y / period) * Math.PI * 2);
  const xA = SIZE / 2 - s * amp;
  const xB = SIZE / 2 + s * amp;
  for (let w = 0; w < strandW; w++) {
    for (let x = xA - w / 2; x < xA + w / 2 + 1; x++) blend(Math.round(x), y, 166, 227, 161, 235);      // mint
    for (let x = xB - w / 2; x < xB + w / 2 + 1; x++) blend(Math.round(x), y, 137, 180, 250, 235);      // blue
  }
  // rung segments (draw every 6 rows)
  if (Math.floor(y) % 6 === 0) {
    for (let x = Math.round(xA); x < Math.round(xB); x++) blend(x, y, 148, 226, 213, 90);
  }
  // "bases": small dots along strands
  if (Math.floor(y) % 28 === 0) {
    const bx = SIZE / 2 + s * amp * 1.18;
    for (let k = 0; k < 38; k++) {
      const ang = (k / 37) * Math.PI * 2;
      blend(Math.round(bx + Math.cos(ang) * 13), y + Math.round(Math.sin(ang) * 13), 250, 179, 135, 120);
    }
  }
}

// Bottom "HCR · PROBE DESIGNER" label band
const bandY = SIZE - 230;
for (let y = bandY; y < SIZE - M - 20; y++) {
  for (let x = M + 60; x < SIZE - M - 60; x++) {
    if (inRoundedRect(x, y, M + 60, bandY, SIZE - M - 60, SIZE - M - 20, 54)) {
      blend(x, y, 15, 17, 26, 230);
    }
  }
}
// Simple horizontal "sequence bars" in the label area
const barYs = [SIZE - 210, SIZE - 196, SIZE - 182];
const barLens = [0.62, 0.9, 0.48];
const barColors = [[166, 227, 161], [137, 180, 250], [250, 179, 135]];
const barX = M + 100;
for (let bi = 0; bi < barYs.length; bi++) {
  const len = Math.round((SIZE - 2 * M - 200) * barLens[bi]);
  for (let y = barYs[bi]; y < barYs[bi] + 10; y++)
    for (let x = barX; x < barX + len; x++) {
      if (x >= M + 80 && x <= SIZE - M - 80) blend(x, y, ...barColors[bi], 255);
    }
}

fs.mkdirSync(path.dirname(outFile), { recursive: true });
fs.writeFileSync(outFile, encodePNG(SIZE, SIZE, px));
console.log(`Wrote ${outFile} (${SIZE}×${SIZE} PNG, ${fs.statSync(outFile).size} bytes)`);