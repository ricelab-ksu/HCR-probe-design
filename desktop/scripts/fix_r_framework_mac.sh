#!/bin/bash
# ============================================================================
# fix_r_framework_mac.sh — relocate a macOS R.framework so it runs from any
# location (no /Library/Frameworks install required).
#
# Background:
#   The CRAN R.framework ships with dylib load commands hardcoded to
#   /Library/Frameworks/R.framework/Versions/<ver>/Resources/lib/...
#   Those absolute paths break when the framework is bundled inside an app.
#
#   `install_name_tool` can rewrite the paths, but R is code-signed, so the
#   edits invalidate the signature. macOS then refuses to run the modified
#   binaries. Solution: rewrite every R framework binary's load commands to
#   @loader_path-relative forms and then re-sign everything ad hoc
#   (codesign -s -), which Apple allows and which clears the invalidated
#   signature.
#
#   @loader_path resolution is computed per-binary, so binaries in nested
#   directories (e.g. library/stats/libs/stats.so) correctly point back to
#   ../lib/ at any depth.
#
# Usage:
#   fix_r_framework_mac.sh /path/to/R.framework
# ============================================================================
set -euo pipefail

FW="${1:?usage: $0 /path/to/R.framework}"
RES="$FW/Resources"

if [ ! -d "$RES" ]; then
  echo "ERROR: R.framework Resources directory not found at: $RES" >&2
  exit 1
fi

OLD_PREFIX="/Library/Frameworks/R.framework/Versions"
patch_count=0

# ---------------------------------------------------------------------------
# Support helpers
# ---------------------------------------------------------------------------
rel_to() { # rel_to <from_dir> <to_dir>  -> relative path
  python3 - "$1" "$2" <<'PY'
import os, sys
print(os.path.relpath(sys.argv[2], sys.argv[1]))
PY
}

echo "Relocating R.framework: $FW"

# ---------------------------------------------------------------------------
# 1. Rewrite dylib install-name ids and every R-framework dependency path
# ---------------------------------------------------------------------------
while IFS= read -r -d '' f; do
  case "$f" in
    *.dSYM*|*.TBD*) continue ;;
  esac

  # The dylib's own install id (first line of `otool -D`)
  idline=""
  if [[ "$f" == *.dylib ]]; then
    idline=$(otool -D "$f" 2>/dev/null | sed -n 2p || true)
  fi

  args=()

  # --- id: /Library/Frameworks/.../lib/libFoo.dylib -> @loader_path/<rel>/libFoo.dylib
  if [[ -n "$idline" && "$idline" == "$OLD_PREFIX"* ]]; then
    rel=$(rel_to "$(dirname "$f")" "$RES/lib")
    args+=( -id "@loader_path/${rel}/$(basename "$idline")" )
    patch_count=$((patch_count + 1))
  fi

  # --- dependencies: /Library/Frameworks/.../lib/libBar.dylib -> @loader_path/<rel>/libBar.dylib
  while read -r dep; do
    [[ -z "$dep" ]] && continue
    if [[ "$dep" == "$OLD_PREFIX"* ]]; then
      base=$(basename "$dep")
      if [ -e "$RES/lib/$base" ] || [ -L "$RES/lib/$base" ]; then
        rel=$(rel_to "$(dirname "$f")" "$RES/lib")
        args+=( -change "$dep" "@loader_path/${rel}/${base}" )
        patch_count=$((patch_count + 1))
      else
        echo "  WARN: skipping dependency not in bundle lib/: $dep" >&2
      fi
    fi
  done < <(otool -L "$f" 2>/dev/null | sed -n '2,$p' | awk '{print $1}')

  if [ ${#args[@]} -gt 0 ]; then
    install_name_tool "${args[@]}" "$f" >/dev/null 2>&1 \
      || echo "  FAIL (install_name_tool): $f" >&2
  fi
done < <(find "$RES" -type f \( -name '*.dylib' -o -name '*.so' \) ! -path '*/.dSYM/*' ! -name '*.TBD' -print0 2>/dev/null)

# exec/R binary's own dependency on libR.dylib must be rewritten too.
for exe in "$RES/bin/exec/R" "$RES/bin/exec/Rscript"; do
  if [ -f "$exe" ]; then
    rel=$(rel_to "$(dirname "$exe")" "$RES/lib")
    while read -r dep; do
      [[ -z "$dep" ]] && continue
      if [[ "$dep" == "$OLD_PREFIX"* ]]; then
        base=$(basename "$dep")
        if [ -e "$RES/lib/$base" ] || [ -L "$RES/lib/$base" ]; then
          install_name_tool -change "$dep" "@loader_path/${rel}/${base}" "$exe" >/dev/null 2>&1 \
            || echo "  FAIL (install_name_tool exec): $exe" >&2
          patch_count=$((patch_count + 1))
        fi
      fi
    done < <(otool -L "$exe" 2>/dev/null | sed -n '2,$p' | awk '{print $1}')
  fi
done

echo "Patched $patch_count load-command paths."

# ---------------------------------------------------------------------------
# 2. Re-sign everything ad hoc so macOS accepts the rewritten binaries.
#    (.so/.dylib/executables, plus plain executables in bin/ and tools/)
# ---------------------------------------------------------------------------
sign_count=0
while IFS= read -r -d '' f; do
  case "$f" in
    *.dSYM*|*.TBD*) continue ;;
  esac
  if [[ "$f" == *".so" || "$f" == *".dylib" || -x "$f" ]]; then
    codesign --force --sign - "$f" >/dev/null 2>&1 || true
    sign_count=$((sign_count + 1))
  fi
done < <(find "$FW" -type f \( -name '*.so' -o -name '*.dylib' -o -perm +111 \) ! -path '*/.dSYM/*' -print0 2>/dev/null)

echo "Re-signed $sign_count binaries (ad hoc)."

# ---------------------------------------------------------------------------
# 3. Sanity check: no /Library/Frameworks load commands should remain.
# ---------------------------------------------------------------------------
leftover=$(find "$RES" -type f \( -name '*.dylib' -o -name '*.so' \) ! -path '*/.dSYM/*' -print0 2>/dev/null | \
  xargs -0 otool -L 2>/dev/null | grep -c "/Library/Frameworks/R\.framework" || true)
echo "Remaining /Library/Frameworks R references: $leftover"

if [ "$leftover" -gt 0 ]; then
  echo "WARNING: some binaries still reference the system R.framework." >&2
else
  echo "OK: R.framework is fully relocatable."
fi