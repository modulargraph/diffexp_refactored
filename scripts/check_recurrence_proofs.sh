#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../docs/proofs/lean"
LEAN_BIN="${LEAN_BIN:-lean}"
"$LEAN_BIN" --version
proof_tmp=$(mktemp -d "${TMPDIR:-/tmp}/diffexp-proofs.XXXXXX")
trap 'rm -rf "$proof_tmp"' EXIT
export LEAN_PATH="$proof_tmp${LEAN_PATH:+:$LEAN_PATH}"
for module in Contracts FiniteQuotient EnclosurePolicy UncertaintyWindow; do
  "$LEAN_BIN" -o "$proof_tmp/$module.olean" "$module.lean"
done
