#!/usr/bin/env bash
# Fetch and verify the arXiv ancillary files used by the Henn nonplanar
# five-point canonical-system example.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
destination=${HENN_NONPLANAR_DATA:-"$repo_root/Examples/OriginalDiffExp/Data/HennNonplanar"}
base_url=https://arxiv.org/src/1812.11160v2/anc

mkdir -p "$destination"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp2-henn-data.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

files=(
  dlogBasisXB.txt
  XB_Atilde.txt
  XB_Boundary_values_X0.txt
  XB_Boundary_values_X1.txt
)
hashes=(
  daf048c1b0b6e6462f12ae8f7c9a9033c8435024d903071c50766b81696a3d36
  3d997ab1eae2a6ca1e2489f3e733080b5b20950bc7440555b76230a5d2350c3d
  132611651f9c2fc5de983965bf170328ea0a54cef0f8cb4a3f65d87182c34a96
  7883537af39b181f84266934b9440df3f67a5eb8f67b576388d7aa29524d9b92
)

hash_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    echo "shasum or sha256sum is required" >&2
    return 2
  fi
}

for index in "${!files[@]}"; do
  file=${files[$index]}
  expected=${hashes[$index]}
  echo "Fetching $file"
  curl --fail --location --silent --show-error \
    "$base_url/$file" --output "$scratch/$file"
  actual=$(hash_file "$scratch/$file")
  if [[ "$actual" != "$expected" ]]; then
    echo "SHA-256 mismatch for $file" >&2
    echo "expected: $expected" >&2
    echo "actual:   $actual" >&2
    exit 1
  fi
  mv "$scratch/$file" "$destination/$file"
done

echo "Verified Henn nonplanar data in $destination"
