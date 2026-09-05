#!/usr/bin/env bash
# Fetch and verify the arXiv ancillary files used by the four canonical
# planar one-mass pentagon systems.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
destination=${PLANAR_ONE_MASS_DATA:-"$repo_root/examples/original/Data/PlanarOneMass"}
base_url=https://arxiv.org/src/2005.04195v1/anc

mkdir -p "$destination"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp-planar-one-mass-data.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

files=(
  alphabet.m
  1loop/diffEq-1loop.m
  1loop/numIntegrals-1loop.m
  zmz/diffEq-zmz.m
  zmz/numIntegrals-zmz.m
  mzz/diffEq-mzz.m
  mzz/numIntegrals-mzz.m
  zzz/diffEq-zzz.m
  zzz/numIntegrals-zzz.m
  zzz/pureBasis-zzz.m
)
hashes=(
  7d685bcd1028b55b974039fe99acd5d1ec0e712af75223212c27040e04038e86
  8436b42770cea4de9ed8ff3e95acf1897af180081a3907af672938fbf8605016
  e772dd2e9402004f40ea72b1bd8ffca5200a2b91f37bd7ed6eb8bd4b8c982642
  7345c3d1fd8e6b248668a6ce072d90866971e73d4663b07d3b43cd6c98e3f018
  abd31762bef73caa217f500422a8b23ccc97ea4e222faeddfba6cd7b5af585c8
  bb768d48123765fce8622709f896ec3cf954c068a53f2fa89b6325cc0f333a1f
  e2424eb7cdc0d994d72bc0031bcc938f2ea562b2727ecbfa798655547652c74f
  43f617d7929ee5f6aaef9d808f76ea2248c3961b0677a3632f72f6c83a9a29b1
  184dfc6a105164776d1406ae0b5e7e3b9f2447091b02479604613fe065b0df66
  68368e419471556e319ec6323580d216af730c372e6b5b55611b40e4f2522369
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
  mkdir -p "$scratch/$(dirname "$file")" "$destination/$(dirname "$file")"
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

echo "Verified planar one-mass data in $destination"
