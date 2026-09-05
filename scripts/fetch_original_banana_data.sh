#!/usr/bin/env bash
# Fetch and verify the differential-equation matrices used by the original
# DiffExp Banana.nb example.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
destination=${ORIGINAL_BANANA_DATA:-"$repo_root/examples/original/Data/Banana"}
revision=784c8229bf92369a03f011a48e161522c8c54bbd
base_url="https://gitlab.com/hiddingm/diffexp/-/raw/$revision/Examples"

mkdir -p "$destination/EqualMass" "$destination/UnequalMass"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/diffexp-banana-data.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

files=(
  EqualMass/dt_0.m
  EqualMass/dt_1.m
  UnequalMass/dmm1_0.m
  UnequalMass/dmm1_1.m
  UnequalMass/dmm2_0.m
  UnequalMass/dmm2_1.m
  UnequalMass/dmm3_0.m
  UnequalMass/dmm3_1.m
  UnequalMass/dmm4_0.m
  UnequalMass/dmm4_1.m
  UnequalMass/dpsq_0.m
  UnequalMass/dpsq_1.m
)
sources=(
  Banana_EqualMass_Matrices/dt_0.m
  Banana_EqualMass_Matrices/dt_1.m
  Banana_Matrices/dmm1_0.m
  Banana_Matrices/dmm1_1.m
  Banana_Matrices/dmm2_0.m
  Banana_Matrices/dmm2_1.m
  Banana_Matrices/dmm3_0.m
  Banana_Matrices/dmm3_1.m
  Banana_Matrices/dmm4_0.m
  Banana_Matrices/dmm4_1.m
  Banana_Matrices/dpsq_0.m
  Banana_Matrices/dpsq_1.m
)
hashes=(
  d87a9039ab70226a0afbb4c4c2c37f2fbfde6c638c30a2771087734ae982ef1b
  c4c597adeda5872d532d3edb82a11bdfbd04456f39e6bd1d5fd538390f96fb71
  a053ab7ee079fd98cf315e1305c46fbdb6d30af6a2854915ed359ad62588be2f
  cf6e7ac6c538869dca55ad18e9cfea2cdeea04ef2164798f3a3efc6d59690a67
  88f010001b717903fbeceb44765691033fc6742f9d480a2fd6846b0077e32c84
  04d71c28c9cb5c5643503aee26eec6486c187eb260b5e0e0aa703f9dcd6299e8
  32bfc4a3c8b89991bbaec57d254e2627cb6dbb7be14116152972aeefae73601c
  1e17d878e87445626c3b8109b5df45c2b851ee9b448cb924b760b94442171e48
  5d6d7613337eb7cefdc212d67f78ff0f98cdd8e4c7853283fd882429537e7b1a
  07079a8a7405176618a3df5edac2e9247c0dce2ff6216ee443fe0593d97cda4f
  22c28c28cb1d3fbc0eb1aec7e3a28756c23f1262f03261f0a6af8207263ca7cd
  56147a32a26068bc84e9f2381399434f46aab4041138b0789e9c9492312ba9c7
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
  source=${sources[$index]}
  expected=${hashes[$index]}
  mkdir -p "$scratch/$(dirname "$file")"
  echo "Fetching $source"
  curl --fail --location --silent --show-error \
    "$base_url/$source" --output "$scratch/$file"
  actual=$(hash_file "$scratch/$file")
  if [[ "$actual" != "$expected" ]]; then
    echo "SHA-256 mismatch for $source" >&2
    echo "expected: $expected" >&2
    echo "actual:   $actual" >&2
    exit 1
  fi
  mv "$scratch/$file" "$destination/$file"
done

echo "Verified original DiffExp Banana data in $destination"
