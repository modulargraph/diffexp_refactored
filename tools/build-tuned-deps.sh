#!/bin/sh
# Builds GMP and FLINT for the CPU of this machine into PREFIX, for DiffExp to
# link against. Distribution packages of GMP are generic builds; a host build
# selects the multiply-and-carry code of the actual processor, and FLINT is
# compiled with -march=native.
#
#   tools/build-tuned-deps.sh PREFIX [JOBS]
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=PREFIX
#   PKG_CONFIG_PATH=PREFIX/lib/pkgconfig cmake --build build
#
# Needs a C compiler, make, m4, curl, tar and xz. MPFR is taken from the
# system. The result is specific to this CPU model: do not copy it elsewhere.
set -eu
PREFIX=${1:?usage: build-tuned-deps.sh PREFIX [JOBS]}
JOBS=${2:-2}
GMP=${GMP_VERSION:-6.3.0}
FLINT=${FLINT_VERSION:-3.4.0}
mkdir -p "$PREFIX/src"
PREFIX=$(cd "$PREFIX" && pwd)
cd "$PREFIX/src"
fetch() { # fetch FILE URL...
  file=$1; shift
  [ -s "$file" ] && return 0
  for url in "$@"; do
    if curl -fsSL --retry 2 -o "$file.part" "$url"; then mv "$file.part" "$file"; return 0; fi
  done
  echo "could not download $file" >&2; exit 1
}
if [ ! -f "$PREFIX/lib/libgmp.a" ] && [ ! -f "$PREFIX/lib/libgmp.so" ] && [ ! -f "$PREFIX/lib/libgmp.dylib" ]; then
  fetch "gmp-$GMP.tar.xz" "https://ftp.gnu.org/gnu/gmp/gmp-$GMP.tar.xz" "https://gmplib.org/download/gmp/gmp-$GMP.tar.xz" \
    "https://ftpmirror.gnu.org/gmp/gmp-$GMP.tar.xz" "https://archive.ubuntu.com/ubuntu/pool/main/g/gmp/gmp_$GMP+dfsg.orig.tar.xz"
  rm -rf "gmp-$GMP"; mkdir "gmp-$GMP"; tar -xf "gmp-$GMP.tar.xz" -C "gmp-$GMP" --strip-components=1
  # A source package without its documentation still configures doc/.
  if [ ! -f "gmp-$GMP/doc/Makefile.in" ]; then
    mkdir -p "gmp-$GMP/doc"
    printf 'all install install-exec install-data install-strip uninstall check installcheck installdirs clean mostlyclean distclean maintainer-clean dvi pdf ps info html tags ctags distdir:\n' > "gmp-$GMP/doc/Makefile.in"
  fi
  (cd "gmp-$GMP" && ./configure --prefix="$PREFIX" --enable-cxx --disable-static >/dev/null && make -j"$JOBS" >/dev/null && make install >/dev/null)
fi
if [ ! -f "$PREFIX/lib/pkgconfig/flint.pc" ]; then
  fetch "flint-$FLINT.tar.gz" "https://github.com/flintlib/flint/releases/download/v$FLINT/flint-$FLINT.tar.gz" \
    "https://flintlib.org/download/flint-$FLINT.tar.gz"
  rm -rf "flint-$FLINT"; mkdir "flint-$FLINT"; tar -xf "flint-$FLINT.tar.gz" -C "flint-$FLINT" --strip-components=1
  (cd "flint-$FLINT" && ./configure --prefix="$PREFIX" --with-gmp="$PREFIX" --enable-arch=native --disable-static >/dev/null &&
    make -j"$JOBS" >/dev/null && make install >/dev/null)
fi
echo "GMP $GMP and FLINT $FLINT for this CPU are in $PREFIX. Configure DiffExp with:"
echo "  PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$PREFIX"
