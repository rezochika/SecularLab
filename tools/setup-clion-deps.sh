#!/usr/bin/env bash
# Called by setup-clion-deps.ps1 using an existing Git for Windows Bash.
set -euo pipefail
if test "$#" -ne 1; then
  echo 'Usage: setup-clion-deps.sh <existing-CLion-installation-folder>' >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
clion=$(cygpath -u "$1")
mingw="$clion/bin/mingw"
export PATH="$mingw/bin:$PATH"
deps="$repo/build/deps"
prefix="$deps/mingw"
source="$deps/src/gmp-6.3.0"
gmp_hash=a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
boost_hash=5e93d582aff26868d581a52ae78c7d8edf3f3064742c6e77901a1f18a437eea9
for tool in gcc g++ ar objdump curl tar sha256sum cygpath; do
  command -v "$tool" >/dev/null || { echo "Required existing tool missing: $tool" >&2; exit 1; }
done
for file in "$mingw/bin/gcc.exe" "$mingw/bin/g++.exe" "$mingw/bin/ar.exe" \
    "$mingw/include/gmp.h" "$mingw/lib/libgmp.a"; do
  test -f "$file" || { echo "Missing CLion bundled dependency: $file" >&2; exit 1; }
done
target=$(g++ -dumpmachine)
test "$target" = x86_64-w64-mingw32 || { echo "Unsupported compiler target: $target" >&2; exit 1; }
version=$(g++ -dumpfullversion)
mkdir -p "$deps/downloads" "$deps/src" "$prefix/include" "$prefix/lib"

# Check both the header ABI and linked library before downloading/building.
cat > "$deps/check-bundled-gmp.cpp" <<'CPP'
#include <gmp.h>
#include <cstring>
#include <iostream>
#include <type_traits>
#if __GNU_MP_VERSION != 6 || __GNU_MP_VERSION_MINOR != 3 || __GNU_MP_VERSION_PATCHLEVEL != 0
#error This setup supports bundled GMP 6.3.0 only
#endif
#if __GMP_LIBGMP_DLL != 0 || GMP_LIMB_BITS != 64 || GMP_NAIL_BITS != 0
#error This setup requires static GMP with 64-bit limbs and no nails
#endif
static_assert(sizeof(void*) == 8 && sizeof(long) == 4);
static_assert(sizeof(mp_limb_t) == 8 && sizeof(mp_size_t) == 4 && sizeof(mp_exp_t) == 4);
static_assert(std::is_same_v<mp_limb_t, unsigned long long>);
int main() {
  if (std::strcmp(gmp_version, "6.3.0") != 0 || mp_bits_per_limb != 64) return 1;
  mpz_t z;
  mpz_init_set_ui(z, 1);
  mpz_mul_2exp(z, z, 200);
  const bool ok = mpz_sizeinbase(z, 2) == 201 && mpz_tstbit(z, 200) == 1;
  mpz_clear(z);
  if (!ok) return 2;
  std::cout << "Bundled GMP 6.3.0/static/64-bit limb ABI: PASS\n";
}
CPP
g++ -std=c++17 -O2 -I"$mingw/include" "$deps/check-bundled-gmp.cpp" \
  "$mingw/lib/libgmp.a" -o "$deps/check-bundled-gmp.exe"
"$deps/check-bundled-gmp.exe"
imports=$(objdump -p "$deps/check-bundled-gmp.exe")
grep -qi 'DLL Name: msvcrt.dll' <<< "$imports" || {
  echo 'This setup requires the existing MSVCRT CLion MinGW runtime.' >&2; exit 1;
}
if grep -Eqi 'DLL Name:.*(gmp|ucrt|api-ms-win-crt)' <<< "$imports"; then
  echo 'Unexpected GMP DLL or UCRT runtime: this setup requires static GMP with MSVCRT.' >&2
  exit 1
fi

fetch() {
  local url=$1 archive=$2 expected=$3
  if test ! -f "$archive"; then
    curl --proto '=https' --tlsv1.2 --fail --location --retry 2 --output "$archive.part" "$url"
    printf '%s  %s\n' "$expected" "$archive.part" | sha256sum --check
    mv "$archive.part" "$archive"
  fi
  printf '%s  %s\n' "$expected" "$archive" | sha256sum --check
}
fetch https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz "$deps/downloads/gmp-6.3.0.tar.xz" "$gmp_hash"
fetch https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.gz "$deps/downloads/boost_1_90_0.tar.gz" "$boost_hash"
if test ! -f "$source/.secular-extracted-$gmp_hash"; then
  tar -xf "$deps/downloads/gmp-6.3.0.tar.xz" -C "$deps/src"
  touch "$source/.secular-extracted-$gmp_hash"
fi
if test ! -f "$prefix/include/boost/.secular-extracted-$boost_hash"; then
  tar -xf "$deps/downloads/boost_1_90_0.tar.gz" -C "$prefix/include" --strip-components 1 boost_1_90_0/boost
  touch "$prefix/include/boost/.secular-extracted-$boost_hash"
fi
cp "$mingw/include/gmp.h" "$prefix/include/gmp.h"
cp "$mingw/lib/libgmp.a" "$prefix/lib/libgmp.a"
cp "$source/gmpxx.h" "$prefix/include/gmpxx.h"
mkdir -p "$prefix/share/licenses/gmp" "$prefix/share/licenses/boost"
cp "$source"/COPYING* "$prefix/share/licenses/gmp/"
tar -xf "$deps/downloads/boost_1_90_0.tar.gz" -C "$prefix/share/licenses/boost" \
  --strip-components 1 boost_1_90_0/LICENSE_1_0.txt

# Separate configuration caches when the existing compiler or bundled GMP changes.
fingerprint=$( { printf '%s\n%s\n' "$target" "$version"; sha256sum "$mingw/include/gmp.h" "$mingw/lib/libgmp.a"; } | sha256sum)
work="$deps/gmpxx-gcc-$version-${fingerprint:0:12}"
mkdir -p "$work"
cd "$work"
if test ! -f configure.complete; then
  printf 'Configuring GMP wrapper prerequisites; log: %s\n' "$(cygpath -m "$work/configure-output.log")"
  # GCC 15 defaults to C23; GMP 6.3's configure tests require the older C mode.
  # Generic CPU selection also omits unused timing-tool assembly and its m4 need.
  # Configure embeds srcdir in C #include strings. Keep it relative so native
  # Windows GCC can resolve it without Git Bash's argument path conversion.
  if ! CFLAGS='-O2 -std=gnu11' ../src/gmp-6.3.0/configure --cache-file=config.cache \
      --disable-assembly --disable-shared --enable-static --enable-cxx \
      --host=none-w64-mingw32 --build=x86_64-w64-mingw32 \
      --prefix="$prefix" > configure-output.log 2>&1; then
    tail -n 30 configure-output.log >&2
    exit 1
  fi
  touch configure.complete
fi
# Only the C++ wrapper is built. Use the preflight-checked bundled header and
# its long-long limb ABI instead of generic configure's LLP64 unsigned long.
sed -i 's/^#define SIZEOF_MP_LIMB_T .*/#define SIZEOF_MP_LIMB_T 8/' config.h
grep -q '^#define SIZEOF_MP_LIMB_T 8$' config.h
cp "$prefix/include/gmp.h" gmp.h
printf 'Building GMP C++ wrapper with existing GCC %s\n' "$version"
for gen in gen-fib gen-fac gen-sieve gen-bases; do
  gcc -O2 -std=gnu11 "$source/$gen.c" -o "$gen.exe"
done
./gen-fib.exe header 64 0 > fib_table.h
./gen-fac.exe 64 0 > fac_table.h
./gen-sieve.exe 64 > sieve_table.h
./gen-bases.exe header 64 0 > mp_bases.h
names=(isfuns ismpf ismpq ismpz ismpznw limits osdoprnti osfuns osmpf osmpq osmpz)
objects=()
for name in "${names[@]}"; do
  g++ -std=c++17 -O2 -D__GMP_WITHIN_GMPXX -I"$prefix/include" -I"$work" -I"$source" \
    -c "$source/cxx/$name.cc" -o "$name.o"
  objects+=("$name.o")
done
ar rcs "$work/libgmpxx.a" "${objects[@]}"
cp "$work/libgmpxx.a" "$prefix/lib/libgmpxx.a.pending"
mv -f "$prefix/lib/libgmpxx.a.pending" "$prefix/lib/libgmpxx.a"

cat > "$deps/check-clion-deps.cpp" <<'CPP'
#include <gmpxx.h>
#include <boost/multiprecision/gmp.hpp>
#include <boost/version.hpp>
#include <iostream>
#include <sstream>
static_assert(BOOST_VERSION == 109000 && GMP_LIMB_BITS == 64);
int main() {
  mpz_class z = 1;
  z <<= 200;
  if (z.get_str() != "1606938044258990275541962092341162602522202993782792835301376") return 1;
  mpq_class q(18, 24);
  q.canonicalize();
  if (q != mpq_class(3, 4)) return 2;
  std::stringstream text;
  text << z << ' ' << q << ' ' << mpf_class("1.25");
  mpz_class zz; mpq_class qq; mpf_class ff;
  text >> zz >> qq >> ff;
  if (!text || zz != z || qq != q || ff != mpf_class("1.25")) return 3;
  using boost::multiprecision::mpf_float_100;
  const mpf_float_100 x("1.25");
  if (x*x != mpf_float_100("1.5625")) return 4;
  std::cout << "GMP C++ streams and Boost 1.90 GMP backend: PASS\n";
}
CPP
g++ -std=c++20 -O2 -I"$prefix/include" "$deps/check-clion-deps.cpp" \
  -L"$prefix/lib" -lgmpxx -lgmp -o "$deps/check-clion-deps.exe"
"$deps/check-clion-deps.exe"
printf 'Dependency prefix ready: %s\n' "$(cygpath -m "$prefix")"
