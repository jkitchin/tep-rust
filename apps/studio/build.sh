#!/usr/bin/env bash
#
# Build TEP Studio into apps/studio/dist.
#
# Three steps and no bundler: compile crates/tepsim-wasm to wasm32, generate the
# wasm-bindgen glue beside it, copy this directory over the top. The output is a
# directory of static files that any host which can serve files can serve, which
# is the whole deployment story (GitHub Pages, a Hugging Face Static Space) and
# the reason this app avoids SharedArrayBuffer and the COOP/COEP headers neither
# of them can set.
#
# Requirements:
#
#   rustup target add wasm32-unknown-unknown
#   cargo install wasm-bindgen-cli --version <the one in Cargo.lock>
#
# The wasm-bindgen-cli version must match the wasm-bindgen dependency exactly.
# It refuses to run otherwise, which is correct: the glue and the module share a
# private ABI that changes between releases. This script reads the required
# version out of Cargo.lock rather than hard-coding it, so a dependency bump
# produces a clear error here instead of a mystery at run time.
#
# wasm-opt is used when it is on PATH and skipped with a note when it is not.
# Measured on this module it is worth about 9 percent of the raw size and
# essentially nothing compressed (87,298 to 79,759 bytes raw; 35,359 to 35,284
# gzipped), because gzip already finds most of what -Oz does. It is an optimizer
# and not a correctness step, and the determinism self-check on the page runs
# either way, so a missing wasm-opt is a note rather than a failure.
#
# The feature flags are not optional. rustc 1.97 emits bulk-memory instructions
# (`memory.copy`), non-trapping float conversions, and sign-extension operators
# (`i32.extend16_s`), and a wasm-opt that is not told to expect them rejects the
# module outright at validation. Trunk's own bundled binaryen hits exactly this,
# and so does Ubuntu's: `apt install binaryen` on noble gives version 108, from
# 2021, which failed the Pages build on sign-extension until `--enable-sign-ext`
# was added here. A binaryen new enough to enable these by default does not need
# the flags and is not harmed by them.
#
# Usage:
#   apps/studio/build.sh              release, optimized
#   apps/studio/build.sh --dev        debug profile, much faster to compile
#   apps/studio/build.sh --no-opt     skip wasm-opt even if it is installed

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${here}/../.." && pwd)"
dist="${here}/dist"

profile="release-wasm"
profile_dir="release-wasm"
run_opt=1

for arg in "$@"; do
  case "${arg}" in
    --dev) profile="dev"; profile_dir="debug" ;;
    --no-opt) run_opt=0 ;;
    *) echo "unknown option: ${arg}" >&2; exit 2 ;;
  esac
done

# The wasm-bindgen version Cargo.lock pins, so the CLI check below compares
# against the real requirement rather than against a number typed here.
wanted="$(
  awk '/^name = "wasm-bindgen"$/ { found = 1; next }
       found && /^version = / { gsub(/[",]/, "", $3); print $3; exit }' \
    "${root}/Cargo.lock"
)"

if ! command -v wasm-bindgen >/dev/null 2>&1; then
  echo "wasm-bindgen is not on PATH. Install it with:" >&2
  echo "  cargo install wasm-bindgen-cli --version ${wanted} --locked" >&2
  exit 1
fi

have="$(wasm-bindgen --version | awk '{print $2}')"
if [ "${have}" != "${wanted}" ]; then
  echo "wasm-bindgen CLI is ${have}, Cargo.lock wants ${wanted}." >&2
  echo "  cargo install wasm-bindgen-cli --version ${wanted} --locked" >&2
  exit 1
fi

echo "==> cargo build -p tepsim-wasm --target wasm32-unknown-unknown (${profile})"
if [ "${profile}" = "dev" ]; then
  cargo build --manifest-path "${root}/Cargo.toml" -p tepsim-wasm \
    --target wasm32-unknown-unknown
else
  cargo build --manifest-path "${root}/Cargo.toml" -p tepsim-wasm \
    --target wasm32-unknown-unknown --profile "${profile}"
fi

wasm="${root}/target/wasm32-unknown-unknown/${profile_dir}/tepsim_wasm.wasm"

rm -rf "${dist}"
mkdir -p "${dist}/js"

echo "==> wasm-bindgen --target web"
wasm-bindgen --target web --out-dir "${dist}/js" --no-typescript "${wasm}"

if [ "${run_opt}" = "1" ] && command -v wasm-opt >/dev/null 2>&1; then
  echo "==> wasm-opt -Oz"
  wasm-opt -Oz --enable-bulk-memory --enable-nontrapping-float-to-int \
    --enable-sign-ext \
    "${dist}/js/tepsim_wasm_bg.wasm" -o "${dist}/js/tepsim_wasm_bg.wasm"
else
  if [ "${run_opt}" = "1" ]; then
    echo "==> wasm-opt not found, skipping (install binaryen for a smaller module)"
  fi
fi

echo "==> copying the page"
cp "${here}/index.html" "${dist}/index.html"
cp "${here}/studio.css" "${dist}/studio.css"
cp "${here}"/js/*.js "${dist}/js/"

# Attribution travels with the artifact, not merely with the repository: the
# upstream NCSA licence's conditions survive into binary distributions, and a
# deployed dist/ is a binary distribution. PLAN.org, "Documentation".
for file in LICENSE LICENSE-NCSA NOTICE.md; do
  if [ -f "${root}/${file}" ]; then cp "${root}/${file}" "${dist}/${file}"; fi
done

echo
echo "==> sizes"
"${here}/measure.sh" "${dist}"

echo
echo "Serve it:"
echo "  python3 -m http.server 8000 --directory ${dist}"
echo "  open http://localhost:8000/"
