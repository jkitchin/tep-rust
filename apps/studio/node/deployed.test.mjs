// The module that actually gets deployed still computes the right numbers.
//
// `apps/studio/build.sh` puts the wasm through wasm-bindgen and then, when
// binaryen is installed, through `wasm-opt -Oz`. Both rewrite the module. This
// asserts that neither changed a single bit of arithmetic, by instantiating the
// built artifact and asking it for the Tier 9 self-check digest, which
// `crates/tepsim/src/tier9.rs` commits as a constant and which the native
// build, three wasm profiles and the browser transport path all agree on.
//
// It is not a duplicate of `cargo xtask tier9`. That checks the module cargo
// produces; this checks the one that ships, after two more tools have touched
// it. `-Oz` is an optimiser, and an optimiser that changed a floating-point
// result would be a serious bug in it and an even more serious one for a
// project whose whole claim is bit-level agreement with Fortran.
//
// Skipped, not failed, when `dist/` has not been built: a checkout that has
// never run `build.sh` has nothing to check, and `npm test` should still work.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const wasmPath = join(here, '..', 'dist', 'js', 'tepsim_wasm_bg.wasm');

/// `tepsim::tier9::CASES[0]`: one hour of the fault-free plant, closed loop,
/// under the faithful Euler integrator, from `teprob.f`'s own seed.
const EXPECTED = 'c8a26889992f1719';

test('the built module reproduces the Tier 9 self-check digest', async (t) => {
  if (!existsSync(wasmPath)) {
    t.skip('dist/ has not been built; run apps/studio/build.sh');
    return;
  }

  const module = await WebAssembly.compile(readFileSync(wasmPath));

  // Every wasm-bindgen import is stubbed with a function that throws. Nothing
  // on the determinism path should cross into JavaScript, and a stub that
  // throws turns "it quietly did" into a loud failure rather than a silent
  // pass. The import module name comes from the module itself so that a
  // wasm-bindgen rename does not turn this into a false failure.
  const imports = WebAssembly.Module.imports(module);
  const namespaces = {};
  for (const { module: ns, name } of imports) {
    namespaces[ns] ??= {};
    namespaces[ns][name] = () => {
      throw new Error(`the digest path called into JavaScript: ${ns}.${name}`);
    };
  }

  const instance = await WebAssembly.instantiate(module, namespaces);
  const digest = instance.exports.tepsim_wasm_self_check_digest();
  const hex = BigInt.asUintN(64, digest).toString(16).padStart(16, '0');

  assert.equal(
    hex,
    EXPECTED,
    `the deployed module computes ${hex}, not ${EXPECTED}: wasm-bindgen or ` +
      `wasm-opt changed the arithmetic, or the model changed and the Tier 9 ` +
      `table was not re-baselined`,
  );
  console.log(`deployed module digest ${hex}, ${imports.length} JS import(s) stubbed`);
});
