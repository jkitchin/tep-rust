// Static checks on the page itself.
//
// Rendering needs a browser and these tests do not pretend otherwise. What they
// catch is the class of breakage that a browser reports as a blank page and a
// line in a console nobody is reading: an element id renamed in the HTML but
// not in the script, a stylesheet the build forgot to copy, a `Math.random`
// that wandered into a file whose whole purpose is to be reproducible.
//
// Everything here reads the built `dist`, not the sources, so it also proves
// the build assembled what it claimed to.

import test from "node:test";
import assert from "node:assert/strict";
import { readFile, readdir } from "node:fs/promises";
import path from "node:path";

import { DIST } from "./harness.mjs";

const html = await readFile(path.join(DIST, "index.html"), "utf8");
// A lint that trips on its own documentation is useless, and these files
// explain at length why they do not use `Math.random` or `SharedArrayBuffer`.
// So the forbidden-pattern checks read code, not prose.
//
// Block comments go, and so do whole lines that are comments. Trailing
// comments after code are left alone, which is a deliberate limitation rather
// than an oversight: stripping those needs a tokenizer to avoid eating the
// `//` in a URL string, and the files here put comments on their own lines.
function stripComments(source) {
  return source
    .replace(/\/\*[\s\S]*?\*\//g, "")
    .split("\n")
    .filter((line) => !/^\s*\/\//.test(line))
    .join("\n");
}

const scripts = Object.fromEntries(
  await Promise.all(
    (await readdir(path.join(DIST, "js")))
      .filter((f) => f.endsWith(".js") && !f.startsWith("tepsim_wasm"))
      .map(async (f) => [f, await readFile(path.join(DIST, "js", f), "utf8")]),
  ),
);

test("every element the app looks up exists in the page", () => {
  const ids = new Set([...html.matchAll(/\bid="([^"]+)"/g)].map((m) => m[1]));
  const looked = new Set(
    [...scripts["app.js"].matchAll(/\$\("([^"]+)"\)/g)].map((m) => m[1]),
  );
  assert.ok(looked.size > 20, `only ${looked.size} lookups found, regex is wrong`);
  for (const id of looked) {
    assert.ok(ids.has(id), `app.js looks up #${id}, which is not in index.html`);
  }
});

test("the page's own assets were built", async () => {
  const refs = [
    ...[...html.matchAll(/<script[^>]+src="([^"]+)"/g)].map((m) => m[1]),
    ...[...html.matchAll(/<link[^>]+href="([^"]+)"/g)].map((m) => m[1]),
  ];
  assert.ok(refs.length >= 2, "expected at least a stylesheet and a script");
  for (const ref of refs) {
    assert.ok(ref.startsWith("./"), `${ref} is not a relative path, so it is not static`);
    await readFile(path.join(DIST, ref)); // throws if the build missed it
  }
});

test("the module and its glue are both in dist", async () => {
  await readFile(path.join(DIST, "js", "tepsim_wasm.js"));
  await readFile(path.join(DIST, "js", "tepsim_wasm_bg.wasm"));
});

test("nothing on the page can reach the simulation with a clock or a coin", () => {
  // Determinism is a hard invariant: a run is a pure function of its scenario.
  // `performance.now` is allowed and is the exception that proves the rule, so
  // it is named here rather than waved through by a looser pattern.
  const forbidden = [
    [/Math\.random/, "Math.random"],
    [/new Date\b/, "new Date"],
    [/Date\.now/, "Date.now"],
    [/crypto\.getRandomValues/, "crypto.getRandomValues"],
  ];
  for (const [name, source] of Object.entries(scripts)) {
    const code = stripComments(source);
    for (const [pattern, what] of forbidden) {
      assert.ok(
        !pattern.test(code),
        `${name} contains ${what}; a run must be a pure function of its scenario`,
      );
    }
  }
});

test("the page ships no third-party JavaScript", () => {
  // The byte budget in PLAN.org survives because nothing is vendored. If that
  // ever changes it should change deliberately, with the licence carried into
  // dist, not by a copy-paste.
  const external = [...html.matchAll(/(?:src|href)="(https?:)?\/\/[^"]+"/g)];
  assert.deepEqual(
    external.map((m) => m[0]),
    [],
    "index.html loads something off-origin; a static host cannot vouch for it",
  );
});

test("the worker is a module and is loaded relative to its own script", () => {
  // `new Worker` resolves a bare relative string against the document, not
  // against the importing module, so the `import.meta.url` form is what makes
  // the js/ subdirectory work at all. This is a one-line mistake that only
  // shows up as a 404 in a browser console.
  assert.match(
    scripts["app.js"],
    /new Worker\(\s*new URL\("\.\/worker\.js", import\.meta\.url\)/,
  );
  assert.match(scripts["app.js"], /type:\s*"module"/);
});

test("the worker transfers chunk buffers rather than cloning them", () => {
  // Without the transfer list every chunk is copied a second time by the
  // structured clone algorithm. It is invisible, it is the dominant cost of a
  // fast run, and it is one array literal.
  assert.match(scripts["worker.js"], /\[values\.buffer\]/);
});

test("SharedArrayBuffer is not used anywhere", () => {
  // It needs COOP and COEP response headers, which neither GitHub Pages nor a
  // Hugging Face Static Space can set. Using it would trade free static hosting
  // for a copy this app is not big enough to notice.
  for (const [name, source] of Object.entries(scripts)) {
    assert.ok(
      !/SharedArrayBuffer/.test(stripComments(source)),
      `${name} uses SharedArrayBuffer`,
    );
  }
});
