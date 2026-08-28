// A run, in a URL fragment.
//
// The point of this project is that a run is a pure function of its scenario:
// the same description gives bit-identical output on x86-64, aarch64 and wasm.
// That claim is only useful if the description travels. A link is the cheapest
// vehicle there is, and it is what lets an instructor hand out an exact
// scenario rather than a data file, and what lets a bug report say "here" and
// mean it.
//
// # The scenario is one token, and this file cannot read it
//
// Everything a run's output depends on arrives here already serialised, as the
// single string `Scenario.text` produces, and leaves the same way. This file
// does not know what a seed is.
//
// That is the whole design, and it is a correction. Until B-0054a the fields
// were enumerated by hand in three places: the key table here, `startRequest`
// in `app.js`, and `buildScenario` in `worker.js`. A field added to `Scenario`
// reached none of them and nothing failed: the link still opened, the run still
// ran, and it was a different run from the one the link claimed. Now there is
// one serialisation, it is versioned, and `crates/tepsim/tests/scenario_text.rs`
// fails to compile when a field is added to `Scenario` without being handled.
//
// What is left here is the view: which channels are plotted, how big a chunk
// is, and how fast to pace. None of those reach the simulation, which is why
// they are still spelled out, and why a mistyped one can safely fall back.
//
// # The fragment, not the query string
//
// Everything after `#` stays in the browser and is never sent to the server.
// GitHub Pages and a Hugging Face Static Space serve one static page, so a
// query string would be an unused round trip carrying a scenario into somebody
// else's logs. The fragment also changes without a navigation, so the address
// bar can track the panel as it is edited.
//
// The scenario text is written verbatim rather than percent-encoded. It is made
// of `A-Z a-z 0-9 - . _ ~ ; = , :` and `+` only, every one of which a fragment
// carries as itself, and
// `the_text_needs_no_percent_encoding_in_a_url_fragment` in
// `crates/tepsim/tests/scenario_text.rs` is what keeps that true.
//
// # Defaults are omitted
//
// A link for the baseline is `#`. Only parts that differ from the defaults are
// written, so the common link is short and a reader can see at a glance what
// was changed.
//
// # Nothing here trusts its input
//
// A fragment is user input from an untrusted link. The view fields are range
// checked and fall back rather than throwing, because a mistyped link should
// open the baseline plant and not a blank page. The scenario token is handed
// back as it was found; the bindings parse it, and they are strict on purpose,
// so a link naming a field this build does not have is refused with a message
// instead of quietly opening a different run.

/** Short keys, so a link with a channel list stays readable. */
const KEYS = {
  scenario: "q",
  channels: "p",
  chunkSamples: "k",
  speedMultiple: "x",
};

function finiteOr(text, fallback, { min = -Infinity, max = Infinity } = {}) {
  const value = Number(text);
  if (!Number.isFinite(value) || value < min || value > max) return fallback;
  return value;
}

/** A comma-separated list of integers, filtered to a range and de-duplicated. */
function intListOr(text, fallback, min, max) {
  if (text === undefined) return fallback;
  if (text === "") return [];
  const seen = new Set();
  for (const part of text.split(",")) {
    const value = Number(part);
    if (!Number.isInteger(value) || value < min || value > max) continue;
    seen.add(value);
  }
  return [...seen].sort((a, b) => a - b);
}

// `decodeURIComponent` throws on a malformed escape such as a lone `%`, and a
// link is user input. The whole file promises not to throw, so a value that
// cannot be decoded is passed through as it was written and then fails its own
// range check like any other nonsense.
function safeDecode(text) {
  try {
    return decodeURIComponent(text);
  } catch {
    return text;
  }
}

/**
 * Render the parts of `link` that differ from `defaults` as a fragment.
 *
 * The leading `#` is not included, so the caller decides whether an empty
 * result means `#` or an untouched address bar.
 *
 * @param {{scenario: string, channels: number[], chunkSamples: number,
 *          speedMultiple: number}} link
 * @param {object} defaults the same shape, for the baseline
 * @returns {string}
 */
export function encodeLink(link, defaults) {
  const parts = [];
  const put = (key, value) => parts.push(`${key}=${value}`);

  // Verbatim: the scenario text is fragment-safe by construction, and
  // percent-encoding it would turn a readable link into noise.
  if (link.scenario !== defaults.scenario) put(KEYS.scenario, link.scenario);

  const channels = [...link.channels].sort((a, b) => a - b).join(",");
  if (channels !== [...defaults.channels].sort((a, b) => a - b).join(",")) {
    put(KEYS.channels, channels);
  }
  if (link.chunkSamples !== defaults.chunkSamples) {
    put(KEYS.chunkSamples, String(link.chunkSamples));
  }
  if (link.speedMultiple !== defaults.speedMultiple) {
    put(KEYS.speedMultiple, String(link.speedMultiple));
  }

  return parts.join("&");
}

/**
 * Parse a fragment back, falling back to `defaults` field by field. Never
 * throws.
 *
 * The scenario comes back as the text it was written as, not as a parsed
 * object: parsing it needs the bindings, and this file is deliberately free of
 * them so that it can be tested without a wasm build.
 *
 * @param {string} fragment with or without a leading `#`
 * @param {object} defaults
 * @param {number} channelCount how many columns a packed row has
 * @returns {{scenario: string, channels: number[], chunkSamples: number,
 *            speedMultiple: number}}
 */
export function decodeLink(fragment, defaults, channelCount = 54) {
  const link = {
    scenario: defaults.scenario,
    channels: [...defaults.channels],
    chunkSamples: defaults.chunkSamples,
    speedMultiple: defaults.speedMultiple,
  };
  const text = fragment.startsWith("#") ? fragment.slice(1) : fragment;
  if (text === "") return link;

  const raw = new Map();
  for (const pair of text.split("&")) {
    if (pair === "") continue;
    const eq = pair.indexOf("=");
    if (eq < 0) continue;
    // Only the first `=` splits the pair. The scenario text is full of them.
    raw.set(pair.slice(0, eq), safeDecode(pair.slice(eq + 1)));
  }

  const scenario = raw.get(KEYS.scenario);
  if (scenario !== undefined && scenario !== "") link.scenario = scenario;

  // Column 0 is the time and is not plottable.
  link.channels = intListOr(raw.get(KEYS.channels), link.channels, 1, channelCount - 1);

  const chunk = finiteOr(raw.get(KEYS.chunkSamples), defaults.chunkSamples, { min: 1 });
  link.chunkSamples = Math.max(1, Math.round(chunk));

  // Zero is a legal speed and means "as fast as this machine goes".
  link.speedMultiple = finiteOr(raw.get(KEYS.speedMultiple), defaults.speedMultiple, {
    min: 0,
  });

  return link;
}
